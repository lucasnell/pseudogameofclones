#include <vector>

// From http://gallery.rcpp.org/articles/using-bigmemory-with-rcpparmadillo/

// /*
//  To enable the functionality provided by Armadillo's various macros,
//  simply include them before you include the RcppArmadillo headers.
//  */
// // #define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// // [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;

// The following header file provides the definitions for the BigMatrix object
#include <bigmemory/BigMatrix.h>

#include <progress.hpp>  // for the progress bar


// /*
//  Bigmemory now accesses Boost headers from the BH package,
//  we need to make sure we do so as well in this Rcpp::depends comment.
//  Boost headers can generate some warning with the default compilation
//  options for R.  To suppress these, we can enable C++11 mode which gets
//  us 'long long' types.
//  */
//
// [[Rcpp::plugins(cpp11)]]




class OneSummary {
public:
    std::vector<double> group_id;
    double total;
    double n;
    OneSummary() : group_id(), total(), n() {}
    OneSummary(const arma::rowvec& current_row, const arma::uvec& group_cols)
        : group_id(), total(0), n(1) {
        group_id.reserve(group_cols.n_elem);
        for (arma::uword i = 0; i < group_cols.n_elem; i++) {
            group_id.push_back(current_row[group_cols[i]]);
        }
    };
    OneSummary(const arma::rowvec& current_row, const arma::uvec& group_cols,
               const arma::uword& summ_col) : group_id(), total(current_row[summ_col]), n(1) {
        group_id.reserve(group_cols.n_elem);
        for (arma::uword i = 0; i < group_cols.n_elem; i++) {
            group_id.push_back(current_row[group_cols[i]]);
        }
    }
    void operator+=(const OneSummary& other) {
        if (other.group_id.size() != group_id.size()) {
            stop("OneSummary sizes don't match");
        }
        total += other.total;
        n++;
        return;
    }
    // This allows you to use whatever function you want to combine summaries
    void combine(const OneSummary& other, double (*f)(double, double)) {
        if (other.group_id.size() != group_id.size()) {
            stop("OneSummary sizes don't match");
        }
        double new_total = f(total, other.total);
        total = new_total;
        n++;
        return;
    }
    bool operator==(const OneSummary& other) const {
        if (other.group_id.size() != group_id.size()) {
            stop("OneSummary sizes don't match");
        }
        bool out = true;
        for (arma::uword i = 0; i < group_id.size(); i++) {
            if (other.group_id[i] != group_id[i]) {
                out = false;
                break;
            }
        }
        return out;
    }
};



class MeanObj {
public:
        std::vector<OneSummary> summaries;

        void add_row(const arma::rowvec& current_row, const arma::uvec& group_cols,
                     const arma::uword& summ_col, const bool& zeros) {
            OneSummary summ;
            if (zeros) {
                // (Sets total field to zero by default)
                summ = OneSummary(current_row, group_cols);
                if (current_row[summ_col] == 0) summ.total = 1;
            } else summ = OneSummary(current_row, group_cols, summ_col);
            auto iter = std::find(summaries.begin(), summaries.end(), summ);
            if (iter != summaries.end()) {
                (*iter) += summ;
            } else {
                summaries.push_back(summ);
            }
            return;
        }

        void fill_matrix(arma::mat& outmat) {

            if (summaries.size() == 0) {
                outmat.set_size(0, 0);
                return;
            }

            outmat.set_size(summaries.size(), summaries[0].group_id.size() + 1);
            for (arma::uword i = 0; i < outmat.n_rows; i++) {
                outmat(i,0) = summaries[i].total / summaries[i].n;
                for (arma::uword j = 1; j < outmat.n_cols; j++) {
                    outmat(i,j) = summaries[i].group_id[j-1];
                }
            }

            return;
        }

};




//' Mean by groups of columns of simulated aphid abundances.
//'
//' Note: It assumes the input matrix has the following columns (in order):
//' plant, line, date, N, pool, and rep.
//' It also assumes that the matrix is ordered by pool, rep, plant, line, then date.
//'
//'
//' @param pBigMat `big.matrix` object pointer (the `@address` slot!).
//' @param group_cols Vector of columns along which means will be calculated.
//' @param zeros Boolean for whether to compute the proportion of zeros rather than
//'     the mean. Defaults to `FALSE`.
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat grouped_mean(SEXP pBigMat, const bool& by_plant, const bool& by_date,
                       const bool& zeros = false,
                       const bool& show_progress = true) {

    // First we tell Rcpp that the object we've been given is an external
    // pointer.
    XPtr<BigMatrix> xpMat(pBigMat);

    // The 4th column should be N
    arma::uword summ_col = 3;
    std::vector<arma::uword> group_cols_ = {4, 5};
    if (by_plant) group_cols_.push_back(0);
    group_cols_.push_back(1);
    if (by_date) group_cols_.push_back(2);
    if (group_cols_.size() == (xpMat->ncol() - 1)) {
        stop("You're not actually summarizing anything.");
    }
    arma::uvec group_cols(group_cols_);

    // The actual data for the matrix is stored in the matrix() field.
    // This is just a pointer to an array, which is laid out in memory in
    // the column major format. Armadillo matrices are also stored in column
    // major format. We can therefore use the advanced `arma::mat` constructor
    // with `copy_aux_mem` set to `false` to effectively "cast" this memory
    // to an object RcppArmadillo understands.
    //
    // Note that this is an 'unsafe' cast, since we're telling armadillo
    // to use existing memory, rather than to create a new matrix. So we need
    // to be careful that the memory we're telling it to use has the correct
    // dimensions!
    //
    arma::uword type = xpMat->matrix_type();
    if (type != 8) stop("Input matrix is not of type double");

    const arma::Mat<double> M((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(),
                              false);

    if (M.n_rows == 0) stop("empty matrix");

    MeanObj means;

    arma::uword n_ticks = M.n_rows / 1000;
    Progress p(n_ticks, show_progress);

    means.add_row(M.row(0), group_cols, summ_col, zeros);
    for (arma::uword i = 1; i < M.n_rows; i++) {
        Rcpp::checkUserInterrupt();
        means.add_row(M.row(i), group_cols, summ_col, zeros);
        if (i % 1000 == 0) p.increment();
    }

    arma::mat out_matrix;
    means.fill_matrix(out_matrix);

    return out_matrix;

}
