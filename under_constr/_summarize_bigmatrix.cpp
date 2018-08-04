#include <vector>

// http://gallery.rcpp.org/articles/using-bigmemory-with-rcpparmadillo/

/*
 To enable the functionality provided by Armadillo's various macros,
 simply include them before you include the RcppArmadillo headers.
 */
// #define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory, RcppProgress)]]

using namespace Rcpp;

// The following header file provides the definitions for the BigMatrix object
#include <bigmemory/BigMatrix.h>
#include <progress.hpp>  // for the progress bar

/*
 Bigmemory now accesses Boost headers from the BH package,
 we need to make sure we do so as well in this Rcpp::depends comment.
 Boost headers can generate some warning with the default compilation
 options for R.  To suppress these, we can enable C++11 mode which gets
 us 'long long' types.
 */

// [[Rcpp::plugins(cpp11, openmp)]]

#ifdef _OPENMP
#include <omp.h>
#endif



/*
 Because I know that the matrix is ordered by pool, rep, plant, line, then date, I can
 create a tree, with the latter categories inside the earlier ones.
 */

struct Group {

    double value;
    arma::uword start;
    arma::uword size;

    Group() : value(), size() {}
    Group(const double& value_) : value(value_), size(1) {}
    Group(const double& value_, const arma::uword& size_)
        : value(value_), size(size_) {}
};

template <typename T>
struct GroupWVec : public Group {

    std::vector<T> inner;

    GroupWVec() : Group(), inner() {}
    GroupWVec(const double& value_) : Group(value_), inner() {}
    GroupWVec(const double& value_, const arma::uword& size_)
        : Group(value_, size_), inner() {}
    T operator[](const arma::uword& idx) const {
        return inner[idx];
    }
    T& operator[](const arma::uword& idx) {
        return inner[idx];
    }
};

typedef Group Line;
typedef GroupWVec<Line> Plant;
typedef GroupWVec<Plant> Rep;
typedef GroupWVec<Rep> Pool;

struct GroupTree {
    std::vector<Pool> pools;
    GroupTree();
    GroupTree(const arma::uword& n_pools) : pools(n_pools) {}
    Pool operator[](const arma::uword& idx) const {
        return pools[idx];
    }
    Pool& operator[](const arma::uword& idx) {
        return pools[idx];
    }
};




// pool_sizes is the # lines present in the pool
// ordered by pool, rep, plant, line, then date.
// [[Rcpp::export]]
SEXP make_group_tree(SEXP pBigMat,
                     const std::vector<arma::uword>& pool_sizes,
                     const arma::uword& n_reps,
                     const arma::uword& n_plants,
                     const arma::uword& max_t) {

    arma::uword n_dates = max_t + 1;

    // First we tell Rcpp that the object we've been given is an external
    // pointer.
    XPtr<BigMatrix> xpMat(pBigMat);
    arma::uword type = xpMat->matrix_type();
    if (type != 8) stop("Input matrix is not of type double");
    const arma::Mat<double> M((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(),
                              false);
    if (M.n_rows == 0) stop("empty matrix");

    arma::uword n_pools = pool_sizes.size();

    XPtr<GroupTree> tree_ptr(new GroupTree(n_pools));
    GroupTree& tree(*tree_ptr);

    arma::uword M_ind = 0;

    // M is ordered by pool, rep, plant, line, then date
    // These objects are grouped as GroupTree > Pool > Rep > Plant > Line

    for (arma::uword i = 0; i < n_pools; i++) {

        Rcpp::checkUserInterrupt();

        if (M_ind >= M.n_rows) {
            Rcout << M_ind << std::endl;
            stop("M_ind too large");
        }

        const arma::uword& n_lines(pool_sizes[i]);

        Pool& pool(tree[i]);
        pool.value = M(M_ind, 4);
        pool.start = M_ind;
        pool.size = n_reps * n_plants * n_lines * n_dates;
        pool.inner.resize(n_reps);

        for (arma::uword j = 0; j < n_reps; j++) {

            Rep& rep(pool[j]);
            rep.value = M(M_ind, 5);
            rep.start = M_ind;
            rep.size = n_plants * n_lines * n_dates;
            rep.inner.resize(n_plants);

            for (arma::uword k = 0; k < n_plants; k++) {

                Plant& plant(rep[k]);
                plant.value = M(M_ind, 0);
                plant.start = M_ind;
                plant.size = n_lines * n_dates;
                plant.inner.resize(n_lines);

                for (arma::uword l = 0; l < n_lines; l++) {

                    Line& line(plant[l]);
                    line.value = M(M_ind, 1);
                    line.start = M_ind;
                    line.size = n_dates;

                    M_ind += n_dates;

                }

            }

        }

    }


    return tree_ptr;
}


// // [[Rcpp::export]]
// arma::mat grouped_mean_tree(SEXP pBigMat,
//                             SEXP pTree,
//                             const bool& by_plant,
//                             const bool& by_date,
//                             const bool& show_progress = true) {
//
//     if (by_plant && by_date) stop("You're not actually summarizing anything.");
//
//     XPtr<BigMatrix> xpMat(pBigMat);
//     arma::uword type = xpMat->matrix_type();
//     if (type != 8) stop("Input matrix is not of type double");
//     const arma::Mat<double> M((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(),
//                               false);
//     if (M.n_rows == 0) stop("empty matrix");
//
//     XPtr<GroupTree> tree_ptr(pTree);
//     GroupTree& tree(*tree_ptr);
//
//     // The 4th column should be N
//     // Column orders: plant, line, date, N, pool, rep
//     arma::uword summ_col = 3;
//     std::vector<arma::uword> group_cols_ = {4, 5};
//     if (by_plant) group_cols_.push_back(0);
//     group_cols_.push_back(1);
//     if (by_date) group_cols_.push_back(2);
//
//     arma::uvec group_cols(group_cols_);
//
//
//
//
// }
