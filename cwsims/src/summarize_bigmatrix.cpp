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
#include <bigmemory/MatrixAccessor.hpp>

#include <progress.hpp>  // for the progress bar

#include <cereal/types/base_class.hpp> // for the cereal::base_class
#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>






/*
 =======================================================================================
 =======================================================================================

 Making, saving, loading tree to describe groups

 =======================================================================================
 =======================================================================================
 */

/*
 Because I know that the matrix is ordered by pool, rep, line, plant, then date, I can
 create a tree, with the latter categories inside the earlier ones.
 */

typedef uint_fast32_t uint32;

struct Group {

    double value;
    uint32 begin;
    uint32 n_rows;

    Group() : value(), n_rows() {}
    Group(const double& value_) : value(value_), n_rows(1) {}
    Group(const double& value_, const uint32& n_rows_)
        : value(value_), n_rows(n_rows_) {}

    bool operator<(const Group& other) {
        return value < other.value;
    }
    bool operator>(const Group& other) {
        return value > other.value;
    }
    bool operator==(const Group& other) {
        return value == other.value;
    }

    template <class Archive>
    void serialize(Archive & arch) {
        arch(value, begin, n_rows);
        return;
    }
};

template <typename T>
struct GroupWVec : public Group {

    std::vector<T> inner;

    GroupWVec() : Group(), inner() {}
    GroupWVec(const double& value_) : Group(value_), inner() {}
    GroupWVec(const double& value_, const uint32& n_rows_)
        : Group(value_, n_rows_), inner() {}
    T operator[](const uint32& idx) const {
        return inner[idx];
    }
    T& operator[](const uint32& idx) {
        return inner[idx];
    }

    template <class Archive>
    void serialize(Archive & arch) {
        // We pass this cast to the base type for each base type we
        // need to serialize.  Do this instead of calling serialize functions
        // directly
        arch(cereal::base_class<Group>(this), inner);
    }
};

typedef Group Plant;
typedef GroupWVec<Plant> Line;
typedef GroupWVec<Line> Rep;
typedef GroupWVec<Rep> Pool;

struct GroupTree {
    std::vector<Pool> pools;
    GroupTree() {}
    GroupTree(const uint32& n_pools) : pools(n_pools) {}
    Pool operator[](const uint32& idx) const {
        return pools[idx];
    }
    Pool& operator[](const uint32& idx) {
        return pools[idx];
    }
    uint32 size() const noexcept {
        return pools.size();
    }

    template <class Archive>
    void save(Archive& arch) const {
        arch(pools);
    }
    template <class Archive>
    void load(Archive& arch) {
        arch(pools);
    }
};




// pool_sizes is the # lines present in the pool
// ordered by pool, rep, line, plant, then date.
// [[Rcpp::export]]
SEXP make_group_tree(SEXP pBigMat,
                     const std::vector<uint32>& pool_sizes,
                     const uint32& n_reps,
                     const uint32& n_plants,
                     const uint32& max_t) {

    uint32 n_dates = max_t + 1;

    // First we tell Rcpp that the object we've been given is an external
    // pointer.
    XPtr<BigMatrix> xpMat(pBigMat);
    uint32 type = xpMat->matrix_type();
    if (type != 8) stop("Input matrix is not of type double");
    const arma::Mat<double> M((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(),
                              false);
    if (M.n_rows == 0) stop("empty matrix");

    uint32 n_pools = pool_sizes.size();

    XPtr<GroupTree> tree_ptr(new GroupTree(n_pools));
    GroupTree& tree(*tree_ptr);

    uint32 M_ind = 0;

    // M is ordered by pool, rep, line, plant, then date
    // These objects are grouped as GroupTree > Pool > Rep > Line > Plant

    for (uint32 i = 0; i < n_pools; i++) {

        Rcpp::checkUserInterrupt();

        const uint32& n_lines(pool_sizes[i]);

        Pool& pool(tree[i]);
        pool.value = M(M_ind, 4);
        pool.begin = M_ind;
        pool.n_rows = n_reps * n_lines * n_plants * n_dates;
        pool.inner.resize(n_reps);

        for (uint32 j = 0; j < n_reps; j++) {

            Rep& rep(pool[j]);
            rep.value = M(M_ind, 5);
            rep.begin = M_ind;
            rep.n_rows = n_lines * n_plants * n_dates;
            rep.inner.resize(n_lines);

            for (uint32 k = 0; k < n_lines; k++) {

                Line& line(rep[k]);
                line.value = M(M_ind, 1);
                line.begin = M_ind;
                line.n_rows = n_plants * n_dates;
                line.inner.resize(n_plants);

                for (uint32 l = 0; l < n_plants; l++) {

                    Plant& plant(line[l]);
                    plant.value = M(M_ind, 0);
                    plant.begin = M_ind;
                    plant.n_rows = n_dates;

                    M_ind += n_dates;

                }

            }

        }

    }


    return tree_ptr;
}


// [[Rcpp::export]]
int save_group_tree(SEXP pTree, const std::string& filename) {

    XPtr<GroupTree> tree_xptr(pTree);
    GroupTree& tree(*tree_xptr);

    std::ofstream os(filename.c_str(), std::ios::binary);
    cereal::BinaryOutputArchive archive(os);

    archive(tree);

    return 0;
}


// [[Rcpp::export]]
SEXP load_group_tree(const std::string& filename) {

    std::ifstream os(filename.c_str(), std::ios::binary);
    cereal::BinaryInputArchive archive(os);

    XPtr<GroupTree> tree_xptr(new GroupTree());
    GroupTree& tree(*tree_xptr);

    archive(tree);

    return tree_xptr;
}




// [[Rcpp::export]]
std::vector<double> view_tree_values(SEXP pTree, const std::string& level) {

    XPtr<GroupTree> tree_xptr(pTree);
    GroupTree& tree(*tree_xptr);

    std::vector<double> out;

    if (level == "pool") {
        out.reserve(tree.pools.size());
        for (const Pool& pool : tree.pools) {
            out.push_back(pool.value);
        }
    } else if (level == "rep") {
        out.reserve(tree.pools.size() * tree.pools[0].inner.size());
        for (const Pool& pool : tree.pools) {
            for (const Rep& rep : pool.inner) {
                out.push_back(rep.value);
            }
        }
    } else stop("Not programmed");

    return out;
}





/*
 =====================================================================================
 =====================================================================================

 Mean using group tree

 =====================================================================================
 =====================================================================================
 */

double sum_M(const arma::mat& M, const uint32& begin, const uint32& n_rows,
             const uint32& summ_col) {
    return arma::accu(M(arma::span(begin, begin + n_rows - 1), arma::span(summ_col)));
}
double zeros_M(const arma::mat& M, const uint32& begin, const uint32& n_rows,
               const uint32& summ_col) {
    double n_zeros = 0;
    for (uint32 j = begin; j < begin + n_rows; j++) {
        if (M(j, summ_col) == 0) n_zeros++;
    }
    return n_zeros;
}


double sum_M_by_date(const arma::mat& M, const uint32& j,
                     const uint32& summ_col, const Line& line) {
    double sum_ = 0;
    for (const Plant& plant : line.inner) {
        sum_ += M(plant.begin + j, summ_col);
    }
    return sum_;
}
double zeros_M_by_date(const arma::mat& M, const uint32& j,
                       const uint32& summ_col, const Line& line) {
    double n_zeros = 0;
    for (const Plant& plant : line.inner) {
        if (M(plant.begin + j, summ_col) == 0) n_zeros++;
    }
    return n_zeros;
}




// [[Rcpp::export]]
arma::mat grouped_mean(SEXP pBigMat,
                       SEXP pTree,
                       const bool& by_plant,
                       const bool& by_date,
                       const bool& zeros = false,
                       const bool& show_progress = true) {

    if (by_plant && by_date) stop("You're not actually summarizing anything.");

    XPtr<BigMatrix> xpMat(pBigMat);
    uint32 type = xpMat->matrix_type();
    if (type != 8) stop("Input matrix is not of type double");
    const arma::Mat<double> M((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(),
                              false);
    if (M.n_rows == 0) stop("empty matrix");

    XPtr<GroupTree> tree_ptr(pTree);
    GroupTree& tree(*tree_ptr);

    // The 4th column should be N
    // Column orders: plant, line, date, N, pool, rep
    uint32 summ_col = 3;

    // Construct output matrix:
    uint32 n_rows = 0;
    uint32 n_cols = 4;
    for (const Pool& pool : tree.pools) {
        for (const Rep& rep : pool.inner) {
            n_rows += rep.inner.size();
        }
    }
    if (by_plant) {
        const Line& line(tree[0][0][0]);
        n_rows *= line.inner.size();
        n_cols++;
    }
    if (by_date) {
        const Plant& plant(tree[0][0][0][0]);
        n_rows *= plant.n_rows;
        n_cols++;
    }

    arma::mat out(n_rows, n_cols);

    Progress p(n_rows, show_progress);

    // Index for which row in out you're referring to:
    uint32 i = 0;

    /*
    ---------
    Fill output matrix:
    ---------
    */

    if (!by_plant && !by_date) {

        std::function<double(const arma::mat&,
                             const uint32&,
                             const uint32&,
                             const uint32&)> sumf;
        if (zeros) {
            sumf = zeros_M;
        } else sumf = sum_M;

        arma::rowvec out_row(n_cols);
        for (const Pool& pool : tree.pools) {
            Rcpp::checkUserInterrupt();
            out_row(0) = pool.value;
            for (const Rep& rep : pool.inner) {
                out_row(1) = rep.value;
                for (const Line& line : rep.inner) {
                    out_row(2) = line.value;
                    double sum_ = sumf(M, line.begin, line.n_rows, summ_col);
                    double mean_ = sum_ / line.n_rows;
                    out_row.tail(1).fill(mean_);
                    out.row(i) = out_row;
                    i++;
                    p.increment();
                }
            }
        }

    } else if (by_plant && !by_date) {

        std::function<double(const arma::mat&,
                             const uint32&,
                             const uint32&,
                             const uint32&)> sumf;
        if (zeros) {
            sumf = zeros_M;
        } else sumf = sum_M;

        arma::rowvec out_row(n_cols);
        for (const Pool& pool : tree.pools) {
            Rcpp::checkUserInterrupt();
            out_row(0) = pool.value;
            for (const Rep& rep : pool.inner) {
                out_row(1) = rep.value;
                for (const Line& line : rep.inner) {
                    out_row(2) = line.value;
                    for (const Plant& plant : line.inner) {
                        out_row(3) = plant.value;
                        double sum_ = sumf(M, plant.begin, plant.n_rows, summ_col);
                        double mean_ = sum_ / plant.n_rows;
                        out_row.tail(1).fill(mean_);
                        out.row(i) = out_row;
                        i++;
                        p.increment();
                    }
                }
            }
        }

    } else if (!by_plant && by_date) {

        std::function<double(const arma::mat&,
                             const uint32&,
                             const uint32&,
                             const Line&)> sumf;
        if (zeros) {
            sumf = zeros_M_by_date;
        } else sumf = sum_M_by_date;

        for (const Pool& pool : tree.pools) {
            Rcpp::checkUserInterrupt();
            for (const Rep& rep : pool.inner) {
                for (const Line& line : rep.inner) {
                    double n_plants = line.inner.size();
                    uint32 n_dates = line[0].n_rows;
                    for (uint32 j = 0; j < n_dates; j++) {
                        // Fill values:
                        out(i+j, 0) = pool.value;
                        out(i+j, 1) = rep.value;
                        out(i+j, 2) = line.value;
                        out(i+j, 3) = j;
                        out(i+j, 4) = sumf(M, j, summ_col, line);
                        out(i+j, 4) /= n_plants;
                    }
                    i += n_dates;
                    p.increment(n_dates);
                }
            }
        }

    } else stop("This should never happen");

    return out;
}


