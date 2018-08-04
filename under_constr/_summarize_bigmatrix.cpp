#include <vector>
#include <string>
#include <fstream>

// http://gallery.rcpp.org/articles/using-bigmemory-with-rcpparmadillo/

/*
 To enable the functionality provided by Armadillo's various macros,
 simply include them before you include the RcppArmadillo headers.
 */
// #define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory, RcppProgress, Rcereal)]]

using namespace Rcpp;

// The following header file provides the definitions for the BigMatrix object
#include <bigmemory/BigMatrix.h>
#include <progress.hpp>  // for the progress bar

#include <cereal/types/base_class.hpp> // for the cereal::base_class
#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>


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

typedef uint_fast32_t uint32;

struct Group {

    double value;
    uint32 start;
    uint32 size;

    Group() : value(), size() {}
    Group(const double& value_) : value(value_), size(1) {}
    Group(const double& value_, const uint32& size_)
        : value(value_), size(size_) {}

    template <class Archive>
    void serialize(Archive & arch) {
        arch(value, start, size);
        return;
    }
};

template <typename T>
struct GroupWVec : public Group {

    std::vector<T> inner;

    GroupWVec() : Group(), inner() {}
    GroupWVec(const double& value_) : Group(value_), inner() {}
    GroupWVec(const double& value_, const uint32& size_)
        : Group(value_, size_), inner() {}
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

typedef Group Line;
typedef GroupWVec<Line> Plant;
typedef GroupWVec<Plant> Rep;
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
// ordered by pool, rep, plant, line, then date.
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

    // M is ordered by pool, rep, plant, line, then date
    // These objects are grouped as GroupTree > Pool > Rep > Plant > Line

    for (uint32 i = 0; i < n_pools; i++) {

        Rcpp::checkUserInterrupt();

        if (M_ind >= M.n_rows) {
            Rcout << M_ind << std::endl;
            stop("M_ind too large");
        }

        const uint32& n_lines(pool_sizes[i]);

        Pool& pool(tree[i]);
        pool.value = M(M_ind, 4);
        pool.start = M_ind;
        pool.size = n_reps * n_plants * n_lines * n_dates;
        pool.inner.resize(n_reps);

        for (uint32 j = 0; j < n_reps; j++) {

            Rep& rep(pool[j]);
            rep.value = M(M_ind, 5);
            rep.start = M_ind;
            rep.size = n_plants * n_lines * n_dates;
            rep.inner.resize(n_plants);

            for (uint32 k = 0; k < n_plants; k++) {

                Plant& plant(rep[k]);
                plant.value = M(M_ind, 0);
                plant.start = M_ind;
                plant.size = n_lines * n_dates;
                plant.inner.resize(n_lines);

                for (uint32 l = 0; l < n_lines; l++) {

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
//     uint32 type = xpMat->matrix_type();
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
//     uint32 summ_col = 3;
//     std::vector<uint32> group_cols_ = {4, 5};
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
