
#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <random>  // normal_distribution
#include <pcg/pcg_random.hpp>   // pcg prng

#ifndef RCPPTHREAD_OVERRIDE_COUT
#define RCPPTHREAD_OVERRIDE_COUT 1    // std::cout override
#endif
#ifndef RCPPTHREAD_OVERRIDE_CERR
#define RCPPTHREAD_OVERRIDE_CERR 1    // std::cerr override
#endif
// #ifndef RCPPTHREAD_OVERRIDE_THREAD
// #define RCPPTHREAD_OVERRIDE_THREAD 1  // std::thread override
// #endif
#include <RcppThread.h>         // multithreading


#include "pseudogameofclones_types.hpp"  // integer types
#include "pcg.hpp"              // runif_01, seed_rng functions
#include "util.hpp"      // thread_check



using namespace Rcpp;


// Convert from 1D to 2D:
void to_2d(uint32& x, uint32& y, const uint32& k,
           const int& x_size, const int& y_size) {
    x = k - y_size * (k / y_size);
    y = k / y_size;
    return;
}
// Convert from 2D to 1D:
void to_1d(uint32& k, const uint32& x, const uint32& y,
           const int& x_size, const int& y_size) {
    k = (y * x_size + x);
    return;
}



//'
//' @noRd
//'
//' @export
//'
//[[Rcpp::export]]
DataFrame target_type_sims2(int x_size,
                            int y_size,
                            const arma::vec& probs,
                            const double& k,
                            const double& sigma2) {

    if (x_size < 2) stop("x_size must be >= 2");
    if (y_size < 2) stop("y_size must be >= 2");
    if (probs.n_elem < 2) stop("length(probs) cannot be < 2");
    if (arma::any(probs <= 0)) stop("probs cannot contain values <= 0");
    if (arma::accu(probs) > 1) stop("sum(probs) cannot > 1");
    if (k < 1) stop("k must be > 1");
    if (sigma2 < 0) stop("sigma2 must be > 0");

    // I'm reducing these by 1 bc I want to sample from 1 to floor(x_size-1)
    // and floor(y_size-1). See description in docs above for why.
    x_size--;
    y_size--;

    // Random number generator:
    pcg32 eng;
    seed_pcg(eng);

    // Normal distribution:
    std::normal_distribution<double> norm_distr(0, 1);


    uint32 n = x_size * y_size;



    // Create covariance matrix with spatial weighting:
    uint32 xi, yi, xj, yj;
    double d, dx, dy;
    arma::mat dmat(n, n, arma::fill::none);
    for (uint32 i = 0; i < n-1; i++) {
        to_2d(xi, yi, i, x_size, y_size);
        for (uint32 j = i+1; j < n; j++) {
            to_2d(xj, yj, j, x_size, y_size);
            dx = static_cast<double>(xi - xj) * static_cast<double>(xi - xj);
            dy = static_cast<double>(yi - yj) * static_cast<double>(yi - yj);
            d = std::sqrt(dx + dy);
            dmat(i,j) = 1 / std::pow(d, k);
            dmat(j,i) = dmat(i,j);
        }
    }
    dmat.diag().fill(sigma2);

    arma::mat iD;
    try {
        iD = arma::chol(dmat);
    } catch(const std::runtime_error& re) {
        std::string err_msg = static_cast<std::string>(re.what());
        if (err_msg == "chol(): decomposition failed") {
            std::string err_msg_out = "Choleski decomposition failed. ";
            err_msg_out += "Increasing  `sigma` can remedy this.";
            throw(Rcpp::exception(err_msg_out.c_str(), false));
        } else {
            std::string err_msg_out = "Runtime error: \n" + err_msg;
            throw(Rcpp::exception(err_msg_out.c_str(), false));
        }
    } catch(const std::exception& ex) {
        std::string err_msg = static_cast<std::string>(ex.what());
        stop("Error occurred: \n" + err_msg);
    } catch(...) {
        stop("Unknown failure occurred.");
    }

    iD = iD.t();

    arma::vec N(n, arma::fill::none);
    for (uint32 i = 0; i < n; i++) N(i) = norm_distr(eng);

    N = iD * N;

    // // Potentially stabler (but slower) way of generating `N`
    // // (as is done in MASS::mvrnorm):
    // arma::vec eigval;
    // arma::mat eigvec;
    // eig_sym(eigval, eigvec, dmat);
    // double thresh = -1e-6 * std::abs(eigval(0U));
    // if (!arma::all(eigval >= thresh)) {
    //     std::string err_msg_out = "Distance matrix is not positive definite. ";
    //     err_msg_out += "Increasing  `sigma` can remedy this.";
    //     stop(err_msg_out.c_str());
    // }
    // arma::mat ev_mat(n, n, arma::fill::zeros);
    // arma::vec N(n, arma::fill::none);
    // for (uint32 i = 0; i < n; i++) {
    //     if (eigval(i) > 0) ev_mat(i,i) = std::sqrt(eigval(i));
    //     N(i) = norm_distr(eng);
    // }
    // N = eigvec * ev_mat * N;


    // use `pnorm` (done manually below) for normal (~N(0,diag(dmat))) to
    // convert from normal distributed data to uniform, then use cumulative
    // probabilities of `props` to generate categorical data.
    std::vector<uint32> type;
    type.reserve(n);
    arma::vec cs_probs = arma::cumsum(probs);
    double u;
    uint32 t;
    uint32 nt = probs.n_elem;
    for (uint32 i = 0; i < n; i++) {
        u = 0.5 * erfc(-1 * N(i) / (std::sqrt(dmat(i,i)) * M_SQRT2));
        t = 0;
        while (t < nt && cs_probs(t) < u) t++;
        type.push_back(t+1U);
    }

    // x and y for each type:
    std::vector<int> x;
    std::vector<int> y;
    x.reserve(n);
    y.reserve(n);
    for (uint32 i = 0; i < n; i++) {
        to_2d(xi, yi, i, x_size, y_size);
        x.push_back(xi);
        y.push_back(yi);
    }


    DataFrame out = DataFrame::create(_["x"] = x, _["y"] = y, _["type"] = type);

    out.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});

    return out;

}
