# ifndef __PSEUDOGAMEOFCLONES_MATH_H
# define __PSEUDOGAMEOFCLONES_MATH_H


#include <RcppArmadillo.h>
#include <cmath>
#include <random>

#include <pcg/pcg_random.hpp>   // pcg prng
#include "pseudogameofclones_types.hpp"
#include "pcg.hpp"              // runif_ab fxn


using namespace Rcpp;








/*
 =====================================================================================
 =====================================================================================
 Logit and inverse logit
 =====================================================================================
 =====================================================================================
 */




inline void logit__(const double& p, double& out) {
    out = std::log(p / (1-p));
    return;
}
inline void inv_logit__(const double& a, double& out) {
    out = 1 / (1 + std::exp(-a));
    return;
}



/*
 =====================================================================================
 =====================================================================================
 Leslie matrix
 =====================================================================================
 =====================================================================================
 */


inline void leslie_matrix__(const arma::uvec& instar_days,
                            const double& surv_juv,
                            const arma::vec& surv_adult,
                            const arma::vec& repro,
                            arma::mat& out) {

    arma::vec tmp;
    uint32 juv_time = arma::accu(instar_days.head(instar_days.n_elem - 1));
    // adult time is last age with >0 survival to next stage, plus one
    // (obviously any stage having 0 survival before this could be a problem,
    //  depending on whether you're trying to test something...)
    uint32 adult_time = arma::as_scalar(arma::find(surv_adult, 1, "last")) + 2;
    uint32 n_stages = adult_time + juv_time;
    // Age-specific survivals
    tmp = arma::vec(n_stages - 1);
    tmp.head(juv_time).fill(surv_juv);
    tmp.tail(adult_time-1) = surv_adult.head(adult_time-1);
    out = arma::diagmat(tmp, -1);
    // Age-specific fecundities
    uint32 n_adult_repos = std::min(repro.n_elem, adult_time);
    out(0, arma::span(juv_time, juv_time + n_adult_repos - 1)) =
        repro.head(n_adult_repos).t();

    return;
}







/*
 =====================================================================================
 =====================================================================================
 Stable age distribution
 =====================================================================================
 =====================================================================================
 */

// This computes the "stable age distribution" from the Leslie matrix, which is
// the proportion of different classes that is required for the population to grow
// exponentially
inline void sad_leslie__(const arma::mat& L, arma::vec& out) {

    arma::cx_vec r_cx;
    arma::cx_mat SAD;

    arma::eig_gen(r_cx, SAD, L);

    arma::vec r = arma::abs(r_cx);

    double rmax = arma::max(r);

    arma::cx_mat SADdist = SAD.cols(arma::find(r == rmax));
    arma::cx_double all_SAD = arma::accu(SADdist);
    SADdist /= all_SAD;
    SADdist.resize(SADdist.n_elem, 1);

    out = arma::real(SADdist);

    return;
}





#endif

