# ifndef __CLONEWARS_MATH_H
# define __CLONEWARS_MATH_H


#include <RcppArmadillo.h>
#include <cmath>
#include <random>

#include <pcg/pcg_random.hpp>   // pcg prng
#include "clonewars_types.hpp"
#include "pcg.hpp"              // runif_ fxns


using namespace Rcpp;



/*
 =====================================================================================
 =====================================================================================
 Beta and truncated normal distributions
 =====================================================================================
 =====================================================================================
 */





/*
 Normal distribution truncated above zero.
 Used for generating `K` bc we never want it to be < 0.
 */
class trunc_normal_distribution {

    double mu;
    double sigma;
    double a_bar;
    double p;

public:

    trunc_normal_distribution() : mu(0), sigma(1), a_bar((0 - mu) / sigma), p(1) {}
    trunc_normal_distribution(const double& mu_, const double& sigma_)
        : mu(mu_), sigma(sigma_), a_bar((0 - mu) / sigma),
          p(R::pnorm5(a_bar, 0, 1, 1, 0)) {}
    trunc_normal_distribution(const trunc_normal_distribution& other)
        : mu(other.mu), sigma(other.sigma), a_bar((0 - mu) / sigma),
          p(R::pnorm5(a_bar, 0, 1, 1, 0)) {}

    trunc_normal_distribution& operator=(const trunc_normal_distribution& other) {
        mu = other.mu;
        sigma = other.sigma;
        a_bar = (0 - mu) / sigma;
        p = R::pnorm5(a_bar, 0, 1, 1, 0);
        return *this;
    }

    double operator()(pcg32& eng) {

        double u = runif_ab(eng, p, 1);

        double x = R::qnorm5(u, 0, 1, 1, 0);
        x = x * sigma + mu;

        return x;
    }


};



/*
 Custom Beta distibution class
 Idea for combining two gammas from https://stackoverflow.com/a/10359049
 Used for generating plant-death-mortality growth-rate modifiers bc they should be
 between 0 and 1.
 */
class beta_distribution {

    std::gamma_distribution<double> X;
    std::gamma_distribution<double> Y;
    double x;
    double y;
    double z;

public:

    beta_distribution() : X(1, 1), Y(1, 1) {}
    beta_distribution(const double& shape1, const double& shape2)
        : X(shape1, 1), Y(shape2, 1) {}
    beta_distribution(const beta_distribution& other)
        : X(other.X), Y(other.Y) {}

    beta_distribution& operator=(const beta_distribution& other) {
        X = other.X;
        Y = other.Y;
        return *this;
    }

    double operator()(pcg32& eng) {

        x = X(eng);
        y = Y(eng);

        z = x / (x + y);

        return z;
    }


};














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

