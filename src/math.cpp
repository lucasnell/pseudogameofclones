#include <RcppArmadillo.h>
#include <cmath>

#include "math.hpp"
#include "clonewars_types.hpp"

using namespace Rcpp;

/*
 =====================================================================================
 =====================================================================================
 Logit and inverse logit
 =====================================================================================
 =====================================================================================
 */


//' Logit and inverse logit functions.
//'
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector logit(NumericVector p) {
    NumericVector out(p.size());
    for (uint32 i = 0; i < p.size(); i++) {
        logit__(p[i], out[i]);
    }
    return out;
}
//' @describeIn logit
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector inv_logit(NumericVector a){
    NumericVector out(a.size());
    for (uint32 i = 0; i < a.size(); i++) {
        inv_logit__(a[i], out[i]);
    }
    return out;
}



/*
 =====================================================================================
 =====================================================================================
 Leslie matrix
 =====================================================================================
 =====================================================================================
 */


//' Create Leslie matrix from aphid info
//'
//' @param instar_days Integer vector of the number of stages (days) per aphid instar.
//' @param surv_juv Single numeric of daily juvenile survival.
//' @param surv_adult Numeric vector of aphid adult survivals by stage.
//' @param repro Numeric vector of aphid reproductive rates by stage.
//'
//'
//' @export
//[[Rcpp::export]]
NumericMatrix leslie_matrix(IntegerVector instar_days, const double& surv_juv,
                            NumericVector surv_adult, NumericVector repro) {
    arma::mat LL;
    leslie_matrix__(as<arma::uvec>(instar_days), surv_juv,
                    as<arma::vec>(surv_adult),
                    as<arma::vec>(repro),
                    LL);
    return wrap(LL);
}




/*
 Carrying capacity for patch.
 It depends on the Leslie matrix for each line's apterous aphids.
 I'm assuming apterous ones drive the carrying capacity, rather than alates, because
 they should be much more numerous.
 */
//' Calculate carrying capacity from aphid and patch info
//'
//' @param apterous Leslie matrix for apterous aphids of this line
//' @param alates Leslie matrix for alates of this line
//' @param alate_prop The proportion of new aphids that are alates
//' @param disp_prop The proportion of alates that disperse away from the patch.
//' @param disp_mort Mortality rate for dispersing alates.
//' @param disp_start Index for the age at which alates disperse (0-indexed).
//' @param K The "carrying capacity" for this patch.
//'
//'
//' @export
//[[Rcpp::export]]
double carrying_capacity(const arma::mat& apterous,
                         const arma::mat& alates,
                         const double& alate_prop,
                         const double& disp_prop,
                         const double& disp_mort,
                         const uint32& disp_start,
                         const double& K) {

    arma::mat L;

    combine_leslies(L, apterous, alates, alate_prop, disp_prop,
                    disp_mort, disp_start);

    arma::cx_vec eigval = arma::eig_gen( L );
    double ev = eigval.max().real();
    double cc = (ev - 1) * K;

    return cc;
}




/*
 =====================================================================================
 =====================================================================================
 Stable age distribution
 =====================================================================================
 =====================================================================================
 */



//' Compute the "stable age distribution" from a Leslie matrix
//'
//' @param leslie Leslie matrix of the population.
//'
//' @export
//'
//' @return Vector with the stable age distribution given the Leslie matrix.
//'
//[[Rcpp::export]]
NumericVector sad_leslie(NumericMatrix leslie) {

    arma::mat L = as<arma::mat>(leslie);
    arma::vec out;

    sad_leslie__(L, out);

    return wrap(out);
}
