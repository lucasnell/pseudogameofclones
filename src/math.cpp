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
