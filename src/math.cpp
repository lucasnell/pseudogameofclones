#include <RcppArmadillo.h>
#include <cmath>

#include "clonewars_types.hpp"

using namespace Rcpp;

/*
 =====================================================================================
 =====================================================================================
 Logit and inverse logit
 =====================================================================================
 =====================================================================================
 */

/*
 --------------
 C++ versions
 --------------
 */

void logit__(const arma::vec& p, arma::vec& out) {
    out = arma::log(p / (1-p));
    return;
}
void inv_logit__(const arma::vec& a, arma::vec& out){
    out = 1 / (1 + arma::exp(-a));
    return;
}

/*
 --------------
 R versions
 --------------
*/

//' Logit and inverse logit functions.
//'
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector logit(NumericVector p) {
    arma::vec out;
    logit__(as<arma::vec>(p), out);
    return wrap(out);
}
//' @describeIn logit
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector inv_logit(NumericVector a){
    arma::vec out;
    inv_logit__(as<arma::vec>(a), out);
    return wrap(out);
}



/*
 =====================================================================================
 =====================================================================================
 Leslie matrix
 =====================================================================================
 =====================================================================================
 */

/*
 --------------
 C++ versions
 --------------
 */

void leslie_matrix__(const arma::uvec& instar_days, const double& surv_juv,
                          const arma::vec& surv_adult, const arma::vec& repro,
                          arma::mat& out) {
    uint32 n_stages = arma::accu(instar_days);
    arma::vec tmp;
    uint32 juv_time = arma::accu(instar_days(arma::span(0, (instar_days.n_elem - 2))));
    // Age-specific survivals
    tmp = arma::vec(n_stages - 1);
    tmp.head(juv_time).fill(surv_juv);
    tmp.tail(n_stages-juv_time-1) = surv_adult(arma::span(0,(n_stages-juv_time-2)));
    out = arma::diagmat(tmp, -1);
    // Age-specific fecundities
    out(0, arma::span(juv_time, juv_time + instar_days(instar_days.n_elem - 1) - 1)) =
        repro(arma::span(0, instar_days(instar_days.n_elem - 1) - 1)).t();

    return;
}



/*
 --------------
 R versions
 --------------
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

// This computes the "stable age distribution" from the Leslie matrix, which is
// the proportion of different classes that is required for the population to grow
// exponentially
void sad_leslie__(const arma::mat& L, arma::vec& out) {

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
