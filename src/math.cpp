#include <RcppArmadillo.h>
#include <cmath>

#include "clonewars_types.hpp"
#include "math.hpp"

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

arma::vec logit__(const arma::vec& p) {
    arma::vec out = arma::log(p / (1-p));
    return out;
}
arma::vec inv_logit__(const arma::vec& a){
    arma::vec out = 1 / (1 + arma::exp(-a));
    return out;
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
    arma::vec out = logit__(as<arma::vec>(p));
    return wrap(out);
}
//' @describeIn logit
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector inv_logit(NumericVector a){
    arma::vec out = inv_logit__(as<arma::vec>(a));
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

arma::mat leslie_matrix__(const arma::uvec& instar_days, const double& surv_juv,
                          const arma::vec& surv_adult, const arma::vec& repro) {
    uint32 n_stages = arma::accu(instar_days);
    arma::vec tmp;
    arma::mat LL;
    uint32 juv_time = arma::accu(instar_days(arma::span(0, (instar_days.n_elem - 2))));
    // Age-specific survivals
    tmp = arma::vec(n_stages - 1);
    tmp.head(juv_time).fill(surv_juv);
    tmp.tail(n_stages-juv_time-1) = surv_adult(arma::span(0,(n_stages-juv_time-2)));
    LL = arma::diagmat(tmp, -1);
    // Age-specific fecundities
    LL(0, arma::span(juv_time, juv_time + instar_days(instar_days.n_elem - 1) - 1)) =
        repro(arma::span(0, instar_days(instar_days.n_elem - 1) - 1)).t();

    return LL;
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

    arma::mat LL = leslie_matrix__(as<arma::uvec>(instar_days), surv_juv,
                                   as<arma::vec>(surv_adult),
                                   as<arma::vec>(repro));

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
// Used for constructing X_0 for a aphid_const class
arma::vec leslie_sad__(const arma::mat& L) {

    arma::cx_vec r_cx;
    arma::cx_mat SAD;

    arma::eig_gen(r_cx, SAD, L);

    arma::vec r = arma::abs(r_cx);

    double rmax = arma::max(r);

    arma::cx_mat SADdist = SAD.cols(arma::find(r == rmax));
    arma::cx_double all_SAD = arma::accu(SADdist);
    SADdist /= all_SAD;
    SADdist.resize(SADdist.n_elem, 1);

    arma::vec SADdist_Re = arma::real(SADdist);

    return SADdist_Re;
}
