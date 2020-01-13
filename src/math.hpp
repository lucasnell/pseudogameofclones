# ifndef _CLONEWARS_MATH_H
# define _CLONEWARS_MATH_H


#include <RcppArmadillo.h>
#include "clonewars_types.hpp"


using namespace Rcpp;

void logit__(const arma::vec& p, arma::vec& out);
void inv_logit__(const arma::vec& a, arma::vec& out);

void leslie_matrix__(const arma::uvec& instar_days, const double& surv_juv,
                     const arma::vec& surv_adult, const arma::vec& repro,
                     arma::mat& out);

void leslie_sad__(const arma::mat& L, arma::vec& out);

#endif
