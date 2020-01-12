# ifndef _CLONEWARS_MATH_H
# define _CLONEWARS_MATH_H


#include <RcppArmadillo.h>
#include "clonewars_types.hpp"


using namespace Rcpp;

arma::vec logit__(const arma::vec& p);
arma::vec inv_logit__(const arma::vec& a);

arma::mat leslie_matrix__(const arma::uvec& instar_days, const double& surv_juv,
                          const arma::vec& surv_adult, const arma::vec& repro);

arma::vec leslie_sad__(const arma::mat& L);

#endif
