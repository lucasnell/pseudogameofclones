
#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <cmath>                // log, exp
#include <random>               // normal distribution
#include <cstdint>              // integer types
#include <algorithm>            // find
#include <pcg/pcg_random.hpp>   // pcg prng
#include "pseudogameofclones_types.hpp"  // integer types
#include "wasps.hpp"             // wasp types



using namespace Rcpp;




//[[Rcpp::export]]
SEXP make_wasps_ptr(const arma::vec& rel_attack,
                    const double& a,
                    const double& k,
                    const double& h,
                    const double& sex_ratio,
                    const double& s_y,
                    const double& sigma_y,
                    const double& mummy_smooth,
                    const uint32& mummy_dev_time) {

    // this field will be adjusted later if demographic error is desired
    bool demog_error = false;
    // these will also potentially vary by field so are not included here
    uint32 delay = 0;
    double Y_0 = 0;
    arma::vec mummy_Y_0(mummy_dev_time, arma::fill::zeros);

    XPtr<WaspPop> wasps_xptr(
            new WaspPop(rel_attack, a, k, h, Y_0, delay, sex_ratio, s_y,
                        sigma_y, demog_error, mummy_Y_0, mummy_smooth), true);

    return wasps_xptr;

}
