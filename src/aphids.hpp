# ifndef __CLONEWARS_APHIDS_H
# define __CLONEWARS_APHIDS_H

#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <cmath>                // log, exp
#include <random>               // normal distribution
#include <cstdint>              // integer types
#include <algorithm>            // find
#include <pcg/pcg_random.hpp>   // pcg prng
#include "math.hpp"             // leslie_matrix and leslie_sad
#include "pcg.hpp"              // runif_ fxns
#include "clonewars_types.hpp"  // integer types



using namespace Rcpp;



// Aphid population
class AphidPop {

    arma::mat leslie_;       // Leslie matrix with survival and reproduction
    arma::vec X_0_;          // initial aphid abundances by stage
    double K_;               // aphid density dependence
    uint32 n_aphid_stages_;    // number of aphid stages (i.e., days)

public:
    // Changing through time
    arma::vec X_t;                // Aphid density at time t
    arma::vec X_t1;               // Aphid density at time t+1

    /*
     Constructors
     */
    // Starting with "stable age distribution" with a given total density
    AphidPop(const double& K,
             const arma::uvec& instar_days,
             const double& surv_juv,
             const arma::vec& surv_adult,
             const arma::vec& repro,
             const double& aphid_density_0)
        : leslie_(),
          X_0_(),
          K_(K),
          n_aphid_stages_(arma::accu(instar_days)),
          X_t(), X_t1() {

        leslie_matrix__(instar_days, surv_juv, surv_adult, repro, leslie_);

        leslie_sad__(leslie_, X_0_);
        X_0_ *= aphid_density_0;

        X_t = X_0_;
        X_t1 = X_0_;

    };
    // Starting with starting densities of each stage directly given:
    AphidPop(const double& K,
             const arma::uvec& instar_days,
             const double& surv_juv,
             const arma::vec& surv_adult,
             const arma::vec& repro,
             const arma::vec& aphid_density_0)
        : leslie_(),
          X_0_(),
          K_(K),
          n_aphid_stages_(arma::accu(instar_days)),
          X_t(), X_t1() {

        leslie_matrix__(instar_days, surv_juv, surv_adult, repro, leslie_);

        X_0_ = aphid_density_0;

        X_t = X_0_;
        X_t1 = X_0_;

    };

    /*
     Total (non-parasitized) aphids
     (Using X at time t+1 bc these get calculated before X get iterated.)
     */
    inline double total_aphids() const {
        return arma::accu(X_t1);
    }


    // Returning references to private members:
    const arma::mat& leslie() const {return leslie_;}
    const arma::vec& X_0() const {return X_0_;}
    const double& K() const {return K_;}
    const uint32& n_aphid_stages() const {return n_aphid_stages_;}


};





#endif
