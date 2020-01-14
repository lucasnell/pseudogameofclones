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



// Generic aphid "type" population: either apterous or alate for a particular clonal line
class AphidTypePop {

protected:

    arma::mat leslie_;       // Leslie matrix with survival and reproduction
    arma::vec X_0_;          // initial aphid abundances by stage
    double K_;               // aphid density dependence
    uint32 n_aphid_stages_;  // number of aphid stages (i.e., days)

public:
    // Changing through time
    arma::vec X_t;                // Aphid density at time t
    arma::vec X_t1;               // Aphid density at time t+1

    /*
     Constructors
     */
    // Starting with "stable age distribution" with a given total density
    AphidTypePop(const double& K,
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
    AphidTypePop(const double& K,
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


// Aphid "type" population for apterous of a particular clonal line
class ApterousPop : public AphidTypePop {

    double alate_rate_;  // rate at which alates are produced

public:

    ApterousPop(const double& K,
                const arma::uvec& instar_days,
                const double& surv_juv,
                const arma::vec& surv_adult,
                const arma::vec& repro,
                const double& aphid_density_0,
                const double& alate_rate)
        : AphidTypePop(K, instar_days, surv_juv, surv_adult, repro, aphid_density_0),
          alate_rate_(alate_rate) {};
    ApterousPop(const double& K,
                const arma::uvec& instar_days,
                const double& surv_juv,
                const arma::vec& surv_adult,
                const arma::vec& repro,
                const arma::vec& aphid_density_0,
                const double& alate_rate)
        : AphidTypePop(K, instar_days, surv_juv, surv_adult, repro, aphid_density_0),
          alate_rate_(alate_rate) {};


    const double& alate_rate() const {return alate_rate_;}

};
// Aphid "type" population for alates of a particular clonal line
class AlatePop : public AphidTypePop {

    double disp_rate_;      // rate at which they leave focal plant
    double disp_mort_;      // mortality of dispersers
    uint32 disp_start_;     // stage in which dispersal starts

public:

    AlatePop(const double& K,
                     const arma::uvec& instar_days,
                     const double& surv_juv,
                     const arma::vec& surv_adult,
                     const arma::vec& repro,
                     const double& aphid_density_0,
                     const double& disp_rate,
                     const double& disp_mort,
                     const uint32& disp_start)
    : AphidTypePop(K, instar_days, surv_juv, surv_adult, repro, aphid_density_0),
      disp_rate_(disp_rate),
      disp_mort_(disp_mort),
      disp_start_(disp_start) {};
    AlatePop(const double& K,
                     const arma::uvec& instar_days,
                     const double& surv_juv,
                     const arma::vec& surv_adult,
                     const arma::vec& repro,
                     const arma::vec& aphid_density_0,
                     const double& disp_rate,
                     const double& disp_mort,
                     const uint32& disp_start)
    : AphidTypePop(K, instar_days, surv_juv, surv_adult, repro, aphid_density_0),
      disp_rate_(disp_rate),
      disp_mort_(disp_mort),
      disp_start_(disp_start) {};


    const double& disp_rate() const {return disp_rate_;}
    const double& disp_mort() const {return disp_mort_;}
    const uint32& disp_start() const {return disp_start_;}

};




// Aphid population: both alates and apterous for one clonal line on a patch
class AphidPop {

    double sigma_;         // environmental standard deviation for aphids
    double rho_;           // environmental correlation among instars
    double demog_mult_;    // multiplier for demographic stochasticity

public:
    ApterousPop apterous;
    AlatePop alates;

    /*
     Constructors.
     Make sure that all input vectors are of length 2!
     */
    // Starting with "stable age distribution" with a given total density
    AphidPop(const double& sigma,
             const double& rho,
             const double& demog_mult,
             const std::vector<double>& K,
             const std::vector<arma::uvec>& instar_days,
             const std::vector<double>& surv_juv,
             const std::vector<arma::vec>& surv_adult,
             const std::vector<arma::vec>& repro,
             const std::vector<double>& aphid_density_0,
             const double& alate_rate,
             const double& disp_rate,
             const double& disp_mort,
             const uint32& disp_start)
        : sigma_(sigma),
          rho_(rho),
          demog_mult_(demog_mult),
          apterous(K[0], instar_days[0], surv_juv[0], surv_adult[0], repro[0],
                   aphid_density_0[0], alate_rate),
          alates(K[1], instar_days[1], surv_juv[1], surv_adult[1], repro[1],
                 aphid_density_0[1], disp_rate, disp_mort, disp_start) {};
    // Starting with starting densities of each stage directly given:
    AphidPop(const double& sigma,
             const double& rho,
             const double& demog_mult,
             const std::vector<double>& K,
             const std::vector<arma::uvec>& instar_days,
             const std::vector<double>& surv_juv,
             const std::vector<arma::vec>& surv_adult,
             const std::vector<arma::vec>& repro,
             const std::vector<arma::vec>& aphid_density_0,
             const double& alate_rate,
             const double& disp_rate,
             const double& disp_mort,
             const uint32& disp_start)
        : sigma_(sigma),
          rho_(rho),
          demog_mult_(demog_mult),
          apterous(K[0], instar_days[0], surv_juv[0], surv_adult[0], repro[0],
                   aphid_density_0[0], alate_rate),
          alates(K[1], instar_days[1], surv_juv[1], surv_adult[1], repro[1],
                 aphid_density_0[1], disp_rate, disp_mort, disp_start) {};



    /*
     Total (non-parasitized) aphids
     (Using X at time t+1 bc these get calculated before X get iterated.)
     */
    inline double total_aphids() const {
        return apterous.total_aphids() + alates.total_aphids();
    }


    const double& sigma() const {return sigma_;}
    const double& rho() const {return rho_;}
    const double& demog_mult() const {return demog_mult_;}

};




#endif
