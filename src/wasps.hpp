# ifndef __CLONEWARS_WASPS_H
# define __CLONEWARS_WASPS_H

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




// Wasp population
class WaspPop {

    arma::vec Y_0_;          // initial wasp abundances by stage
    double sex_ratio_;       // proportion of female wasps
    double K_y_;             // parasitized aphid density dependence
    double s_y_;             // parasitoid adult daily survival
    arma::uvec mum_days_;    // number of days per mummy stage (aphid alive & dead)
    uint32 n_wasp_stages_;   // number of wasp stages (i.e., days)

public:

    // Changing through time
    arma::vec Y_t;                // Wasp density at time t
    arma::vec Y_t1;               // Wasp density at time t+1

    // Constructor
    WaspPop(const arma::uvec& mum_days,
            const double& wasp_density_0,
            const double& sex_ratio,
            const double& K_y,
            const double& s_y)
        : Y_0_(),
          sex_ratio_(sex_ratio),
          K_y_(K_y),
          s_y_(s_y),
          mum_days_(mum_days),
          n_wasp_stages_(arma::accu(mum_days) + 1),
          Y_t(),
          Y_t1() {

        Y_0_ = arma::zeros<arma::vec>(arma::accu(mum_days) + 1);
        Y_0_.back() = wasp_density_0;

        Y_t = Y_0_;
        Y_t1 = Y_0_;

    };

    // (Using Y at time t+1 bc these get calculated before Y get iterated.)
    // Total living, but parasitized aphids
    double total_living_paras() const {
        return arma::accu(Y_t1.head(mum_days_.front()));
    }
    // Total adult wasps
    double total_adult_wasps() const {
        return Y_t1.back();
    }

    // Returning references to private members:
    const arma::vec& Y_0() const { return Y_0_;}
    const double& sex_ratio() const { return sex_ratio_;}
    const double& K_y() const { return K_y_;}
    const double& s_y() const { return s_y_;}
    const arma::uvec& mum_days() const { return mum_days_;}
    const uint32& n_wasp_stages() const { return n_wasp_stages_;}


};






// Wasp attack
class WaspAttack {

    arma::vec rel_attack_;   // relative wasp attack rates by aphid stage
    double a_;               // overall parasitoid attack rate
    double k_;               // aggregation parameter of the nbinom distribution
    double h_;               // parasitoid attack rate handling time
    arma::vec attack_surv_;  // survival rates of singly & multiply attacked aphids

public:
    // Changing through time
    arma::vec A;                  // attack probabilities at time t

    // Constructor
    WaspAttack(const arma::vec& rel_attack,
               const double& a,
               const double& k,
               const double& h,
               const arma::vec& attack_surv)
        : rel_attack_(rel_attack),
          a_(a),
          k_(k),
          h_(h),
          attack_surv_(attack_surv),
          A(arma::zeros<arma::vec>(rel_attack.n_elem)) {};

    /*
     Update attack probabilities
     Equation 6 from Meisner et al. (2014)
     Note: rel_attack is equivalent to p_i
    */
    void iterate_A(const double& Y_m, const double& x) {

        arma::vec mm = (a_ * rel_attack_ * Y_m) / (h_ * x + 1);
        arma::vec AA = (1 + mm / k_);
        if (attack_surv_.n_elem < 2 || arma::accu(attack_surv_) == 0) {
            A = arma::pow(AA, -k_);
        } else {
            A = arma::pow(AA, -k_) +
                attack_surv_(0) * mm % arma::pow(AA, -k_-1) +
                attack_surv_(1) * (1-(arma::pow(AA, -k_) + mm % arma::pow(AA, -k_-1)));
        }
        return;
    }


    // Returning references to private members:
    const arma::vec& rel_attack() const { return rel_attack_; }
    const double& a() const { return a_; }
    const double& k() const { return k_; }
    const double& h() const { return h_; }
    const arma::vec& attack_surv() const { return attack_surv_; }


};




#endif
