# ifndef __CLONEWARS_WASPS_H
# define __CLONEWARS_WASPS_H


#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <cmath>                // log, exp
#include <random>               // normal distribution
#include <cstdint>              // integer types
#include <algorithm>            // find
#include <pcg/pcg_random.hpp>   // pcg prng
#include "clonewars_types.hpp"  // integer types



using namespace Rcpp;




// Wasp attack info
class WaspAttack {

public:

    arma::vec rel_attack;    // relative wasp attack rates by aphid stage
    double a;                // overall parasitoid attack rate
    double k;                // aggregation parameter of the nbinom distribution
    double h;                // parasitoid attack rate handling time


    // Constructors
    WaspAttack()
        : rel_attack(), a(), k(), h() {};
    WaspAttack(const arma::vec& rel_attack_,
               const double& a_,
               const double& k_,
               const double& h_)
        : rel_attack(rel_attack_),
          a(a_),
          k(k_),
          h(h_) {};
    WaspAttack(const WaspAttack& other)
        : rel_attack(other.rel_attack),
          a(other.a),
          k(other.k),
          h(other.h) {};
    WaspAttack& operator=(const WaspAttack& other) {
        rel_attack = other.rel_attack;
        a = other.a;
        k = other.k;
        h = other.h;
        return *this;
    }

    /*
     Compute attack probabilities
     Equation 6 from Meisner et al. (2014)
     Note: rel_attack is equivalent to p_i
     */
    arma::vec A(const double& Y_m,
                const double& x,
                const arma::vec& attack_surv) const {

        arma::vec A_ = (a * rel_attack * Y_m) / (h * x + 1);
        arma::vec AA = (1 + A_ / k);
        if (attack_surv.n_elem < 2 || arma::accu(attack_surv) == 0) {
            A_ = arma::pow(AA, -k);
        } else {
            A_ = arma::pow(AA, -k) +
                attack_surv(0) * A_ % arma::pow(AA, -k-1) +
                attack_surv(1) * (1-(arma::pow(AA, -k) + A_ % arma::pow(AA, -k-1)));
        }
        return A_;
    }



};





// Mummy population
class MummyPop {

public:

    arma::vec Y_0;          // initial mummy densities
    /*
     Proportion of mummies that will NOT take exactly 3 days to develop.
     As this value approaches 2/3, it will provide greater smoothing of
     wasp numbers through time.
     */
    double smooth;


    // Changing through time
    arma::vec Y;            // Mummy density

    // Constructors
    MummyPop() : Y_0(4, arma::fill::zeros), Y(4, arma::fill::zeros) {};
    MummyPop(const arma::vec& Y_0_, const double& smooth_)
        : Y_0(arma::join_vert(arma::vec(1, arma::fill::zeros), Y_0_)),
          smooth(smooth_),
          Y(arma::join_vert(arma::vec(1, arma::fill::zeros), Y_0_)) {};
    MummyPop(const MummyPop& other) : Y_0(other.Y_0), smooth(other.smooth),
    Y(other.Y) {};
    MummyPop& operator=(const MummyPop& other) {
        Y_0 = other.Y_0;
        smooth = other.smooth;
        Y = other.Y;
        return *this;
    }


    // Update # mummies
    // `nm` is # aphids that are newly "mummified"
    void update(const double& pred_rate,
                const double& nm) {

        // Go backwards through stages to avoid conflicts...
        for (uint32 i = Y.n_elem-1; i > 0; i--) {
            Y(i) = (1 - pred_rate) * Y(i-1);
        }

        // Add newly mummified over three days:
        if (smooth > 0) {
            Y(0) = nm * smooth / 2;
            Y(1) += (nm * (1 - smooth));
            Y(2) += (nm * smooth / 2);
        } else Y(1) += nm;

        return;

    }

    // Clearing a plant kills all mummies
    inline void clear() {
        Y.fill(0);
    }
    // Clearing part of plant
    inline void clear(const double& surv) {
        Y *= surv;
    }


};



// Adult wasp population
class WaspPop {

public:

    WaspAttack attack;      // info for attack rates
    double Y_0;             // initial adult wasp density
    uint32 delay;           // when to add initial wasps
    double sex_ratio;       // proportion of female wasps
    double s_y;             // parasitoid adult daily survival

    // for process error:
    std::normal_distribution<double> norm_distr;
    double sigma_y;

    // Changing through time
    double Y;               // Wasp density
    double x;               // Total number of unparasitized aphids


    // Constructors
    WaspPop()
        : attack(), Y_0(), delay(), sex_ratio(), s_y(), norm_distr(),
          sigma_y(), Y(), x() {};
    WaspPop(const arma::vec& rel_attack_,
            const double& a_,
            const double& k_,
            const double& h_,
            const double& Y_0_,
            const uint32& delay_,
            const double& sex_ratio_,
            const double& s_y_,
            const double& sigma_y_)
        : attack(rel_attack_, a_, k_, h_),
          Y_0(Y_0_),
          delay(delay_),
          sex_ratio(sex_ratio_),
          s_y(s_y_),
          norm_distr(0, 1),
          sigma_y(sigma_y_),
          Y((delay_ == 0) ? Y_0_ : 0.0),
          x(0) {};
    WaspPop(const WaspPop& other)
        : attack(other.attack),
          Y_0(other.Y_0),
          delay(other.delay),
          sex_ratio(other.sex_ratio),
          s_y(other.s_y),
          norm_distr(other.norm_distr),
          sigma_y(other.sigma_y),
          Y(other.Y),
          x(other.x) {};
    WaspPop& operator=(const WaspPop& other) {
        attack = other.attack;
        Y_0 = other.Y_0;
        delay = other.delay;
        sex_ratio = other.sex_ratio;
        s_y = other.s_y;
        norm_distr = other.norm_distr;
        sigma_y = other.sigma_y;
        Y = other.Y;
        x = other.x;
        return *this;
    }


    // Return attack matrix
    arma::vec A(const arma::vec& attack_surv) const {
        arma::vec A_ = attack.A(Y, x, attack_surv);
        return A_;
    }

    // Check for whether to add initial wasps:
    void add_Y_0_check(const uint32& t) {
        if (t == delay) Y += Y_0;
        return;
    }

    /*
     Update # adult wasps
     `old_mums` is # mummies at time t that are in the last mummy stage
     before merging
     */
    void update(const double& old_mums,
                pcg32& eng) {
        double max_Y = old_mums + Y;
        if (max_Y == 0) return;
        Y *= s_y;
        Y += (sex_ratio * old_mums);
        Y *= std::exp(norm_distr(eng) * sigma_y);
        if (Y > max_Y) Y = max_Y; // make sure it doesn't exceed what's possible
        return;
    }
    // Same as above, but with no stochasticity
    void update(const double& old_mums) {
        if (old_mums + Y == 0) return;
        Y *= s_y;
        Y += (sex_ratio * old_mums);
        return;
    }

};









#endif
