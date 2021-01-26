# ifndef __CLONEWARS_WASPS_H
# define __CLONEWARS_WASPS_H


#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <cmath>                // log, exp
#include <random>               // normal distribution
#include <cstdint>              // integer types
#include <algorithm>            // find
#include <pcg/pcg_random.hpp>   // pcg prng
#include "pcg.hpp"              // runif_ fxns
#include "clonewars_types.hpp"  // integer types



using namespace Rcpp;




// Wasp attack info
class WaspAttack {

    arma::vec rel_attack;    // relative wasp attack rates by aphid stage
    double a;                // overall parasitoid attack rate
    double k;                // aggregation parameter of the nbinom distribution
    double h;                // parasitoid attack rate handling time

public:

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





// Wasp population (+ mummies)
class WaspPop {

    WaspAttack attack;      // info for attack rates
    arma::vec Y_0;          // initial mummy and adult wasp densities
    double sex_ratio;       // proportion of female wasps
    double s_y;             // parasitoid adult daily survival

    std::normal_distribution<double> norm_distr;    // for process error

public:

    // Changing through time
    arma::vec Y;            // Wasp and mummy density

    // Constructors
    WaspPop() : attack(), Y_0(), sex_ratio(), s_y(), Y(), norm_distr() {};
    WaspPop(const arma::vec& rel_attack_,
            const double& a_,
            const double& k_,
            const double& h_,
            const arma::vec& Y_0_,
            const double& sex_ratio_,
            const double& s_y_,
            const double& sigma_y)
        : attack(rel_attack_, a_, k_, h_),
          Y_0(Y_0_),
          sex_ratio(sex_ratio_),
          s_y(s_y_),
          Y(Y_0_),
          norm_distr(0, sigma_y) {};
    WaspPop(const WaspPop& other)
        : attack(other.attack),
          Y_0(other.Y_0),
          sex_ratio(other.sex_ratio),
          s_y(other.s_y),
          Y(other.Y),
          norm_distr(other.norm_distr) {};
    WaspPop& operator=(const WaspPop& other) {
        attack = other.attack;
        Y_0 = other.Y_0;
        sex_ratio = other.sex_ratio;
        s_y = other.s_y;
        Y = other.Y;
        norm_distr = other.norm_distr;
        return *this;
    }

    // Total # mummies
    double total_mummies() const {
        return arma::accu(Y.head(Y.n_elem - 1));
    }
    // Total adult wasps
    double total_wasps() const {
        return Y.back();
    }


    arma::vec A(const double& x,
                const arma::vec& attack_surv) const {

        arma::vec A_ = attack.A(Y.back(), x, attack_surv);

        return A_;

    }

    // Update # mummies and adult wasps
    // `nm` is # aphids that are newly "mummified"
    void update(const double& pred_rate,
                const double& nm,
                pcg32& eng) {

        update(pred_rate, nm);

        Y.back() *= std::exp(norm_distr(eng));

        return;

    }
    // Same as above, but with no stochasticity
    void update(const double& pred_rate,
                const double& nm) {

        // Go backwards through stages to avoid conflicts...
        // adult wasps
        Y.back() = s_y * Y.back() + sex_ratio * Y(Y.n_elem-2);
        // mummies
        for (uint32 t = 0; t < Y.n_elem-2; t++) Y(t+1) = pred_rate * Y(t);
        // Newly mummified:
        Y.front() = nm;

        return;

    }

    // Clearing a patch only kills mummies
    inline void clear() {
        Y.head(Y.n_elem-1).fill(0);
    }


};









#endif
