# ifndef __PSEUDOGAMEOFCLONES_WASPS_H
# define __PSEUDOGAMEOFCLONES_WASPS_H


#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <cmath>                // log, exp
#include <random>               // normal distribution
#include <cstdint>              // integer types
#include <algorithm>            // find
#include <pcg/pcg_random.hpp>   // pcg prng
#include "pseudogameofclones_types.hpp"  // integer types



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
     Compute probabilities of aphids surviving wasp attacks and of
     aphids being successfully mummified.
     These don't necessarily add to 1 bc of the possibility of mutual mortality,
     which is especially common for superparasitism.

     Survivals are equation 6 from Meisner et al. (2014), where
     `rel_attack` is equivalent to p_i
     */
    void A_mats(const double& Y_m,
                const double& x,
                arma::vec& A_surv,
                const arma::vec& attack_surv) const {

        uint32 n_max_attacks = attack_surv.n_elem;
        uint32 n_stages = rel_attack.n_elem;

        if (A_surv.n_elem != n_stages) A_surv.set_size(n_stages);

        if (Y_m == 0) {
            A_surv.ones();
            return;
        }

        // Mean of negative binomial distribution:
        arma::vec A_bar = rel_attack * (a * Y_m) / (h * x + 1);

        // Probabilities of being attacked by stage and number of attacks,
        // where the last probability is the prob of being attacked
        // **at least** `n_max_attacks` times:
        arma::mat attack_probs(n_stages, n_max_attacks+1U, arma::fill::none);

        // Probability of being attacked zero times:
        attack_probs.col(0) = arma::pow((1 + A_bar / k), -k);

        // Now we (optionally) expand to being attacked >0 times.
        // The last term will be 1 - (sum of all other probs) because it refers
        // to the prob of being attacked **at least** `n_max_attacks` times.
        if (n_max_attacks == 1) {
            attack_probs.col(1) = 1 - attack_probs.col(0);
        } else if (n_max_attacks > 1) {
            // Two terms that get used for all:
            arma::vec Aa = 1 + k / A_bar;
            arma::vec Ab = 1 + A_bar / k;
            double prod = 1;
            // Will sum all other cols into this vector in the loop below:
            arma::vec A_prob_sums = attack_probs.col(0);
            // Do for all in `attack_surv` and add to `A_`:
            for (uint32 i = 1; i < n_max_attacks; i++) {
                double nt = i; // number of times attacked
                prod *= ((k - 1) / nt + 1);
                attack_probs.col(i) = arma::pow(Aa, -nt) % arma::pow(Ab, -k);
                attack_probs.col(i) *= prod;
                // To avoid an extra loop (see below):
                A_prob_sums += attack_probs.col(i);
            }
            // Now fill in the last column which is 1 - sum(other cols):
            attack_probs.col(n_max_attacks) = (1 - A_prob_sums);
        }

        /*
         Converting attack probabilities to survivals and successful
         mummifications
         */
        for (uint32 j = 0; j < n_stages; j++) {
            A_surv(j) = attack_probs(j, 0);
            for (uint32 i = 0; i < n_max_attacks; i++) {
                A_surv(j) += attack_probs(j, i+1) * attack_surv(i);
            }
        }


        return;
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

    // Clearing a field kills all mummies
    inline void clear() {
        Y.fill(0);
    }
    // Clearing part of field
    inline void clear(const double& surv) {
        Y *= surv;
    }


};



// Adult wasp and mummy population
class WaspPop {

public:

    MummyPop mummies;       // mummy population
    WaspAttack attack;      // info for attack rates
    double Y_0;             // initial adult wasp density
    uint32 delay;           // when to add initial wasps
    double sex_ratio;       // proportion of female wasps
    double s_y;             // parasitoid adult daily survival

    // for process error:
    std::normal_distribution<double> norm_distr;
    double sigma_y;
    // for demographic error:
    std::binomial_distribution<uint32> binom_distr;
    bool demog_error;       // whether to include demographic stochasticity

    // Changing through time
    double Y;               // Wasp density
    double x;               // Total number of unparasitized aphids


    // Constructors
    WaspPop()
        : mummies(), attack(), Y_0(), delay(), sex_ratio(), s_y(), norm_distr(),
          sigma_y(), binom_distr(), demog_error(), Y(), x() {};
    WaspPop(const arma::vec& rel_attack_,
            const double& a_,
            const double& k_,
            const double& h_,
            const double& Y_0_,
            const uint32& delay_,
            const double& sex_ratio_,
            const double& s_y_,
            const double& sigma_y_,
            const bool& demog_error_,
            const arma::vec& mummy_Y_0_,
            const double& mummy_smooth_)
        : mummies(MummyPop(mummy_Y_0_, mummy_smooth_)),
          attack(rel_attack_, a_, k_, h_),
          Y_0(Y_0_),
          delay(delay_),
          sex_ratio(sex_ratio_),
          s_y(s_y_),
          norm_distr(0, 1),
          sigma_y(sigma_y_),
          binom_distr(),
          demog_error(demog_error_),
          Y((delay_ == 0) ? Y_0_ : 0.0),
          x(0) {};
    WaspPop(const WaspPop& other)
        : mummies(other.mummies),
          attack(other.attack),
          Y_0(other.Y_0),
          delay(other.delay),
          sex_ratio(other.sex_ratio),
          s_y(other.s_y),
          norm_distr(other.norm_distr),
          sigma_y(other.sigma_y),
          binom_distr(other.binom_distr),
          demog_error(other.demog_error),
          Y(other.Y),
          x(other.x) {};
    WaspPop& operator=(const WaspPop& other) {
        mummies = other.mummies;
        attack = other.attack;
        Y_0 = other.Y_0;
        delay = other.delay;
        sex_ratio = other.sex_ratio;
        s_y = other.s_y;
        norm_distr = other.norm_distr;
        sigma_y = other.sigma_y;
        binom_distr = other.binom_distr;
        demog_error = other.demog_error;
        Y = other.Y;
        x = other.x;
        return *this;
    }


    // Adjust error parameters, delay, and starting abundances
    void adjust_error_Y0(const bool& demog_error_,
                         const bool& environ_error_,
                         const uint32& delay_,
                         const double& adult_Y0,
                         double mummy_Y0) {

        demog_error = demog_error_;
        delay = delay_;
        if (!environ_error_) sigma_y = 0;
        if (environ_error_ && sigma_y == 0) {
            stop("Cannot have environmental error with sigma_y == 0");
        }

        Y_0 = adult_Y0;
        double mum_y0_sum = arma::accu(mummies.Y_0);
        if (mum_y0_sum != 1 && mum_y0_sum > 0) mummy_Y0 /= mum_y0_sum;
        mummies.Y_0 *= mummy_Y0;
        // refresh starting conditions:
        if (delay_ == 0) {
            Y = Y_0;
        } else Y = 0;
        mummies.Y = mummies.Y_0;
        return;
    }


    // Fills matrix for Pr(aphids survive)
    void A_mats(arma::vec& A_surv,
                const arma::vec& attack_surv) const {
        attack.A_mats(Y, x, A_surv, attack_surv);
        return;
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
        if (demog_error) {
            // Number of adult females that survives is binomial with `s_y`
            // as probability and number of adult females as number of trials.
            uint32 n = static_cast<uint32>(std::round(Y));
            std::binomial_distribution<uint32>::param_type prms(n, s_y);
            binom_distr.param(prms);
            Y = static_cast<double>(binom_distr(eng));
            // Number of mummies that are female is binomial with `sex_ratio`
            // as probability and number of old mummies as number of trials.
            n = static_cast<uint32>(std::round(old_mums));
            prms = std::binomial_distribution<uint32>::param_type(n, sex_ratio);
            binom_distr.param(prms);
            Y += static_cast<double>(binom_distr(eng));
        } else {
            Y *= s_y;
            Y += (sex_ratio * old_mums);
        }
        if (sigma_y > 0) {
            Y *= std::exp(norm_distr(eng) * sigma_y);
            // make sure it doesn't exceed what's possible:
            if (Y > max_Y) Y = max_Y;
        }
        return;
    }


};









#endif
