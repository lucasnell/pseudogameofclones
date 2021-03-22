# ifndef __CLONEWARS_APHIDS_H
# define __CLONEWARS_APHIDS_H

#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "clonewars_types.hpp"  // integer types
#include "wasps.hpp"            // wasp classes
#include "math.hpp"             // inv_logit__



using namespace Rcpp;


// This is necessary for dispersal methods and ApterousPop::alate_prop
class OnePatch;
// Necessary here to declare friendship
class AphidPop;



// Generic aphid "type" population: either apterous or alate for a particular clonal line
class AphidTypePop {

protected:

    arma::mat leslie_;       // Leslie matrix with survival and reproduction
    arma::vec X_0_;          // initial aphid abundances by stage


public:
    // Changing through time
    arma::vec X;                // Aphid density

    /*
     Constructors
     */
    AphidTypePop() : leslie_(), X_0_(), X() {};
    AphidTypePop(const arma::mat& leslie_mat,
                 const arma::vec& aphid_density_0)
        : leslie_(leslie_mat),
          X_0_(aphid_density_0),
          X(aphid_density_0) {};

    AphidTypePop(const AphidTypePop& other)
        : leslie_(other.leslie_),
          X_0_(other.X_0_),
          X(other.X) {};

    AphidTypePop& operator=(const AphidTypePop& other) {
        leslie_ = other.leslie_;
        X_0_ = other.X_0_;
        X = other.X;
        return *this;
    }


    /*
     Total aphids
     */
    inline double total_aphids() const {
        return arma::accu(X);
    }

    // Kill all aphids
    inline void clear() {
        X.fill(0);
        return;
    }
    // Kill some of aphids
    inline void clear(const double& surv) {
        X *= surv;
        return;
    }

    // Add process error:
    void process_error(const double& z,
                       const double& sigma_x,
                       const double& rho,
                       const double& demog_mult,
                       std::normal_distribution<double>& norm_distr,
                       pcg32& eng);

    // Returning references to private members:
    const arma::mat& leslie() const {return leslie_;}
    const arma::vec& X_0() const {return X_0_;}


};


// Aphid "type" population for apterous of a particular clonal line
class ApterousPop : public AphidTypePop {

    friend class AphidPop;

    // Parameters for logit(Pr(alates)) ~ b0 + b1 * N
    MEMBER(double, alate_b0)
    MEMBER(double, alate_b1)

public:

    ApterousPop() : AphidTypePop(), alate_b0_(0), alate_b1_(0) {};
    ApterousPop(const arma::mat& leslie_mat,
                const arma::vec& aphid_density_0,
                const double& alate_b0,
                const double& alate_b1)
        : AphidTypePop(leslie_mat, aphid_density_0),
          alate_b0_(alate_b0),
          alate_b1_(alate_b1){};

    ApterousPop(const ApterousPop& other)
        : AphidTypePop(other),
          alate_b0_(other.alate_b0_),
          alate_b1_(other.alate_b1_){};

    ApterousPop& operator=(const ApterousPop& other) {
        AphidTypePop::operator=(other);
        alate_b0_ = other.alate_b0_;
        alate_b1_ = other.alate_b1_;
        return *this;
    }


    // logit(Pr(alates)) ~ b0 + b1 * z, where `z` is # aphids (all lines)
    double alate_prop(const OnePatch* patch) const;

};
// Aphid "type" population for alates of a particular clonal line
class AlatePop : public AphidTypePop {

    friend class AphidPop;

    MEMBER(double, disp_rate)      // rate at which they leave focal plant
    MEMBER(double, disp_mort)      // mortality of dispersers
    MEMBER(uint32, disp_start)     // index for stage in which dispersal starts

public:

    AlatePop() : AphidTypePop(), disp_rate_(0), disp_mort_(0), disp_start_(0) {};
    AlatePop(const arma::mat& leslie_mat,
             const arma::vec& aphid_density_0,
             const double& disp_rate,
             const double& disp_mort,
             const uint32& disp_start)
        : AphidTypePop(leslie_mat, aphid_density_0),
          disp_rate_(disp_rate),
          disp_mort_(disp_mort),
          disp_start_(disp_start) {};

    AlatePop(const AlatePop& other)
        : AphidTypePop(other),
          disp_rate_(other.disp_rate_),
          disp_mort_(other.disp_mort_),
          disp_start_(other.disp_start_) {};

    AlatePop& operator=(const AlatePop& other) {
        AphidTypePop::operator=(other);
        disp_rate_ = other.disp_rate_;
        disp_mort_ = other.disp_mort_;
        disp_start_ = other.disp_start_;
        return *this;
    }



};

// Aphid "type" population for parasitized (but alive) aphids
class ParasitizedPop : public AphidTypePop {

    friend class AphidPop;

protected:

     arma::vec s;    // vector of survival rates of parasitized aphids by day

public:

    ParasitizedPop() : AphidTypePop(), s() {};
    ParasitizedPop(const arma::mat& leslie_mat,
                   const uint32& living_days)
        : AphidTypePop(arma::mat(), arma::vec(living_days, arma::fill::zeros)),
          s(arma::diagvec(leslie_mat, -1)) {
            s.resize(living_days);
    };
    ParasitizedPop(const ParasitizedPop& other)
        : AphidTypePop(other),
          s(other.s) {};

    ParasitizedPop& operator=(const ParasitizedPop& other) {
        AphidTypePop::operator=(other);
        s = other.s;
        return *this;
    }



};




// Aphid population: both alates and apterous for one clonal line on a patch
class AphidPop {

    double sigma_x;          // environmental standard deviation for aphids
    double rho;            // environmental correlation among instars
    double demog_mult;     // multiplier for demographic stochasticity
    /*
     Vector of length 2 with survival rates of singly & multiply attacked
     aphids, respectively:
     */
    arma::vec attack_surv;
    // for process error:
    mutable std::normal_distribution<double> norm_distr =
        std::normal_distribution<double>(0, 1);
    // samples total dispersers:
    mutable std::poisson_distribution<uint32> pois_distr =
        std::poisson_distribution<uint32>(1);
    // samples dead dispersers and alates:
    mutable std::binomial_distribution<uint32> bino_distr =
        std::binomial_distribution<uint32>(1, 0.1);

    // Process error for all stages, plus checks so that they don't exceed
    // what's possible
    void process_error(const arma::vec& apterous_Xt,
                       const arma::vec& alates_Xt,
                       const arma::vec& paras_Xt,
                       const double& z,
                       pcg32& eng);




public:
    std::string aphid_name;    // unique identifying name for this aphid line
    ApterousPop apterous;
    AlatePop alates;
    ParasitizedPop paras;
    bool extinct;

    /*
     Constructors.
     */
    AphidPop()
        : sigma_x(0), rho(0), demog_mult(0), attack_surv(2, arma::fill::zeros),
          aphid_name(""), apterous(), alates(), paras(), extinct(false) {};

    // Make sure `leslie_mat` has 3 slices and `aphid_density_0` has two columns!
    AphidPop(const std::string& aphid_name_,
             const double& sigma_x_,
             const double& rho_,
             const double& demog_mult_,
             const arma::vec& attack_surv_,
             const arma::cube& leslie_mat,
             const arma::mat& aphid_density_0,
             const double& alate_b0,
             const double& alate_b1,
             const double& disp_rate,
             const double& disp_mort,
             const uint32& disp_start,
             const uint32& living_days)
        : sigma_x(sigma_x_),
          rho(rho_),
          demog_mult(demog_mult_),
          attack_surv(attack_surv_),
          aphid_name(aphid_name_),
          apterous(leslie_mat.slice(0), aphid_density_0.col(0), alate_b0, alate_b1),
          alates(leslie_mat.slice(1), aphid_density_0.col(1), disp_rate, disp_mort,
                 disp_start),
          paras(leslie_mat.slice(2), living_days),
          extinct(false) {};

    AphidPop(const AphidPop& other)
        : sigma_x(other.sigma_x),
          rho(other.rho),
          demog_mult(other.demog_mult),
          attack_surv(other.attack_surv),
          norm_distr(other.norm_distr),
          pois_distr(other.pois_distr),
          bino_distr(other.bino_distr),
          aphid_name(other.aphid_name),
          apterous(other.apterous),
          alates(other.alates),
          paras(other.paras),
          extinct(other.extinct) {};

    AphidPop& operator=(const AphidPop& other) {

        sigma_x = other.sigma_x;
        rho = other.rho;
        demog_mult = other.demog_mult;
        attack_surv = other.attack_surv;
        norm_distr = other.norm_distr;
        pois_distr = other.pois_distr;
        bino_distr = other.bino_distr;
        aphid_name = other.aphid_name;
        apterous = other.apterous;
        alates = other.alates;
        paras = other.paras;
        extinct = other.extinct;

        return *this;

    };


    /*
     Total aphids
     */
    inline double total_aphids() const {
        double ta = apterous.total_aphids() + alates.total_aphids() +
            paras.total_aphids();
        return ta;
    }

    /*
     Returns vector of abundances of adults that would be moved between cages,
     given that `disp_prop` is the proportion of winged adults that will be
     moved.
     */
    arma::vec remove_dispersers(const double& disp_prop) {
        arma::vec D = alates.X;
        uint32 ds = alates.disp_start();
        D.head(ds).fill(0);
        /*
         Uncomment below if you want to assume that only young adult alates
         are moved.
         This could be because older ones are already settled on plants, and
         so that the same alates aren't being moved back and forth.
         */
        // for (uint32 i = (ds + 7); i < D.n_elem; i++) D(i) = 0;
        D *= disp_prop;
        alates.X -= D;
        return D;
    }


    // Kill all aphids
    inline void clear() {
        apterous.clear();
        alates.clear();
        paras.clear();
        extinct = true;
        return;
    }
    // Kill some aphids
    inline void clear(const double& surv) {
        apterous.clear(surv);
        alates.clear(surv);
        paras.clear(surv);
        return;
    }


    /*
     Calculate dispersal of this line to all other patches.
     Emigration doesn't necessarily == immigration due to disperser mortality.
    */
    void calc_dispersal(const OnePatch* patch,
                        arma::mat& emigrants,
                        arma::mat& immigrants,
                        pcg32& eng) const;
    // Same, but with no stochasticity
    void calc_dispersal(const OnePatch* patch,
                        arma::mat& emigrants,
                        arma::mat& immigrants) const;

    // Update new aphid abundances, return the # newly mummified aphids
    double update(const OnePatch* patch,
                  const WaspPop* wasps,
                  const arma::vec& emigrants,
                  const arma::vec& immigrants,
                  pcg32& eng);
    // Same as above, but no randomness in alate production:
    double update(const OnePatch* patch,
                  const WaspPop* wasps,
                  const arma::vec& emigrants,
                  const arma::vec& immigrants);

};




#endif
