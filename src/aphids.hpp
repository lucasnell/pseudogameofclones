# ifndef __CLONEWARS_APHIDS_H
# define __CLONEWARS_APHIDS_H

#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "math.hpp"             // leslie_matrix and leslie_sad
#include "clonewars_types.hpp"  // integer types



using namespace Rcpp;


// This is necessary for dispersal methods
class OnePatch;
// Necessary here to declare friendship
class AphidPop;



// Generic aphid "type" population: either apterous or alate for a particular clonal line
class AphidTypePop {

protected:

    arma::mat leslie_;       // Leslie matrix with survival and reproduction
    arma::vec X_0_;          // initial aphid abundances by stage
    uint32 n_aphid_stages_;  // number of aphid stages (i.e., days)


public:
    // Changing through time
    arma::vec X_t;                // Aphid density at time t
    arma::vec X_t1;               // Aphid density at time t+1

    /*
     Constructors
     */
    AphidTypePop()
        : leslie_(),
          X_0_(),
          n_aphid_stages_(),
          X_t(),
          X_t1() {};
    // Starting with "stable age distribution" with a given total density
    AphidTypePop(const arma::uvec& instar_days,
                 const double& surv_juv,
                 const arma::vec& surv_adult,
                 const arma::vec& repro,
                 const double& aphid_density_0)
        : leslie_(),
          X_0_(),
          n_aphid_stages_(arma::accu(instar_days)),
          X_t(), X_t1() {

        leslie_matrix__(instar_days, surv_juv, surv_adult, repro, leslie_);

        leslie_sad__(leslie_, X_0_);
        X_0_ *= aphid_density_0;

        X_t = X_0_;
        X_t1 = X_0_;

    };
    // Starting with starting densities of each stage directly given:
    AphidTypePop(const arma::uvec& instar_days,
                 const double& surv_juv,
                 const arma::vec& surv_adult,
                 const arma::vec& repro,
                 const arma::vec& aphid_density_0)
        : leslie_(),
          X_0_(),
          n_aphid_stages_(arma::accu(instar_days)),
          X_t(), X_t1() {

        leslie_matrix__(instar_days, surv_juv, surv_adult, repro, leslie_);

        X_0_ = aphid_density_0;

        X_t = X_0_;
        X_t1 = X_0_;

    };


    AphidTypePop(const AphidTypePop& other)
        : leslie_(other.leslie_),
          X_0_(other.X_0_),
          n_aphid_stages_(other.n_aphid_stages_),
          X_t(other.X_t),
          X_t1(other.X_t1) {};

    AphidTypePop& operator=(const AphidTypePop& other) {

        leslie_ = other.leslie_;
        X_0_ = other.X_0_;
        n_aphid_stages_ = other.n_aphid_stages_;
        X_t = other.X_t;
        X_t1 = other.X_t1;

        return *this;
    }


    /*
     Total (non-parasitized) aphids
     (Using X at time t+1 bc these get calculated before X get iterated.)
     */
    inline double total_aphids() const {
        return arma::accu(X_t1);
    }

    // Kill all aphids
    inline void clear() {
        X_t1.fill(0);
        return;
    }

    // Add process error:
    void process_error(const double& z,
                       const double& sigma,
                       const double& rho,
                       const double& demog_mult,
                       std::normal_distribution<double>& norm_distr,
                       pcg32& eng);

    // Returning references to private members:
    const arma::mat& leslie() const {return leslie_;}
    const arma::vec& X_0() const {return X_0_;}
    const uint32& n_aphid_stages() const {return n_aphid_stages_;}


};


// Aphid "type" population for apterous of a particular clonal line
class ApterousPop : public AphidTypePop {

    friend class AphidPop;

    double alate_prop_;  // proportion of new offspring that are alates

public:

    ApterousPop() : AphidTypePop(), alate_prop_(0) {};
    ApterousPop(const arma::uvec& instar_days,
                const double& surv_juv,
                const arma::vec& surv_adult,
                const arma::vec& repro,
                const double& aphid_density_0,
                const double& alate_prop)
        : AphidTypePop(instar_days, surv_juv, surv_adult, repro, aphid_density_0),
          alate_prop_(alate_prop) {};
    ApterousPop(const arma::uvec& instar_days,
                const double& surv_juv,
                const arma::vec& surv_adult,
                const arma::vec& repro,
                const arma::vec& aphid_density_0,
                const double& alate_prop)
        : AphidTypePop(instar_days, surv_juv, surv_adult, repro, aphid_density_0),
          alate_prop_(alate_prop) {};

    ApterousPop(const ApterousPop& other)
        : AphidTypePop(other),
          alate_prop_(other.alate_prop_) {};

    ApterousPop& operator=(const ApterousPop& other) {
        AphidTypePop::operator=(other);
        alate_prop_ = other.alate_prop_;
        return *this;
    }


    const double& alate_prop() const {return alate_prop_;}

};
// Aphid "type" population for alates of a particular clonal line
class AlatePop : public AphidTypePop {

    friend class AphidPop;

    double disp_rate_;      // rate at which they leave focal plant
    double disp_mort_;      // mortality of dispersers
    uint32 disp_start_;     // index for stage in which dispersal starts

public:

    AlatePop() : AphidTypePop(), disp_rate_(0), disp_mort_(0), disp_start_(0) {};
    AlatePop(const arma::uvec& instar_days,
             const double& surv_juv,
             const arma::vec& surv_adult,
             const arma::vec& repro,
             const double& aphid_density_0,
             const double& disp_rate,
             const double& disp_mort)
        : AphidTypePop(instar_days, surv_juv, surv_adult, repro, aphid_density_0),
          disp_rate_(disp_rate),
          disp_mort_(disp_mort),
          disp_start_(arma::accu(instar_days.head(instar_days.n_elem - 1)) - 1) {};
    AlatePop(const arma::uvec& instar_days,
             const double& surv_juv,
             const arma::vec& surv_adult,
             const arma::vec& repro,
             const arma::vec& aphid_density_0,
             const double& disp_rate,
             const double& disp_mort)
        : AphidTypePop(instar_days, surv_juv, surv_adult, repro, aphid_density_0),
          disp_rate_(disp_rate),
          disp_mort_(disp_mort),
          disp_start_(arma::accu(instar_days.head(instar_days.n_elem - 1)) - 1) {};

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


    const double& disp_rate() const {return disp_rate_;}
    const double& disp_mort() const {return disp_mort_;}
    const uint32& disp_start() const {return disp_start_;}

};




// Aphid population: both alates and apterous for one clonal line on a patch
class AphidPop {

    double sigma_;         // environmental standard deviation for aphids
    double rho_;           // environmental correlation among instars
    double demog_mult_;    // multiplier for demographic stochasticity
    mutable std::normal_distribution<double> norm_distr;    // for process error
    mutable std::poisson_distribution<uint32> pois_distr;   // samples total dispersers
    mutable std::binomial_distribution<uint32> bino_distr;  // samples dead dispersers


public:
    std::string aphid_name;    // unique identifying name for this aphid line
    ApterousPop apterous;
    AlatePop alates;

    /*
     Constructors.
     Make sure that all input vectors are of length 2!
     */
    AphidPop()
        : sigma_(0), rho_(0), demog_mult_(0), norm_distr(0,1), pois_distr(1),
          bino_distr(1, 0.1), aphid_name(""), apterous(), alates() {};
    // Starting with "stable age distribution" with a given total density
    AphidPop(const std::string& aphid_name_,
             const double& sigma,
             const double& rho,
             const double& demog_mult,
             const std::vector<arma::uvec>& instar_days,
             const std::vector<double>& surv_juv,
             const std::vector<arma::vec>& surv_adult,
             const std::vector<arma::vec>& repro,
             const std::vector<double>& aphid_density_0,
             const double& alate_prop,
             const double& disp_rate,
             const double& disp_mort)
        : aphid_name(aphid_name_),
          sigma_(sigma),
          rho_(rho),
          demog_mult_(demog_mult),
          norm_distr(0,1),
          pois_distr(1),
          bino_distr(1, 0.1),
          apterous(instar_days[0], surv_juv[0], surv_adult[0], repro[0],
                   aphid_density_0[0], alate_prop),
          alates(instar_days[1], surv_juv[1], surv_adult[1], repro[1],
                 aphid_density_0[1], disp_rate, disp_mort) {};
    // Starting with starting densities of each stage directly given:
    AphidPop(const std::string& aphid_name_,
             const double& sigma,
             const double& rho,
             const double& demog_mult,
             const std::vector<arma::uvec>& instar_days,
             const std::vector<double>& surv_juv,
             const std::vector<arma::vec>& surv_adult,
             const std::vector<arma::vec>& repro,
             const std::vector<arma::vec>& aphid_density_0,
             const double& alate_prop,
             const double& disp_rate,
             const double& disp_mort)
        : aphid_name(aphid_name_),
          sigma_(sigma),
          rho_(rho),
          demog_mult_(demog_mult),
          norm_distr(0,1),
          pois_distr(1),
          bino_distr(1, 0.1),
          apterous(instar_days[0], surv_juv[0], surv_adult[0], repro[0],
                   aphid_density_0[0], alate_prop),
          alates(instar_days[1], surv_juv[1], surv_adult[1], repro[1],
                 aphid_density_0[1], disp_rate, disp_mort) {};



    /*
     Total (non-parasitized) aphids
     (Using X at time t+1 bc these get calculated before X get iterated.)
     */
    inline double total_aphids() const {
        return apterous.total_aphids() + alates.total_aphids();
    }

    // Kill all aphids
    inline void clear() {
        apterous.clear();
        alates.clear();
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

    // Update new aphid abundances
    void update_pop(const double& z,
                    const double& S,
                    const double& pred_rate,
                    const arma::vec& emigrants,
                    const arma::vec& immigrants,
                    pcg32& eng);
    // Same as above, but no randomness in alate production:
    void update_pop(const double& z,
                    const double& S,
                    const double& pred_rate,
                    const arma::vec& emigrants,
                    const arma::vec& immigrants);

    const double& sigma() const {return sigma_;}
    const double& rho() const {return rho_;}
    const double& demog_mult() const {return demog_mult_;}

};




#endif
