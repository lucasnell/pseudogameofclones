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


public:
    // Changing through time
    arma::vec X_t;                // Aphid density at time t
    arma::vec X_t1;               // Aphid density at time t+1

    /*
     Constructors
     */
    AphidTypePop() : leslie_(), X_0_(), X_t(), X_t1() {};
    AphidTypePop(const arma::mat& leslie_mat,
                 const arma::vec& aphid_density_0)
        : leslie_(leslie_mat),
          X_0_(aphid_density_0),
          X_t(aphid_density_0),
          X_t1(aphid_density_0) {};

    AphidTypePop(const AphidTypePop& other)
        : leslie_(other.leslie_),
          X_0_(other.X_0_),
          X_t(other.X_t),
          X_t1(other.X_t1) {};

    AphidTypePop& operator=(const AphidTypePop& other) {

        leslie_ = other.leslie_;
        X_0_ = other.X_0_;
        X_t = other.X_t;
        X_t1 = other.X_t1;

        return *this;
    }


    /*
     Total aphids
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
    ApterousPop(const arma::mat& leslie_mat,
                const arma::vec& aphid_density_0,
                const double& alate_prop)
        : AphidTypePop(leslie_mat, aphid_density_0),
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
    bool extinct;

    /*
     Constructors.
     */
    AphidPop()
        : sigma_(0), rho_(0), demog_mult_(0), norm_distr(0,1), pois_distr(1),
          bino_distr(1, 0.1), aphid_name(""), apterous(), alates() {};

    // Make sure `leslie_mat` has 2 slices and `aphid_density_0` has two columns!
    AphidPop(const std::string& aphid_name_,
             const double& sigma,
             const double& rho,
             const double& demog_mult,
             const arma::cube& leslie_mat,
             const arma::mat& aphid_density_0,
             const double& alate_prop,
             const double& disp_rate,
             const double& disp_mort,
             const uint32& disp_start)
        : aphid_name(aphid_name_),
          sigma_(sigma),
          rho_(rho),
          demog_mult_(demog_mult),
          norm_distr(0,1),
          pois_distr(1),
          bino_distr(1, 0.1),
          apterous(leslie_mat.slice(0), aphid_density_0.col(0), alate_prop),
          alates(leslie_mat.slice(1), aphid_density_0.col(1), disp_rate, disp_mort,
                 disp_start),
          extinct(false) {};


    /*
     Total aphids
     */
    inline double total_aphids() const {
        if (extinct) return 0;
        return apterous.total_aphids() + alates.total_aphids();
    }

    // Kill all aphids
    inline void clear() {
        apterous.clear();
        alates.clear();
        extinct = true;
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
    void update_pop(const OnePatch* patch,
                    const arma::vec& emigrants,
                    const arma::vec& immigrants,
                    pcg32& eng);
    // Same as above, but no randomness in alate production:
    void update_pop(const OnePatch* patch,
                    const arma::vec& emigrants,
                    const arma::vec& immigrants);

    const double& sigma() const {return sigma_;}
    const double& rho() const {return rho_;}
    const double& demog_mult() const {return demog_mult_;}

};




#endif
