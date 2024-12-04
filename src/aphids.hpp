# ifndef __PSEUDOGAMEOFCLONES_APHIDS_H
# define __PSEUDOGAMEOFCLONES_APHIDS_H

#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "pseudogameofclones_types.hpp"  // integer types
#include "wasps.hpp"            // wasp classes



using namespace Rcpp;


// Necessary for ApterousPop::alate_prop
class NewOneField;
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
    void process_error(const arma::vec& Xt,
                       const double& sigma_x,
                       const double& rho,
                       const bool& demog_error,
                       const double& aphids_sum,
                       std::normal_distribution<double>& norm_distr,
                       pcg32& eng);

    // Returning references to private members:
    const arma::mat& leslie() const {return leslie_;}
    const arma::vec& X_0() const {return X_0_;}


};


// Aphid "type" population for apterous of a particular clonal line
class ApterousPop : public AphidTypePop {

    friend class AphidPop;

public:

    // Parameters for logit(Pr(alates)) ~ b0 + b1 * N
    double alate_b0_;
    double alate_b1_;


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
    double alate_prop(const NewOneField* field) const;

};
// Aphid "type" population for alates of a particular clonal line
class AlatePop : public AphidTypePop {

    friend class AphidPop;

    // index for stage in which across-field dispersal starts
    uint32 field_disp_start_;


public:

    AlatePop() : AphidTypePop(), field_disp_start_(0) {};
    AlatePop(const arma::mat& leslie_mat,
             const arma::vec& aphid_density_0,
             const uint32& field_disp_start)
        : AphidTypePop(leslie_mat, aphid_density_0),
          field_disp_start_(field_disp_start) {};

    AlatePop(const AlatePop& other)
        : AphidTypePop(other),
          field_disp_start_(other.field_disp_start_) {};

    AlatePop& operator=(const AlatePop& other) {
        AphidTypePop::operator=(other);
        field_disp_start_ = other.field_disp_start_;
        return *this;
    }

    uint32 const& field_disp_start() const { return field_disp_start_; }


};

// Aphid "type" population for parasitized (but alive) aphids
class ParasitizedPop : public AphidTypePop {

    friend class AphidPop;

protected:

    arma::vec s;        // vector of survival rates of parasitized aphids by day

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






// Aphid population: both alates and apterous for one clonal line in a field
class AphidPop {

    double sigma_x;          // environmental standard deviation for aphids
    double rho;              // environmental correlation among instars
    bool demog_error;        // whether to include demographic stochasticity

    uint32 adult_age;           // age at which aphids are adults
                                // (assumed == `field_disp_start`)

    /*
     Vector of length >=2 with survival probabilities of singly & multiply
     attacked aphids:
     */
    arma::vec attack_surv;

    // for process error:
    mutable std::normal_distribution<double> norm_distr =
        std::normal_distribution<double>(0, 1);

    // Process error for all stages, plus checks so that they don't exceed
    // what's possible
    void process_error(const arma::vec& apterous_Xt,
                       const arma::vec& alates_Xt,
                       const arma::vec& paras_Xt,
                       pcg32& eng) {

        // Total aphids for this line:
        double aphids_sum = arma::accu(apterous_Xt) + arma::accu(alates_Xt) +
            arma::accu(paras_Xt);

        apterous.process_error(apterous_Xt, sigma_x, rho, demog_error, aphids_sum,
                               norm_distr, eng);
        alates.process_error(alates_Xt, sigma_x, rho, demog_error, aphids_sum,
                             norm_distr, eng);
        paras.process_error(paras_Xt, sigma_x, rho, demog_error, aphids_sum,
                            norm_distr, eng);

        return;
    }





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
        : sigma_x(0), rho(0), demog_error(false), adult_age(0),
          attack_surv(2, arma::fill::zeros),
          aphid_name(""), apterous(), alates(), paras(), extinct(false) {};

    // Make sure `leslie_mat` has 3 slices and `aphid_density_0` has two columns!
    AphidPop(const std::string& aphid_name_,
             const double& sigma_x_,
             const double& rho_,
             const bool& demog_error_,
             const arma::vec& attack_surv_,
             const arma::cube& leslie_mat,
             const arma::mat& aphid_density_0,
             const double& alate_b0,
             const double& alate_b1,
             const uint32& field_disp_start,
             const uint32& living_days)
        : sigma_x(sigma_x_),
          rho(rho_),
          demog_error(demog_error_),
          adult_age(field_disp_start),
          attack_surv(attack_surv_),
          aphid_name(aphid_name_),
          apterous(leslie_mat.slice(0), aphid_density_0.col(0), alate_b0, alate_b1),
          alates(leslie_mat.slice(1), aphid_density_0.col(1), field_disp_start),
          paras(leslie_mat.slice(2), living_days),
          extinct(false) {};

    AphidPop(const AphidPop& other)
        : sigma_x(other.sigma_x),
          rho(other.rho),
          demog_error(other.demog_error),
          adult_age(other.adult_age),
          attack_surv(other.attack_surv),
          norm_distr(other.norm_distr),
          aphid_name(other.aphid_name),
          apterous(other.apterous),
          alates(other.alates),
          paras(other.paras),
          extinct(other.extinct) {};

    AphidPop& operator=(const AphidPop& other) {

        sigma_x = other.sigma_x;
        rho = other.rho;
        demog_error = other.demog_error;
        adult_age = other.adult_age;
        attack_surv = other.attack_surv;
        norm_distr = other.norm_distr;
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
     Returns vector of abundances of adults that would be moved between fields,
     given that `disp_prop` is the proportion of winged adults that will be
     moved.
     */
    arma::vec remove_field_dispersers(const double& disp_prop) {
        arma::vec D = alates.X;
        uint32 ds = alates.field_disp_start();
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


    // Update new aphid abundances, return the # newly mummified aphids
    double update(const NewOneField* field,
                  const WaspPop* wasps,
                  pcg32& eng);

    /*
     Adult and juvenile numbers for alates and apterous:
     */
    inline double total_adult_apterous() const {
        return arma::accu(apterous.X.tail(apterous.X.n_elem - adult_age));
    }
    inline double total_juven_apterous() const {
        return arma::accu(apterous.X.head(adult_age));
    }
    inline double total_adult_alates() const {
        return arma::accu(alates.X.tail(alates.X.n_elem - adult_age));
    }
    inline double total_juven_alates() const {
        return arma::accu(alates.X.head(adult_age));
    }


};




#endif
