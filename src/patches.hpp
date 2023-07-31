# ifndef __GAMEOFCLONES_PATCHES_H
# define __GAMEOFCLONES_PATCHES_H


#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "gameofclones_types.hpp"  // integer types
#include "math.hpp"             // distributions
#include "aphids.hpp"           // aphid classes
#include "wasps.hpp"            // wasp classes



using namespace Rcpp;

// Necessary here to declare friendship
class AllFields;



// Info for a single perturbation (it happens to all plants within a field)
struct PerturbInfo {

    // When to perturb.
    uint32 time;
    // Where to perturb.
    uint32 field;
    // Number to multiply abundance by.
    double multiplier;
    /*
     Who to perturb.
     Indices less than the # aphids points to the particular
     aphid line, equal to the # aphids is mummies, and greater than
     the # aphids is adult wasps.
     */
    uint32 index;

    PerturbInfo() : time(), field(), multiplier(), index() {}
    PerturbInfo(const uint32& when, const uint32& where,
                const uint32& who, const double& how)
        : time(when), field(where), multiplier(how), index(who) {}
    PerturbInfo& operator=(const PerturbInfo& other) {
        time = other.time;
        field = other.field;
        multiplier = other.multiplier;
        index = other.index;
        return *this;
    }

};




// Info for clearing plants
template <typename T>
struct PlantClearingInfo {

    uint32 ind;
    double N;
    bool wilted;
    T thresh_info; // age or # aphids

    PlantClearingInfo(const uint32& ind_,
                      const double& N_,
                      const bool& wilted_,
                      const T& thresh_info_)
        : ind(ind_), N(N_), wilted(wilted_), thresh_info(thresh_info_) {}
    PlantClearingInfo(const PlantClearingInfo<T>& other)
        : ind(other.ind), N(other.N), thresh_info(other.thresh_info) {}

    // This returns true when `other` should be cleared first
    bool operator<(const PlantClearingInfo<T>& other) const {
        if (other.wilted && !wilted) return true;
        if (!other.wilted && wilted) return false;
        return thresh_info < other.thresh_info;
    }
    // This returns true when `this` should be cleared first
    bool operator>(const PlantClearingInfo<T>& other) const {
        if (other.wilted && !wilted) return false;
        if (!other.wilted && wilted) return true;
        return thresh_info > other.thresh_info;
    }

};


/*
One plant or a number of plants so close that aphids freely disperse
 across them.
*/
class OnePlant {


    /*
     Adjust for potential extinction or re-colonization:
     */
    void extinct_colonize(const uint32& i);


    // Updates total # living aphids in this plant (`z`), and
    // checks to see if plant should be wilted.
    void update_z_wilted();

    MEMBER(bool,wilted)


public:

    std::vector<AphidPop> aphids;   // aphids on this plant
    MummyPop mummies;               // mummies on this plant
    bool empty;                     // whether no aphids are on this plant
    double pred_rate;               // predation on aphids and mummies
    double K;                       // unparasitized aphid carrying capacity
    double K_y;                     // parasitized aphid carrying capacity
    double z = 0;                   // sum of all living aphids at time t
    double S = 0;                   // effect of density dependence on aphids
    double S_y = 0;                 // effect of dd on parasitized aphids
    uint32 n_plants;                // total # plants
    uint32 this_j;                  // index for this plant
    uint32 age = 0;                 // age of this plant
    double wilted_prop;             // prop. of CC that kills plant
    double wilted_mort;             // growth-rate modifier due to wilting
    double extinct_N;               // threshold for calling something extinct
    double max_mum_density;         // maximum mummy density (ignored if zero)



    OnePlant()
        : wilted_(false), aphids(), mummies(), empty(true), pred_rate(0),
          K(0), K_y(1), n_plants(1), this_j(0),
          wilted_prop(1), wilted_mort(1), extinct_N() {};

    /*
     In `aphid_density_0` below, rows are aphid stages, columns are types
     (alate vs apterous), and slices are aphid lines.
     In `leslie_mat` below, items in vector are aphid lines, slices are
     alate/apterous/parasitized.
     */
    OnePlant(const double& sigma_x,
             const double& rho,
             const bool& aphid_demog_error,
             const arma::mat& attack_surv_,
             const double& K_,
             const double& K_y_,
             const double& wilted_prop_,
             const double& wilted_mort_,
             const std::vector<std::string>& aphid_name,
             const std::vector<arma::cube>& leslie_mat,
             const arma::cube& aphid_density_0,
             const std::vector<double>& alate_b0,
             const std::vector<double>& alate_b1,
             const std::vector<double>& alate_plant_disp_p,
             const std::vector<double>& disp_mort,
             const std::vector<uint32>& disp_start,
             const std::vector<uint32>& living_days,
             const double& pred_rate_,
             const uint32& n_plants_,
             const uint32& this_j_,
             const double& extinct_N_,
             const arma::vec& mum_density_0,
             const double& mum_smooth,
             const double& max_mum_density_)
        : wilted_(false),
          aphids(),
          mummies(mum_density_0, mum_smooth),
          empty(true),
          pred_rate(pred_rate_),
          K(K_),
          K_y(K_y_),
          n_plants(n_plants_),
          this_j(this_j_),
          wilted_prop(wilted_prop_),
          wilted_mort(wilted_mort_),
          extinct_N(extinct_N_),
          max_mum_density(max_mum_density_) {

        uint32 n_lines = aphid_name.size();

        aphids.reserve(n_lines);

        for (uint32 i = 0; i < n_lines; i++) {
            AphidPop ap(aphid_name[i], sigma_x, rho, aphid_demog_error,
                        attack_surv_.col(i),
                        leslie_mat[i], aphid_density_0.slice(i),
                        alate_b0[i], alate_b1[i], alate_plant_disp_p[i],
                        disp_mort[i], disp_start[i], living_days[i]);
            aphids.push_back(ap);
            double N = aphids.back().total_aphids();
            if (N < extinct_N) {
                aphids.back().clear();
            } else if (empty) empty = false;
        }

    };



    OnePlant(const OnePlant& other)
        : wilted_(false), aphids(other.aphids), mummies(other.mummies),
          empty(other.empty), pred_rate(other.pred_rate),
          K(other.K), K_y(other.K_y), z(other.z),
          S(other.S), S_y(other.S_y), n_plants(other.n_plants),
          this_j(other.this_j), age(other.age), wilted_prop(other.wilted_prop),
          wilted_mort(other.wilted_mort), extinct_N(other.extinct_N),
          max_mum_density(other.max_mum_density) {};

    OnePlant& operator=(const OnePlant& other) {
        wilted_ = other.wilted_;
        aphids = other.aphids;
        mummies = other.mummies;
        empty = other.empty;
        pred_rate = other.pred_rate;
        K = other.K;
        K_y = other.K_y;
        z = other.z;
        S = other.S;
        S_y = other.S_y;
        n_plants = other.n_plants;
        this_j = other.this_j;
        age = other.age;
        wilted_prop = other.wilted_prop;
        wilted_mort = other.wilted_mort;
        extinct_N = other.extinct_N;
        max_mum_density = other.max_mum_density;
        return *this;
    }


    inline uint32 size() const noexcept {
        return aphids.size();
    }

    AphidPop& operator[](const uint32& idx) {
        return aphids[idx];
    }
    const AphidPop& operator[](const uint32& idx) const {
        return aphids[idx];
    }




    /*
     Clear to no aphids or mummies
     */
    void clear(const double& K_,
               const double& K_y_,
               const double& wilted_mort_) {
        for (AphidPop& ap : aphids) ap.clear();
        mummies.clear();
        empty = true;
        wilted_ = false;
        age = 0;
        K = K_;
        K_y = K_y_;
        wilted_mort = wilted_mort_;
        return;
    }
    /*
     Clear some of the aphids and mummies
     */
    void clear(const double& K_,
               const double& K_y_,
               const double& wilted_mort_,
               const double& surv) {
        empty = true;
        for (AphidPop& ap : aphids) {
            ap.clear(surv);
            double N = ap.total_aphids();
            if (N < extinct_N) {
                ap.clear();
            } else {
                empty = false;
                ap.extinct = false;
            }
        }
        mummies.clear(surv);
        if (arma::accu(mummies.Y) < extinct_N) mummies.Y.fill(0);
        wilted_ = false;
        age = 0;
        K = K_;
        K_y = K_y_;
        wilted_mort = wilted_mort_;
        return;
    }

    // Total (living) aphids on plant
    inline double total_aphids() const {
        double ta = 0;
        for (const AphidPop& ap : aphids) {
            ta += ap.total_aphids();
        }
        return ta;
    }
    // Total UNparasitized aphids on plant
    inline double total_unpar_aphids() const {
        double ta = 0;
        for (const AphidPop& ap : aphids) {
            ta += ap.apterous.total_aphids() + ap.alates.total_aphids();
        }
        return ta;
    }
    // Total mummies on plant
    inline double total_mummies() const {
        double tm = arma::accu(mummies.Y);
        return tm;
    }

    /*
     Add dispersal info to `emigrants` and `immigrants` cubes.
     In these cubes, rows are aphid stages, columns are plants,
     and slices are aphid lines.
    */
    void calc_dispersal(arma::cube& emigrants,
                        arma::cube& immigrants,
                        pcg32& eng) const {
        for (uint32 i = 0; i < aphids.size(); i++) {
            aphids[i].calc_dispersal(this, emigrants.slice(i),
                                     immigrants.slice(i), eng);
        }
        return;
    }

    // Same thing as above, but overloaded for not including dispersal stoch.
    void calc_dispersal(arma::cube& emigrants,
                        arma::cube& immigrants) const {
        for (uint32 i = 0; i < aphids.size(); i++) {
            aphids[i].calc_dispersal(this, emigrants.slice(i),
                                     immigrants.slice(i));
        }
        return;
    }


    /*
     Iterate one time step, after calculating dispersal numbers
     */
    void update(const arma::cube& emigrants,
                const arma::cube& immigrants,
                const WaspPop* wasps,
                pcg32& eng);


};




/*
 One field of plants.

 For the constructors below, all vector arguments should have a length equal
 to the number of aphid lines, except for `K`, `aphid_density_0`, and
 `pred_rate_`;
 these should have a length equal to the number of plants.
 Each item in `aphid_density_0` should have a length equal to the number of
 aphid lines.

 */

class OneField {

    friend class AllFields;

    // For new Ks if desired:
    trunc_normal_distribution tnorm_distr;

    // For new wilted_morts if desired:
    beta_distribution beta_distr;

    double mean_K_;                 // mean of distribution of `K` for plants
    double sd_K_;                   // sd of distribution of `K` for plants

    // Helps calculate carrying capacity for parasitized aphids
    // (`K_y = K * K_y_mult`):
    double K_y_mult;

    double shape1_wilted_mort_;      // shape1 of distr of `wilted_mort`
    double shape2_wilted_mort_;      // shape2 of distr of `wilted_mort`

    double extinct_N;               // used here for the wasps
    bool constant_wasps;            // keep wasp abundance constant?


    // Set K and K_y
    void set_K(double& K, double& K_y, pcg32& eng) {
        if (sd_K_ == 0) {
            K = mean_K_;
        } else {
            K = tnorm_distr(eng);
        }
        K_y = K * K_y_mult;
        return;
    }
    // Set plant-wilted mortality parameter
    void set_wilted_mort(double& wilted_mort, pcg32& eng) {
        if (shape2_wilted_mort_ == 0) {
            wilted_mort = shape1_wilted_mort_;
        } else {
            wilted_mort = beta_distr(eng);
        }
        return;
    }

    inline void set_wasp_info(double& old_mums) {
        wasps.x = 0;
        old_mums = 0;
        for (OnePlant& p : plants) {
            wasps.x += p.total_unpar_aphids();
            old_mums += p.mummies.Y.back();
        }
        return;
    }


    // Do the actual clearing of plants while avoiding extinction
    template <typename T>
    void do_clearing(std::vector<PlantClearingInfo<T>>& clear_plants,
                     double& remaining,
                     std::vector<bool>& wilted,
                     const double& clear_surv,
                     pcg32& eng);


public:

    std::vector<OnePlant> plants;
    WaspPop wasps;
    arma::cube emigrants;
    arma::cube immigrants;


    OneField()
        : tnorm_distr(), beta_distr(), mean_K_(), sd_K_(), K_y_mult(),
          shape1_wilted_mort_(), shape2_wilted_mort_(), extinct_N(),
          constant_wasps(),
          plants(), wasps(), emigrants(), immigrants() {};

    /*
     In `aphid_density_0` below, rows are aphid stages, columns are types
     (alate vs apterous), and slices are aphid lines.
     In `leslie_mat` below, slices are aphid lines.
     */
    OneField(const double& sigma_x,
               const double& sigma_y,
               const double& rho,
               const bool& aphid_demog_error,
               const bool& wasp_demog_error,
               const double& mean_K,
               const double& sd_K,
               const double& K_y_mult_,
               const double& wilted_prop,
               const double& shape1_wilted_mort,
               const double& shape2_wilted_mort,
               const arma::mat& attack_surv_,
               const std::vector<std::string>& aphid_name,
               const std::vector<arma::cube>& leslie_mat,
               const std::vector<arma::cube>& aphid_density_0,
               const std::vector<double>& alate_b0,
               const std::vector<double>& alate_b1,
               const std::vector<double>& alate_plant_disp_p,
               const std::vector<double>& disp_mort,
               const std::vector<uint32>& disp_start,
               const std::vector<uint32>& living_days,
               const double& pred_rate,
               const double& extinct_N_,
               const arma::mat& mum_density_0,
               const double& mum_smooth_,
               const double& max_mum_density_,
               const arma::vec& rel_attack_,
               const double& a_,
               const double& k_,
               const double& h_,
               const double& wasp_density_0_,
               const uint32& wasp_delay_,
               const double& sex_ratio_,
               const double& s_y_,
               const bool& constant_wasps_,
               pcg32& eng)
        : tnorm_distr(),
          beta_distr(),
          mean_K_(mean_K),
          sd_K_(sd_K),
          K_y_mult(K_y_mult_),
          shape1_wilted_mort_(shape1_wilted_mort),
          shape2_wilted_mort_(shape2_wilted_mort),
          extinct_N(extinct_N_),
          constant_wasps(constant_wasps_),
          plants(),
          wasps(rel_attack_, a_, k_, h_, wasp_density_0_, wasp_delay_,
                sex_ratio_, s_y_, sigma_y, wasp_demog_error),
          emigrants(),
          immigrants() {


        /*
         None of these hyperparameters can be <= 0, so if they're <= then
         that makes the associated parameter whose distribution it describes
         (e.g., K for `sd_K_`) not have stochasticity.
         This means that the first hyperparameter provided will be the value
         used for ALL of that parameter WITHOUT BEING TRANSFORMED.
         */
        if (sd_K_ <= 0) {
            sd_K_ = 0;
        } else {
            tnorm_distr = trunc_normal_distribution(mean_K_, sd_K_);
        }
        if (shape2_wilted_mort_ <= 0) {
            shape2_wilted_mort_ = 0;
        } else {
            beta_distr = beta_distribution(shape1_wilted_mort_,
                                           shape2_wilted_mort_);
        }

        uint32 n_plants = aphid_density_0.size();

        uint32 n_lines = aphid_name.size();
        uint32 n_stages = leslie_mat.front().n_rows;

        double K, K_y, wilted_mort;
        plants.reserve(n_plants);
        for (uint32 j = 0; j < n_plants; j++) {
            set_K(K, K_y, eng);
            set_wilted_mort(wilted_mort, eng);
            OnePlant ap(sigma_x, rho, aphid_demog_error, attack_surv_,
                        K, K_y, wilted_prop, wilted_mort,
                        aphid_name, leslie_mat,
                        aphid_density_0[j], alate_b0, alate_b1,
                        alate_plant_disp_p, disp_mort,
                        disp_start, living_days, pred_rate,
                        n_plants, j, extinct_N_,
                        mum_density_0.col(j), mum_smooth_, max_mum_density_);
            plants.push_back(ap);
        }

        emigrants = arma::zeros<arma::cube>(n_stages, n_plants, n_lines);
        immigrants = arma::zeros<arma::cube>(n_stages, n_plants, n_lines);

    }

    OneField(const OneField& other)
        : tnorm_distr(other.tnorm_distr),
          beta_distr(other.beta_distr),
          mean_K_(other.mean_K_),
          sd_K_(other.sd_K_),
          K_y_mult(other.K_y_mult),
          shape1_wilted_mort_(other.shape1_wilted_mort_),
          shape2_wilted_mort_(other.shape2_wilted_mort_),
          extinct_N(other.extinct_N),
          constant_wasps(other.constant_wasps),
          plants(other.plants),
          wasps(other.wasps),
          emigrants(other.emigrants),
          immigrants(other.immigrants) {};

    OneField& operator=(const OneField& other) {
        tnorm_distr = other.tnorm_distr;
        beta_distr = other.beta_distr;
        mean_K_ = other.mean_K_;
        sd_K_ = other.sd_K_;
        K_y_mult = other.K_y_mult;
        shape1_wilted_mort_ = other.shape1_wilted_mort_;
        shape2_wilted_mort_ = other.shape2_wilted_mort_;
        extinct_N = other.extinct_N;
        constant_wasps = other.constant_wasps;
        plants = other.plants;
        wasps = other.wasps;
        emigrants = other.emigrants;
        immigrants = other.immigrants;
        return *this;
    };


    inline uint32 size() const noexcept {
        return plants.size();
    }

    OnePlant& operator[](const uint32& idx) {
        return plants[idx];
    }
    const OnePlant& operator[](const uint32& idx) const {
        return plants[idx];
    }

    // Remove dispersers from this field:
    arma::mat remove_dispersers(const double& disp_prop) {

        uint32 n_lines = plants[0].aphids.size();
        uint32 n_stages = plants[0].aphids[0].alates.X.n_elem;

        arma::mat D = arma::zeros<arma::mat>(n_stages, n_lines);

        arma::vec Di;
        for (OnePlant& p : plants) {
            for (uint32 i = 0; i < n_lines; i++) {
                Di = p.aphids[i].remove_dispersers(disp_prop);
                D.col(i) += Di;
            }
        }

        return D;

    }

    // Add dispersers from another field:
    void add_dispersers(const arma::mat& D) {

        uint32 n_lines = plants[0].aphids.size();

        double n_plants = static_cast<double>(plants.size());

        for (uint32 i = 0; i < n_lines; i++) {
            arma::vec DD = D.col(i) / n_plants;
            for (OnePlant& p : plants) {
                p.aphids[i].alates.X += DD;
            }
        }

        return;

    }

    // Calculate among-plant dispersal:
    inline void calc_dispersal(pcg32& eng) {

        // Dispersal from previous generation
        emigrants.fill(0);
        immigrants.fill(0);
        for (const OnePlant& p : plants) {
            p.calc_dispersal(emigrants, immigrants, eng);
        }
        return;
    }
    inline void calc_dispersal() {
        emigrants.fill(0);
        immigrants.fill(0);
        for (const OnePlant& p : plants) {
            p.calc_dispersal(emigrants, immigrants);
        }
        return;
    }

    inline void update(pcg32& eng) {
        /*
         Once `calc_dispersal` has updated inside `emigrants` and `immigrants`,
         we can update the populations using those dispersal numbers.
         */
        // First update info for wasps before iterating:
        double old_mums;
        set_wasp_info(old_mums);
        // Then we can update aphids and mummies:
        for (OnePlant& p : plants) {
            p.update(emigrants, immigrants, &wasps, eng);
        }
        // Lastly update adult wasps (if constant_wasps = false):
        if (!constant_wasps) {
            wasps.update(old_mums, eng);
            if (wasps.Y < extinct_N) wasps.Y = 0;
        }
        return;
    }


    // Clear plants by either a maximum age or total abundance
    void clear_plants(const uint32& max_age,
                       const double& clear_surv,
                       pcg32& eng);
    void clear_plants(const double& max_N,
                       const double& clear_surv,
                       pcg32& eng);


    // Total (living) aphids in patch
    inline double total_aphids() const {
        double ta = 0;
        for (const OnePlant& op : plants) {
            ta += op.total_aphids();
        }
        return ta;
    }


};










// ============================================================================
// ============================================================================
// ============================================================================
// ============================================================================






class AllFields {

    pcg32 eng;   // RNG

    uint32 max_age;
    double max_N;


public:

    std::vector<OneField> fields;

    double clear_surv;
    double alate_field_disp_p;
    double wasp_disp_m0;
    double wasp_disp_m1;
    std::vector<double> wasp_field_attract;
    bool disp_error;
    double extinct_N;

    // Total stages for all wasps, mummies, and aphids. Used for output.
    uint32 total_stages;

    AllFields()
        : eng(), max_age(), max_N(), fields(), clear_surv(),
          alate_field_disp_p(), wasp_disp_m0(), wasp_disp_m1(),
          wasp_field_attract(),
          disp_error(), extinct_N(), total_stages() {};

    AllFields(const AllFields& other)
        : eng(other.eng),
          max_age(other.max_age), max_N(other.max_N),
          fields(other.fields),
          clear_surv(other.clear_surv),
          alate_field_disp_p(other.alate_field_disp_p),
          wasp_disp_m0(other.wasp_disp_m0),
          wasp_disp_m1(other.wasp_disp_m1),
          wasp_field_attract(other.wasp_field_attract),
          disp_error(other.disp_error),
          extinct_N(other.extinct_N),
          total_stages(other.total_stages){};

    AllFields& operator=(const AllFields& other) {
        eng = other.eng;
        max_age = other.max_age;
        max_N = other.max_N;
        fields = other.fields;
        clear_surv = other.clear_surv;
        alate_field_disp_p = other.alate_field_disp_p;
        wasp_disp_m0 = other.wasp_disp_m0;
        wasp_disp_m1 = other.wasp_disp_m1;
        wasp_field_attract = other.wasp_field_attract;
        disp_error = other.disp_error;
        extinct_N = other.extinct_N;
        total_stages = other.total_stages;
        return *this;
    };


    AllFields(const uint32& n_fields,
              const uint32& max_age_,
              const double& max_N_,
              const std::vector<uint32>& wasp_delay,
              const double& sigma_x,
              const double& sigma_y,
              const double& rho,
              const bool& aphid_demog_error,
              const bool& wasp_demog_error,
              const double& mean_K,
              const double& sd_K,
              const std::vector<double>& K_y_mult,
              const double& wilted_prop,
              const double& shape1_wilted_mort,
              const double& shape2_wilted_mort,
              const arma::mat& attack_surv,
              const std::vector<std::string>& aphid_name,
              const std::vector<arma::cube>& leslie_mat,
              const std::vector<arma::cube>& aphid_density_0,
              const std::vector<double>& alate_b0,
              const std::vector<double>& alate_b1,
              const std::vector<double>& alate_plant_disp_p,
              const std::vector<double>& disp_mort,
              const std::vector<uint32>& disp_start,
              const bool& disp_error_,
              const std::vector<uint32>& living_days,
              const std::vector<double>& pred_rate,
              const double& extinct_N_,
              const arma::mat& mum_density_0,
              const double& mum_smooth,
              const double& max_mum_density,
              const arma::vec& rel_attack,
              const double& a,
              const double& k,
              const double& h,
              const std::vector<double>& wasp_density_0,
              const double& sex_ratio,
              const std::vector<double>& s_y,
              const std::vector<bool>& constant_wasps,
              const double& clear_surv_,
              const double& alate_field_disp_p_,
              const double& wasp_disp_m0_,
              const double& wasp_disp_m1_,
              const std::vector<double>& wasp_field_attract_,
              const std::vector<uint64>& seeds)
        : eng(), max_age(max_age_), max_N(max_N_),
          fields(),
          clear_surv(clear_surv_),
          alate_field_disp_p(alate_field_disp_p_),
          wasp_disp_m0(wasp_disp_m0_),
          wasp_disp_m1(wasp_disp_m1_),
          wasp_field_attract(wasp_field_attract_),
          disp_error(disp_error_),
          extinct_N(extinct_N_),
          total_stages(0) {

        seed_pcg(eng, seeds);

        //' Make sure `wasp_field_attract` sums to 1 (negative values
        //' and a sum <= 0 are already checked for in sim_gameofclones_cpp):
        double wfa_sum = std::accumulate(wasp_field_attract.begin(),
                                         wasp_field_attract.end(), 0.0);
        for (double& x : wasp_field_attract) x /= wfa_sum;

        fields.reserve(n_fields);
        for (uint32 i = 0; i < n_fields; i++) {
            fields.push_back(
                OneField(sigma_x, sigma_y, rho, aphid_demog_error, wasp_demog_error,
                         mean_K, sd_K, K_y_mult[i],
                         wilted_prop, shape1_wilted_mort, shape2_wilted_mort,
                         attack_surv, aphid_name, leslie_mat, aphid_density_0,
                         alate_b0, alate_b1, alate_plant_disp_p,
                         disp_mort, disp_start, living_days, pred_rate[i],
                         extinct_N, mum_density_0, mum_smooth, max_mum_density,
                         rel_attack, a, k, h,
                         wasp_density_0[i], wasp_delay[i],
                         sex_ratio,
                         s_y[i], constant_wasps[i], eng));
            total_stages += 1;
            const OneField& field(fields[i]);
            for (const OnePlant& plant : field.plants) {
                total_stages += plant.mummies.Y.n_elem;
                for (const AphidPop& aphid : plant.aphids) {
                    total_stages += aphid.apterous.X.n_elem;
                    total_stages += aphid.alates.X.n_elem;
                    total_stages += aphid.paras.X.n_elem;
                }
            }
        }

    };

    inline uint32 size() const noexcept {
        return fields.size();
    }

    OneField& operator[](const uint32& idx) {
        return fields[idx];
    }
    const OneField& operator[](const uint32& idx) const {
        return fields[idx];
    }

    uint32 n_plants() const {
        if (fields.size() == 0) return 0;
        return fields[0].plants.size();
    }
    uint32 n_lines() const {
        if (n_plants() == 0) return 0;
        return fields[0].plants[0].aphids.size();
    }

    void reseed(const std::vector<uint64>& seeds) {
        seed_pcg(eng, seeds);
        return;
    }


    // Alate dispersal across fields:
    void across_field_disp_alates() {
        if (fields.size() > 1 && alate_field_disp_p > 0) {
            arma::mat D = fields[0].remove_dispersers(alate_field_disp_p);
            for (uint32 i = 1; i < fields.size(); i++) {
                D += fields[i].remove_dispersers(alate_field_disp_p);
            }
            D /= static_cast<double>(fields.size());
            for (OneField& c : fields) c.add_dispersers(D);
        }
        return;
    }

    // Wasp dispersal across fields:
    void across_field_disp_wasps() {
        if (fields.size() > 1 && wasp_disp_m0 > 0) {
            double from_wasp_pool = 0;
            if (wasp_disp_m1 != 0) {
                double p_out, z;
                for (OneField& field : fields) {
                    z = field.total_aphids();
                    p_out = wasp_disp_m0 * std::exp(-wasp_disp_m1 * z);
                    from_wasp_pool += (field.wasps.Y * p_out);
                    field.wasps.Y *= (1 - p_out);
                }
            } else {
                for (OneField& field : fields) {
                    from_wasp_pool += (field.wasps.Y * wasp_disp_m0);
                    field.wasps.Y *= (1 - wasp_disp_m0);
                }
            }
            for (uint32 i = 0; i < fields.size(); i++) {
                OneField& field(fields[i]);
                field.wasps.Y += (from_wasp_pool * wasp_field_attract[i]);
            }
        }
        return;
    }

    void do_perturb(std::deque<PerturbInfo>& perturbs, const uint32& t);


    // Update for one time step.
    // Returns true if all fields/plants are empty
    bool update(const uint32& t,
                std::deque<PerturbInfo>& perturbs,
                std::deque<uint32>& check_for_clear);


    // Clear plants by age or total N thresholds (or neither)
    void clear_plants(const uint32& t,
                      std::deque<uint32>& check_for_clear);

    List to_list() const;

    // It's assumed this vector is in the same order as the outputs are above!
    // Make sure this happens from the R side.
    void from_vector(std::vector<double>& N);

    // Set new parameters
    void set_new_pars(const double& K_,
                      const std::vector<double>& alate_b0_,
                      const std::vector<double>& alate_b1_,
                      const double& alate_field_disp_p_,
                      const std::vector<double>& K_y_mult_,
                      const std::vector<double>& s_y_,
                      const double& a_,
                      const double& k_,
                      const double& h_,
                      const double& wasp_disp_m0_,
                      const double& wasp_disp_m1_,
                      const std::vector<double>& wasp_field_attract_,
                      const double& mum_smooth_,
                      const std::vector<double>& pred_rate_,
                      const uint32& max_plant_age_,
                      const double& clear_surv_);



};




#endif
