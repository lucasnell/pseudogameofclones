# ifndef __CLONEWARS_PATCHES_H
# define __CLONEWARS_PATCHES_H


#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "clonewars_types.hpp"  // integer types
#include "math.hpp"             // distributions
#include "aphids.hpp"           // aphid classes
#include "wasps.hpp"            // wasp classes
#include "pcg.hpp"              // runif_ fxns



using namespace Rcpp;





// Info for a single perturbation (it happens to all patches)
struct PerturbInfo {

    // When to perturb.
    uint32 time;
    // Number to multiply abundance by.
    double multiplier;
    /*
     Who to perturb.
     Indices less than the # aphids points to the particular
     aphid line, equal to the # aphids is mummies, and greater than
     the # aphids is adult wasps.
     */
    uint32 index;

    PerturbInfo() : time(), multiplier(), index() {}
    PerturbInfo(const uint32& when, const uint32& who, const double& how)
        : time(when), multiplier(how), index(who) {}
    PerturbInfo& operator=(const PerturbInfo& other) {
        time = other.time;
        multiplier = other.multiplier;
        index = other.index;
        return *this;
    }

};




// Info for clearing patches
template <typename T>
struct PatchClearingInfo {

    uint32 ind;
    double N;
    bool wilted;
    T thresh_info; // age or # aphids

    PatchClearingInfo(const uint32& ind_,
                      const double& N_,
                      const bool& wilted_,
                      const T& thresh_info_)
        : ind(ind_), N(N_), wilted(wilted_), thresh_info(thresh_info_) {}
    PatchClearingInfo(const PatchClearingInfo<T>& other)
        : ind(other.ind), N(other.N), thresh_info(other.thresh_info) {}

    // This returns true when `other` should be cleared first
    bool operator<(const PatchClearingInfo<T>& other) const {
        if (other.wilted && !wilted) return true;
        if (!other.wilted && wilted) return false;
        return thresh_info < other.thresh_info;
    }
    // This returns true when `this` should be cleared first
    bool operator>(const PatchClearingInfo<T>& other) const {
        if (other.wilted && !wilted) return false;
        if (!other.wilted && wilted) return true;
        return thresh_info > other.thresh_info;
    }

};





/*
One patch of continuous habitat.
*/
class OnePatch {


    /*
     Carrying capacity for patch.
     It depends on the Leslie matrix for each line's apterous and alates.
     I'm not including parasitized aphids because they shouldn't be too numerous.
     */
    double carrying_capacity() const;


    /*
     Adjust for potential extinction or re-colonization:
     */
    void extinct_colonize(const uint32& i);


    // Updates total # living aphids in this patch (`z`), and
    // checks to see if plant should be wilted.
    void update_z_wilted();

    MEMBER(bool,wilted)


public:

    std::vector<AphidPop> aphids;   // aphids in this patch
    MummyPop mummies;               // mummies in this patch
    bool empty;                     // whether no aphids are on this patch
    double pred_rate;               // predation on aphids
    double K;                       // unparasitized aphid carrying capacity
    double K_y;                     // parasitized aphid carrying capacity
    double z = 0;                   // sum of all living aphids at time t
    double S = 0;                   // effect of density dependence on aphids
    double S_y = 0;                 // effect of dd on parasitized aphids
    uint32 n_patches;               // total # patches
    uint32 this_j;                  // index for this patch
    uint32 age = 0;                 // age of this patch
    double death_prop;              // proportion of carrying capacity that kills plant
    double death_mort;              // growth-rate modifier once plants start dying
    double extinct_N;               // threshold for calling an aphid line extinct



    OnePatch()
        : wilted_(false), aphids(), mummies(), empty(true), pred_rate(0),
          K(0), K_y(1),
          n_patches(1), this_j(0), death_prop(1), death_mort(1), extinct_N() {};

    /*
     In `aphid_density_0` below, rows are aphid stages, columns are types (alate vs
     apterous), and slices are aphid lines.
     In `leslie_mat` below, items in vector are aphid lines, slices are
     alate/apterous/parasitized.
     */
    OnePatch(const double& sigma_x,
             const double& rho,
             const double& demog_mult,
             const arma::mat& attack_surv_,
             const double& K_,
             const double& K_y_,
             const double& death_prop_,
             const double& death_mort_,
             const std::vector<std::string>& aphid_name,
             const std::vector<arma::cube>& leslie_mat,
             const arma::cube& aphid_density_0,
             const std::vector<double>& alate_b0,
             const std::vector<double>& alate_b1,
             const std::vector<double>& disp_rate,
             const std::vector<double>& disp_mort,
             const std::vector<uint32>& disp_start,
             const std::vector<uint32>& living_days,
             const double& pred_rate_,
             const uint32& n_patches_,
             const uint32& this_j_,
             const double& extinct_N_,
             const arma::vec& mum_density_0)
        : wilted_(false),
          aphids(),
          mummies(mum_density_0),
          empty(true),
          pred_rate(pred_rate_),
          K(K_),
          K_y(K_y_),
          n_patches(n_patches_),
          this_j(this_j_),
          death_prop(death_prop_),
          death_mort(death_mort_),
          extinct_N(extinct_N_) {

        uint32 n_lines = aphid_name.size();

        aphids.reserve(n_lines);

        for (uint32 i = 0; i < n_lines; i++) {
            AphidPop ap(aphid_name[i], sigma_x, rho, demog_mult,
                        attack_surv_.col(i),
                        leslie_mat[i], aphid_density_0.slice(i),
                        alate_b0[i], alate_b1[i], disp_rate[i], disp_mort[i],
                        disp_start[i], living_days[i]);
            aphids.push_back(ap);
            double N = aphids.back().total_aphids();
            if (N < extinct_N) {
                aphids.back().clear();
            } else if (empty) empty = false;
        }

    };



    OnePatch(const OnePatch& other)
        : wilted_(false), aphids(other.aphids), mummies(other.mummies),
          empty(other.empty), pred_rate(other.pred_rate),
          K(other.K), K_y(other.K_y), z(other.z),
          S(other.S), S_y(other.S_y), n_patches(other.n_patches),
          this_j(other.this_j), age(other.age), death_prop(other.death_prop),
          death_mort(other.death_mort), extinct_N(other.extinct_N) {};

    OnePatch& operator=(const OnePatch& other) {
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
        n_patches = other.n_patches;
        this_j = other.this_j;
        age = other.age;
        death_prop = other.death_prop;
        death_mort = other.death_mort;
        extinct_N = other.extinct_N;
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
               const double& death_mort_) {
        for (AphidPop& ap : aphids) ap.clear();
        mummies.clear();
        empty = true;
        wilted_ = false;
        age = 0;
        K = K_;
        K_y = K_y_;
        death_mort = death_mort_;
        return;
    }
    /*
     Clear some of the aphids and mummies
     */
    void clear(const double& K_,
               const double& K_y_,
               const double& death_mort_,
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
        death_mort = death_mort_;
        return;
    }

    // Total (living) aphids on patch
    inline double total_aphids() const {
        double ta = 0;
        for (const AphidPop& ap : aphids) {
            ta += ap.total_aphids();
        }
        return ta;
    }
    // Total UNparasitized aphids on patch
    inline double total_unpar_aphids() const {
        double ta = 0;
        for (const AphidPop& ap : aphids) {
            ta += ap.apterous.total_aphids() + ap.alates.total_aphids();
        }
        return ta;
    }
    // Total mummies on patch
    inline double total_mummies() const {
        double tm = arma::accu(mummies.Y);
        return tm;
    }

    /*
     Add dispersal info to `emigrants` and `immigrants` cubes.
     In these cubes, rows are aphid stages, columns are patches,
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

    // Same thing as above, but overloaded for not including dispersal stochasticity
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
    // Same but minus stochasticity
    void update(const arma::cube& emigrants,
                const arma::cube& immigrants,
                const WaspPop* wasps);


};




/*
 Class for all patches.

 For the constructors below, all vector arguments should have a length equal to the
 number of aphid lines, except for `K`, `aphid_density_0`, and `pred_rate_`;
 these should have a length equal to the number of patches.
 Each item in `aphid_density_0` should have a length equal to the number of aphid lines.

 */

class AllPatches {

    // For new Ks if desired:
    trunc_normal_distribution tnorm_distr;

    // For new death_morts if desired:
    beta_distribution beta_distr;

    double mean_K_;                 // mean of distribution of `K` for plants
    double sd_K_;                   // sd of distribution of `K` for plants

    // Helps calculate carrying capacity for parasitized wasps
    // (`K_y = K * K_y_mult`):
    double K_y_mult;

    double shape1_death_mort_;      // shape1 of distribution of `death_mort` for plants
    double shape2_death_mort_;      // shape2 of distribution of `death_mort` for plants

    double extinct_N;               // used here for the wasps


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
    // Set plant-death mortality parameter
    void set_death_mort(double& death_mort, pcg32& eng) {
        if (shape2_death_mort_ == 0) {
            death_mort = shape1_death_mort_;
        } else {
            death_mort = beta_distr(eng);
        }
        return;
    }

    inline void set_wasp_info(double& old_mums) {
        wasps.x = 0;
        old_mums = 0;
        for (OnePatch& p : patches) {
            wasps.x += p.total_unpar_aphids();
            old_mums += p.mummies.Y.back();
        }
        return;
    }


    // Do the actual clearing of patches while avoiding extinction
    template <typename T>
    void do_clearing(std::vector<PatchClearingInfo<T>>& clear_patches,
                     double& remaining,
                     std::vector<bool>& wilted,
                     const double& clear_surv,
                     pcg32& eng);


public:

    std::vector<OnePatch> patches;
    WaspPop wasps;
    arma::cube emigrants;
    arma::cube immigrants;


    AllPatches()
        : tnorm_distr(), beta_distr(), mean_K_(), sd_K_(), K_y_mult(),
          shape1_death_mort_(), shape2_death_mort_(), extinct_N(),
          patches(), wasps(), emigrants(), immigrants() {};

    /*
     In `aphid_density_0` below, rows are aphid stages, columns are types (alate vs
     apterous), and slices are aphid lines.
     In `leslie_mat` below, slices are aphid lines.
     */
    AllPatches(const double& sigma_x,
               const double& sigma_y,
               const double& rho,
               const double& demog_mult,
               const double& mean_K,
               const double& sd_K,
               const double& K_y_mult_,
               const double& death_prop,
               const double& shape1_death_mort,
               const double& shape2_death_mort,
               const arma::mat& attack_surv_,
               const std::vector<std::string>& aphid_name,
               const std::vector<arma::cube>& leslie_mat,
               const std::vector<arma::cube>& aphid_density_0,
               const std::vector<double>& alate_b0,
               const std::vector<double>& alate_b1,
               const std::vector<double>& disp_rate,
               const std::vector<double>& disp_mort,
               const std::vector<uint32>& disp_start,
               const std::vector<uint32>& living_days,
               const std::vector<double>& pred_rate,
               const double& extinct_N_,
               const arma::mat& mum_density_0,
               const arma::vec& rel_attack_,
               const double& a_,
               const double& k_,
               const double& h_,
               const double& wasp_density_0_,
               const double& sex_ratio_,
               const double& s_y_,
               pcg32& eng)
        : tnorm_distr(),
          beta_distr(),
          mean_K_(mean_K),
          sd_K_(sd_K),
          K_y_mult(K_y_mult_),
          shape1_death_mort_(shape1_death_mort),
          shape2_death_mort_(shape2_death_mort),
          extinct_N(extinct_N_),
          patches(),
          wasps(rel_attack_, a_, k_, h_, wasp_density_0_,
                sex_ratio_, s_y_, sigma_y),
          emigrants(),
          immigrants() {


        /*
         None of these hyperparameters can be <= 0, so if they're <= then that makes
         the associated parameter whose distribution it describes (e.g., K for `sd_K_`)
         not have stochasticity.
         This means that the first hyperparameter provided will be the value
         used for ALL of that parameter WITHOUT BEING TRANSFORMED.
         */
        if (sd_K_ <= 0) {
            sd_K_ = 0;
        } else {
            tnorm_distr = trunc_normal_distribution(mean_K_, sd_K_);
        }
        if (shape2_death_mort_ <= 0) {
            shape2_death_mort_ = 0;
        } else {
            beta_distr = beta_distribution(shape1_death_mort_,
                                           shape2_death_mort_);
        }

        uint32 n_patches = aphid_density_0.size();
        // (We know pred_rate.size() == n_patches bc it's check inside sim_clonewars_cpp)

        uint32 n_lines = aphid_name.size();
        uint32 n_stages = leslie_mat.front().n_rows;

        double K, K_y, death_mort;
        patches.reserve(n_patches);
        for (uint32 j = 0; j < n_patches; j++) {
            set_K(K, K_y, eng);
            set_death_mort(death_mort, eng);
            OnePatch ap(sigma_x, rho, demog_mult, attack_surv_,
                        K, K_y, death_prop, death_mort,
                        aphid_name, leslie_mat,
                        aphid_density_0[j], alate_b0, alate_b1, disp_rate, disp_mort,
                        disp_start, living_days, pred_rate[j], n_patches, j, extinct_N_,
                        mum_density_0.col(j));
            patches.push_back(ap);
        }

        emigrants = arma::zeros<arma::cube>(n_stages, n_patches, n_lines);
        immigrants = arma::zeros<arma::cube>(n_stages, n_patches, n_lines);

    }

    AllPatches(const AllPatches& other)
        : tnorm_distr(other.tnorm_distr),
          beta_distr(other.beta_distr),
          mean_K_(other.mean_K_),
          sd_K_(other.sd_K_),
          K_y_mult(other.K_y_mult),
          shape1_death_mort_(other.shape1_death_mort_),
          shape2_death_mort_(other.shape2_death_mort_),
          extinct_N(other.extinct_N),
          patches(other.patches),
          wasps(other.wasps),
          emigrants(other.emigrants),
          immigrants(other.immigrants) {};

    AllPatches& operator=(const AllPatches& other) {
        tnorm_distr = other.tnorm_distr;
        beta_distr = other.beta_distr;
        mean_K_ = other.mean_K_;
        sd_K_ = other.sd_K_;
        K_y_mult = other.K_y_mult;
        shape1_death_mort_ = other.shape1_death_mort_;
        shape2_death_mort_ = other.shape2_death_mort_;
        extinct_N = other.extinct_N;
        patches = other.patches;
        wasps = other.wasps;
        emigrants = other.emigrants;
        immigrants = other.immigrants;
        return *this;
    };


    inline uint32 size() const noexcept {
        return patches.size();
    }

    OnePatch& operator[](const uint32& idx) {
        return patches[idx];
    }
    const OnePatch& operator[](const uint32& idx) const {
        return patches[idx];
    }

    inline void calc_dispersal(pcg32& eng) {

        // Dispersal from previous generation
        emigrants.fill(0);
        immigrants.fill(0);
        for (const OnePatch& p : patches) {
            p.calc_dispersal(emigrants, immigrants, eng);
        }
        return;
    }
    inline void calc_dispersal() {
        emigrants.fill(0);
        immigrants.fill(0);
        for (const OnePatch& p : patches) {
            p.calc_dispersal(emigrants, immigrants);
        }
        return;
    }

    inline void update(pcg32& eng) {
        /*
         Once `calc_dispersal` has updated inside `emigrants` and `immigrants`, we
         can update the populations using those dispersal numbers.
         */
        // First update info for wasps before iterating:
        double old_mums;
        set_wasp_info(old_mums);
        // Then we can update aphids and mummies:
        for (OnePatch& p : patches) {
            p.update(emigrants, immigrants, &wasps, eng);
        }
        // Lastly update adult wasps:
        wasps.update(old_mums, eng);
        if (wasps.Y < extinct_N) wasps.Y = 0;
        return;
    }
    inline void update() {
        double old_mums;
        set_wasp_info(old_mums);
        for (OnePatch& p : patches) {
            p.update(emigrants, immigrants, &wasps);
        }
        wasps.update(old_mums);
        if (wasps.Y < extinct_N) wasps.Y = 0;
        return;
    }


    // Clear patches by either a maximum age or total abundance
    void clear_patches(const uint32& max_age,
                       const double& clear_surv,
                       pcg32& eng);
    void clear_patches(const double& max_N,
                       const double& clear_surv,
                       pcg32& eng);


};





#endif
