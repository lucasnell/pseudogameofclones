# ifndef __CLONEWARS_PATCHES_H
# define __CLONEWARS_PATCHES_H


#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "clonewars_types.hpp"  // integer types
#include "math.hpp"             // distributions
#include "aphids.hpp"           // aphid classes
#include "pcg.hpp"              // runif_ fxns



using namespace Rcpp;





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
     Adjust for potential extinction or re-colonization:
     */
    void extinct_colonize(const uint32& i);



public:

    std::vector<AphidPop> aphids;   // aphids in this patch
    bool empty;                     // boolean for whether no aphids are on this patch
    double pred_rate;               // predation on aphids
    double K;                       // aphid carrying capacity
    double z = 0;                   // sum of all aphids at time t
    double S = 0;                   // effect of density dependence, based on z and K
    uint32 n_patches;               // total # patches
    uint32 this_j;                  // index for this patch
    uint32 age = 0;                 // age of this patch
    uint32 death_age;               // age at which plant starts to die
    double death_mort;              // growth-rate modifier once plants start dying
    double extinct_N;               // threshold for calling an aphid line extinct


    OnePatch()
        : aphids(), empty(true), pred_rate(0), K(0),
          n_patches(1), this_j(0), death_age(0), death_mort(1), extinct_N() {};

    /*
     In `aphid_density_0` below, rows are aphid stages, columns are types (alate vs
     apterous), and slices are aphid lines.
     In `leslie_mat` below, items in vector are aphid lines, slices are alate/apterous.
     */
    OnePatch(const double& sigma,
             const double& rho,
             const double& demog_mult,
             const double& K_,
             const uint32& death_age_,
             const double& death_mort_,
             const std::vector<std::string>& aphid_name,
             const std::vector<arma::cube>& leslie_mat,
             const arma::cube& aphid_density_0,
             const std::vector<double>& alate_prop,
             const std::vector<double>& disp_rate,
             const std::vector<double>& disp_mort,
             const std::vector<uint32>& disp_start,
             const double& pred_rate_,
             const uint32& n_patches_,
             const uint32& this_j_,
             const double& extinct_N_)
        : aphids(),
          empty(true),
          pred_rate(pred_rate_),
          K(K_),
          n_patches(n_patches_),
          this_j(this_j_),
          death_age(death_age_),
          death_mort(death_mort_),
          extinct_N(extinct_N_) {

        uint32 n_lines = aphid_name.size();

        aphids.reserve(n_lines);

        for (uint32 i = 0; i < n_lines; i++) {
            AphidPop ap(aphid_name[i], sigma, rho, demog_mult,
                        leslie_mat[i], aphid_density_0.slice(i),
                        alate_prop[i], disp_rate[i], disp_mort[i],
                        disp_start[i]);
            aphids.push_back(ap);
            double N = aphids.back().total_aphids();
            if (N < extinct_N) {
                aphids.back().clear();
            } else if (empty) empty = false;
        }

    };



    OnePatch(const OnePatch& other)
        : aphids(other.aphids), empty(other.empty), pred_rate(other.pred_rate),
          K(other.K), z(other.z), S(other.S), n_patches(other.n_patches),
          this_j(other.this_j), age(other.age), death_age(other.death_age),
          death_mort(other.death_mort), extinct_N(other.extinct_N) {};

    OnePatch& operator=(const OnePatch& other) {
        aphids = other.aphids;
        empty = other.empty;
        pred_rate = other.pred_rate;
        K = other.K;
        z = other.z;
        S = other.S;
        n_patches = other.n_patches;
        this_j = other.this_j;
        age = other.age;
        death_age = other.death_age;
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

    const bool wilted() const {
        return age > death_age;
    }

    /*
     Clear to no aphids
     */
    void clear(const uint32& death_age_,
               const double& K_,
               const double& death_mort_) {
        for (AphidPop& ap : aphids) ap.clear();
        empty = true;
        age = 0;
        death_age = death_age_;
        K = K_;
        death_mort = death_mort_;
        return;
    }

    // Total aphids on patch
    inline double total_aphids() const {
        double ta = 0;
        for (const AphidPop& ap : aphids) {
            ta += ap.total_aphids();
        }
        return ta;
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
            aphids[i].calc_dispersal(this, emigrants.slice(i), immigrants.slice(i), eng);
        }
        return;
    }

    // Same thing as above, but overloaded for not including dispersal stochasticity
    void calc_dispersal(arma::cube& emigrants,
                        arma::cube& immigrants) const {
        for (uint32 i = 0; i < aphids.size(); i++) {
            aphids[i].calc_dispersal(this, emigrants.slice(i), immigrants.slice(i));
        }
        return;
    }


    /*
     Iterate one time step, after calculating dispersal numbers
     */
    void update_pops(const arma::cube& emigrants,
                     const arma::cube& immigrants,
                     pcg32& eng);
    // Same but minus stochasticity
    void update_pops(const arma::cube& emigrants,
                     const arma::cube& immigrants);


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

    // For new death_ages if desired:
    std::lognormal_distribution<double> lnorm_distr;

    // For new death_morts if desired:
    beta_distribution beta_distr;

    double mean_K_;                 // mean of distribution of `K` for plants
    double sd_K_;                   // sd of distribution of `K` for plants

    double meanlog_death_age_;      // meanlog of distribution of `death_age` for plants
    double sdlog_death_age_;        // sdlog of distribution of `death_age` for plants

    double shape1_death_mort_;      // shape1 of distribution of `death_mort` for plants
    double shape2_death_mort_;      // shape2 of distribution of `death_mort` for plants

    double get_K(pcg32& eng) {
        double K;
        if (sd_K_ == 0) {
            K = mean_K_;
        } else {
            K = tnorm_distr(eng);
        }
        return K;
    }
    uint32 get_death_age(pcg32& eng) {
        uint32 death_age;
        if (sdlog_death_age_ == 0) {
            death_age = static_cast<uint32>(std::round(meanlog_death_age_));
        } else {
            death_age = static_cast<uint32>(std::round(lnorm_distr(eng)));
        }
        return death_age;
    }
    double get_death_mort(pcg32& eng) {
        double death_mort;
        if (shape2_death_mort_ == 0) {
            death_mort = shape1_death_mort_;
        } else {
            death_mort = beta_distr(eng);
        }
        return death_mort;
    }


    // Do the actual clearing of patches while avoiding extinction
    template <typename T>
    void do_clearing(std::vector<PatchClearingInfo<T>>& clear_patches,
                     double& remaining,
                     std::vector<bool>& wilted,
                     pcg32& eng);


public:

    std::vector<OnePatch> patches;
    arma::cube emigrants;
    arma::cube immigrants;


    AllPatches()
        : mean_K_(), sd_K_(),
          meanlog_death_age_(), sdlog_death_age_(),
          shape1_death_mort_(), shape2_death_mort_(),
          patches(), emigrants(), immigrants() {};

    /*
     In `aphid_density_0` below, rows are aphid stages, columns are types (alate vs
     apterous), and slices are aphid lines.
     In `leslie_mat` below, slices are aphid lines.
     */
    AllPatches(const double& sigma,
               const double& rho,
               const double& demog_mult,
               const double& mean_K,
               const double& sd_K,
               const double& meanlog_death_age,
               const double& sdlog_death_age,
               const double& shape1_death_mort,
               const double& shape2_death_mort,
               const std::vector<std::string>& aphid_name,
               const std::vector<arma::cube>& leslie_mat,
               const std::vector<arma::cube>& aphid_density_0,
               const std::vector<double>& alate_prop,
               const std::vector<double>& disp_rate,
               const std::vector<double>& disp_mort,
               const std::vector<uint32>& disp_start,
               const std::vector<double>& pred_rate,
               const double& extinct_N,
               pcg32& eng)
        : mean_K_(mean_K),
          sd_K_(sd_K),
          meanlog_death_age_(meanlog_death_age),
          sdlog_death_age_(sdlog_death_age),
          shape1_death_mort_(shape1_death_mort),
          shape2_death_mort_(shape2_death_mort),
          patches(),
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
        if (sdlog_death_age_ <= 0) {
            sdlog_death_age_ = 0;
        } else {
            lnorm_distr = std::lognormal_distribution<double>(meanlog_death_age_,
                                                              sdlog_death_age_);
        }
        if (shape2_death_mort_ <= 0) {
            shape2_death_mort_ = 0;
        } else {
            beta_distr = beta_distribution(shape1_death_mort_,
                                           shape2_death_mort_);
        }

        uint32 n_patches = std::min(aphid_density_0.size(), pred_rate.size());
        uint32 n_lines = aphid_name.size();
        uint32 n_stages = leslie_mat.front().n_rows;

        patches.reserve(n_patches);
        for (uint32 j = 0; j < n_patches; j++) {
            double K = get_K(eng);
            uint32 death_age = get_death_age(eng);
            double death_mort = get_death_mort(eng);
            OnePatch ap(sigma, rho, demog_mult, K, death_age, death_mort,
                        aphid_name, leslie_mat,
                        aphid_density_0[j], alate_prop, disp_rate, disp_mort,
                        disp_start, pred_rate[j], n_patches, j, extinct_N);
            patches.push_back(ap);
        }

        emigrants = arma::zeros<arma::cube>(n_stages, n_patches, n_lines);
        immigrants = arma::zeros<arma::cube>(n_stages, n_patches, n_lines);

    }

    AllPatches(const AllPatches& other)
        : mean_K_(other.mean_K_),
          sd_K_(other.sd_K_),
          meanlog_death_age_(other.meanlog_death_age_),
          sdlog_death_age_(other.sdlog_death_age_),
          shape1_death_mort_(other.shape1_death_mort_),
          shape2_death_mort_(other.shape2_death_mort_),
          patches(other.patches),
          emigrants(other.emigrants),
          immigrants(other.immigrants) {};

    AllPatches& operator=(const AllPatches& other) {
        mean_K_ = other.mean_K_;
        sd_K_ = other.sd_K_;
        meanlog_death_age_ = other.meanlog_death_age_;
        sdlog_death_age_ = other.sdlog_death_age_;
        shape1_death_mort_ = other.shape1_death_mort_;
        shape2_death_mort_ = other.shape2_death_mort_;
        patches = other.patches;
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

    inline void update_pops(pcg32& eng) {
        /*
         Once `calc_dispersal` has updated inside `emigrants` and `immigrants`, we
         can update the populations using those dispersal numbers.
         */
        for (OnePatch& p : patches) {
            p.update_pops(emigrants, immigrants, eng);
        }
        return;
    }
    inline void update_pops() {
        for (OnePatch& p : patches) {
            p.update_pops(emigrants, immigrants);
        }
        return;
    }


    // Clear patches by either a maximum age or total abundance
    void clear_patches(const uint32& max_age,
                       pcg32& eng);
    void clear_patches(const double& max_N,
                       pcg32& eng);


};





#endif
