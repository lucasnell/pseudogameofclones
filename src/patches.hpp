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


/*
 Combine Leslie matrices for apterous and alates into one matrix
 */
inline void combine_leslies(arma::mat& L,
                            const arma::mat& apterous,
                            const arma::mat& alates,
                            const double& alate_prop,
                            const double& disp_rate,
                            const double& disp_mort,
                            const uint32& disp_start) {

    // Should be the rows and columns for both matrices:
    arma::uword n = apterous.n_rows;

    L = arma::zeros<arma::mat>(2 * n, 2 * n);

    L(arma::span(0, n - 1), arma::span(0, n - 1)) = apterous;
    L(arma::span(n, 2 * n - 1), arma::span(n, 2 * n - 1)) = alates;

    L.row(0).head(n) *= (1 - alate_prop);
    L.row(n).head(n) = apterous.row(0) * alate_prop;

    L.row(0).tail(n) = alates.row(0);
    L.row(n).tail(n).fill(0);


    for (uint32 i = disp_start; i < (n-1); i++) {
        L(n+i+1, n+i) *= (1 - disp_mort * disp_rate);
    }

    return;
}






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


    // Updates total # unparasitized aphids (`x`), # living aphids (`z`), and
    // checks to see if plant should be wilted.
    void update_x_z_wilted();

    MEMBER(bool,wilted)


public:

    std::vector<AphidPop> aphids;   // aphids in this patch
    WaspPop wasps;                  // wasps in this patch
    bool empty;                     // whether no aphids/wasps are on this patch
    double pred_rate;               // predation on aphids
    double K;                       // unparasitized aphid carrying capacity
    double K_y;                     // parasitized aphid carrying capacity
    double z = 0;                   // sum of all living aphids at time t
    double x = 0;                   // Sum of unparasitized aphids at time t
    double S = 0;                   // effect of density dependence on aphids
    double S_y = 0;                 // effect of dd on parasitized aphids
    uint32 n_patches;               // total # patches
    uint32 this_j;                  // index for this patch
    uint32 age = 0;                 // age of this patch
    double death_prop;              // proportion of carrying capacity that kills plant
    double death_mort;              // growth-rate modifier once plants start dying
    double extinct_N;               // threshold for calling an aphid line extinct



    OnePatch()
        : wilted_(false), aphids(), wasps(), empty(true), pred_rate(0), K(0), K_y(1),
          n_patches(1), this_j(0), death_prop(1), death_mort(1), extinct_N() {};

    /*
     In `aphid_density_0` below, rows are aphid stages, columns are types (alate vs
     apterous), and slices are aphid lines.
     In `leslie_mat` below, items in vector are aphid lines, slices are alate/apterous.
     */
    OnePatch(const double& sigma,
             const double& rho,
             const double& demog_mult,
             const std::vector<arma::vec>& attack_surv_,
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
             const double& pred_rate_,
             const uint32& n_patches_,
             const uint32& this_j_,
             const double& extinct_N_,
             const arma::vec& rel_attack_,
             const double& a_,
             const double& k_,
             const double& h_,
             const arma::uvec& mum_days_,
             const double& wasp_density_0_,
             const double& sex_ratio_,
             const double& K_y_,
             const double& s_y_)
        : wilted_(false),
          aphids(),
          wasps(rel_attack_, a_, k_, h_, mum_days_,
                wasp_density_0_, sex_ratio_, K_y_, s_y_),
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
            AphidPop ap(aphid_name[i], sigma, rho, demog_mult, attack_surv_[i],
                        leslie_mat[i], aphid_density_0.slice(i),
                        alate_b0[i], alate_b1[i], disp_rate[i], disp_mort[i],
                        disp_start[i]);
            aphids.push_back(ap);
            double N = aphids.back().total_aphids();
            if (N < extinct_N) {
                aphids.back().clear();
            } else if (empty) empty = false;
        }

    };



    OnePatch(const OnePatch& other)
        : wilted_(false), aphids(other.aphids), wasps(other.wasps), empty(other.empty), pred_rate(other.pred_rate),
          K(other.K), K_y(other.K_y), z(other.z), x(other.x),
          S(other.S), S_y(other.S_y), n_patches(other.n_patches),
          this_j(other.this_j), age(other.age), death_prop(other.death_prop),
          death_mort(other.death_mort), extinct_N(other.extinct_N) {};

    OnePatch& operator=(const OnePatch& other) {
        wilted_ = other.wilted_;
        aphids = other.aphids;
        wasps = other.wasps;
        empty = other.empty;
        pred_rate = other.pred_rate;
        K = other.K;
        K_y = other.K_y;
        z = other.z;
        x = other.x;
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
     Clear to no aphids or wasps
     */
    void clear(const double& K_,
               const double& K_y_,
               const double& death_mort_) {
        for (AphidPop& ap : aphids) ap.clear();
        wasps.clear();  // this kills mummies only
        empty = true;
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
                pcg32& eng);
    // Same but minus stochasticity
    void update(const arma::cube& emigrants,
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

    // For new death_morts if desired:
    beta_distribution beta_distr;

    double mean_K_;                 // mean of distribution of `K` for plants
    double sd_K_;                   // sd of distribution of `K` for plants

    // Helps calculate carrying capacity for parasitized wasps
    // (`K_y = K * K_y_mult`):
    double K_y_mult;

    double shape1_death_mort_;      // shape1 of distribution of `death_mort` for plants
    double shape2_death_mort_;      // shape2 of distribution of `death_mort` for plants

    // wasp dispersal rate
    double wasp_disp;


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
        : tnorm_distr(), beta_distr(), mean_K_(), sd_K_(), K_y_mult(),
          shape1_death_mort_(), shape2_death_mort_(), wasp_disp(),
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
               const double& K_y_mult_,
               const double& death_prop,
               const double& shape1_death_mort,
               const double& shape2_death_mort,
               const std::vector<arma::vec>& attack_surv_,
               const std::vector<std::string>& aphid_name,
               const std::vector<arma::cube>& leslie_mat,
               const std::vector<arma::cube>& aphid_density_0,
               const std::vector<double>& alate_b0,
               const std::vector<double>& alate_b1,
               const std::vector<double>& disp_rate,
               const std::vector<double>& disp_mort,
               const std::vector<uint32>& disp_start,
               const double& wasp_disp_,
               const std::vector<double>& pred_rate,
               const double& extinct_N,
               const arma::vec& rel_attack_,
               const double& a_,
               const double& k_,
               const double& h_,
               const arma::uvec& mum_days_,
               const double& wasp_density_0_,
               const double& sex_ratio_,
               const double& K_y_,
               const double& s_y_,
               pcg32& eng)
        : tnorm_distr(),
          beta_distr(),
          mean_K_(mean_K),
          sd_K_(sd_K),
          K_y_mult(K_y_mult_),
          shape1_death_mort_(shape1_death_mort),
          shape2_death_mort_(shape2_death_mort),
          wasp_disp(wasp_disp_),
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
            OnePatch ap(sigma, rho, demog_mult, attack_surv_, K, K_y, death_prop, death_mort,
                        aphid_name, leslie_mat,
                        aphid_density_0[j], alate_b0, alate_b1, disp_rate, disp_mort,
                        disp_start, pred_rate[j], n_patches, j, extinct_N,
                        rel_attack_, a_, k_, h_, mum_days_,
                        wasp_density_0_, sex_ratio_, K_y_, s_y_);
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
          patches(other.patches),
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

    inline void update(pcg32& eng) {
        /*
         Once `calc_dispersal` has updated inside `emigrants` and `immigrants`, we
         can update the populations using those dispersal numbers.
         */
        double mean_wasps = 0;
        for (OnePatch& p : patches) {
            p.update(emigrants, immigrants, eng);
            mean_wasps += p.wasps.total_wasps();
        }
        mean_wasps /= static_cast<double>(patches.size());
        // dispersal of adult wasps comes last:
        for (OnePatch& p : patches) {
            p.wasps.Y.back() *= (1 - wasp_disp);
            p.wasps.Y.back() += wasp_disp * mean_wasps;
        }
        return;
    }
    inline void update() {
        double mean_wasps = 0;
        for (OnePatch& p : patches) {
            p.update(emigrants, immigrants);
            mean_wasps += p.wasps.total_wasps();
        }
        mean_wasps /= static_cast<double>(patches.size());
        // dispersal of adult wasps comes last:
        for (OnePatch& p : patches) {
            p.wasps.Y.back() *= (1 - wasp_disp);
            p.wasps.Y.back() += wasp_disp * mean_wasps;
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
