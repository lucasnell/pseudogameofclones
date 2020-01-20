# ifndef __CLONEWARS_PATCHES_H
# define __CLONEWARS_PATCHES_H


#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "clonewars_types.hpp"  // integer types
#include "aphids.hpp"           // aphid classes
#include "pcg.hpp"              // runif_ fxns



using namespace Rcpp;







/*
One patch of continuous habitat.
*/
class OnePatch {


    /*
     Adjust for potential extinction or re-colonization:
     */
    void extinct_colonize(const uint32& i) {
        AphidPop& ap(aphids[i]);
        double N = ap.total_aphids();
        if (N >= extinct_N) {
            // Newly colonized:
            if (empty) empty = false;
            if (ap.extinct) ap.extinct = false;
        }
        // Newly extinct:
        if (N < extinct_N && !ap.extinct) {
            ap.clear();
        }
        return;
    }

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
    double extinct_N;               // threshold for calling an aphid line extinct


    OnePatch()
        : aphids(), empty(true), pred_rate(0), K(0),
          n_patches(1), this_j(0), extinct_N() {};

    /*
     In `aphid_density_0` below, rows are aphid stages, columns are types (alate vs
     apterous), and slices are aphid lines.
     In `leslie_mat` below, slices are aphid lines.
     */
    OnePatch(const double& sigma,
             const double& rho,
             const double& demog_mult,
             const double& K_,
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
          this_j(other.this_j), age(other.age), extinct_N(other.extinct_N) {};

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
     Clear to no aphids
     */
    void clear() {
        for (AphidPop& ap : aphids) ap.clear();
        empty = true;
        age = 0;
        return;
    }
    // Same, then reset `K`
    void clear(const double& K_) {
        clear();
        K = K_;
        return;
    }

    // Total aphids on patch
    inline double total_aphids() const {
        double ta = 0;
        if (!empty) {
            for (const AphidPop& ap : aphids) {
                ta += ap.total_aphids();
            }
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
                     pcg32& eng) {

        z = 0;
        for (const AphidPop& ap : aphids) z += ap.total_aphids();

        S = 1 / (1 + z / K);

        empty = true;

        for (uint32 i = 0; i < aphids.size(); i++) {

            // Update population, including process error and dispersal:
            aphids[i].update_pop(this,
                                 emigrants.slice(i).col(this_j),
                                 immigrants.slice(i).col(this_j),
                                 eng);

            // Adjust for potential extinction or re-colonization:
            extinct_colonize(i);

        }

        age++;

        return;

    }
    // Same but minus stochasticity
    void update_pops(const arma::cube& emigrants,
                     const arma::cube& immigrants) {

        z = 0;
        for (const AphidPop& ap : aphids) z += ap.total_aphids();

        S = 1 / (1 + z / K);

        empty = true;

        for (uint32 i = 0; i < aphids.size(); i++) {
            aphids[i].update_pop(this,
                                 emigrants.slice(i).col(this_j),
                                 immigrants.slice(i).col(this_j));
            extinct_colonize(i);
        }

        age++;

        return;

    }


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
    std::normal_distribution<double> norm_distr = std::normal_distribution<double>(0,1);
    double mu_K_;       // mean of normal distribution of `K` for plants
    double sigma_K_;    // sd of normal distribution of `K` for plants
    bool K_error;       // whether to include stochasticity in `K`



public:

    std::vector<OnePatch> patches;
    arma::cube emigrants;
    arma::cube immigrants;


    AllPatches()
        : mu_K_(), sigma_K_(), patches(), emigrants(), immigrants() {};

    /*
     In `aphid_density_0` below, rows are aphid stages, columns are types (alate vs
     apterous), and slices are aphid lines.
     In `leslie_mat` below, slices are aphid lines.
     */
    AllPatches(const double& sigma,
               const double& rho,
               const double& demog_mult,
               const double& mu_K,
               const double& sigma_K,
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
        : mu_K_(mu_K), sigma_K_(sigma_K), K_error(true), patches(),
          emigrants(), immigrants() {

        if (sigma_K_ <= 0) {
            sigma_K_ = 0;
            K_error = false;
        }
        uint32 n_patches = std::min(aphid_density_0.size(), pred_rate.size());
        uint32 n_lines = aphid_name.size();
        uint32 n_stages = leslie_mat.front().n_rows;

        patches.reserve(n_patches);
        for (uint32 j = 0; j < n_patches; j++) {
            double K = mu_K_;
            if (K_error) K += sigma_K_ * norm_distr(eng);
            OnePatch ap(sigma, rho, demog_mult, K, aphid_name, leslie_mat,
                        aphid_density_0[j], alate_prop, disp_rate, disp_mort,
                        disp_start, pred_rate[j], n_patches, j, extinct_N);
            patches.push_back(ap);
        }

        emigrants.set_size(n_stages, n_patches, n_lines);
        immigrants.set_size(n_stages, n_patches, n_lines);

    }

    AllPatches(const AllPatches& other)
        : mu_K_(other.mu_K_),
          sigma_K_(other.sigma_K_),
          patches(other.patches),
          emigrants(other.emigrants),
          immigrants(other.immigrants) {};

    AllPatches& operator=(const AllPatches& other) {
        mu_K_ = other.mu_K_;
        sigma_K_ = other.sigma_K_;
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
    inline void clear_patches(const uint32& max_age,
                              pcg32& eng) {
        if (K_error) {
            for (OnePatch& p : patches) {
                if (p.age > max_age) {
                    double K_= mu_K_ + sigma_K_ * norm_distr(eng);
                    p.clear(K_);
                }
            }
        } else {
            for (OnePatch& p : patches) {
                if (p.age > max_age) p.clear();
            }
        }

        return;
    }
    inline void clear_patches(const double& max_N,
                              pcg32& eng) {
        if (K_error) {
            for (OnePatch& p : patches) {
                double N = p.total_aphids();
                if (N > max_N) {
                    double K_= mu_K_ + sigma_K_ * norm_distr(eng);
                    p.clear(K_);
                }
            }
        } else {
            for (OnePatch& p : patches) {
                double N = p.total_aphids();
                if (N > max_N) p.clear();
            }
        }
        return;
    }


};





#endif
