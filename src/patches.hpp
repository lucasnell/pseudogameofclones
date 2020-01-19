# ifndef __CLONEWARS_PATCHES_H
# define __CLONEWARS_PATCHES_H


#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "clonewars_types.hpp"  // integer types
#include "aphids.hpp"           // aphid classes
#include "math.hpp"             // leslie_matrix and leslie_sad
#include "pcg.hpp"              // runif_ fxns



using namespace Rcpp;







/*
One patch of continuous habitat.
*/
class OnePatch {

    // aphid intraspecific density dependence
    inline double S_fun() const {
        return 1 / (1 + z / K);
    }

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
        : aphids(), extinct(), empty(true), pred_rate(0), K(0),
          n_patches(1), this_j(0) {};

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
        return this;
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
            for (AphidPop& ap : aphids) {
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

        S = S_fun(z);

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

        S = S_fun(z);

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

public:

    std::vector<OnePatch> patches;
    arma::cube emigrants;
    arma::cube immigrants;


    AllPatches() : patches(), emigrants(), immigrants() {};

    /*
     In `aphid_density_0` below, rows are aphid stages, columns are types (alate vs
     apterous), and slices are aphid lines.
     In `leslie_mat` below, slices are aphid lines.
     */
    AllPatches(const double& sigma,
               const double& rho,
               const double& demog_mult,
               const std::vector<double>& K,
               const std::vector<std::string>& aphid_name,
               const std::vector<arma::cube>& leslie_mat,
               const std::vector<arma::cube>& aphid_density_0,
               const std::vector<double>& alate_prop,
               const std::vector<double>& disp_rate,
               const std::vector<double>& disp_mort,
               const std::vector<uint32>& disp_start,
               const std::vector<double>& pred_rate,
               const double& extinct_N) {

        uint32 n_patches = std::min(aphid_density_0.size(), pred_rate.size());
        uint32 n_lines = aphid_name.size();
        uint32 n_stages = leslie_mat.front().n_rows;

        patches.reserve(n_patches);
        for (uint32 j = 0; j < n_patches; j++) {
            OnePatch ap(sigma, rho, demog_mult, K[j], aphid_name, leslie_mat,
                        aphid_density_0[j], alate_prop, disp_rate, disp_mort,
                        disp_start, pred_rate[j], n_patches, j, extinct_N);
            patches.push_back(ap);
        }

        emigrants.set_size(n_stages, n_patches, n_lines);
        immigrants.set_size(n_stages, n_patches, n_lines);

    }

    AllPatches(const AllPatches& other)
        : patches(other.patches),
          emigrants(other.emigrants),
          immigrants(other.immigrants) {};

    AllPatches& operator=(const AllPatches& other) {
        patches = other.patches;
        emigrants = other.emigrants;
        immigrants = other.immigrants;
        return this;
    };



    void iterate(pcg32& eng) {
        // Dispersal from previous generation
        emigrants.fill(0);
        immigrants.fill(0);
        for (const OnePatch& p : patches) {
            p.calc_dispersal(emigrants, immigrants, eng);
        }
        /*
         Once that's updated inside `emigrants` and `immigrants`, we can update the
         populations using those dispersal numbers.
         */
        for (OnePatch& p : patches) {
            p.update_pops(emigrants, immigrants, eng);
        }
        return;
    }
    void iterate() {
        emigrants.fill(0);
        immigrants.fill(0);
        for (const OnePatch& p : patches) {
            p.calc_dispersal(emigrants, immigrants);
        }
        for (OnePatch& p : patches) {
            p.update_pops(emigrants, immigrants);
        }
        return;
    }


    // Clear patches by either a maximum age or total abundance
    void clear_patches(const uint32& max_age) {
        for (OnePatch& p : patches) {
            if (p.age > max_age) p.clear();
        }
        return;
    }
    void clear_patches(const double& max_N) {
        for (OnePatch& p : patches) {
            double N = p.total_aphids();
            if (N > max_N) p.clear();
        }
        return;
    }
    // Same, but reset `K` when they're cleared
    void clear_patches(const uint32& max_age,
                       std::normal_distribution<double>& distr,
                       pcg32& eng) {
        for (OnePatch& p : patches) {
            if (p.age > max_age) {
                double K_= distr(eng);
                p.clear(K_);
            }
        }
        return;
    }
    void clear_patches(const double& max_N,
                       std::normal_distribution<double>& distr,
                       pcg32& eng) {
        for (OnePatch& p : patches) {
            double N = p.total_aphids();
            if (N > max_N) {
                double K_= distr(eng);
                p.clear(K_);
            }
        }
        return;
    }

};





#endif
