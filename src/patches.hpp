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
    inline double S() const {
        return 1 / (1 + z / K);
    }

    /*
     Adjust for potential extinction or re-colonization:
     */
    void extinct_colonize(const uint32& i, const double& extinct_N) {
        double N = aphids[i].total_aphids();
        if (N > extinct_N && empty) empty = false;
        // Newly extinct:
        if (N < extinct_N && !extinct[i]) {
            aphids[i].clear();
            extinct[i] = true;
        }
        // Newly colonized:
        if (N > extinct_N && extinct[i]) {
            extinct[i] = false;
        }
        return;
    }

public:

    std::vector<AphidPop> aphids;   // aphids in this patch
    std::vector<bool> extinct;      // keeping track of line extinctions
    bool empty;                     // boolean for whether no aphids are on this patch
    double pred_rate;               // predation on aphids
    double z;                       // sum of all aphids at time t
    double K;                       // aphid carrying capacity
    uint32 n_patches;               // total # patches
    uint32 this_j;                  // index for this patch
    uint32 age = 0;                 // age of this patch


    OnePatch()
        : aphids(), extinct(), empty(true), pred_rate(0), z(0),
          n_patches(1), this_j(0) {};
    // Starting all aphids with "stable age distribution" with a given total density
    OnePatch(const std::vector<std::string>& aphid_name_,
             const double& sigma,
             const double& rho,
             const double& demog_mult,
             const double& K_,
             const std::vector<std::vector<arma::uvec>>& instar_days,
             const std::vector<std::vector<double>>& surv_juv,
             const std::vector<std::vector<arma::vec>>& surv_adult,
             const std::vector<std::vector<arma::vec>>& repro,
             const std::vector<std::vector<double>>& aphid_density_0,
             const std::vector<double>& alate_prop,
             const std::vector<double>& disp_rate,
             const std::vector<double>& disp_mort,
             const std::vector<uint32>& disp_start,
             const double& pred_rate_,
             const uint32& n_patches_,
             const uint32& this_j_)
        : aphids(),
          extinct(aphid_name_.size(), false),
          empty(true),
          pred_rate(pred_rate_),
          z(0),
          K(K_),
          n_patches(n_patches_),
          this_j(this_j_) {

        aphids.reserve(n_lines);

        for (uint32 i = 0; i < n_lines; i++) {
            AphidPop ap(aphid_name_[i], sigma, rho, demog_mult,
                        instar_days[i], surv_juv[i], surv_adult[i], repro[i],
                        aphid_density_0[i], alate_prop[i], disp_rate[i], disp_mort[i],
                        disp_start[i]);
            aphids.push_back(ap);
            const double& zi(aphid_density_0[i]);
            z += zi;
            if (zi <= 0) extinct[i] = true;
        }

    };

    // Starting each line with starting densities of each stage directly given:
    OnePatch(const std::vector<std::string>& aphid_name_,
             const double& sigma,
             const double& rho,
             const double& demog_mult,
             const double& K_,
             const std::vector<std::vector<arma::uvec>>& instar_days,
             const std::vector<std::vector<double>>& surv_juv,
             const std::vector<std::vector<arma::vec>>& surv_adult,
             const std::vector<std::vector<arma::vec>>& repro,
             const std::vector<std::vector<arma::vec>>& aphid_density_0,
             const std::vector<double>& alate_prop,
             const std::vector<double>& disp_rate,
             const std::vector<double>& disp_mort,
             const std::vector<uint32>& disp_start,
             const double& pred_rate_,
             const uint32& n_patches_,
             const uint32& this_j_)
        : aphids(),
          extinct(aphid_name_.size(), false),
          empty(true),
          pred_rate(pred_rate_),
          z(0),
          K(K_),
          n_patches(n_patches_),
          this_j(this_j_) {

        uint32 n_lines = aphid_name_.size();

        aphids.reserve(n_lines);

        for (uint32 i = 0; i < n_lines; i++) {
            AphidPop ap(aphid_name_[i], sigma, rho, demog_mult,
                        instar_days[i], surv_juv[i], surv_adult[i], repro[i],
                        aphid_density_0[i], alate_prop[i], disp_rate[i], disp_mort[i],
                        disp_start[i]);
            aphids.push_back(ap);
            double zi = arma::accu(aphid_density_0[i]);
            z += zi;
            if (zi <= 0) extinct[i] = true;
        }

    };


    OnePatch(const OnePatch& other)
        : aphids(other.aphids), extinct(other.extinct), empty(other.empty),
          pred_rate(other.pred_rate), z(other.z), K(other.K), n_patches(other.n_patches),
          this_j(other.this_j), age(other.age) {};

    OnePatch& operator=(const OnePatch& other) {
        aphids = other.aphids;
        extinct = other.extinct;
        empty = other.empty;
        pred_rate = other.pred_rate;
        z = other.z;
        K = other.K;
        n_patches = other.n_patches;
        this_j = other.this_j;
        age = other.age;
        return this;
    }


    /*
     Clear to no aphids
     */
    void clear() {
        for (uint32 i = 0; i < aphids.size(); i++) {
            aphids[i].clear();
            extinct[i] = true;
        }
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
                     const double& extinct_N,
                     pcg32& eng) {

        z = 0;
        for (const AphidPop& ap : aphids) z += ap.total_aphids();

        double S_ = S(z);

        empty = true;

        for (uint32 i = 0; i < aphids.size(); i++) {

            // Update population, including process error and dispersal:
            aphids[i].update_pop(z, S_, pred_rate,
                                 emigrants.slice(i).col(this_j),
                                 immigrants.slice(i).col(this_j),
                                 eng);

            // Adjust for potential extinction or re-colonization:
            extinct_colonize(i, extinct_N);

        }

        age++;

        return;

    }
    // Same but minus stochasticity
    void update_pops(const arma::cube& emigrants,
                     const arma::cube& immigrants,
                     const double& extinct_N) {

        z = 0;
        for (const AphidPop& ap : aphids) z += ap.total_aphids();

        double S_ = S(z);

        empty = true;

        for (uint32 i = 0; i < aphids.size(); i++) {
            aphids[i].update_pop(z, S_, pred_rate,
                                 emigrants.slice(i).col(this_j),
                                 immigrants.slice(i).col(this_j));
            extinct_colonize(i, extinct_N);
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
    // Starting all aphids with "stable age distribution" with a given total density
    AllPatches(const std::vector<std::string>& aphid_name_,
               const double& sigma,
               const double& rho,
               const double& demog_mult,
               const std::vector<double>& K,
               const std::vector<std::vector<arma::uvec>>& instar_days,
               const std::vector<std::vector<double>>& surv_juv,
               const std::vector<std::vector<arma::vec>>& surv_adult,
               const std::vector<std::vector<arma::vec>>& repro,
               const std::vector<std::vector<std::vector<double>>>& aphid_density_0,
               const std::vector<double>& alate_prop,
               const std::vector<double>& disp_rate,
               const std::vector<double>& disp_mort,
               const std::vector<uint32>& disp_start,
               const std::vector<double>& pred_rate_) {

        uint32 n_patches = std::min(aphid_density_0.size(), pred_rate_.size());
        uint32 n_lines = aphid_name_.size();
        uint32 n_stages = arma::accu(instar_days.front().front());

        patches.reserve(n_patches);
        for (uint32 j = 0; j < n_patches; j++) {
            OnePatch ap(aphid_name_, sigma, rho, demog_mult, K[j], instar_days,
                        surv_juv, surv_adult, repro, aphid_density_0[j],
                        alate_prop, disp_rate, disp_mort, disp_start,
                        pred_rate_[j], n_patches, j);
            patches.push_back(ap);
        }

        emigrants.set_size(n_stages, n_patches, n_lines);
        immigrants.set_size(n_stages, n_patches, n_lines);

    }
    // Starting each line with starting densities of each stage directly given:
    AllPatches(const std::vector<std::string>& aphid_name_,
               const double& sigma,
               const double& rho,
               const double& demog_mult,
               const std::vector<double>& K,
               const std::vector<std::vector<arma::uvec>>& instar_days,
               const std::vector<std::vector<double>>& surv_juv,
               const std::vector<std::vector<arma::vec>>& surv_adult,
               const std::vector<std::vector<arma::vec>>& repro,
               const std::vector<std::vector<std::vector<arma::vec>>>& aphid_density_0,
               const std::vector<double>& alate_prop,
               const std::vector<double>& disp_rate,
               const std::vector<double>& disp_mort,
               const std::vector<uint32>& disp_start,
               const std::vector<double>& pred_rate_) {

        uint32 n_patches = std::min(aphid_density_0.size(), pred_rate_.size());

        patches.reserve(n_patches);
        for (uint32 j = 0; j < n_patches; j++) {
            OnePatch ap(aphid_name_, sigma, rho, demog_mult, K[j], instar_days,
                        surv_juv, surv_adult, repro, aphid_density_0[j],
                        alate_prop, disp_rate, disp_mort, disp_start,
                        pred_rate_[j], n_patches, j);
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



    void iterate(const double& extinct_N, pcg32& eng) {
        emigrants.fill(0);
        immigrants.fill(0);
        for (const OnePatch& p : patches) {
            p.calc_dispersal(emigrants, immigrants, eng);
        }
        for (OnePatch& p : patches) {
            p.update_pops(emigrants, immigrants, extinct_N, eng);
        }
        return;
    }
    void iterate(const double& extinct_N) {
        emigrants.fill(0);
        immigrants.fill(0);
        for (const OnePatch& p : patches) {
            p.calc_dispersal(emigrants, immigrants);
        }
        for (OnePatch& p : patches) {
            p.update_pops(emigrants, immigrants, extinct_N);
        }
        return;
    }


};




#endif
