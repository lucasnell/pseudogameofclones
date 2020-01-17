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

public:

    std::vector<AphidPop> aphids;   // aphids in this patch
    std::vector<bool> extinct;      // keeping track of line extinctions
    bool empty;                     // boolean for whether no aphids are on this patch
    double pred_rate;               // predation on aphids
    double z;                       // Sum of all aphids at time t
    uint32 n_patches;               // Total # patches
    uint32 this_j;                  // Index for this patch


    OnePatch() : aphids(), extinct(), empty(true), pred_rate(0), z(0) {};
    // Starting all aphids with "stable age distribution" with a given total density
    OnePatch(const std::vector<std::string>& aphid_name_,
             const std::vector<double>& sigma,
             const std::vector<double>& rho,
             const std::vector<double>& demog_mult,
             const std::vector<std::vector<double>>& K,
             const std::vector<std::vector<arma::uvec>>& instar_days,
             const std::vector<std::vector<double>>& surv_juv,
             const std::vector<std::vector<arma::vec>>& surv_adult,
             const std::vector<std::vector<arma::vec>>& repro,
             const std::vector<std::vector<double>>& aphid_density_0,
             const std::vector<double>& alate_rate,
             const std::vector<double>& disp_rate,
             const std::vector<double>& disp_mort,
             const std::vector<uint32>& disp_start,
             const double& pred_rate_)
        : aphids(),
          extinct(aphid_name_.size(), false),
          empty(true),
          pred_rate(pred_rate_),
          z(0) {

        aphids.reserve(n_lines);

        for (uint32 i = 0; i < n_lines; i++) {
            AphidPop ap(aphid_name_[i], sigma[i], rho[i], demog_mult[i], K[i],
                        instar_days[i], surv_juv[i], surv_adult[i], repro[i],
                        aphid_density_0[i], alate_rate[i], disp_rate[i], disp_mort[i],
                        disp_start[i]);
            aphids.push_back(ap);
            const double& zi(aphid_density_0[i]);
            z += zi;
            if (zi <= 0) extinct[i] = true;
        }

    };

    // Starting each line with starting densities of each stage directly given:
    OnePatch(const std::vector<std::string>& aphid_name_,
             const std::vector<double>& sigma,
             const std::vector<double>& rho,
             const std::vector<double>& demog_mult,
             const std::vector<std::vector<double>>& K,
             const std::vector<std::vector<arma::uvec>>& instar_days,
             const std::vector<std::vector<double>>& surv_juv,
             const std::vector<std::vector<arma::vec>>& surv_adult,
             const std::vector<std::vector<arma::vec>>& repro,
             const std::vector<std::vector<arma::vec>>& aphid_density_0,
             const std::vector<double>& alate_rate,
             const std::vector<double>& disp_rate,
             const std::vector<double>& disp_mort,
             const std::vector<uint32>& disp_start,
             const double& pred_rate_)
        : aphids(),
          extinct(aphid_name_.size(), false),
          empty(true),
          pred_rate(pred_rate_),
          z(0) {

        uint32 n_lines = aphid_name_.size();

        aphids.reserve(n_lines);

        for (uint32 i = 0; i < n_lines; i++) {
            AphidPop ap(aphid_name_[i], sigma[i], rho[i], demog_mult[i], K[i],
                        instar_days[i], surv_juv[i], surv_adult[i], repro[i],
                        aphid_density_0[i], alate_rate[i], disp_rate[i], disp_mort[i],
                        disp_start[i]);
            aphids.push_back(ap);
            double zi = arma::accu(aphid_density_0[i]);
            z += zi;
            if (zi <= 0) extinct[i] = true;
        }

    };

    /*
     Clear to no aphids
     */
    void clear() {
        for (uint32 i = 0; i < aphids.size(); i++) {
            aphids[i].clear();
            extinct[i] = true;
        }
        empty = true;
        return;
    }

    /*
     Emigration of one line from this patch to all other patches:
    */
    arma::vec emigration(const uint32& j,
                         const uint32& n_patches,
                         pcg32& eng);

    /*
     Same thing as above, but overloaded for not including dispersal stochasticity
    */
    arma::vec emigration(const uint32& j,
                         const uint32& n_patches);


    /*
     Iterate one time step
     */
    void iterate(const arma::mat& emigrants,
                 const arma::mat& immigrants,
                 const double& extinct_N,
                 pcg32& eng);

};




#endif
