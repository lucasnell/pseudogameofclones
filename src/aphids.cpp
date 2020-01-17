
#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "clonewars_types.hpp"  // integer types
#include "aphids.hpp"           // aphid classes
#include "patches.hpp"          // patch classes
#include "math.hpp"             // leslie_matrix and leslie_sad
#include "pcg.hpp"              // runif_ fxns



using namespace Rcpp;


/*
 Emigration and immigration of this line to all other patches.

 (`this_j` is the index for the patch this line's population is on)

 These do not necessarily match up due to mortality of dispersers.

 For both `emigrants` and `immigrants` matrices, rows are aphid stages and columns
 are patches.
 */

void AphidPop::dispersal(const OnePatch* patch,
                         arma::mat& emigrants,
                         arma::mat& immigrants,
                         pcg32& eng) const {

    const uint32& this_j(patch->this_j);
    const uint32& n_patches(patch->n_patches);

    // Do this in higher-level function:
    // if (emigrants.n_rows != alates.n_aphid_stages() || emigrants.n_cols != n_patches) {
    //     emigrants.set_size(alates.n_aphid_stages(), n_patches);
    // }
    // if (immigrants.n_rows != alates.n_aphid_stages() || immigrants.n_cols != n_patches) {
    //     immigrants.set_size(alates.n_aphid_stages(), n_patches);
    // }
    //
    // emigrants.fill(0);
    // immigrants.fill(0);


    if (n_patches == 1 || alates.disp_rate() <= 0) return;

    // Abundance for alates. (Only adult alates can disperse.)
    const arma::vec& X_disp(alates.X_t);

    arma::rowvec n_leaving(n_patches);
    arma::rowvec n_leaving_alive(n_patches);

    // Sample dispersal for each dispersing stage:
    for (uint32 i = alates.disp_start(); i < X_disp.n_elem; i++) {

        if (X_disp(i) < 1) continue;

        /*
         Calculate emigration, or the # aphids that leave the patch:
         */
        // Calculate lambda. (We're splitting dispersers over # of other patches.)
        double lambda_ = alates.disp_rate() * X_disp(i) /
            static_cast<double>(n_patches - 1);
        pois_distr.param(std::poisson_distribution<uint32>::param_type(lambda_));
        for (uint32 j = 0; j < n_patches; j++) {
            if (j == this_j) {
                n_leaving(j) = 0;
            } else n_leaving(j) = pois_distr(eng);
        }
        // Making absolutely sure that dispersal never exceeds the number possible:
        double total_emigrants = arma::accu(n_leaving);
        if (total_emigrants > X_disp(i)) {
            double extras = total_emigrants - X_disp(i);
            std::vector<uint32> extra_inds;
            extra_inds.reserve(n_patches);
            for (uint32 j = 0; j < n_patches; j++) {
                if (n_leaving(j) > 0) extra_inds.push_back(j);
            }
            while (extras > 0) {
                uint32 rnd = runif_01(eng) * extra_inds.size();
                n_leaving(extra_inds[rnd])--;
                if (n_leaving(extra_inds[rnd]) == 0) {
                    extra_inds.erase(extra_inds.begin() + rnd);
                }
                extras--;
                total_emigrants--;
            }
        }

        emigrants(i, this_j) = total_emigrants;

        /*
         Calculate immigration, or the number leaving that stay alive to get to
         another patch.
         For the focal line and patch, stage i, to patch j, the upper limit to this is
         `emigrants(i,j)`.
         */
        if (alates.disp_mort() <= 0) {
            // If mortality is zero, then all emigrants survive:
            immigrants.row(i) = n_leaving;
        } else if (alates.disp_mort() < 1) {
            // If mortality is >0 and <1, then we have to sample for the # that survive:
            for (uint32 j = 0; j < n_patches; j++) {
                if (j == this_j || n_leaving(j) == 0) continue;
                // Sample number leaving that stay alive to immigrate:
                bino_distr.param(std::binomial_distribution<uint32>::param_type(
                        n_leaving(j), 1 - alates.disp_mort()));
                immigrants(i,j) += static_cast<double>(bino_distr(eng));
            }
        }

    }


    return;
}




// Same as above, but with no stochasticity
void AphidPop::dispersal(const OnePatch* patch,
                         arma::mat& emigrants,
                         arma::mat& immigrants) const {

    const uint32& this_j(patch->this_j);
    const uint32& n_patches(patch->n_patches);

    if (n_patches == 1 || alates.disp_rate() <= 0) return;

    // Abundance for alates. (Only adult alates can disperse.)
    const arma::vec& X_disp(alates.X_t);

    arma::rowvec n_leaving(n_patches);
    arma::rowvec n_leaving_alive(n_patches);

    // Dispersal for each dispersing stage:
    for (uint32 i = alates.disp_start(); i < X_disp.n_elem; i++) {

        if (X_disp(i) == 0) continue;

        double lambda_ = alates.disp_rate() * X_disp(i) /
            static_cast<double>(n_patches - 1);
        n_leaving.fill(lambda_);
        n_leaving(this_j, this_j) = 0;
        emigrants(i, this_j) = arma::accu(n_leaving);

        if (alates.disp_mort() <= 0) {
            immigrants.row(i) = n_leaving;
        } else if (alates.disp_mort() < 1) {
            immigrants.row(i) = n_leaving * (1 - alates.disp_mort());
        }

    }


    return;
}


