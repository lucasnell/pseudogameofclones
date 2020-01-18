
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





// Add process error
void AphidTypePop::process_error(const double& z,
                                 const double& sigma,
                                 const double& rho,
                                 const double& demog_mult,
                                 std::normal_distribution<double>& norm_distr,
                                 pcg32& eng) {

    if (demog_mult == 0 || sigma == 0) return;

    uint32 n_stages = X_t1.n_elem;

    arma::mat Se(n_stages, n_stages, arma::fill::zeros);

    Se = (sigma*sigma + demog_mult * std::min(0.5, 1 / std::abs(1 + z))) *
        (rho * arma::mat(n_stages,n_stages,arma::fill::ones) +
        (1-rho) * arma::mat(n_stages,n_stages,arma::fill::eye));

    // chol doesn't work with zeros on diagonal
    arma::uvec non_zero = arma::find(Se.diag() > 0);

    /*
     Cholesky decomposition of Se so output has correct variance-covariance matrix
       "a vector of independent normal random variables,
       when multiplied by the transpose of the Cholesky deposition of [Se] will
       have covariance matrix equal to [Se]."
     */
    arma::mat chol_decomp = arma::chol(Se(non_zero,non_zero)).t();

    // Random numbers from distribution N(0,1)
    arma::vec E(non_zero.n_elem);
    for (uint32 i = 0; i < E.n_elem; i++) E(i) = rnorm_distr(eng);

    // Making each element of E have correct variance-covariance matrix
    E = chol_decomp * E;

    // Plugging in errors into the X[t+1] vector
    for (uint32 i = 0; i < non_zero.n_elem; i++) {
        X_t1(non_zero(i)) *= std::exp(E(i));
    }


    /*
     Because we used normal distributions to approximate demographic and environmental
     stochasticity, it is possible for aphids and parasitoids to
     "spontaneously appear" when the estimate of e(t) is large. To disallow this
     possibility, the number of aphids and parasitized aphids in a given age class
     on day t was not allowed to exceed the number in the preceding age class on
     day t – 1.
    */
    for (uint32 i = 1; i < X_t1.n_elem; i++) {
        if (X_t1(i) > X_t(i-1)) X_t1(i) = X_t(i-1);
    }

    return;

}











/*
 Emigration and immigration of this line to all other patches.

 (`this_j` is the index for the patch this line's population is on)

 These do not necessarily match up due to mortality of dispersers.

 For both `emigrants` and `immigrants` matrices, rows are aphid stages and columns
 are patches.
 */

void AphidPop::calc_dispersal(const OnePatch* patch,
                              arma::mat& emigrants,
                              arma::mat& immigrants,
                              pcg32& eng) const {

    const uint32& this_j(patch->this_j);
    const uint32& n_patches(patch->n_patches);

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
void AphidPop::calc_dispersal(const OnePatch* patch,
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







void AphidPop::update_pop(const double& z,
                          const double& pred_rate,
                          const arma::vec& emigrants,
                          const arma::vec& immigrants,
                          pcg32& eng) {

    // Basic updates for each:
    apterous.X_t = apterous.X_t1;
    apterous.X_t1 = (pred_rate * S(z)) % (apterous.leslie_ * apterous.X_t);
    alates.X_t = alates.X_t1;
    alates.X_t1 = (pred_rate * S(z)) % (alates.leslie_ * alates.X_t);

    // Process error
    apterous.process_error(z, sigma_, rho_, demog_mult_, norm_distr, eng);
    alates.process_error(z, sigma_, rho_, demog_mult_, norm_distr, eng);

    // Sample for # offspring from apterous aphids that are alates:
    double new_alates = 0;
    if (apterous.alate_prop_ > 0 && apterous.X_t1.front() > 0) {
        double lambda_ = apterous.alate_prop_ * apterous.X_t1.front();
        pois_distr.param(std::poisson_distribution<uint32>::param_type(lambda_));
        new_alates = static_cast<double>(pois_distr(eng));
    }

    /*
     All alate offspring are assumed to be apterous,
     so the only way to get new alates is from apterous aphids.
     */
    apterous.X_t1.front() -= new_alates;
    apterous.X_t1.front() += alates.X_t1.front();
    alates.X_t1.front() = new_alates;

    // Finally add immigrants and subtract emigrants
    alates.X_t1 += immigrants;
    alates.X_t1 -= emigrants;

    return;
}

// Same as above, but no randomness in alate production:
void AphidPop::update_pop(const double& z,
                          const double& pred_rate,
                          const arma::vec& emigrants,
                          const arma::vec& immigrants) {

    // Basic updates for each:
    apterous.X_t = apterous.X_t1;
    apterous.X_t1 = (pred_rate * S(z)) % (apterous.leslie_ * apterous.X_t);
    alates.X_t = alates.X_t1;
    alates.X_t1 = (pred_rate * S(z)) % (alates.leslie_ * alates.X_t);
    // # offspring from apterous aphids that are alates:
    double new_alates = apterous.alate_prop_ * apterous.X_t1.front();

    /*
     All alate offspring are assumed to be apterous,
     so the only way to get new alates is from apterous aphids.
     */
    apterous.X_t1.front() -= new_alates;
    apterous.X_t1.front() += alates.X_t1.front();
    alates.X_t1.front() = new_alates;

    // Finally add immigrants and subtract emigrants
    alates.X_t1 += immigrants;
    alates.X_t1 -= emigrants;

    return;
}

