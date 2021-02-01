
#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "clonewars_types.hpp"  // integer types
#include "aphids.hpp"           // aphid classes
#include "wasps.hpp"            // wasp classes
#include "patches.hpp"          // patch classes
#include "pcg.hpp"              // runif_ fxns



using namespace Rcpp;





// Add process error
void AphidTypePop::process_error(const double& z,
                                 const double& sigma_x,
                                 const double& rho,
                                 const double& demog_mult,
                                 std::normal_distribution<double>& norm_distr,
                                 pcg32& eng) {

    if (demog_mult == 0 || sigma_x == 0) return;

    // If `X` is at time t+1, this is at time t. It's used at the very bottom of this fxn.
    arma::vec X_t = X;

    uint32 n_stages = X.n_elem;

    arma::mat Se(n_stages, n_stages, arma::fill::zeros);

    Se = (sigma_x*sigma_x + demog_mult * std::min(0.5, 1 / std::abs(1 + z))) *
        (rho * arma::mat(n_stages,n_stages,arma::fill::ones) +
        (1-rho) * arma::mat(n_stages,n_stages,arma::fill::eye));

    // chol doesn't work with zeros on diagonal
    arma::uvec non_zero = arma::find(Se.diag() > 0);

    /*
     Cholesky decomposition of Se so output has correct variance-covariance
     matrix:
       "a vector of independent normal random variables,
       when multiplied by the transpose of the Cholesky deposition of [Se] will
       have covariance matrix equal to [Se]."
     */
    arma::mat chol_decomp = arma::chol(Se(non_zero,non_zero)).t();

    // Random numbers from distribution N(0,1)
    arma::vec E(non_zero.n_elem);
    for (uint32 i = 0; i < E.n_elem; i++) E(i) = norm_distr(eng);

    // Making each element of E have correct variance-covariance matrix
    E = chol_decomp * E;

    // Plugging in errors into the X[t+1] vector
    for (uint32 i = 0; i < non_zero.n_elem; i++) {
        X(non_zero(i)) *= std::exp(E(i));
    }


    /*
     Because we used normal distributions to approximate demographic and environmental
     stochasticity, it is possible for aphids and parasitoids to
     "spontaneously appear" when the estimate of e(t) is large. To disallow this
     possibility, the number of aphids and parasitized aphids in a given age class
     on day t was not allowed to exceed the number in the preceding age class on
     day t – 1.
    */
    for (uint32 i = 1; i < X.n_elem; i++) {
        if (X(i) > X_t(i-1)) X(i) = X_t(i-1);
    }

    return;

}







// logit(Pr(alates)) ~ b0 + b1 * z, where `z` is # aphids (all lines)
double ApterousPop::alate_prop(const OnePatch* patch) const {

    const double lap = alate_b0_ + alate_b1_ * patch->z;
    double ap;
    inv_logit__(lap, ap);
    return ap;

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

    if (arma::accu(alates.X) == 0 || n_patches == 1 || alates.disp_rate() <= 0) return;

    // Abundance for alates. (Only adult alates can disperse.)
    const arma::vec& X_disp(alates.X);

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
            immigrants.row(i) += n_leaving;
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

    if (arma::accu(alates.X) == 0 || n_patches == 1 || alates.disp_rate() <= 0) return;


    // Abundance for alates. (Only adult alates can disperse.)
    const arma::vec& X_disp(alates.X);

    arma::rowvec n_leaving(n_patches);
    arma::rowvec n_leaving_alive(n_patches);

    // Dispersal for each dispersing stage:
    for (uint32 i = alates.disp_start(); i < X_disp.n_elem; i++) {

        if (X_disp(i) == 0) continue;

        double lambda_ = alates.disp_rate() * X_disp(i) /
            static_cast<double>(n_patches - 1);
        n_leaving.fill(lambda_);
        n_leaving(this_j) = 0;
        emigrants(i, this_j) = arma::accu(n_leaving);

        if (alates.disp_mort() <= 0) {
            immigrants.row(i) += n_leaving;
        } else if (alates.disp_mort() < 1) {
            immigrants.row(i) += n_leaving * (1 - alates.disp_mort());
        }

    }


    return;
}




/*
 Update living aphids (both parasitized and un-parasitized), then
 return the number of newly mummified aphids,
 to be added to the wasps.
 */

double AphidPop::update(const OnePatch* patch,
                        const WaspPop* wasps,
                        const arma::vec& emigrants,
                        const arma::vec& immigrants,
                        pcg32& eng) {


    // First subtract emigrants and add immigrants:
    alates.X -= emigrants;
    alates.X += immigrants;

    double nm = 0; // newly mummified

    if (arma::accu(alates.X) > 0 || arma::accu(apterous.X) > 0) {

        const double& z(patch->z);
        const double& S(patch->S);
        const double& S_y(patch->S_y);
        arma::vec A = wasps->A(attack_surv);
        double pred_surv = 1 - patch->pred_rate;


        // Basic updates for non-parasitized aphids:
        arma::mat LX_apt = apterous.leslie_ * apterous.X;
        arma::mat LX_ala = alates.leslie_ * alates.X;
        apterous.X = pred_surv * S * A % LX_apt;
        alates.X = pred_surv * S * A % LX_ala;

        double np = 0; // newly parasitized
        np += pred_surv * S_y * arma::as_scalar((1 - A).t() * LX_apt);
        np += pred_surv * S_y * arma::as_scalar((1 - A).t() * LX_ala);

        nm += pred_surv * paras.X.back();  // newly mummified

        // alive but parasitized
        for (uint32 i = 1; i < paras.X.n_elem; i++) {
            paras.X(i) = pred_surv * paras.s(i) * S_y * paras.X(i-1);
        }
        paras.X.front() = np;


        // Process error
        apterous.process_error(z, sigma_x, rho, demog_mult, norm_distr, eng);
        alates.process_error(z, sigma_x, rho, demog_mult, norm_distr, eng);
        paras.process_error(z, sigma_x, rho, demog_mult, norm_distr, eng);


        // Sample for # offspring from apterous aphids that are alates:
        double new_alates = 0;
        double alate_prop = apterous.alate_prop(patch);
        if (alate_prop > 0 && apterous.X.front() > 0) {
            double lambda_ = alate_prop * apterous.X.front();
            pois_distr.param(std::poisson_distribution<uint32>::param_type(lambda_));
            new_alates = static_cast<double>(pois_distr(eng));
            if (new_alates > apterous.X.front()) new_alates = apterous.X.front();
        }

        /*
         All alate offspring are assumed to be apterous,
         so the only way to get new alates is from apterous aphids.
         */
        apterous.X.front() -= new_alates;
        apterous.X.front() += alates.X.front(); // <-- we assume alates make apterous
        alates.X.front() = new_alates;

    }

    return nm;
}

// Same as above, but no randomness in alate production:
double AphidPop::update(const OnePatch* patch,
                        const WaspPop* wasps,
                        const arma::vec& emigrants,
                        const arma::vec& immigrants) {

    // First subtract emigrants and add immigrants:
    alates.X -= emigrants;
    alates.X += immigrants;

    double nm = 0; // newly mummified

    if (arma::accu(alates.X) > 0 || arma::accu(apterous.X) > 0) {

        const double& S(patch->S);
        const double& S_y(patch->S_y);
        arma::vec A = wasps->A(attack_surv);
        double pred_surv = 1 - patch->pred_rate;

        // Basic updates for unparasitized aphids:
        arma::mat LX_apt = apterous.leslie_ * apterous.X;
        arma::mat LX_ala = alates.leslie_ * alates.X;
        apterous.X = (pred_surv * S * A) % LX_apt;
        alates.X = (pred_surv * S * A) % LX_ala;

        double np = 0; // newly parasitized
        np += pred_surv * S_y * arma::as_scalar((1 - A).t() * LX_apt);
        np += pred_surv * S_y * arma::as_scalar((1 - A).t() * LX_ala);

        nm += pred_surv * paras.X.back();  // newly mummified

        // alive but parasitized
        for (uint32 i = 1; i < paras.X.n_elem; i++) {
            paras.X(i) = pred_surv * paras.s(i) * S_y * paras.X(i-1);
        }
        paras.X.front() = np;

        // # offspring from apterous aphids that are alates:
        double new_alates = apterous.alate_prop(patch);
        new_alates *= apterous.X.front();


        /*
         All alate offspring are assumed to be apterous,
         so the only way to get new alates is from apterous aphids.
         */
        apterous.X.front() -= new_alates;
        apterous.X.front() += alates.X.front(); // <-- we assume alates make apterous
        alates.X.front() = new_alates;

    }


    return nm;
}


