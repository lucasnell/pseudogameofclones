
#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "pseudogameofclones_types.hpp"  // integer types
#include "aphids.hpp"           // aphid classes
#include "wasps.hpp"            // wasp classes
#include "patches.hpp"          // field and plant classes
#include "pcg.hpp"              // runif_01 fxn
#include "math.hpp"             // inv_logit__ fxn



using namespace Rcpp;





// Add process error
void AphidTypePop::process_error(const arma::vec& Xt,
                                 const double& sigma_x,
                                 const double& rho,
                                 const bool& demog_error,
                                 const double& aphids_sum,
                                 std::normal_distribution<double>& norm_distr,
                                 pcg32& eng) {

    if (!demog_error && sigma_x == 0) return;

    uint32 n_stages = X.n_elem;

    // Variance for all process error:
    double sigma2 = sigma_x * sigma_x;
    if (demog_error) sigma2 += std::min(0.5, 1 / (1 + aphids_sum));

    arma::mat Se = sigma2 *
        (rho * arma::mat(n_stages, n_stages, arma::fill::ones) +
        (1 - rho) * arma::mat(n_stages, n_stages, arma::fill::eye));

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
    for (uint32 i = 1; i < n_stages; i++) {
        if (X(i) > Xt(i-1)) X(i) = Xt(i-1);
    }

    return;

}



// logit(Pr(alates)) ~ b0 + b1 * z, where `z` is # aphids (all lines)
double ApterousPop::alate_prop(const OnePlant* plant) const {
    const double lap = alate_b0_ + alate_b1_ * plant->z;
    double ap;
    inv_logit__(lap, ap);
    return ap;
}









/*
 Update living aphids (both parasitized and un-parasitized), then
 return the number of newly mummified aphids,
 to be added to the wasps.
 */

double AphidPop::update(const OnePlant* plant,
                        const WaspPop* wasps,
                        pcg32& eng) {

    double nm = 0; // newly mummified

    if (total_aphids() > 0) {

        const double& S(plant->S);
        const double& S_y(plant->S_y);

        arma::vec A_surv;
        wasps->A_mats(A_surv, attack_surv);
        arma::vec A_mumm = 1 - A_surv;

        // making adult alates not able to be parasitized:
        arma::vec A_surv_ala = A_surv;
        arma::vec A_mumm_ala = A_mumm;
        for (uint32 i = alates.field_disp_start_; i < alates.leslie_.n_cols; i++) {
            A_surv_ala[i] = 1;
            A_mumm_ala[i] = 0;
        }

        double pred_surv = 1 - plant->pred_rate;

        // Starting abundances (used in `process_error`):
        arma::vec apterous_Xt = apterous.X;
        arma::vec alates_Xt = alates.X;
        arma::vec paras_Xt = paras.X;

        // Basic updates for non-parasitized aphids:
        arma::mat LX_apt = apterous.leslie_ * apterous.X;
        arma::mat LX_ala = alates.leslie_ * alates.X;
        apterous.X = (pred_surv * S * A_surv) % LX_apt;
        alates.X = (pred_surv * S * A_surv_ala) % LX_ala;

        double np = 0; // newly parasitized
        np += pred_surv * S_y * arma::as_scalar(A_mumm.t() * LX_apt);
        np += pred_surv * S_y * arma::as_scalar(A_mumm_ala.t() * LX_ala);

        nm += pred_surv * paras.X.back();  // newly mummified

        // alive but parasitized
        if (paras.X.n_elem > 1) {
            for (uint32 i = paras.X.n_elem - 1; i > 0; i--) {
                paras.X(i) = pred_surv * paras.s(i) * S_y * paras.X(i-1);
            }
        }
        paras.X.front() = np;


        // Process error
        if (demog_error || sigma_x > 0) {
            process_error(apterous_Xt, alates_Xt, paras_Xt, eng);
        }

        // # offspring from apterous aphids that are alates:
        double alate_prop = apterous.alate_prop(plant);
        double new_alates = apterous.X.front() * alate_prop;

        /*
         All offspring from alates are assumed to be apterous,
         so the only way to get new alates is from apterous aphids.
         */
        apterous.X.front() -= new_alates;
        apterous.X.front() += alates.X.front(); // <-- we assume alates make apterous
        alates.X.front() = new_alates;

    }

    return nm;
}



