
#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "gameofclones_types.hpp"  // integer types
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
                                 std::normal_distribution<double>& norm_distr,
                                 pcg32& eng) {

    if (!demog_error && sigma_x == 0) return;

    uint32 n_stages = X.n_elem;

    // Standard deviations by stage (all the same if no demographic error):
    arma::vec stdevs(n_stages, arma::fill::value(sigma_x));
    if (demog_error) {
        //' To add demographic error, convert to variances, combine,
        //' then convert back to stdev:
        stdevs *= sigma_x;
        for (uint32 i = 0; i < n_stages; i++) {
            stdevs(i) += std::min(0.5, 1 / (1 + Xt(i)));
            stdevs(i) = std::sqrt(stdevs(i));
        }
    }

    arma::mat Se = (rho * arma::mat(n_stages, n_stages, arma::fill::ones) +
        (1 - rho) * arma::mat(n_stages, n_stages, arma::fill::eye));
    // Convert from matrix of correlations to variances and covariances:
    for (uint32 i = 0; i < n_stages; i++) {
        for (uint32 j = 0; j < n_stages; j++) {
            Se(i,j) *= stdevs(i);
            Se(i,j) *= stdevs(j);
        }
    }

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
 Emigration and immigration of this line to all other plants.
 This function is run individually on `apterous`, `alates`, and `paras` aphids
 inside the `AphidPop::cal_plant_dispersal` function.

 (`this_j` is the index for the plant this line's population is on)

 These do not necessarily match up due to mortality of dispersers.

 For both `emigrants` and `immigrants` matrices, rows are aphid stages and columns
 are plants.
 */
void AphidPop::one_calc_plant_dispersal__(const arma::vec& Xt,
                                          const OnePlant* plant,
                                          arma::mat& emigrants,
                                          arma::mat& immigrants) const {

    if (total_aphids() == 0 || plant->n_plants == 1 ||
        aphid_plant_disp_p <= 0) return;

    const uint32& this_j(plant->this_j);
    const uint32& n_plants(plant->n_plants);
    double n_plants_dbl = static_cast<double>(n_plants);

    double n_leaving;
    double n_arriving;

    // Dispersal for each dispersing stage:
    for (uint32 i = plant_disp_start; i < Xt.n_elem; i++) {

        if (Xt(i) == 0) continue;

        n_leaving = aphid_plant_disp_p * Xt(i);
        emigrants(i, this_j) = n_leaving;

        n_arriving = (n_leaving / n_plants_dbl) * (1 - plant_disp_mort);
        immigrants.row(i) += n_arriving;

    }


    return;
}




/*
 Calculate dispersal of this line to all other plants.
 Emigration doesn't necessarily == immigration due to disperser mortality.
 Note: keep the definition of this in the *.cpp file bc it won't
 compile otherwise.
 */
void AphidPop::calc_plant_dispersal(const OnePlant* plant,
                                    AphidPlantDisps& emigrants,
                                    AphidPlantDisps& immigrants) const {

    one_calc_plant_dispersal__(apterous.X, plant, emigrants.apterous,
                               immigrants.apterous);
    one_calc_plant_dispersal__(alates.X, plant, emigrants.alates,
                               immigrants.alates);
    one_calc_plant_dispersal__(paras.X, plant, emigrants.paras,
                               immigrants.paras);

    return;

}







/*
 Update living aphids (both parasitized and un-parasitized), then
 return the number of newly mummified aphids,
 to be added to the wasps.
 */

double AphidPop::update(const OnePlant* plant,
                        const WaspPop* wasps,
                        const AphidPlantDisps& emigrants,
                        const AphidPlantDisps& immigrants,
                        pcg32& eng) {


    // First subtract emigrants and add immigrants:
    const uint32& this_j(plant->this_j);
    apterous.X -= emigrants.apterous.col(this_j);
    apterous.X += immigrants.apterous.col(this_j);
    alates.X -= emigrants.alates.col(this_j);
    alates.X += immigrants.alates.col(this_j);
    paras.X -= emigrants.paras.col(this_j);
    paras.X += immigrants.paras.col(this_j);

    double nm = 0; // newly mummified

    if (total_aphids() > 0) {

        const double& S(plant->S);
        const double& S_y(plant->S_y);
        arma::vec A_apt = wasps->A(attack_surv);
        // making adult alates not able to be parasitized:
        arma::vec A_ala = A_apt;
        for (uint32 i = alates.field_disp_start_; i < alates.leslie_.n_cols; i++) {
            A_ala[i] = 1;
        }

        double pred_surv = 1 - plant->pred_rate;

        // Starting abundances (used in `process_error`):
        arma::vec apterous_Xt = apterous.X;
        arma::vec alates_Xt = alates.X;
        arma::vec paras_Xt = paras.X;


        // Basic updates for non-parasitized aphids:
        arma::mat LX_apt = apterous.leslie_ * apterous.X;
        arma::mat LX_ala = alates.leslie_ * alates.X;
        apterous.X = (pred_surv * S * A_apt) % LX_apt;
        alates.X = (pred_surv * S * A_ala) % LX_ala;

        double np = 0; // newly parasitized
        np += pred_surv * S_y * arma::as_scalar((1 - A_apt).t() * LX_apt);
        np += pred_surv * S_y * arma::as_scalar((1 - A_ala).t() * LX_ala);

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
         I previously sampled for alate vs apterous offspring, but this
         isn't necessary with proper demographic stochasticity, plus
         the implementation below would work strangely with the continuous
         densities we use.
         */
        // // Sample for # offspring from apterous aphids that are alates:
        // double new_alates = 0;
        // double alate_prop = apterous.alate_prop(plant);
        // if (alate_prop > 0 && apterous.X.front() > 0) {
        //     double lambda_ = alate_prop * apterous.X.front();
        //     pois_distr.param(std::poisson_distribution<uint32>::param_type(lambda_));
        //     new_alates = static_cast<double>(pois_distr(eng));
        //     if (new_alates > apterous.X.front()) new_alates = apterous.X.front();
        // }

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



