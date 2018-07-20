#include <RcppArmadillo.h>
#include <random> // C++11 distributions
#include <pcg/pcg_random.hpp> // pcg PRNG
#include <progress.hpp>  // for the progress bar

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cwsims_types.hpp"
#include "pcg.hpp"

using namespace Rcpp;



//' Sample the number of days PAST plant death.
//'
//' @noRd
//'
inline sint32 sample_plant_days(const double& days_mean,
                                const double& days_sd,
                                std::normal_distribution<double>& normal_rng,
                                pcg32& eng) {
    double days_ = -1 * std::round(days_mean + normal_rng(eng) * days_sd);
    if (days_ > -1) days_ = -1;
    sint32 days = static_cast<sint32>(days_);
    // This offsets adding one to this (i.e., using `++`) before the first iteration:
    days--;
    return days;
}



//' Update Z, plant_days, and D_lambdas
//'
//' @noRd
//'
void update_Z_pd_L(arma::vec& Z, arma::ivec& plant_days, arma::mat& D_lambdas,
                   const uint32& t, const arma::cube& N_out,
                   const arma::rowvec& A,
                   const arma::rowvec& D_0, const arma::rowvec& D_1) {
    for (uint32 i = 0; i < Z.n_elem; i++) {
        Z(i) = 0;
        // Update # plant days PAST death for time t+1:
        plant_days(i)++;
        for (uint32 j = 0; j < D_lambdas.n_cols; j++) {
            if (N_out(i,j,t) == 0) {
                D_lambdas(i, j) = 0;
                continue;
            }
            D_lambdas(i, j) = std::exp(D_0(j) + D_1(j) * std::log(N_out(i,j,t)));
            Z[i] += (A(j) * N_out(i,j,t));
        }
    }
    return;
}


/*
 Update dispersal with no Poisson overdispersion
 */
void update_dispersal(arma::mat& emigrants, arma::mat& immigrants,
                      std::poisson_distribution<uint32>& poisson_rng,
                      pcg32& eng,
                      const arma::mat& D_lambdas) {
    emigrants.fill(0);
    immigrants.fill(0);
    uint32 n_lines = emigrants.n_cols;
    uint32 n_plants = emigrants.n_rows;
    for (uint32 j = 0; j < n_lines; j++) {
        uint32 dispersed_;
        for (uint32 from_i = 0; from_i < n_plants; from_i++) {
            for (uint32 to_i = 0; to_i < n_plants; to_i++) {
                if (from_i == to_i) continue;
                if (D_lambdas(from_i, j) == 0) continue;
                // Calculate lambda and reset distribution to it:
                double lambda_ = D_lambdas(from_i, j) /
                    (static_cast<double>(n_plants) - 1.0);
                poisson_rng.param(
                    std::poisson_distribution<uint32>::param_type(lambda_));
                // Draw from distribution
                dispersed_ = poisson_rng(eng);
                emigrants(from_i, j) += dispersed_;
                immigrants(to_i, j) += dispersed_;
            }
        }
    }
    return;
}


/*
 Overloaded version for overdispersion:
 */
void update_dispersal(arma::mat& emigrants, arma::mat& immigrants,
                      std::poisson_distribution<uint32>& poisson_rng,
                      std::normal_distribution<double>& normal_rng,
                      pcg32& eng,
                      const arma::mat& D_lambdas,
                      const double& overdispersion_sd) {
    emigrants.fill(0);
    immigrants.fill(0);
    uint32 n_lines = emigrants.n_cols;
    uint32 n_plants = emigrants.n_rows;
    for (uint32 j = 0; j < n_lines; j++) {
        uint32 dispersed_;
        for (uint32 from_i = 0; from_i < n_plants; from_i++) {
            for (uint32 to_i = 0; to_i < n_plants; to_i++) {
                if (from_i == to_i) continue;
                if (D_lambdas(from_i, j) == 0) continue;
                // Calculate lambda and reset distribution to it:
                double lambda_ = D_lambdas(from_i, j) /
                    (static_cast<double>(n_plants) - 1.0);
                poisson_rng.param(
                    std::poisson_distribution<uint32>::param_type(lambda_));
                // Draw from distribution
                dispersed_ = poisson_rng(eng);
                /* Add overdispersion: */
                dispersed_ += std::round(normal_rng(eng) * overdispersion_sd);
                emigrants(from_i, j) += dispersed_;
                immigrants(to_i, j) += dispersed_;
            }
        }
    }
    return;

}


//' Simulate one cage.
//'
//'
//' @noRd
//'
void sim_cage(const arma::mat& N_0, const uint32& max_t,
              const arma::rowvec& R, const arma::rowvec& A,
              const arma::rowvec& D_0, const arma::rowvec& D_1,
              const double& process_error,
              const double& plant_death_age_mean,
              const double& plant_death_age_sd,
              const std::vector<uint32>& repl_times,
              const sint32& repl_age,
              const arma::mat& log_morts,
              pcg32& eng, arma::cube& N_out) {

    uint32 n_plants = N_0.n_rows;
    uint32 n_lines = N_0.n_cols;

    std::poisson_distribution<uint32> poisson_rng(1.0);
    std::normal_distribution<double> normal_rng(0.0, 1.0);

    // Fill with initial values:
    N_out.slice(0) = N_0;

    // Matrix to keep track of extinctions:
    arma::mat extinct(n_plants, n_lines, arma::fill::zeros);
    // Vector to keep track of # days PAST death for plants
    // (these will start out negative)
    arma::ivec plant_days(n_plants);
    for (uint i = 0; i < n_plants; i++) {
        plant_days(i) = sample_plant_days(plant_death_age_mean, plant_death_age_sd,
                                          normal_rng, eng);
    }
    // Index to keep track of position in `repl_times`
    uint32 repl_ind = 0;
    // Matrices keeping track of numbers of dispersed aphids:
    arma::mat emigrants(n_plants, n_lines, arma::fill::zeros);
    arma::mat immigrants(n_plants, n_lines, arma::fill::zeros);
    // Dispersal lambdas (for Poisson distr.) for all lines and plants
    arma::mat D_lambdas(n_plants, n_lines);
    // Summed density dependences * N_t among all lines for each plant
    arma::vec Z(n_plants);

    for (uint32 t = 0; t < max_t; t++) {

        // Update Z, plant_days, and D_lambdas:
        update_Z_pd_L(Z, plant_days, D_lambdas, t, N_out, A, D_0, D_1);
        // Generate numbers of dispersed aphids:
        update_dispersal(emigrants, immigrants, poisson_rng, eng, D_lambdas);

        // Start out new N based on previous time step
        N_out.slice(t+1) = N_out.slice(t);

        /*
         Go back through to simulate at time t+1
        */
        for (uint32 i = 0; i < n_plants; i++) {

            for (uint32 j = 0; j < n_lines; j++) {

                // Calculate the net influx of this aphid line from other plants.
                const double& immigration(immigrants(i, j));
                const double& emigration(emigrants(i, j));
                double net_income = immigration - emigration;

                // If it's extinct and no one's coming in, make 0 and skip the rest:
                if (extinct(i,j) == 1 && net_income <= 0) {
                    N_out(i,j,t+1) = 0;
                    continue;
                }
                // If it's extinct and immigrants are coming in, adjust extinct matrix:
                if (extinct(i,j) == 1 && net_income > 0) extinct(i,j) = 0;

                /*
                Add growth, density dependence, and process error.
                If the plant is past plant-death age, then plant-death-induced
                mortality replaces normal growth and density dependence.
                */
                if (plant_days[i] <= 0) {
                    N_out(i,j,t+1) *= std::exp(R[j] * (1 - Z[i]) +
                        normal_rng(eng) * process_error);
                } else {
                    // (` - 1` below is to convert to C++ indices)
                    N_out(i,j,t+1) *= std::exp(log_morts(plant_days[i] - 1, j) +
                        normal_rng(eng) * process_error);
                }

                // Add dispersal:
                N_out(i,j,t+1) += net_income;

                // Check for extinction:
                if (N_out(i,j,t+1) < 1) {
                    N_out(i,j,t+1) = 0;
                    extinct(i, j) = 1;
                }

            }

        }


        /*
        If it's a replacement time point, disperse all aphids from replaced
        plants to others
        */
        // To make sure this doesn't past bounds:
        if (repl_ind >= repl_times.size()) continue;
        // Now check:
        if (t == repl_times[repl_ind]) {

            // Adjust index for `repl_times`
            repl_ind++;

            /*
             Find the plants to replace and those not to, and update plant_days for
             those that are replaced.
              */
            std::vector<uint32> replaced;
            std::vector<uint32> not_replaced;
            for (uint32 i = 0; i < n_plants; i++) {
                if (plant_days(i) >= repl_age) {
                    replaced.push_back(i);
                    plant_days(i) = sample_plant_days(plant_death_age_mean,
                               plant_death_age_sd, normal_rng, eng);
                } else {
                    not_replaced.push_back(i);
                }
            }

            // If none to replace, then continue
            if (replaced.size() == 0) continue;

            /*
             In the rare event that all plants need replaced, aphids will
             be dispersed evenly across all plants:
             */
            if (not_replaced.size() == 0) {
                // Adjust pools of aphid lines in replaced and non-replaced plants:
                for (uint32 j = 0; j < n_lines; j++) {
                    // Pool of aphids of line j:
                    double aphid_pool_j = 0;
                    // Calculate aphid_pool_j:
                    for (uint32 i = 0; i < n_plants; i++) {
                        aphid_pool_j += N_out(i,j,t+1);
                    }
                    // Now add a portion to all plants:
                    aphid_pool_j /= static_cast<double>(n_plants);
                    // Make it extinct if it's trying to add < 1 aphid:
                    if (aphid_pool_j < 1) {
                        extinct.col(j).fill(1);
                        N_out(arma::span::all, arma::span(j), arma::span(t+1)).fill(0);
                    // Otherwise, fill with `aphid_pool_j`
                    } else {
                        N_out(arma::span::all, arma::span(j), arma::span(t+1)).fill(
                                aphid_pool_j);
                    }
                }
                continue;
            }

            /*
             Otherwise, disperse across all non-replaced plants only:
             */
            // Adjust pools of aphid lines in replaced and non-replaced plants:
            for (uint32 j = 0; j < n_lines; j++) {
                // Pool of aphids of line j to disperse across non-replaced plants:
                double aphid_pool_j = 0;
                // Calculate aphid_pool_j and remove from replaced plant(s):
                for (const uint32& i : replaced) {
                    aphid_pool_j += N_out(i,j,t+1);
                    N_out(i,j,t+1) = 0;
                    extinct(i, j) = 1;
                }
                // Now add a portion to non-replaced plant(s):
                aphid_pool_j /= static_cast<double>(not_replaced.size());
                for (const uint32& i : not_replaced) {
                    // Check that it's not adding < 1 aphid to an extinct plant:
                    if (extinct(i,j) == 1 && aphid_pool_j < 1) continue;
                    // If not, add `aphid_pool_j`:
                    N_out(i,j,t+1) += aphid_pool_j;
                }
            }

        }


    }

    return;
}



//' log(mortality) at a time x after plant death starts.
//'
//' Doing this here keeps from having to do this calculation many times.
//' `max_t` should be more than enough.
//'
//' @noRd
//'
arma::mat make_log_morts(const uint32& max_t,
                         const arma::rowvec& plant_mort_0,
                         const arma::rowvec& plant_mort_1) {

    uint32 n_lines = plant_mort_0.n_elem;
    uint32 n_times = max_t;

    arma::mat log_morts(n_times, n_lines);
    for (uint32 x = 0; x < n_times; x++) {
        arma::rowvec alpha = plant_mort_0 + plant_mort_1 *
            static_cast<double>(x + 1);  // ` + 1` is to go from C++ index to # past death
        // Doing log(inverse_logit(alpha)):
        log_morts.row(x) = alpha - arma::log(arma::exp(alpha) + 1);
    }
    return log_morts;
}





//' Simulate multiple cages.
//'
//'
//' @inheritParams sim_cages
//'
//' @noRd
//'
//'
//[[Rcpp::export]]
std::vector<arma::cube> sim_cages_(const uint32& n_cages,
                                   const arma::mat& N_0,
                                   const uint32& max_t,
                                   const arma::rowvec& R,
                                   const arma::rowvec& A,
                                   const arma::rowvec& D_0,
                                   const arma::rowvec& D_1,
                                   const double& process_error,
                                   const arma::rowvec& plant_mort_0,
                                   const arma::rowvec& plant_mort_1,
                                   const double& plant_death_age_mean,
                                   const double& plant_death_age_sd,
                                   const std::vector<uint32>& repl_times,
                                   const sint32& repl_age,
                                   const uint32& n_cores,
                                   const bool& show_progress) {


    // log(mortality) at a time x after plant death starts, to keep from having to do
    // this calculation many times.
    const arma::mat log_morts = make_log_morts(max_t, plant_mort_0, plant_mort_1);

    const std::vector<std::vector<uint64> > seeds = mc_seeds(n_cores);

    const uint32 n_plants = N_0.n_rows;
    const uint32 n_lines = N_0.n_cols;

    const sint32 repl_age_ = -1 * static_cast<sint32>(repl_age);

    std::string err_msg = "";

    if (R.n_elem != n_lines) {
        err_msg += "ERROR: R.n_elem != n_lines\n";
    }
    if (A.n_elem != n_lines) {
        err_msg += "ERROR: A.n_elem != n_lines\n";
    }
    if (D_0.n_elem != n_lines) {
        err_msg += "ERROR: D_0.n_elem != n_lines\n";
    }
    if (D_1.n_elem != n_lines) {
        err_msg += "ERROR: D_1.n_elem != n_lines\n";
    }
    if (plant_mort_0.n_elem != n_lines) {
        err_msg += "ERROR: plant_mort_0.n_elem != n_lines\n";
    }
    if (plant_mort_1.n_elem != n_lines) {
        err_msg += "ERROR: plant_mort_1.n_elem != n_lines\n";
    }
    if (repl_times.size() == 0) {
        err_msg += "ERROR: repl_times.size() == 0\n";
    }

    if (err_msg.size() > 0) throw(Rcpp::exception(err_msg.c_str(), false));



    Progress p(n_cages, show_progress);

    // Create output object:
    std::vector<arma::cube> reps_out(n_cages, arma::cube(n_plants, n_lines, max_t+1));

#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(n_cores) if(n_cores > 1)
{
#endif

    std::vector<uint64> active_seeds;

    // Write the active seed per core or just write one of the seeds.
#ifdef _OPENMP
    uint32 active_thread = omp_get_thread_num();
    active_seeds = seeds[active_thread];
#else
    active_seeds = seeds[0];
#endif

    // pcg prng
    pcg32 eng = seeded_pcg(active_seeds);

    // Parallelize the Loop
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (uint32 r = 0; r < n_cages; r++) {
        if (p.is_aborted()) continue;
        arma::cube& rep_cube(reps_out[r]);
        sim_cage(N_0, max_t, R, A, D_0, D_1, process_error,
                 plant_death_age_mean, plant_death_age_sd,
                 repl_times, repl_age_, log_morts, eng, rep_cube);
        p.increment();
    }

#ifdef _OPENMP
}
#endif

    return reps_out;

}
