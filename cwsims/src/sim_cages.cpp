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


//' Simulate one cage.
//'
//'
//' @noRd
//'
void sim_cage(const arma::mat& N_0, const uint32& max_t,
              const arma::rowvec& R, const arma::rowvec& A,
              const arma::rowvec& D_0, const arma::rowvec& D_1,
              const double& process_error,
              const uint32& plant_death_age,
              const std::vector<uint32>& repl_times,
              const std::vector<std::vector<uint32> >& repl_plants,
              const arma::vec& log_morts,
              pcg32& eng, arma::cube& N_out,
              Progress& p) {

    uint32 n_plants = N_0.n_rows;
    uint32 n_lines = N_0.n_cols;

    // Fill with initial values:
    N_out.slice(0) = N_0;

    // Matrix to keep track of extinctions:
    arma::mat extinct(n_plants, n_lines, arma::fill::zeros);
    // Vector to keep track of plant ages:
    std::vector<uint32> plant_ages(n_plants, 0);
    // Index to keep track of position in `repl_times`
    uint32 repl_ind = 0;
    // Matrices keeping track of numbers of dispersed aphids:
    arma::mat emigrants(n_plants, n_lines, arma::fill::zeros);
    arma::mat immigrants(n_plants, n_lines, arma::fill::zeros);
    // Dispersal lambdas (for Poisson distr.) for all lines and plants
    arma::mat D_lambdas(n_plants, n_lines);
    // Summed density dependences * N_t among all lines for each plant
    arma::vec Z(n_plants);

    std::poisson_distribution<uint32> poisson_rng(1.0);
    std::normal_distribution<double> normal_rng(0.0, 1.0);

    for (uint32 t = 0; t < max_t; t++) {

        if (p.is_aborted()) return;


        /*
         Go through once to get/update some parameters first:
        */
        for (uint32 i = 0; i < n_plants; i++) {
            Z[i] = 0;
            // Update plant ages to time t+1:
            plant_ages[i]++;
            for (uint32 j = 0; j < n_lines; j++) {
                D_lambdas(i, j) = std::exp(D_0(j) + D_1(j) * std::log(N_out(i,j,t)));
                Z[i] += (A(j) * N_out(i,j,t));
            }
        }
        // Generate numbers of dispersed aphids:
        for (uint32 j = 0; j < n_lines; j++) {
            uint32 dispersed_;
            for (uint32 from_i = 0; from_i < n_plants; from_i++) {
                for (uint32 to_i = 0; to_i < n_plants; to_i++) {
                    if (from_i == to_i) continue;
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

        /*
         Go back through to simulate at time t+1
        */
        for (uint32 i = 0; i < n_plants; i++) {

            for (uint32 j = 0; j < n_lines; j++) {

                // Calculate the net influx of this aphid line from other plants.
                const double& immigration(immigrants(i, j));
                const double& emigration(emigrants(i, j));

                // If it's extinct and no one's coming in, skip the rest
                if (extinct(i,j) == 1 && immigration == 0) continue;

                // Start out new N based on previous time step
                N_out(i,j,t+1) = N_out(i,j,t);
                /*
                Now add growth, density dependence, and process error.
                If the plant is past plant-death age, then plant-death-induced
                mortality replaces normal growth and density dependence.
                */
                if (plant_ages[i] <= plant_death_age) {
                    N_out(i,j,t+1) *= std::exp(R[j] * (1 - Z[i]) +
                        normal_rng(eng) * process_error);
                } else {
                    // (` - 1` below is to convert to C++ indices)
                    uint32 after_death = plant_ages[i] - plant_death_age - 1;
                    N_out(i,j,t+1) *= std::exp(log_morts[after_death] +
                        normal_rng(eng) * process_error);
                }

                // Add dispersal:
                N_out(i,j,t+1) += (immigration - emigration);

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

            const std::vector<uint32>& replaced(repl_plants[repl_ind]);

            uint32 n_replaced = replaced.size();

            // Update replaced plants' ages:
            for (const uint32& i : replaced) plant_ages[i] = 0;

            // Make vector of non-replaced plants:
            std::vector<uint32> non_replaced;
            uint32 n_non_replaced = n_plants - n_replaced;
            non_replaced.reserve(n_non_replaced);
            for (uint32 i = 0; i < n_plants; i++) {
                auto iter = std::find(replaced.begin(), replaced.end(), i);
                if (iter == replaced.end()) non_replaced.push_back(i);
            }

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
                aphid_pool_j /= static_cast<double>(n_non_replaced);
                for (const uint32& i : non_replaced) {
                    // Check that it's not adding < 1 aphid to an extinct plant:
                    if (extinct(i,j) == 1 && aphid_pool_j < 1) continue;
                    // If not, add `aphid_pool_j`:
                    N_out(i,j,t+1) += aphid_pool_j;
                }
            }
            // Now adjust index for `repl_times`
            repl_ind += 1;
        }


    }

    p.increment();

    return;
}



//' log(mortality) at a time x after plant death starts.
//'
//' Doing this here keeps from having to do this calculation many times.
//' max(`max_t`, `max_t - plant_death_age`) should be more than enough.
//'
//' @noRd
//'
arma::vec make_log_morts(const uint32& max_t,
                         const std::vector<double>& plant_mort_coefs,
                         const uint32& plant_death_age) {

    uint32 n_times = 1;
    if (max_t > plant_death_age) n_times = max_t - plant_death_age;

    arma::vec log_morts(n_times);
    for (uint32 x = 0; x < n_times; x++) {
        double alpha = plant_mort_coefs[0] + plant_mort_coefs[1] *
            static_cast<double>(x + 1);  // ` + 1` is to go from C++ index to # past death
        // Doing log(inverse_logit(alpha)):
        log_morts[x] = alpha - std::log(std::exp(alpha) + 1);
    }
    return log_morts;
}




// Get total max over nested vectors of numeric objects
// (or anything that can be compared)
//
template<typename T>
T nested_max2(const std::vector<std::vector<T> >& input) {
    T total_max = input[0][0];
    for (auto iter = input.begin(); iter != input.end(); ++iter) {
        T max_ = *std::max_element((*iter).begin(), (*iter).end());
        if (max_ > total_max) total_max = max_;
    }
    return total_max;
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
                                   const std::vector<double>& plant_mort_coefs,
                                   const uint32& plant_death_age,
                                   const std::vector<uint32>& repl_times,
                                   const std::vector<std::vector<uint32> >& repl_plants,
                                   const uint32& n_cores,
                                   const bool& show_progress) {


    // log(mortality) at a time x after plant death starts, to keep from having to do
    // this calculation many times.
    const arma::vec log_morts = make_log_morts(max_t, plant_mort_coefs, plant_death_age);

    const std::vector<std::vector<uint64> > seeds = mc_seeds(n_cores);

    const uint32 n_plants = N_0.n_rows;
    const uint32 n_lines = N_0.n_cols;

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
    if (plant_mort_coefs.size() != 2) {
        err_msg += "ERROR: plant_mort_coefs.size() != 2\n";
    }
    if (repl_times.size() == 0) {
        err_msg += "ERROR: repl_times.size() == 0\n";
    }
    if (repl_plants.size() != repl_times.size()) {
        err_msg += "ERROR: repl_plants.size() != repl_times.size()\n";
    }
    if (nested_max2<uint32>(repl_plants) >= n_plants) {
        err_msg += "ERROR: repl_plants has value >= n_plants\n";
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
        arma::cube& rep_cube(reps_out[r]);
        sim_cage(N_0, max_t, R, A, D_0, D_1, process_error, plant_death_age,
                 repl_times, repl_plants, log_morts, eng, rep_cube, p);
    }

#ifdef _OPENMP
}
#endif

    return reps_out;

}
