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


inline double inv_logit(const double& alpha) {
    return 1 / (1 + std::exp(-1 * alpha));
}





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



//' Update Z and plant_days
//'
//' @noRd
//'
void update_Z_pd(arma::vec& Z, arma::ivec& plant_days,
                 const uint32& t, const arma::cube& N_out,
                 const arma::rowvec& A) {
    for (uint32 i = 0; i < N_out.n_rows; i++) {
        Z(i) = 0;
        // Update # plant days PAST death for time t+1:
        plant_days(i)++;
        for (uint32 j = 0; j < N_out.n_cols; j++) {
            Z(i) += (A(j) * N_out(i,j,t));
        }
    }
    return;
}






// Negative binomial simulator (I can't get the STL one working)
class nbRNG {
    std::poisson_distribution<arma::sword> poisson_rng;
    std::gamma_distribution<double> gamma_rng;


public:

    nbRNG() : poisson_rng(1.0), gamma_rng(1.0, 1.0) {};

    // Set Gamma distribution based on input theta and mu parameters
    void set_gamma(const double& theta, const double& mu) {

        double p = theta / (theta + mu);
        double g_shape = theta;
        double g_scale = (1 - p) / p;

        gamma_rng.param(
            std::gamma_distribution<double>::param_type(g_shape, g_scale));

        return;
    }

    // This is for after the gamma distribution has been set.
    arma::sword sample(pcg32& eng) {
        double lambda_ = gamma_rng(eng);
        poisson_rng.param(
            std::poisson_distribution<arma::sword>::param_type(lambda_));
        return poisson_rng(eng);
    }

};






/*
 Emigration of one line from one plant to all other plants:
 */
arma::ivec emigration(const uint32& i, const uint32& j, const uint32& t,
                      const arma::cube& N_out, const arma::mat& D_mat,
                      const arma::ivec& plant_days,
                      nbRNG& nb_rng, pcg32& eng) {

    arma::ivec disp_out(N_out.n_rows, arma::fill::zeros);

    if (N_out(i, j, t) == 0) return disp_out;

    const arma::rowvec& D(D_mat.row(j));

    // Basing Pr(dispersal > 0) on the # days past plant death:
    const arma::sword& past_death(plant_days(i));
    double p = inv_logit(D(0) + D(1) * past_death + D(2) * past_death * past_death);
    // Sample ~ U(0,1), and if it's > p, then no dispersal occurs:
    double u = runif_01(eng);
    if (u > p) return disp_out;

    // Else we have to sample from negative binomial:
    // Calculate mu:
    double mu_ = D(3) * N_out(i,j,t);
    // Because we're splitting mu_ over the number of plants:
    double np = static_cast<double>(N_out.n_rows);
    mu_ /= np;
    // Now for the theta parameter:
    double theta_ = D(4) / np;
    // Update inner Gamma distribution based on this:
    nb_rng.set_gamma(theta_, mu_);
    for (uint32 ii = 0; ii < disp_out.n_elem; ii++) {
        if (ii == i) {
            /*
             No point in sampling this combination bc it's dispersal back to
             the same plant, so summing immigrants and emigrants would cancel out anyway.
             */
            disp_out(ii) = 0;
        } else {
            disp_out(ii) = nb_rng.sample(eng);
        }
    }

    return disp_out;
}

/*
 Update dispersal
 */
void update_dispersal(arma::imat& emigrants, arma::imat& immigrants,
                      const uint32& t,
                      const arma::cube& N_out, const arma::mat& D_mat,
                      const arma::ivec& plant_days,
                      nbRNG& nb_rng, pcg32& eng) {

    immigrants.fill(0);
    uint32 n_lines = emigrants.n_cols;
    uint32 n_plants = emigrants.n_rows;

    for (uint32 j = 0; j < n_lines; j++) {
        for (uint32 i = 0; i < n_plants; i++) {
            // Emigrants of line j from plant i to all others:
            arma::ivec e_ij = emigration(i, j, t, N_out, D_mat, plant_days, nb_rng, eng);
            emigrants(i, j) = arma::accu(e_ij);
            immigrants.col(j) += e_ij;
        }
    }

    return;

}







//' Replace plants
//'
//' @noRd
//'
void replace_plants(uint32& repl_ind,
                    arma::ivec& plant_days,
                    arma::cube& N_out,
                    arma::Mat<unsigned short>& extinct,
                    std::normal_distribution<double>& normal_rng,
                    pcg32& eng,
                    const uint32& t,
                    const sint32& repl_age,
                    const double& plant_death_age_mean,
                    const double& plant_death_age_sd) {

    uint32 n_plants = N_out.n_rows;
    uint32 n_lines = N_out.n_cols;

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
    if (replaced.size() == 0) return;

    /*
    In the rare event that all plants need replaced, aphids will
    be dispersed evenly across all new plants:
    */
    if (not_replaced.size() == 0) {
        // Adjust pools of aphid lines in replaced and non-replaced plants:
        for (uint32 j = 0; j < n_lines; j++) {
            // Pool of aphids of line j:
            double aphid_pool_j = arma::accu(N_out(arma::span::all, arma::span(j),
                                                   arma::span(t+1)));
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
        return;
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
            extinct(i,j) = 0;
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
              const arma::mat& D_mat,
              const double& process_error,
              const double& plant_death_age_mean,
              const double& plant_death_age_sd,
              const std::vector<uint32>& repl_times,
              const sint32& repl_age,
              const arma::mat& log_morts,
              const double& extinct_N,
              const bool& by_cage,
              pcg32& eng, arma::cube& N_out) {

    uint32 n_plants = N_0.n_rows;
    uint32 n_lines = N_0.n_cols;

    std::normal_distribution<double> normal_rng(0.0, 1.0);
    nbRNG nb_rng;

    // Fill with initial values:
    N_out.slice(0) = N_0;

    // Matrix to keep track of extinctions:
    arma::Mat<unsigned short> extinct(n_plants, n_lines, arma::fill::zeros);
    // Vector to keep track of # days PAST death for plants
    // (these will start out negative)
    arma::ivec plant_days(n_plants);
    for (uint i = 0; i < n_plants; i++) {
        plant_days(i) = sample_plant_days(plant_death_age_mean, plant_death_age_sd,
                                          normal_rng, eng);
    }
    // Index to keep track of position in `repl_times`
    uint32 repl_ind = 0;
    if (repl_times.front() == 0) repl_ind++;
    // Matrices keeping track of numbers of dispersed aphids:
    arma::imat emigrants(n_plants, n_lines);
    arma::imat immigrants(n_plants, n_lines);
    // Summed total aphids per plant
    arma::vec Ntot(n_plants);
    // Summed density dependences * N_t among all lines for each plant
    arma::vec Z(n_plants);

    for (uint32 t = 0; t < max_t; t++) {

        // Update Z and plant_days:
        update_Z_pd(Z, plant_days, t, N_out, A);
        // Generate numbers of dispersed aphids:
        update_dispersal(emigrants, immigrants, t, N_out, D_mat, plant_days, nb_rng, eng);

        // Start out new N based on previous time step
        N_out.slice(t+1) = N_out.slice(t);

        /*
         Go back through to simulate at time t+1
        */
        for (uint32 i = 0; i < n_plants; i++) {

            for (uint32 j = 0; j < n_lines; j++) {

                // Dispersal of this aphid line to and from other plants.
                double imm = immigrants(i, j);
                double em = emigrants(i, j);

                // If it's extinct and no one's coming in, make 0 and skip the rest:
                if (extinct(i,j) == 1 && (imm - em) <= 0.0) {
                    N_out(i,j,t+1) = 0;
                    continue;
                }
                // If it's extinct and immigrants are coming in, adjust extinct matrix:
                if (extinct(i,j) == 1 && (imm - em) > 0.0) {
                    extinct(i,j) = 0;
                }

                /*
                Add growth, density dependence, and process error.
                If the plant is past plant-death age, then
                 (1) plant-death-induced mortality replaces normal growth and
                    density dependence and
                 (2) immigrants experience this mortality as well.
                */
                if (plant_days[i] <= 0) {
                    N_out(i,j,t+1) *= std::exp(R[j] * (1 - Z[i]) +
                        normal_rng(eng) * process_error);
                } else {
                    // (` - 1` below is to convert to C++ indices)
                    N_out(i,j,t+1) *= std::exp(log_morts(plant_days[i] - 1, j) +
                        normal_rng(eng) * process_error);
                    imm *= std::exp(log_morts(plant_days[i] - 1, j));
                    if (imm < 1) imm = 0;
                }

                // Add dispersal:
                N_out(i,j,t+1) += (imm - em);

                // Check for extinction:
                if (N_out(i,j,t+1) < extinct_N) {
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
        if ((t+1) == repl_times[repl_ind]) {
            replace_plants(repl_ind, plant_days, N_out, extinct, normal_rng, eng, t,
                           repl_age, plant_death_age_mean, plant_death_age_sd);
        }


    }

    // Condense by cage if requested:
    if (by_cage) {
        for (uint32 j = 0; j < n_lines; j++) {
            for (uint32 t = 0; t <= max_t; t++) {
                double cage_sum = arma::accu(N_out(arma::span(), arma::span(j),
                                                   arma::span(t)));
                N_out(0, j, t) = cage_sum;
            }
        }
        N_out.resize(1, n_lines, max_t + 1);
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
                                   const arma::mat& D_mat,
                                   const double& process_error,
                                   const arma::rowvec& plant_mort_0,
                                   const arma::rowvec& plant_mort_1,
                                   const double& plant_death_age_mean,
                                   const double& plant_death_age_sd,
                                   const std::vector<uint32>& repl_times,
                                   const sint32& repl_age,
                                   const double& extinct_N,
                                   const uint32& n_cores,
                                   const bool& by_cage,
                                   const bool& show_progress) {


    // log(mortality) at a time x after plant death starts, to keep from having to do
    // this calculation many times.
    const arma::mat log_morts = make_log_morts(max_t, plant_mort_0, plant_mort_1);

    const std::vector<std::vector<uint64> > seeds = mc_seeds(n_cores);

    const uint32 n_plants = N_0.n_rows;
    const uint32 n_lines = N_0.n_cols;

    const sint32 repl_age_ = -1 * repl_age;

    std::string err_msg = "";

    if (R.n_elem != n_lines) {
        err_msg += "ERROR: R.n_elem != n_lines\n";
    }
    if (A.n_elem != n_lines) {
        err_msg += "ERROR: A.n_elem != n_lines\n";
    }
    if (D_mat.n_rows != n_lines) {
        err_msg += "ERROR: D_mat.n_rows != n_lines\n";
    }
    if (D_mat.n_cols != 5) {
        err_msg += "ERROR: D_mat.n_cols != 5\n";
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
        sim_cage(N_0, max_t, R, A, D_mat, process_error,
                 plant_death_age_mean, plant_death_age_sd,
                 repl_times, repl_age_, log_morts, extinct_N, by_cage, eng, rep_cube);
        p.increment();
    }

#ifdef _OPENMP
}
#endif


    return reps_out;

}
