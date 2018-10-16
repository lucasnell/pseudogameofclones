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



//[[Rcpp::plugins(cpp11)]]



//' Sample the number of days PAST plant death.
//'
//' @noRd
//'
inline sint32 sample_plant_days(const double& days_mean,
                                const double& days_sd,
                                const arma::ivec& plant_death_ages, const uint32& i,
                                std::normal_distribution<double>& normal_rng,
                                pcg32& eng) {
    if (plant_death_ages.n_elem > 0) {
        return -1 * plant_death_ages(i) - 1;
    }
    double days_ = -1 * std::round(days_mean + normal_rng(eng) * days_sd);
    if (days_ > -1) days_ = -1;
    sint32 days = static_cast<sint32>(days_);
    // This offsets adding one to this (i.e., using `++`) before the first iteration:
    days--;
    return days;
}








/*
 Emigration of one line from one plant to all other plants:
 */
arma::vec emigration(const uint32& i, const uint32& j,
                     const arma::mat& Nt, const arma::mat& D_vec,
                     std::poisson_distribution<arma::sword>& poisson_rng, pcg32& eng) {

    arma::vec disp_out(Nt.n_rows, arma::fill::zeros);

    if (Nt(i, j) == 0) return disp_out;

    // Else we have to sample from Poisson:
    // Calculate lambda:
    double lambda_ = D_vec(j) * Nt(i,j);
    // Because we're splitting mu_ over the number of plants:
    double np = static_cast<double>(Nt.n_rows);
    lambda_ /= np;
    // Update Poisson distribution based on this:
    poisson_rng.param(
        std::poisson_distribution<arma::sword>::param_type(lambda_));
    for (uint32 ii = 0; ii < disp_out.n_elem; ii++) {
        if (ii == i) {
            /*
             No point in sampling this combination bc it's dispersal back to
             the same plant, so summing immigrants and emigrants would cancel out anyway.
             */
            disp_out(ii) = 0;
        } else {
            disp_out(ii) = static_cast<double>(poisson_rng(eng));
        }
    }
    // Making absolutely sure that dispersal never exceeds the number possible:
    double total_emigrants = arma::accu(disp_out);
    if (total_emigrants > Nt(i,j)) {
        double extras = total_emigrants - Nt(i,j);
        std::vector<uint32> extra_inds;
        extra_inds.reserve(Nt.n_rows);
        for (uint32 ii = 0; ii < Nt.n_rows; ii++) {
            if (disp_out(ii) > 0) extra_inds.push_back(ii);
        }
        while (extras > 0) {
            uint32 rnd = runif_01(eng) * extra_inds.size();
            disp_out(extra_inds[rnd])--;
            if (disp_out(extra_inds[rnd]) == 0) {
                extra_inds.erase(extra_inds.begin() + rnd);
            }
            extras--;
        }
    }

    return disp_out;
}

/*
 Same thing as above, but overloaded for not including dispersal stochasticity
 */
arma::vec emigration(const uint32& i, const uint32& j,
                      const arma::mat& Nt, const arma::mat& D_vec) {


    if (Nt(i, j) == 0) return arma::vec(Nt.n_rows, arma::fill::zeros);

    // Else we have to calculate E(# dispersed)
    // Calculate lambda:
    double lambda_ = D_vec(j) * Nt(i,j);
    // Because we're splitting mu_ over the number of plants:
    double np = static_cast<double>(Nt.n_rows);
    lambda_ /= np;

    arma::vec disp_out(Nt.n_rows);
    disp_out.fill(lambda_);
    disp_out(i) = 0;

    return disp_out;
}


/*
 Update dispersal
 */
void update_dispersal(arma::mat& emigrants, arma::mat& immigrants,
                      const arma::mat& Nt, const arma::mat& D_vec,
                      std::poisson_distribution<arma::sword>& poisson_rng, pcg32& eng) {

    immigrants.fill(0);
    uint32 n_lines = emigrants.n_cols;
    uint32 n_plants = emigrants.n_rows;

    for (uint32 j = 0; j < n_lines; j++) {
        for (uint32 i = 0; i < n_plants; i++) {
            // Emigrants of line j from plant i to all others:
            arma::vec e_ij = emigration(i, j, Nt, D_vec, poisson_rng, eng);
            emigrants(i, j) = arma::accu(e_ij);
            immigrants.col(j) += e_ij;
        }
    }

    return;

}

/*
 Overloaded for not including dispersal stochasticity
 */
void update_dispersal(arma::mat& emigrants, arma::mat& immigrants,
                      const arma::mat& Nt, const arma::mat& D_vec) {

    immigrants.fill(0);
    uint32 n_lines = emigrants.n_cols;
    uint32 n_plants = emigrants.n_rows;

    for (uint32 j = 0; j < n_lines; j++) {
        for (uint32 i = 0; i < n_plants; i++) {
            // Emigrants of line j from plant i to all others:
            arma::vec e_ij = emigration(i, j, Nt, D_vec);
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
                    arma::mat& Nt1,
                    arma::Mat<unsigned short>& extinct,
                    std::normal_distribution<double>& normal_rng,
                    pcg32& eng,
                    const sint32& repl_age,
                    const double& plant_death_age_mean,
                    const double& plant_death_age_sd,
                    const arma::ivec& plant_death_ages) {

    uint32 n_plants = Nt1.n_rows;
    uint32 n_lines = Nt1.n_cols;

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
            plant_days(i) = sample_plant_days(plant_death_age_mean, plant_death_age_sd,
                       plant_death_ages, i,
                       normal_rng, eng);
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
            double aphid_pool_j = arma::accu(Nt1(arma::span::all, arma::span(j)));
            // Now add a portion to all plants:
            aphid_pool_j /= static_cast<double>(n_plants);
            // Make it extinct if it's trying to add < 1 aphid:
            if (aphid_pool_j < 1) {
                extinct.col(j).fill(1);
                Nt1(arma::span::all, arma::span(j)).fill(0);
            // Otherwise, fill with `aphid_pool_j`
            } else {
                Nt1(arma::span::all, arma::span(j)).fill(aphid_pool_j);
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
            aphid_pool_j += Nt1(i,j);
            Nt1(i,j) = 0;
            extinct(i, j) = 1;
        }
        // Now add a portion to non-replaced plant(s):
        aphid_pool_j /= static_cast<double>(not_replaced.size());
        for (const uint32& i : not_replaced) {
            // Check that it's not adding < 1 aphid to an extinct plant:
            if (extinct(i,j) == 1 && aphid_pool_j < 1) continue;
            // If not, add `aphid_pool_j`:
            Nt1(i,j) += aphid_pool_j;
            extinct(i,j) = 0;
        }
    }
    return;
}



/*
 Originally N_out was arma::cube(n_plants, n_lines, max_t+1)
 Nt and Nt1 are arma::mat(n_plants, n_lines)

 // Mean N for each line across all plants:
 arma::cube N_bycage(max_t+1, n_lines, n_cages);
 // Mean proportion of zeros for each plant across all lines and all plants:
 arma::mat zero_byplant(n_plants, n_cages);
 */
void summarize_output(const arma::mat& Nt1,
                      arma::cube& N_bycage,
                      arma::mat& logN_byplant,
                      arma::mat& zero_byplant,
                      const uint32& r,
                      const uint32& t) {

    N_bycage(arma::span(t+1), arma::span(), arma::span(r)) = arma::sum(Nt1, 0);

    arma::vec plant_sums = arma::sum(Nt1, 1);  // sum all lines for each plant
    for (uint32 i = 0; i < Nt1.n_rows; i++) {
        if (plant_sums(i) == 0) zero_byplant(i, r)++;
        logN_byplant(i, r) += std::log(plant_sums(i) + 1);
    }
    return;
}





//' Simulate one cage.
//'
//'
//' @noRd
//'
void sim_cage(const arma::mat& N_0,
              const arma::rowvec& R,
              const arma::rowvec& A,
              const arma::vec& D_vec,
              const double& process_error,
              const bool& disp_error,
              const double& plant_death_age_mean,
              const double& plant_death_age_sd,
              const arma::ivec& plant_death_ages,
              const std::vector<uint32>& repl_times,
              const sint32& repl_age,
              const arma::mat& log_morts,
              const double& extinct_N,
              const double& repl_threshold,
              pcg32& eng,
              arma::cube& N_bycage,
              arma::mat& logN_byplant,
              arma::mat& zero_byplant,
              const uint32& r) {

    uint32 n_plants = N_0.n_rows;
    uint32 n_lines = N_0.n_cols;
    uint32 max_t = N_bycage.n_rows - 1;

    std::normal_distribution<double> normal_rng(0.0, 1.0);
    std::poisson_distribution<arma::sword> poisson_rng(1.0);

    // Mean density dependence:
    double mean_alpha = arma::mean(A);

    // N at time t:
    arma::mat Nt = N_0;
    // N at time t+1:
    arma::mat Nt1 = Nt;

    // Matrix to keep track of extinctions:
    arma::Mat<unsigned short> extinct(n_plants, n_lines, arma::fill::zeros);
    // Vector to keep track of # days PAST death for plants
    // (these will start out negative)
    arma::ivec plant_days(n_plants);
    for (uint i = 0; i < n_plants; i++) {
        plant_days(i) = sample_plant_days(plant_death_age_mean, plant_death_age_sd,
                   plant_death_ages, i,
                   normal_rng, eng);
    }
// Index to keep track of position in `repl_times`
    uint32 repl_ind = 0;
    if (repl_times.front() == 0) repl_ind++;
    // Matrices keeping track of numbers of dispersed aphids:
    arma::mat emigrants(n_plants, n_lines);
    arma::mat immigrants(n_plants, n_lines);
    // Sum by row (i.e., by plant):
    arma::vec Z = arma::sum(Nt1, 1);

    for (uint32 t = 0; t < max_t; t++) {

        // Update # plant days PAST death for time t+1:
        plant_days++;
        // Generate numbers of dispersed aphids:
        if (disp_error) {
            // Disperal is Poisson process:
            update_dispersal(emigrants, immigrants, Nt, D_vec, poisson_rng, eng);
        } else {
            // Dispersal is entirely deterministic:
            update_dispersal(emigrants, immigrants, Nt, D_vec);
        }

        /*
         Go back through to simulate at time t+1
        */
        for (uint32 i = 0; i < n_plants; i++) {

            /*
             Density dependence for this plant (will be multiplied by each line's alpha):
             */
            double DD_i = arma::as_scalar(A * Nt.row(i).t()) / mean_alpha;

            for (uint32 j = 0; j < n_lines; j++) {

                // Dispersal of this aphid line to and from other plants.
                double imm = immigrants(i, j);
                double em = emigrants(i, j);

                // If it's extinct and no one's coming in, make 0 and skip the rest:
                if (extinct(i,j) == 1 && (imm - em) <= 0.0) {
                    Nt1(i,j) = 0;
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
                    Nt1(i,j) *= std::exp(R[j] * (1 - A[j] * DD_i) +
                        normal_rng(eng) * process_error);
                } else {
                    // (` - 1` below is to convert to C++ indices)
                    Nt1(i,j) *= std::exp(log_morts(plant_days[i] - 1, j) +
                        normal_rng(eng) * process_error);
                    imm *= std::exp(log_morts(plant_days[i] - 1, j));
                    if (imm < 1) imm = 0;
                }

                // Add dispersal:
                Nt1(i,j) += (imm - em);

                // Check for extinction:
                if (Nt1(i,j) < extinct_N) {
                    Nt1(i,j) = 0;
                    extinct(i, j) = 1;
                }

            }

        }


        // /*
        //  If it's a replacement time point, disperse all aphids from replaced
        //  plants to others
        // */
        // // To make sure this doesn't past bounds:
        // if (repl_ind < repl_times.size()) {
        //     // Now check:
        //     if ((t+1) == repl_times[repl_ind]) {
        //         replace_plants(repl_ind, plant_days, Nt1, extinct, normal_rng, eng,
        //                        repl_age, plant_death_age_mean, plant_death_age_sd);
        //     }
        // }

        // Summarize output:
        summarize_output(Nt1, N_bycage, logN_byplant, zero_byplant, r, t);

        // Update rows by plant:
        Z = arma::sum(Nt1, 1);

        for (uint32 i = 0; i < n_plants; i++) {
            if (Z(i) > repl_threshold || plant_days(i) >= 0) {
                Z(i) = 0;
                Nt1.row(i).zeros();
                plant_days(i) = sample_plant_days(plant_death_age_mean,
                           plant_death_age_sd,
                           plant_death_ages, i,
                           normal_rng, eng);
            }

        }

        // Iterate for time t next round
        // This also means I don't have to set Nt1 = Nt at the beginning of each round
        Nt = Nt1;

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
List sim_cages_(const uint32& n_cages,
                const arma::mat& N_0,
                const uint32& max_t,
                const arma::rowvec& R,
                const arma::rowvec& A,
                const arma::vec& D_vec,
                const double& process_error,
                const bool& disp_error,
                const arma::rowvec& plant_mort_0,
                const arma::rowvec& plant_mort_1,
                const double& plant_death_age_mean,
                const double& plant_death_age_sd,
                const arma::ivec& plant_death_ages,
                const std::vector<uint32>& repl_times,
                const sint32& repl_age,
                const double& extinct_N,
                const double& repl_threshold,
                const uint32& n_cores,
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
    if (D_vec.n_elem != n_lines) {
        err_msg += "ERROR: D_vec.n_elem != n_lines\n";
    }
    if (plant_mort_0.n_elem != n_lines) {
        err_msg += "ERROR: plant_mort_0.n_elem != n_lines\n";
    }
    if (plant_mort_1.n_elem != n_lines) {
        err_msg += "ERROR: plant_mort_1.n_elem != n_lines\n";
    }
    if (plant_death_ages.n_elem != 0 && plant_death_ages.n_elem != n_plants) {
        err_msg += "ERROR: plant_death_ages.n_elem != 0 && ";
        err_msg += "plant_death_ages.n_elem != n_plants\n";
    }
    if (arma::any(plant_death_ages <= 0)) {
        err_msg += "ERROR: arma::any(plant_death_ages <= 0)\n";
    }
    if (repl_times.size() == 0) {
        err_msg += "ERROR: repl_times.size() == 0\n";
    }

    if (err_msg.size() > 0) throw(Rcpp::exception(err_msg.c_str(), false));



    Progress p(n_cages, show_progress);

    /*
     --------------
     Create output objects
     --------------
     */
    // Total N for each line across all plants:
    arma::cube N_bycage(max_t+1, n_lines, n_cages);
    // Mean log(N) for each plant across all lines and all plants:
    arma::mat logN_byplant(n_plants, n_cages, arma::fill::zeros);
    // Mean proportion of zeros for each plant across all lines and all plants:
    arma::mat zero_byplant(n_plants, n_cages, arma::fill::zeros);

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
        N_bycage(arma::span(0), arma::span(), arma::span(r)) = arma::sum(N_0, 0);
        sim_cage(N_0, R, A, D_vec, process_error, disp_error,
                 plant_death_age_mean, plant_death_age_sd, plant_death_ages,
                 repl_times, repl_age_, log_morts, extinct_N, repl_threshold,
                 eng,
                 N_bycage, logN_byplant, zero_byplant, r);
        p.increment();
    }

#ifdef _OPENMP
}
#endif

    logN_byplant /= static_cast<double>(max_t);
    zero_byplant /= static_cast<double>(max_t);

    List out = List::create(_["N"] = wrap(N_bycage),
                            _["X"] = wrap(logN_byplant),
                            _["Z"] = wrap(zero_byplant));

    return out;

}
