
/*
 *****************************************************************************
 *****************************************************************************

 This file is for simulation functions that will be available from R.

 *****************************************************************************
 *****************************************************************************
 */

#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <deque>                // deque
#include <pcg/pcg_random.hpp>   // pcg prng
#include <progress.hpp>         // for the progress bar
#ifdef _OPENMP
#include <omp.h>                // OpenMP
#endif


#include "clonewars_types.hpp"  // integer types
#include "aphids.hpp"           // aphid classes
#include "patches.hpp"          // patch classes
#include "pcg.hpp"              // runif_ fxns


//' Check that the number of threads doesn't exceed the number available, and change
//' to 1 if OpenMP isn't enabled.
//'
//' @noRd
//'
inline void thread_check(uint32& n_threads) {

#ifdef _OPENMP
    if (n_threads == 0) n_threads = 1;
    if (n_threads > omp_get_max_threads()) {
        std::string max_threads = std::to_string(omp_get_max_threads());
        std::string err_msg = "\nThe number of requested threads (" +
            std::to_string(n_threads) + ") exceeds the max available on the system (" +
            max_threads + ").";
        stop(err_msg.c_str());
    }
#else
    n_threads = 1;
#endif

    return;
}

// For checking for user interrupts every N iterations:
inline bool interrupt_check(uint32& iters,
                            Progress& prog_bar,
                            const uint32& N = 100) {
    ++iters;
    if (iters > N) {
        if (prog_bar.is_aborted() || prog_bar.check_abort()) return true;
        prog_bar.increment(iters);
        iters = 0;
    }
    return false;
}





// Calculate the number of rows per rep.
uint32 calc_rep_rows(const uint32& max_t,
                     const uint32& save_every,
                     const uint32& n_lines,
                     const uint32& n_patches) {

    uint32 n_saves = (max_t / save_every) + 1;
    if (max_t % save_every > 0) n_saves++;

    uint32 nrows_ = n_patches * n_lines * n_saves;
    nrows_ *= 2;  // for separate alate vs apterous N

    return nrows_;
}


/*


 leslie_matrix__(instar_days, surv_juv, surv_adult, repro, leslie_);
 leslie_sad__(leslie_, X_0_);
 X_0_ *= aphid_density_0;

 OR

 leslie_matrix__(instar_days, surv_juv, surv_adult, repro, leslie_);
 X_0_ = aphid_density_0;
 */

/*
 To calculate dispersal start date
 arma::accu(instar_days.head(instar_days.n_elem - 1)) - 1
 */


struct RepSummary {

    std::vector<uint32> rep;
    std::vector<uint32> time;
    std::vector<uint32> patch;
    std::vector<std::string> line;
    std::vector<std::string> type;
    std::vector<double> N;

    RepSummary() : rep(), time(), patch(), line(), type(), N(), r() {};

    void reserve(const uint32& rep_,
                 const uint32& max_t,
                 const uint32& save_every,
                 const uint32& n_lines,
                 const uint32& n_patches) {
        uint32 n_rows = calc_rep_rows(max_t, save_every, n_lines, n_patches);
        rep.reserve(n_rows);
        time.reserve(n_rows);
        patch.reserve(n_rows);
        line.reserve(n_rows);
        type.reserve(n_rows);
        N.reserve(n_rows);
        r = rep_;
    }

    // This version used when assimilating all reps into the first one
    void reserve(const uint32& n) {
        rep.reserve(n);
        time.reserve(n);
        patch.reserve(n);
        line.reserve(n);
        type.reserve(n);
        N.reserve(n);
        return;
    }

    void push_back(const uint32& t,
                   const AllPatches& patches) {

        for (uint32 j = 0; j < patches.size(); j++) {
            const OnePatch& patch(patches[j]);
            for (uint32 i = 0; i < patch.size(); i++) {
                const AphidPop& aphid(patch[i]);
                append_one__(t, j, aphid.aphid_name, aphid.alates.total_aphids(),
                             aphid.apterous.total_aphids());
            }
        }

        return;
    }


    void clear() {

        time.clear();
        patch.clear();
        line.clear();
        type.clear();
        N.clear();

        // to clear memory:
        time.shrink_to_fit();
        patch.shrink_to_fit();
        line.shrink_to_fit();
        type.shrink_to_fit();
        N.shrink_to_fit();
    }


    void assimilate(RepSummary& other) {

        for (uint32 i = 0; i < other.time.size(); i++) {

            rep.push_back(other.rep[i]);
            time.push_back(other.time[i]);
            patch.push_back(other.patch[i]);
            line.push_back(other.line[i]);
            type.push_back(other.type[i]);
            N.push_back(other.N[i]);

        }

        other.clear();

        return;

    }


private:

    uint32 r;

    inline void append_one__(const uint32& t, const uint32& p, const std::string& l,
                             const double& N_ala, const double& N_apt) {

        rep.push_back(r);
        rep.push_back(r);

        time.push_back(t);
        time.push_back(t);

        patch.push_back(p);
        patch.push_back(p);

        line.push_back(l);
        line.push_back(l);

        type.push_back("alate");
        type.push_back("apterous");

        N.push_back(N_ala);
        N.push_back(N_apt);

        return;
    }
};



/*
 `T` should be uint32 for using age to determine whether a patch gets cleared.
 It should be double for using aphid abundance to determine whether a patch gets cleared.
 */

template <typename T>
RepSummary one_rep__(const T& clear_threshold,
                     const uint32& rep,
                     std::deque<uint32> check_for_clear,
                     const uint32& max_t,
                     const uint32& save_every,
                     const double& mean_K,
                     const double& sd_K,
                     const double& meanlog_death_age,
                     const double& sdlog_death_age,
                     const double& shape1_death_mort,
                     const double& shape2_death_mort,
                     const bool& disp_error,
                     const bool& demog_error,
                     const double& sigma_x,
                     const double& rho,
                     const double& extinct_N,
                     const std::vector<std::string>& aphid_name,
                     const std::vector<arma::cube>& leslie_mat,
                     const std::vector<arma::cube>& aphid_density_0,
                     const std::vector<double>& alate_prop,
                     const std::vector<double>& disp_rate,
                     const std::vector<double>& disp_mort,
                     const std::vector<uint32>& disp_start,
                     const std::vector<double>& pred_rate,
                     Progress& prog_bar,
                     int& status_code,
                     pcg32& eng) {

    double demog_mult = 1;
    if (!demog_error) demog_mult = 0;
    // either demographic or environmental error
    bool process_error = demog_error || (sigma_x > 0);

    uint32 n_lines = aphid_name.size();
    uint32 n_patches = aphid_density_0.size();

    RepSummary summary;

    summary.reserve(rep, max_t, save_every, n_lines, n_patches);

    uint32 iters = 0;

    AllPatches patches(sigma_x, rho, demog_mult, mean_K, sd_K,
                       meanlog_death_age, sdlog_death_age,
                       shape1_death_mort, shape2_death_mort,
                       aphid_name, leslie_mat, aphid_density_0, alate_prop,
                       disp_rate, disp_mort, disp_start, pred_rate, extinct_N, eng);

    summary.push_back(0, patches);

    for (uint32 t = 1; t <= max_t; t++) {

        if (interrupt_check(iters, prog_bar)) {
            status_code = -1;
            return summary;
        }

        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        double N = -1;
        for (uint32 i = 0; i < patches.patches.size(); i++) {
            const OnePatch& p(patches.patches[i]);
            if (N != -1 && p.total_aphids() != N) {
                Rcout << "p.this_j = " << p.this_j << std::endl;
                Rcout << "t = " << t << std::endl;
                stop("\n\nWHY, GOD, WHY?? (#0) \n\n");
            }
            N = p.total_aphids();
        }
        // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


        if (disp_error) {
            patches.calc_dispersal(eng);
        } else patches.calc_dispersal();





        // // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        //
        // patches.emigrants.fill(0);
        // patches.immigrants.fill(0);
        //
        // // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


        // emigrants = arma::zeros<arma::cube>(n_stages, n_patches, n_lines);
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        arma::vec e_sum(n_patches);
        arma::vec i_sum(n_patches);
        for (uint32 i = 0; i < n_patches; i++) {
            e_sum(i) = arma::accu(patches.emigrants.col(i));
            i_sum(i) = arma::accu(patches.immigrants.col(i));
        }
        arma::vec e_unqs = arma::unique(e_sum);
        arma::vec i_unqs = arma::unique(i_sum);
        if (e_unqs.n_elem > 1 || i_unqs.n_elem > 1) {
            Rcout << std::endl << "Emigrants:" << std::endl;
            e_sum.print();
            Rcout << std::endl;
            Rcout << std::endl << "Immigrants:" << std::endl;
            i_sum.print();
            Rcout << std::endl;
            stop("\n\nDuplicate E or I.\n\n");
        }
        // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


        if (process_error) {
            patches.update_pops(eng);
        } else patches.update_pops();

        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        N = -1;
        for (uint32 i = 0; i < patches.patches.size(); i++) {
            const OnePatch& p(patches.patches[i]);
            if (N != -1 && p.total_aphids() != N) {
                Rcout << "p.this_j = " << p.this_j << std::endl;
                Rcout << "t = " << t << std::endl;
                stop("\n\nWHY, GOD, WHY??\n\n");
            }
            N = p.total_aphids();
        }
        // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


        if (t % save_every == 0 || t == max_t) summary.push_back(t, patches);

        // If all patches are empty, then stop this rep.
        // It's important to do this before clearing patches.
        bool all_empty = true;
        for (const OnePatch& p : patches.patches) {
            if (!p.empty) {
                all_empty = false;
                break;
            }
        }
        if (all_empty) break;


        if (!check_for_clear.empty() && t == check_for_clear.front()) {
            check_for_clear.pop_front();
            patches.clear_patches(clear_threshold, eng);
        }
    }


    return summary;

}





//[[Rcpp::export]]
DataFrame sim_clonewars_cpp(const uint32& n_reps,
                            const uint32& max_plant_age,
                            const double& max_N,
                            const std::deque<uint32>& check_for_clear,
                            const uint32& max_t,
                            const uint32& save_every,
                            const double& mean_K,
                            const double& sd_K,
                            const double& meanlog_death_age,
                            const double& sdlog_death_age,
                            const double& shape1_death_mort,
                            const double& shape2_death_mort,
                            const bool& disp_error,
                            const bool& demog_error,
                            const double& sigma_x,
                            const double& rho,
                            const double& extinct_N,
                            const std::vector<std::string>& aphid_name,
                            const std::vector<arma::cube>& leslie_mat,
                            const std::vector<arma::cube>& aphid_density_0,
                            const std::vector<double>& alate_prop,
                            const std::vector<double>& disp_rate,
                            const std::vector<double>& disp_mort,
                            const std::vector<uint32>& disp_start,
                            const std::vector<double>& pred_rate,
                            uint32 n_threads,
                            const bool& show_progress) {

    // Check that # threads isn't too high and change to 1 if not using OpenMP:
    thread_check(n_threads);

    if (n_reps == 0) {
        stop("\nERROR: n_reps == 0\n");
    }


    uint32 n_lines = aphid_name.size();
    if (leslie_mat.size() != n_lines) {
        stop("\nERROR: leslie_mat.size() != n_lines\n");
    }
    if (alate_prop.size() != n_lines) {
        stop("\nERROR: alate_prop.size() != n_lines\n");
    }
    if (disp_rate.size() != n_lines) {
        stop("\nERROR: disp_rate.size() != n_lines\n");
    }
    if (disp_mort.size() != n_lines) {
        stop("\nERROR: disp_mort.size() != n_lines\n");
    }
    if (disp_start.size() != n_lines) {
        stop("\nERROR: disp_start.size() != n_lines\n");
    }

    uint32 n_patches = aphid_density_0.size();
    if (pred_rate.size() != n_patches) {
        stop("\nERROR: pred_rate.size() != n_patches\n");
    }

    if (max_plant_age > 0 && max_N > 0) {
        stop("\nERROR: max_plant_age > 0 && max_N > 0\n");
    }
    if (max_plant_age == 0 && max_N <= 0) {
        stop("\nERROR: max_plant_age == 0 && max_N <= 0\n");
    }


    Progress prog_bar(max_t * n_reps, show_progress);
    std::vector<int> status_codes(n_threads, 0);

    // Generate seeds for random number generators (1 set of seeds per rep)
    const std::vector<std::vector<uint64>> seeds = mt_seeds(n_reps);

    std::vector<RepSummary> summaries(n_reps);


#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(n_threads) if (n_threads > 1)
{
#endif

    // Write the active seed per thread or just write one of the seeds.
#ifdef _OPENMP
    uint32 active_thread = omp_get_thread_num();
#else
    uint32 active_thread = 0;
#endif
    int& status_code(status_codes[active_thread]);

    pcg32 eng;

    // Parallelize the Loop
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (uint32 i = 0; i < n_reps; i++) {
        if (status_code != 0) continue;
        seed_pcg(eng, seeds[i]);
        if (max_plant_age > 0) {
            summaries[i] = one_rep__<uint32>(max_plant_age, i, check_for_clear, max_t,
                                             save_every, mean_K, sd_K,
                                             meanlog_death_age, sdlog_death_age,
                                             shape1_death_mort, shape2_death_mort,
                                             disp_error,
                                             demog_error, sigma_x, rho, extinct_N,
                                             aphid_name, leslie_mat, aphid_density_0,
                                             alate_prop, disp_rate, disp_mort,
                                             disp_start, pred_rate, prog_bar,
                                             status_code, eng);
        } else {
            summaries[i] = one_rep__<double>(max_N, i, check_for_clear, max_t,
                                             save_every, mean_K, sd_K,
                                             meanlog_death_age, sdlog_death_age,
                                             shape1_death_mort, shape2_death_mort,
                                             disp_error,
                                             demog_error, sigma_x, rho, extinct_N,
                                             aphid_name, leslie_mat, aphid_density_0,
                                             alate_prop, disp_rate, disp_mort,
                                             disp_start, pred_rate, prog_bar,
                                             status_code, eng);
        }
    }

#ifdef _OPENMP
}
#endif


    /*
     When # reps > 1, combine all SimSummary objects into the first one.
     This is to make it easier to add them to the output data frame.
     */
    if (n_reps > 1) {
        summaries[0].reserve(summaries[0].N.size() * summaries.size());
        for (uint32 i = 1; i < n_reps; i++) {
            summaries[0].assimilate(summaries[i]);
        }
    }


    DataFrame out = DataFrame::create(_["rep"] = summaries[0].rep,
                                      _["time"] = summaries[0].time,
                                      _["patch"] = summaries[0].patch,
                                      _["line"] = summaries[0].line,
                                      _["type"] = summaries[0].type,
                                      _["N"] = summaries[0].N);

    return out;
}
