
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
void calc_rep_rows(uint32& n_rows,
                   uint32& n_rows_wasps,
                   const uint32& max_t,
                   const uint32& save_every,
                   const uint32& n_lines,
                   const uint32& n_patches) {

    // # time points you'll save:
    n_rows_wasps = (max_t / save_every) + 1;
    if (max_t % save_every > 0) n_rows_wasps++;

    n_rows = n_patches * n_lines * n_rows_wasps;
    n_rows *= 2;  // `*2` for separate alate vs apterous
    n_rows += n_patches * n_rows_wasps;  // for mummies

    return;
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
    std::vector<uint32> wasp_rep;
    std::vector<uint32> wasp_time;
    std::vector<double> wasp_N;

    RepSummary()
        : rep(), time(), patch(), line(), type(), N(),
          wasp_rep(), wasp_time(), wasp_N(), r() {};

    void reserve(const uint32& rep_,
                 const uint32& max_t,
                 const uint32& save_every,
                 const uint32& n_lines,
                 const uint32& n_patches) {
        uint32 n_rows, n_rows_wasps;
        calc_rep_rows(n_rows, n_rows_wasps, max_t, save_every,
                      n_lines, n_patches);
        rep.reserve(n_rows);
        time.reserve(n_rows);
        patch.reserve(n_rows);
        line.reserve(n_rows);
        type.reserve(n_rows);
        N.reserve(n_rows);
        wasp_rep.reserve(n_rows_wasps);
        wasp_time.reserve(n_rows_wasps);
        wasp_N.reserve(n_rows_wasps);
        r = rep_;
    }

    // This version used when assimilating all reps into the first one
    void reserve(const uint32& n, const uint32& nw) {
        rep.reserve(n);
        time.reserve(n);
        patch.reserve(n);
        line.reserve(n);
        type.reserve(n);
        N.reserve(n);
        wasp_rep.reserve(nw);
        wasp_time.reserve(nw);
        wasp_N.reserve(nw); // `nw` is for the wasp vectors
        return;
    }

    void push_back(const uint32& t,
                   const AllPatches& patches) {

        // Everything but wasps:
        for (uint32 j = 0; j < patches.size(); j++) {
            const OnePatch& patch(patches[j]);
            for (uint32 i = 0; i < patch.size(); i++) {
                const AphidPop& aphid(patch[i]);
                append_living_aphids__(t, j, aphid.aphid_name,
                                       aphid.alates.total_aphids(),
                                       aphid.apterous.total_aphids());
            }
            append_mummies__(t, j, patch.total_mummies());
        }

        wasp_rep.push_back(r);
        wasp_time.push_back(t);
        wasp_N.push_back(patches.wasps.Y);

        return;
    }


    void clear() {

        time.clear();
        patch.clear();
        line.clear();
        type.clear();
        N.clear();
        wasp_rep.clear();
        wasp_time.clear();
        wasp_N.clear();

        // to clear memory:
        time.shrink_to_fit();
        patch.shrink_to_fit();
        line.shrink_to_fit();
        type.shrink_to_fit();
        N.shrink_to_fit();
        wasp_rep.shrink_to_fit();
        wasp_time.shrink_to_fit();
        wasp_N.shrink_to_fit();
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

        for (uint32 i = 0; i < other.wasp_time.size(); i++) {

            wasp_rep.push_back(other.wasp_rep[i]);
            wasp_time.push_back(other.wasp_time[i]);
            wasp_N.push_back(other.wasp_N[i]);

        }

        other.clear();

        return;

    }


private:

    uint32 r;

    inline void append_living_aphids__(const uint32& t,
                                       const uint32& p,
                                       const std::string& l,
                                       const double& N_ala,
                                       const double& N_apt) {

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

    inline void append_mummies__(const uint32& t,
                                 const uint32& p,
                                 const double& N_mum) {
        rep.push_back(r);
        time.push_back(t);
        patch.push_back(p);
        line.push_back("");
        type.push_back("mummy");
        N.push_back(N_mum);
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
                     const double& K_y_mult,
                     const double& death_prop,
                     const double& shape1_death_mort,
                     const double& shape2_death_mort,
                     const arma::mat& attack_surv,
                     const bool& disp_error,
                     const bool& demog_error,
                     const double& sigma_x,
                     const double& sigma_y,
                     const double& rho,
                     const double& extinct_N,
                     const std::vector<std::string>& aphid_name,
                     const std::vector<arma::cube>& leslie_mat,
                     const std::vector<arma::cube>& aphid_density_0,
                     const std::vector<double>& alate_b0,
                     const std::vector<double>& alate_b1,
                     const std::vector<double>& disp_rate,
                     const std::vector<double>& disp_mort,
                     const std::vector<uint32>& disp_start,
                     const std::vector<uint32>& living_days,
                     const std::vector<double>& pred_rate,
                     const arma::mat& mum_density_0,
                     const arma::vec& rel_attack,
                     const double& a,
                     const double& k,
                     const double& h,
                     const double& wasp_density_0,
                     const double& sex_ratio,
                     const double& s_y,
                     Progress& prog_bar,
                     int& status_code,
                     pcg32& eng) {

    double demog_mult = 1;
    if (!demog_error) demog_mult = 0;
    // either type of environmental error
    bool process_error = (sigma_x > 0) || (sigma_y > 0);

    uint32 n_lines = aphid_name.size();
    uint32 n_patches = aphid_density_0.size();

    RepSummary summary;

    summary.reserve(rep, max_t, save_every, n_lines, n_patches);

    uint32 iters = 0;

    AllPatches patches(sigma_x, sigma_y, rho, demog_mult, mean_K, sd_K, K_y_mult,
                       death_prop,
                       shape1_death_mort, shape2_death_mort, attack_surv,
                       aphid_name, leslie_mat, aphid_density_0, alate_b0, alate_b1,
                       disp_rate, disp_mort, disp_start, living_days, pred_rate,
                       extinct_N, mum_density_0, rel_attack, a, k, h,
                       wasp_density_0, sex_ratio, s_y, eng);

    summary.push_back(0, patches);

    for (uint32 t = 1; t <= max_t; t++) {

        if (interrupt_check(iters, prog_bar)) {
            status_code = -1;
            return summary;
        }

        if (disp_error) {
            patches.calc_dispersal(eng);
        } else patches.calc_dispersal();

        if (process_error) {
            patches.update(eng);
        } else patches.update();

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

/*
 This template checks for any negative values, and returns an error if
 it finds any.
 It works on std::vector<>, arma::Col<>, arma::Row<>, arma::Mat<>,
 and arma::Cube<>.
 */
//
template <typename T>
void negative_check(const T& input, const std::string& name) {
    for (uint32 i = 0; i < input.size(); i++) {
        if (input[i] < 0) {
            std::string msg = "\nERROR: " + name + " contains negative items\n";
            stop(msg.c_str());
        }
    }
    return;
}
// Same but for one double
void one_negative_check(const double& input, const std::string& name) {
    if (input < 0) {
        std::string msg = "\nERROR: " + name + " is negative\n";
        stop(msg.c_str());
    }
    return;
}


// Same but checks that the numbers are in range [0,1]
template <typename T>
void non_prop_check(const T& input, const std::string& name) {
    for (uint32 i = 0; i < input.size(); i++) {
        if (input[i] < 0 || input[i] > 1) {
            std::string msg = "\nERROR: " + name + " contains items outside " +
                "the range [0,1]\n";
            stop(msg.c_str());
        }
    }
    return;
}
// Same but for one double
void one_non_prop_check(const double& input, const std::string& name) {
    if (input < 0 || input > 1) {
        std::string msg = "\nERROR: " + name + " is outside the range [0,1]\n";
        stop(msg.c_str());
    }
    return;
}


// Same but checks for numbers <= 0
template <typename T>
void positive_check(const T& input, const std::string& name) {
    for (uint32 i = 0; i < input.size(); i++) {
        if (input[i] <= 0) {
            std::string msg = "\nERROR: " + name + " contains items <= 0\n";
            stop(msg.c_str());
        }
    }
    return;
}
// Same but for unsigned integer
void one_positive_check(const uint32& input, const std::string& name) {
    if (input == 0) {
        std::string msg = "\nERROR: " + name + " is zero\n";
        stop(msg.c_str());
    }
    return;
}





void check_args(const uint32& n_reps,
                const uint32& n_lines,
                const uint32& n_patches,
                const uint32& max_plant_age,
                const double& max_N,
                const std::deque<uint32>& check_for_clear,
                const uint32& max_t,
                const uint32& save_every,
                const double& mean_K,
                const double& sd_K,
                const double& K_y_mult,
                const double& death_prop,
                const double& shape1_death_mort,
                const double& shape2_death_mort,
                const arma::mat& attack_surv,
                const bool& disp_error,
                const bool& demog_error,
                const double& sigma_x,
                const double& sigma_y,
                const double& rho,
                const double& extinct_N,
                const std::vector<std::string>& aphid_name,
                const std::vector<arma::cube>& leslie_mat,
                const std::vector<arma::cube>& aphid_density_0,
                const std::vector<double>& alate_b0,
                const std::vector<double>& alate_b1,
                const std::vector<double>& disp_rate,
                const std::vector<double>& disp_mort,
                const std::vector<uint32>& disp_start,
                const std::vector<uint32>& living_days,
                const std::vector<double>& pred_rate,
                const arma::mat& mum_density_0,
                const arma::vec& rel_attack,
                const double& a,
                const double& k,
                const double& h,
                const double& wasp_density_0,
                const double& sex_ratio,
                const double& s_y,
                uint32& n_threads) {


    /*
     ===============================================================
     Checking dimensions
     ===============================================================
     */
    uint32 n_stages = aphid_density_0.front().n_rows;

    if (n_reps == 0) stop("\nERROR: n_reps == 0\n");
    if (n_lines == 0) stop("\nERROR: n_lines == 0\n");
    if (n_patches == 0) stop("\nERROR: n_patches == 0\n");
    if (n_stages == 0) stop("\nERROR: n_stages == 0\n");

    if (leslie_mat.size() != n_lines) {
        stop("\nERROR: leslie_mat.size() != n_lines\n");
    }
    if (alate_b0.size() != n_lines) {
        stop("\nERROR: alate_b0.size() != n_lines\n");
    }
    if (alate_b1.size() != n_lines) {
        stop("\nERROR: alate_b1.size() != n_lines\n");
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
    if (living_days.size() != n_lines) {
        stop("\nERROR: living_days.size() != n_lines\n");
    }
    if (attack_surv.n_cols != n_lines || attack_surv.n_rows != 2) {
        stop(std::string("\nERROR: attack_surv.n_cols != n_lines || ") +
            std::string("attack_surv.n_rows != 2\n"));
    }
    if (leslie_mat.size() != n_lines) {
        stop("\nERROR: leslie_mat.size() != n_lines\n");
    }



    // In `aphid_density_0`, rows are aphid stages, columns are types (alate vs
    // apterous), and slices are aphid lines.
    if (aphid_density_0.size() != n_patches) {
        stop("\nERROR: aphid_density_0.size() != n_patches\n");
    }
    for (uint32 i = 0; i < n_patches; i++) {
        if (aphid_density_0[i].n_rows != n_stages) {
            stop("\nERROR: aphid_density_0[i].n_rows != n_stages\n");
        }
        if (aphid_density_0[i].n_cols != 2) {
            stop("\nERROR: aphid_density_0[i].n_cols != 2\n");
        }
        if (aphid_density_0[i].n_slices != n_lines) {
            stop("\nERROR: aphid_density_0[i].n_slices != n_lines\n");
        }
    }


    // In `leslie_mat`, items in vector are aphid lines, slices are
    // alate/apterous/parasitized.
    if (leslie_mat.size() != n_lines) {
        stop("\nERROR: leslie_mat.size() != n_lines\n");
    }
    for (uint32 i = 0; i < n_lines; i++) {
        if (leslie_mat[i].n_rows != n_stages) {
            stop("\nERROR: leslie_mat[i].n_rows != n_stages\n");
        }
        if (leslie_mat[i].n_cols != n_stages) {
            stop("\nERROR: leslie_mat[i].n_cols != n_stages\n");
        }
        if (leslie_mat[i].n_slices != 3) {
            stop("\nERROR: leslie_mat[i].n_slices != 3\n");
        }
    }


    if (rel_attack.n_elem != n_stages) {
        stop("\nERROR: rel_attack.n_elem != n_stages\n");
    }

    if (pred_rate.size() != n_patches) {
        stop("\nERROR: pred_rate.size() != n_patches\n");
    }
    if (mum_density_0.n_cols != n_patches) {
        stop("\nERROR: mum_density_0.n_cols != n_patches\n");
    }
    if (mum_density_0.n_rows == 0) {
        stop("\nERROR: mum_density_0.n_rows == 0\n");
    }



    /*
     ===============================================================
     Checking values
     ===============================================================
     */

    // Check that # threads isn't too high and change to 1 if not using OpenMP:
    thread_check(n_threads);

    // doubles that must be >= 0
    one_negative_check(max_N, "max_N");
    one_negative_check(mean_K, "mean_K");
    one_negative_check(sd_K, "sd_K");
    one_negative_check(K_y_mult, "K_y_mult");
    one_negative_check(shape1_death_mort, "shape1_death_mort");
    one_negative_check(shape2_death_mort, "shape2_death_mort");
    one_negative_check(sigma_x, "sigma_x");
    one_negative_check(sigma_y, "sigma_y");
    one_negative_check(rho, "rho");
    one_negative_check(extinct_N, "extinct_N");
    one_negative_check(a, "a");
    one_negative_check(k, "k");
    one_negative_check(h, "h");
    one_negative_check(wasp_density_0, "wasp_density_0");

    // objects containing doubles that must be >= 0
    negative_check<arma::mat>(attack_surv, "attack_surv");
    for (uint32 i = 0; i < leslie_mat.size(); i++) {
        negative_check<arma::cube>(leslie_mat[i], "leslie_mat");
    }
    for (uint32 i = 0; i < aphid_density_0.size(); i++) {
        negative_check<arma::cube>(aphid_density_0[i], "aphid_density_0");
    }
    negative_check<std::vector<double>>(disp_rate, "disp_rate");
    negative_check<std::vector<double>>(disp_mort, "disp_mort");
    negative_check<std::vector<double>>(pred_rate, "pred_rate");
    negative_check<arma::mat>(mum_density_0, "mum_density_0");
    negative_check<arma::vec>(rel_attack, "rel_attack");

    // doubles that must be >= 0 and <= 1
    one_non_prop_check(death_prop, "death_prop");
    one_non_prop_check(sex_ratio, "sex_ratio");
    one_non_prop_check(s_y, "s_y");

    // below top rows should be >= 0 and <= 1:
    auto iter = leslie_mat.begin();
    arma::cube tmp((*iter).n_rows - 1, (*iter).n_cols, (*iter).n_slices);
    for (; iter != leslie_mat.end(); iter++) {
        tmp = (*iter)(arma::span(1, iter->n_rows - 1),
               arma::span::all, arma::span::all);
        non_prop_check<arma::cube>(tmp, "leslie_mat");
    }

    // integers that must be > 0
    one_positive_check(max_t, "max_t");
    positive_check<std::vector<uint32>>(living_days, "living_days");



    /*
     ===============================================================
     Other checks
     ===============================================================
     */

    if (max_plant_age > 0 && max_N > 0) {
        stop("\nERROR: max_plant_age > 0 && max_N > 0\n");
    }
    if (max_plant_age == 0 && max_N <= 0) {
        stop("\nERROR: max_plant_age == 0 && max_N <= 0\n");
    }


    return;
}





//[[Rcpp::export]]
List sim_clonewars_cpp(const uint32& n_reps,
                       const uint32& max_plant_age,
                       const double& max_N,
                       const std::deque<uint32>& check_for_clear,
                       const uint32& max_t,
                       const uint32& save_every,
                       const double& mean_K,
                       const double& sd_K,
                       const double& K_y_mult,
                       const double& death_prop,
                       const double& shape1_death_mort,
                       const double& shape2_death_mort,
                       const arma::mat& attack_surv,
                       const bool& disp_error,
                       const bool& demog_error,
                       const double& sigma_x,
                       const double& sigma_y,
                       const double& rho,
                       const double& extinct_N,
                       const std::vector<std::string>& aphid_name,
                       const std::vector<arma::cube>& leslie_mat,
                       const std::vector<arma::cube>& aphid_density_0,
                       const std::vector<double>& alate_b0,
                       const std::vector<double>& alate_b1,
                       const std::vector<double>& disp_rate,
                       const std::vector<double>& disp_mort,
                       const std::vector<uint32>& disp_start,
                       const std::vector<uint32>& living_days,
                       const std::vector<double>& pred_rate,
                       const arma::mat& mum_density_0,
                       const arma::vec& rel_attack,
                       const double& a,
                       const double& k,
                       const double& h,
                       const double& wasp_density_0,
                       const double& sex_ratio,
                       const double& s_y,
                       uint32 n_threads,
                       const bool& show_progress) {

    uint32 n_lines = aphid_name.size();
    uint32 n_patches = aphid_density_0.size();

    check_args(n_reps, n_lines, n_patches,
               max_plant_age, max_N, check_for_clear, max_t, save_every,
               mean_K, sd_K, K_y_mult, death_prop,
               shape1_death_mort, shape2_death_mort,
               attack_surv, disp_error, demog_error,
               sigma_x, sigma_y, rho, extinct_N, aphid_name,
               leslie_mat, aphid_density_0, alate_b0, alate_b1,
               disp_rate, disp_mort, disp_start, living_days,
               pred_rate, mum_density_0, rel_attack, a, k, h,
               wasp_density_0, sex_ratio, s_y, n_threads);


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
                                             save_every, mean_K, sd_K, K_y_mult,
                                             death_prop,
                                             shape1_death_mort, shape2_death_mort,
                                             attack_surv, disp_error,
                                             demog_error, sigma_x, sigma_y, rho, extinct_N,
                                             aphid_name, leslie_mat, aphid_density_0,
                                             alate_b0, alate_b1, disp_rate, disp_mort,
                                             disp_start, living_days, pred_rate,
                                             mum_density_0, rel_attack, a, k,
                                             h, wasp_density_0, sex_ratio, s_y,
                                             prog_bar, status_code, eng);
        } else {
            summaries[i] = one_rep__<double>(max_N, i, check_for_clear, max_t,
                                             save_every, mean_K, sd_K, K_y_mult,
                                             death_prop,
                                             shape1_death_mort, shape2_death_mort,
                                             attack_surv, disp_error,
                                             demog_error, sigma_x, sigma_y, rho, extinct_N,
                                             aphid_name, leslie_mat, aphid_density_0,
                                             alate_b0, alate_b1, disp_rate, disp_mort,
                                             disp_start, living_days, pred_rate,
                                             mum_density_0, rel_attack, a, k,
                                             h, wasp_density_0, sex_ratio, s_y,
                                             prog_bar, status_code, eng);
        }
    }

#ifdef _OPENMP
}
#endif


    /*
     When # reps > 1, combine all RepSummary objects into the first one.
     This is to make it easier to add them to the output data frame.
     */
    RepSummary& summ(summaries.front());
    if (n_reps > 1) {
        summ.reserve(summ.N.size() * summaries.size(),
                     summ.wasp_N.size() * summaries.size());
        for (uint32 i = 1; i < n_reps; i++) {
            summ.assimilate(summaries[i]);
        }
    }

    List out = List::create(_["aphids"] = DataFrame::create(
                                _["rep"] = summ.rep,
                                _["time"] = summ.time,
                                _["patch"] = summ.patch,
                                _["line"] = summ.line,
                                _["type"] = summ.type,
                                _["N"] = summ.N),
                            _["wasps"] = DataFrame::create(
                                _["rep"] = summ.wasp_rep,
                                _["time"] = summ.wasp_time,
                                _["wasps"] = summ.wasp_N));

    return out;
}
