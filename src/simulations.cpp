
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
#include "patches.hpp"          // field and plant classes
#include "pcg.hpp"              // mt_seeds seed_pcg fxns



/*
 Allows you to verify that you're able to use multiple threads.
 */
//[[Rcpp::export]]
bool using_openmp() {
    bool out = false;
#ifdef _OPENMP
    out = true;
#endif
    return out;
}




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
                   const uint32& n_fields,
                   const uint32& n_plants) {

    // # time points you'll save:
    n_rows_wasps = (max_t / save_every) + 1;
    if (max_t % save_every > 0) n_rows_wasps++;
    n_rows_wasps *= n_fields;

    n_rows = n_plants * n_lines * n_rows_wasps;
    n_rows *= 3;  // `*3` for separate alate vs apterous vs parasitized
    n_rows += n_plants * n_rows_wasps;  // for mummies

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
    std::vector<uint32> field;
    std::vector<uint32> plant;
    std::vector<std::string> line;
    std::vector<std::string> type;
    std::vector<double> N;
    std::vector<uint32> wasp_rep;
    std::vector<uint32> wasp_time;
    std::vector<uint32> wasp_field;
    std::vector<double> wasp_N;

    RepSummary()
        : rep(), time(), field(), plant(), line(), type(), N(),
          wasp_rep(), wasp_time(), wasp_field(), wasp_N(), r() {};

    void reserve(const uint32& rep_,
                 const uint32& max_t,
                 const uint32& save_every,
                 const uint32& n_lines,
                 const uint32& n_fields,
                 const uint32& n_plants) {
        uint32 n_rows, n_rows_wasps;
        calc_rep_rows(n_rows, n_rows_wasps, max_t, save_every,
                      n_lines, n_fields, n_plants);
        rep.reserve(n_rows);
        time.reserve(n_rows);
        field.reserve(n_rows);
        plant.reserve(n_rows);
        line.reserve(n_rows);
        type.reserve(n_rows);
        N.reserve(n_rows);
        wasp_rep.reserve(n_rows_wasps);
        wasp_time.reserve(n_rows_wasps);
        wasp_field.reserve(n_rows_wasps);
        wasp_N.reserve(n_rows_wasps);
        r = rep_;
    }

    // This version used when assimilating all reps into the first one
    void reserve(const uint32& n, const uint32& nw) {
        rep.reserve(n);
        time.reserve(n);
        field.reserve(n);
        plant.reserve(n);
        line.reserve(n);
        type.reserve(n);
        N.reserve(n);
        wasp_rep.reserve(nw);
        wasp_time.reserve(nw);
        wasp_field.reserve(nw);
        wasp_N.reserve(nw); // `nw` is for the wasp vectors
        return;
    }


    void push_back(const uint32& t,
                   const std::vector<OneField>& fields) {

        for (uint32 k = 0; k < fields.size(); k++) {

            const OneField& field(fields[k]);

            // Everything but wasps:
            for (uint32 j = 0; j < field.size(); j++) {
                const OnePlant& plant(field[j]);
                for (uint32 i = 0; i < plant.size(); i++) {
                    const AphidPop& aphid(plant[i]);
                    append_living_aphids__(t, k, j, aphid.aphid_name,
                                           aphid.alates.total_aphids(),
                                           aphid.apterous.total_aphids(),
                                           aphid.paras.total_aphids());
                }
                append_mummies__(t, k, j, plant.total_mummies());
            }

            wasp_rep.push_back(r);
            wasp_time.push_back(t);
            wasp_field.push_back(k);
            wasp_N.push_back(field.wasps.Y);
        }

        return;
    }


    void clear() {

        time.clear();
        field.clear();
        plant.clear();
        line.clear();
        type.clear();
        N.clear();
        wasp_rep.clear();
        wasp_time.clear();
        wasp_field.clear();
        wasp_N.clear();

        // to clear memory:
        time.shrink_to_fit();
        field.shrink_to_fit();
        plant.shrink_to_fit();
        line.shrink_to_fit();
        type.shrink_to_fit();
        N.shrink_to_fit();
        wasp_rep.shrink_to_fit();
        wasp_time.shrink_to_fit();
        wasp_field.shrink_to_fit();
        wasp_N.shrink_to_fit();
    }


    void assimilate(RepSummary& other) {

        for (uint32 i = 0; i < other.time.size(); i++) {

            rep.push_back(other.rep[i]);
            time.push_back(other.time[i]);
            field.push_back(other.field[i]);
            plant.push_back(other.plant[i]);
            line.push_back(other.line[i]);
            type.push_back(other.type[i]);
            N.push_back(other.N[i]);

        }

        for (uint32 i = 0; i < other.wasp_time.size(); i++) {

            wasp_rep.push_back(other.wasp_rep[i]);
            wasp_time.push_back(other.wasp_time[i]);
            wasp_field.push_back(other.wasp_field[i]);
            wasp_N.push_back(other.wasp_N[i]);

        }

        other.clear();

        return;

    }


private:

    uint32 r;

    inline void append_living_aphids__(const uint32& t,
                                       const uint32& c,
                                       const uint32& p,
                                       const std::string& l,
                                       const double& N_ala,
                                       const double& N_apt,
                                       const double& N_par) {

        rep.push_back(r);
        rep.push_back(r);
        rep.push_back(r);

        time.push_back(t);
        time.push_back(t);
        time.push_back(t);

        field.push_back(c);
        field.push_back(c);
        field.push_back(c);

        plant.push_back(p);
        plant.push_back(p);
        plant.push_back(p);

        line.push_back(l);
        line.push_back(l);
        line.push_back(l);

        type.push_back("alate");
        type.push_back("apterous");
        type.push_back("parasitized");

        N.push_back(N_ala);
        N.push_back(N_apt);
        N.push_back(N_par);

        return;
    }

    inline void append_mummies__(const uint32& t,
                                 const uint32& c,
                                 const uint32& p,
                                 const double& N_mum) {
        rep.push_back(r);
        time.push_back(t);
        field.push_back(c);
        plant.push_back(p);
        line.push_back("");
        type.push_back("mummy");
        N.push_back(N_mum);
        return;
    }
};




inline void do_perturb(std::deque<PerturbInfo>& perturbs,
                       std::vector<OneField>& fields,
                       const uint32& t,
                       const uint32& n_lines,
                       const double& extinct_N,
                       const bool& do_pop) {

    if (perturbs.empty()) return;

    uint32 i = 0;
    while (i < perturbs.size()) {
        const PerturbInfo& pert(perturbs[i]);
        if (pert.time != t) break;
        double mult = pert.multiplier;
        OneField& field(fields[pert.field]);
        if (pert.index < n_lines) { // aphids (non-parasitized)
            for (OnePlant& p : field.plants) {
                AphidPop& aphids(p.aphids[pert.index]);
                aphids.apterous.clear(mult);
                aphids.alates.clear(mult);
                if (aphids.total_aphids() < extinct_N) aphids.clear();
                if ((p.total_aphids() + p.total_mummies()) == 0) {
                    p.empty = true;
                }
            }
        } else if (pert.index == n_lines) { // mummies + parasitized aphids
            for (OnePlant& p : field.plants) {
                p.mummies.clear(mult);
                if (arma::accu(p.mummies.Y) < extinct_N) p.mummies.Y.fill(0);
                for (AphidPop& aphids : p.aphids) {
                    aphids.paras.clear(mult);
                    if (aphids.total_aphids() < extinct_N) aphids.clear();
                }
                if ((p.total_aphids() + p.total_mummies()) == 0) {
                    p.empty = true;
                }
            }
        } else { // adult wasps
            field.wasps.Y *= mult;
            if (field.wasps.Y < extinct_N) field.wasps.Y = 0;
        }
        i++;
    }

    if (do_pop && i > 0) {
        perturbs.erase(perturbs.begin(), perturbs.begin() + i);
    }

    return;
}





/*
 `T` should be uint32 for using age to determine whether a plant gets cleared.
 It should be double for using aphid abundance to determine whether a plant gets cleared.
 */

template <typename T>
RepSummary one_rep__(const T& clear_threshold,
                     const uint32& rep,
                     const uint32& n_fields,
                     std::deque<uint32> check_for_clear,
                     const double& clear_surv,
                     const uint32& max_t,
                     const uint32& save_every,
                     const double& mean_K,
                     const double& sd_K,
                     const std::vector<double>& K_y_mult,
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
                     const double& alate_disp_prop,
                     const std::vector<double>& disp_rate,
                     const std::vector<double>& disp_mort,
                     const std::vector<uint32>& disp_start,
                     const std::vector<uint32>& living_days,
                     const std::vector<double>& pred_rate,
                     const arma::mat& mum_density_0,
                     const double& mum_smooth,
                     const double& max_mum_density,
                     const arma::vec& rel_attack,
                     const double& a,
                     const double& k,
                     const double& h,
                     const std::vector<double>& wasp_density_0,
                     const std::vector<uint32>& wasp_delay,
                     const double& wasp_dispersal_p,
                     const double& sex_ratio,
                     const std::vector<double>& s_y,
                     const std::vector<bool>& constant_wasps,
                     const std::vector<uint32>& perturb_when,
                     const std::vector<uint32>& perturb_where,
                     const std::vector<uint32>& perturb_who,
                     const std::vector<double>& perturb_how,
                     Progress& prog_bar,
                     int& status_code,
                     pcg32& eng) {

    double demog_mult = 1;
    if (!demog_error) demog_mult = 0;
    // either type of environmental error
    bool process_error = (sigma_x > 0) || (sigma_y > 0);

    uint32 n_lines = aphid_name.size();
    uint32 n_plants = aphid_density_0.size();

    RepSummary summary;

    summary.reserve(rep, max_t, save_every, n_lines, n_fields, n_plants);

    uint32 iters = 0;

    std::vector<OneField> fields;
    fields.reserve(n_fields);
    for (uint32 i = 0; i < n_fields; i++) {
        fields.push_back(
            OneField(sigma_x, sigma_y, rho, demog_mult, mean_K, sd_K, K_y_mult[i],
                    death_prop,
                    shape1_death_mort, shape2_death_mort, attack_surv,
                    aphid_name, leslie_mat, aphid_density_0, alate_b0, alate_b1,
                    disp_rate, disp_mort, disp_start, living_days, pred_rate,
                    extinct_N, mum_density_0, mum_smooth, max_mum_density,
                    rel_attack, a, k, h, 0, sex_ratio,
                    s_y[i], constant_wasps[i], eng));
        if (wasp_delay[i] == 0) fields.back().wasps.Y = wasp_density_0[i];
    }



    std::deque<PerturbInfo> perturbs(perturb_when.size());
    for (uint32 i = 0; i < perturb_when.size(); i++) {
        perturbs[i] = PerturbInfo(perturb_when[i], perturb_where[i],
                                  perturb_who[i], perturb_how[i]);
    }


    summary.push_back(0, fields);

    for (uint32 t = 1; t <= max_t; t++) {

        if (interrupt_check(iters, prog_bar)) {
            status_code = -1;
            return summary;
        }

        // Perturbations
        do_perturb(perturbs, fields, t, aphid_name.size(), extinct_N, true);

        if (n_fields > 1 && alate_disp_prop > 0 &&
            !check_for_clear.empty() && t == check_for_clear.front()) {

            arma::mat D = fields.front().remove_dispersers(alate_disp_prop);
            for (uint32 i = 1; i < fields.size(); i++) {
                D += fields[i].remove_dispersers(alate_disp_prop);
            }
            D /= static_cast<double>(n_fields);
            for (OneField& c : fields) c.add_dispersers(D);

        }

        // used to calculate wasp dispersal if desired:
        double from_wasp_pool = 0;
        if (wasp_dispersal_p > 0 && n_fields > 1) {
            double wasp_Y_bar = 0;
            for (OneField& field : fields) wasp_Y_bar += field.wasps.Y;
            wasp_Y_bar /= static_cast<double>(n_fields);
            from_wasp_pool = wasp_Y_bar * wasp_dispersal_p;
        }

        for (uint32 i = 0; i < n_fields; i++) {

            OneField& field(fields[i]);

            if (disp_error) {
                field.calc_dispersal(eng);
            } else field.calc_dispersal();

            if (process_error) {
                field.update(eng);
            } else field.update();

            // Add from wasp dispersal pool
            if (wasp_dispersal_p > 0 && n_fields > 1) {
                field.wasps.Y *= (1 - wasp_dispersal_p);
                field.wasps.Y += from_wasp_pool;
            }

            if (t == wasp_delay[i]) field.wasps.Y += wasp_density_0[i];

        }

        if (t % save_every == 0 || t == max_t) summary.push_back(t, fields);

        // If all fields are empty, then stop this rep.
        // It's important to do this before clearing plants.
        bool all_empty = true;
        for (uint32 i = 0; i < n_fields; i++) {
            for (const OnePlant& p : fields[i].plants) {
                if (!p.empty) {
                    all_empty = false;
                    break;
                }
            }
            if (!all_empty) break;
        }
        if (all_empty) break;

        if (!check_for_clear.empty() && t == check_for_clear.front()) {
            check_for_clear.pop_front();
            for (uint32 i = 0; i < n_fields; i++) {
                fields[i].clear_plants(clear_threshold, clear_surv, eng);
            }
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
                const uint32& n_fields,
                const uint32& n_plants,
                const uint32& max_plant_age,
                const double& max_N,
                const std::deque<uint32>& check_for_clear,
                const double& clear_surv,
                const uint32& max_t,
                const uint32& save_every,
                const double& mean_K,
                const double& sd_K,
                const std::vector<double>& K_y_mult,
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
                const double& mum_smooth,
                const arma::vec& rel_attack,
                const double& a,
                const double& k,
                const double& h,
                const std::vector<double>& wasp_density_0,
                const std::vector<uint32>& wasp_delay,
                const double& wasp_dispersal_p,
                const double& sex_ratio,
                const std::vector<double>& s_y,
                const std::vector<bool>& constant_wasps,
                const std::vector<uint32>& perturb_when,
                const std::vector<uint32>& perturb_where,
                const std::vector<uint32>& perturb_who,
                const std::vector<double>& perturb_how,
                uint32& n_threads) {


    /*
     ===============================================================
     Checking dimensions
     ===============================================================
     */
    uint32 n_stages = aphid_density_0.front().n_rows;

    if (n_reps == 0) stop("\nERROR: n_reps == 0\n");
    if (n_lines == 0) stop("\nERROR: n_lines == 0\n");
    if (n_fields == 0) stop("\nERROR: n_fields == 0\n");
    if (n_plants == 0) stop("\nERROR: n_plants == 0\n");
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
    if (aphid_density_0.size() != n_plants) {
        stop("\nERROR: aphid_density_0.size() != n_plants\n");
    }
    for (uint32 i = 0; i < n_plants; i++) {
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

    if (K_y_mult.size() != n_fields) {
        stop("\nERROR: K_y_mult.size() != n_fields\n");
    }
    if (wasp_density_0.size() != n_fields) {
        stop("\nERROR: wasp_density_0.size() != n_fields\n");
    }
    if (wasp_delay.size() != n_fields) {
        stop("\nERROR: wasp_delay.size() != n_fields\n");
    }
    if (s_y.size() != n_fields) {
        stop("\nERROR: s_y.size() != n_fields\n");
    }
    if (constant_wasps.size() != n_fields) {
        stop("\nERROR: constant_wasps.size() != n_fields\n");
    }

    if (pred_rate.size() != n_plants) {
        stop("\nERROR: pred_rate.size() != n_plants\n");
    }
    if (mum_density_0.n_cols != n_plants) {
        stop("\nERROR: mum_density_0.n_cols != n_plants\n");
    }
    if (mum_density_0.n_rows == 0) {
        stop("\nERROR: mum_density_0.n_rows == 0\n");
    }

    if (perturb_when.size() != perturb_where.size() ||
        perturb_when.size() != perturb_who.size() ||
        perturb_when.size() != perturb_how.size()) {
        stop("\nERROR: perturb_* args not same size.");
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
    one_negative_check(shape1_death_mort, "shape1_death_mort");
    one_negative_check(shape2_death_mort, "shape2_death_mort");
    one_negative_check(sigma_x, "sigma_x");
    one_negative_check(sigma_y, "sigma_y");
    one_negative_check(rho, "rho");
    one_negative_check(extinct_N, "extinct_N");
    one_negative_check(a, "a");
    one_negative_check(k, "k");
    one_negative_check(h, "h");

    // objects containing doubles that must be >= 0
    negative_check<std::vector<double>>(K_y_mult, "K_y_mult");
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
    negative_check<std::vector<double>>(wasp_density_0, "wasp_density_0");

    // doubles / double vectors that must be >= 0 and <= 1
    one_non_prop_check(mum_smooth, "mum_smooth");
    one_non_prop_check(clear_surv, "clear_surv");
    one_non_prop_check(death_prop, "death_prop");
    one_non_prop_check(sex_ratio, "sex_ratio");
    one_non_prop_check(wasp_dispersal_p, "wasp_dispersal_p");
    non_prop_check<std::vector<double>>(s_y, "s_y");

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
                       const uint32& n_fields,
                       const uint32& max_plant_age,
                       const double& max_N,
                       const std::deque<uint32>& check_for_clear,
                       const double& clear_surv,
                       const uint32& max_t,
                       const uint32& save_every,
                       const double& mean_K,
                       const double& sd_K,
                       const std::vector<double>& K_y_mult,
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
                       const double& alate_disp_prop,
                       const std::vector<double>& disp_rate,
                       const std::vector<double>& disp_mort,
                       const std::vector<uint32>& disp_start,
                       const std::vector<uint32>& living_days,
                       const std::vector<double>& pred_rate,
                       const arma::mat& mum_density_0,
                       const double& mum_smooth,
                       const double& max_mum_density,
                       const arma::vec& rel_attack,
                       const double& a,
                       const double& k,
                       const double& h,
                       const std::vector<double>& wasp_density_0,
                       const std::vector<uint32>& wasp_delay,
                       const double& wasp_dispersal_p,
                       const double& sex_ratio,
                       const std::vector<double>& s_y,
                       const std::vector<bool>& constant_wasps,
                       const std::vector<uint32>& perturb_when,
                       const std::vector<uint32>& perturb_where,
                       const std::vector<uint32>& perturb_who,
                       const std::vector<double>& perturb_how,
                       uint32 n_threads,
                       const bool& show_progress) {

    uint32 n_lines = aphid_name.size();
    uint32 n_plants = aphid_density_0.size();

    check_args(n_reps, n_lines, n_fields, n_plants,
               max_plant_age, max_N, check_for_clear, clear_surv,
               max_t, save_every,
               mean_K, sd_K, K_y_mult, death_prop,
               shape1_death_mort, shape2_death_mort,
               attack_surv, disp_error, demog_error,
               sigma_x, sigma_y, rho, extinct_N, aphid_name,
               leslie_mat, aphid_density_0, alate_b0, alate_b1,
               disp_rate, disp_mort, disp_start, living_days,
               pred_rate, mum_density_0, mum_smooth, rel_attack, a, k, h,
               wasp_density_0, wasp_delay, wasp_dispersal_p,
               sex_ratio, s_y, constant_wasps,
               perturb_when, perturb_where, perturb_who, perturb_how, n_threads);

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
            summaries[i] = one_rep__<uint32>(max_plant_age, i, n_fields, check_for_clear,
                                             clear_surv, max_t,
                                             save_every, mean_K, sd_K, K_y_mult,
                                             death_prop,
                                             shape1_death_mort, shape2_death_mort,
                                             attack_surv, disp_error,
                                             demog_error, sigma_x, sigma_y, rho,
                                             extinct_N,
                                             aphid_name, leslie_mat, aphid_density_0,
                                             alate_b0, alate_b1,
                                             alate_disp_prop,
                                             disp_rate, disp_mort,
                                             disp_start, living_days, pred_rate,
                                             mum_density_0, mum_smooth,
                                             max_mum_density,
                                             rel_attack, a, k,
                                             h, wasp_density_0, wasp_delay,
                                             wasp_dispersal_p,
                                             sex_ratio, s_y, constant_wasps,
                                             perturb_when, perturb_where,
                                             perturb_who, perturb_how,
                                             prog_bar, status_code, eng);
        } else {
            summaries[i] = one_rep__<double>(max_N, i, n_fields, check_for_clear,
                                             clear_surv, max_t,
                                             save_every, mean_K, sd_K, K_y_mult,
                                             death_prop,
                                             shape1_death_mort, shape2_death_mort,
                                             attack_surv, disp_error,
                                             demog_error, sigma_x, sigma_y, rho,
                                             extinct_N,
                                             aphid_name, leslie_mat, aphid_density_0,
                                             alate_b0, alate_b1,
                                             alate_disp_prop,
                                             disp_rate, disp_mort,
                                             disp_start, living_days, pred_rate,
                                             mum_density_0, mum_smooth,
                                             max_mum_density,
                                             rel_attack, a, k,
                                             h, wasp_density_0, wasp_delay,
                                             wasp_dispersal_p,
                                             sex_ratio, s_y, constant_wasps,
                                             perturb_when, perturb_where,
                                             perturb_who, perturb_how,
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
                                _["field"] = summ.field,
                                _["plant"] = summ.plant,
                                _["line"] = summ.line,
                                _["type"] = summ.type,
                                _["N"] = summ.N),
                            _["wasps"] = DataFrame::create(
                                _["rep"] = summ.wasp_rep,
                                _["time"] = summ.wasp_time,
                                _["field"] = summ.wasp_field,
                                _["wasps"] = summ.wasp_N));

    return out;
}
