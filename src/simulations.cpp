
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
#include <utility>              // std::pair, std::make_pair
#include <deque>                // deque
#include <pcg/pcg_random.hpp>   // pcg prng
#include <progress.hpp>         // for the progress bar
#ifdef _OPENMP
#include <omp.h>                // OpenMP
#endif


#include "pseudogameofclones_types.hpp"  // integer types
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
                   const uint32& n_plants,
                   const bool& sep_adults) {

    // Rows for wasps is just the number of time points * number of fields:
    n_rows_wasps = (max_t / save_every) + 1;
    if (max_t % save_every > 0) n_rows_wasps++;
    n_rows_wasps *= n_fields;

    // For aphids, you also include plants, lines, and type (alate, apterous,
    // parasitized). If `sep_adults == TRUE`, then you have separate
    // rows for adult vs juvenile apterous and adult vs juvenile alates.
    n_rows = n_plants * n_lines * n_rows_wasps;
    if (sep_adults) {
        /*
         for (1) juvenile alate vs (2) adult alate vs
         (3) juvenile apterous vs (4) adult apterous vs (5) parasitized
         */
        n_rows *= 5;
    } else n_rows *= 3;  // for (1) alate vs (2) apterous vs (3) parasitized
    n_rows += n_plants * n_rows_wasps;  // for mummies

    return;
}




struct RepSummary {

    bool sep_adults;

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

    // Don't allow default constructor
    RepSummary() = delete; // will never be generated

    RepSummary(const bool& sep_adults_)
        : sep_adults(sep_adults_),
          rep(), time(), field(), plant(), line(), type(), N(),
          wasp_rep(), wasp_time(), wasp_field(), wasp_N(), r() {};

    void reserve(const uint32& rep_,
                 const uint32& max_t,
                 const uint32& save_every,
                 const uint32& n_lines,
                 const uint32& n_fields,
                 const uint32& n_plants) {
        uint32 n_rows, n_rows_wasps;
        calc_rep_rows(n_rows, n_rows_wasps, max_t, save_every,
                      n_lines, n_fields, n_plants, sep_adults);
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
                   const AllFields& fields) {

        for (uint32 k = 0; k < fields.size(); k++) {

            const OneField& field(fields[k]);

            // Everything but wasps:
            for (uint32 j = 0; j < field.size(); j++) {
                const OnePlant& plant(field[j]);
                for (uint32 i = 0; i < plant.size(); i++) {
                    const AphidPop& aphid(plant[i]);
                    if (sep_adults) {
                        append_living_aphids__(t, k, j, aphid.aphid_name,
                                               aphid.total_adult_alates(),
                                               aphid.total_juven_alates(),
                                               aphid.total_adult_apterous(),
                                               aphid.total_juven_apterous(),
                                               aphid.paras.total_aphids());
                    } else {
                        append_living_aphids__(t, k, j, aphid.aphid_name,
                                               aphid.alates.total_aphids(),
                                               aphid.apterous.total_aphids(),
                                               aphid.paras.total_aphids());
                    }
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

        for (uint32 i = 0; i < 3; i++) {
            rep.push_back(r);
            time.push_back(t);
            field.push_back(c);
            plant.push_back(p);
            line.push_back(l);
        }

        type.push_back("alate");
        type.push_back("apterous");
        type.push_back("parasitized");

        N.push_back(N_ala);
        N.push_back(N_apt);
        N.push_back(N_par);

        return;
    }



    inline void append_living_aphids__(const uint32& t,
                                       const uint32& c,
                                       const uint32& p,
                                       const std::string& l,
                                       const double& N_adult_ala,
                                       const double& N_juven_ala,
                                       const double& N_adult_apt,
                                       const double& N_juven_apt,
                                       const double& N_par) {

        for (uint32 i = 0; i < 5; i++) {
            rep.push_back(r);
            time.push_back(t);
            field.push_back(c);
            plant.push_back(p);
            line.push_back(l);
        }

        type.push_back("adult alate");
        type.push_back("juvenile alate");
        type.push_back("adult apterous");
        type.push_back("juvenile apterous");
        type.push_back("parasitized");

        N.push_back(N_adult_ala);
        N.push_back(N_juven_ala);
        N.push_back(N_adult_apt);
        N.push_back(N_juven_apt);
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









void one_rep__(RepSummary& summary,
               AllFields& fields,
               const uint32& rep,
               const uint32& max_t,
               const uint32& save_every,
               const std::vector<uint32>& perturb_when,
               const std::vector<uint32>& perturb_where,
               const std::vector<uint32>& perturb_who,
               const std::vector<double>& perturb_how,
               Progress& prog_bar,
               int& status_code) {

    uint32 n_lines = fields.n_lines();
    uint32 n_plants = fields.n_plants();
    uint32 n_fields = fields.size();

    summary.reserve(rep, max_t, save_every, n_lines, n_fields, n_plants);

    uint32 iters = 0;

    std::deque<PerturbInfo> perturbs(perturb_when.size());
    for (uint32 i = 0; i < perturb_when.size(); i++) {
        perturbs[i] = PerturbInfo(perturb_when[i], perturb_where[i],
                                  perturb_who[i], perturb_how[i]);
    }

    summary.push_back(0, fields);


    for (uint32 t = 1; t <= max_t; t++) {

        if (interrupt_check(iters, prog_bar)) {
            status_code = -1;
            return;
        }

        // Basic updates for perturbations, dispersal, and population dynamics.
        // Returns true if all plants/fields are empty.
        // It's important to do this (& to save output) before clearing plants.
        bool all_empty = fields.update(t, perturbs);

        if (t % save_every == 0 || t == max_t) summary.push_back(t, fields);


        // If all fields are empty, then stop this rep.
        if (all_empty) break;

    }


    return;

}



void one_rep__(RepSummary& summary,
               std::vector<AllStageInfo>& stage_ts,
               AllFields& fields,
               const uint32& rep,
               const uint32& max_t,
               const uint32& save_every,
               const std::vector<uint32>& perturb_when,
               const std::vector<uint32>& perturb_where,
               const std::vector<uint32>& perturb_who,
               const std::vector<double>& perturb_how,
               Progress& prog_bar,
               int& status_code) {

    uint32 n_lines = fields.n_lines();
    uint32 n_plants = fields.n_plants();
    uint32 n_fields = fields.size();

    summary.reserve(rep, max_t, save_every, n_lines, n_fields, n_plants);

    stage_ts.reserve(max_t);

    uint32 iters = 0;

    std::deque<PerturbInfo> perturbs(perturb_when.size());
    for (uint32 i = 0; i < perturb_when.size(); i++) {
        perturbs[i] = PerturbInfo(perturb_when[i], perturb_where[i],
                                  perturb_who[i], perturb_how[i]);
    }

    summary.push_back(0, fields);

    for (uint32 t = 1; t <= max_t; t++) {

        if (interrupt_check(iters, prog_bar)) {
            status_code = -1;
            return;
        }

        // Basic updates for perturbations, dispersal, and population dynamics.
        // Returns true if all plants/fields are empty.
        // It's important to do this (& to save output) before clearing plants.
        bool all_empty = fields.update(t, perturbs);

        if (t % save_every == 0 || t == max_t) summary.push_back(t, fields);

        // If all fields are empty, then stop this rep.
        if (all_empty) break;

        stage_ts.push_back(fields.out_all_info());


    }


    return;

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


// Same but number must be > 0
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
                const uint32& max_t,
                const uint32& save_every,
                const std::vector<double>& K,
                const std::vector<double>& K_y,
                const arma::mat& attack_surv,
                const bool& aphid_demog_error,
                const bool& wasp_demog_error,
                const double& sigma_x,
                const double& sigma_y,
                const double& rho,
                const double& extinct_N,
                const std::vector<std::string>& aphid_name,
                const std::vector<arma::cube>& leslie_mat,
                const std::vector<arma::cube>& aphid_density_0,
                const std::vector<double>& alate_b0,
                const std::vector<double>& alate_b1,
                const std::vector<uint32>& field_disp_start,
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
                const double& wasp_disp_m0,
                const double& wasp_disp_m1,
                const std::vector<double>& wasp_field_attract,
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
    if (field_disp_start.size() != n_lines) {
        stop("\nERROR: field_disp_start.size() != n_lines\n");
    }
    if (living_days.size() != n_lines) {
        stop("\nERROR: living_days.size() != n_lines\n");
    }
    if (attack_surv.n_cols != n_lines || attack_surv.n_rows < 2) {
        stop(std::string("\nERROR: attack_surv.n_cols != n_lines || ") +
            std::string("attack_surv.n_rows < 2\n"));
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

    if (K_y.size() != n_fields) {
        stop("\nERROR: K_y.size() != n_fields\n");
    }
    if (wasp_density_0.size() != n_fields) {
        stop("\nERROR: wasp_density_0.size() != n_fields\n");
    }
    if (wasp_delay.size() != n_fields) {
        stop("\nERROR: wasp_delay.size() != n_fields\n");
    }
    if (wasp_field_attract.size() != n_fields) {
        stop("\nERROR: wasp_field_attract.size() != n_fields\n");
    }
    if (s_y.size() != n_fields) {
        stop("\nERROR: s_y.size() != n_fields\n");
    }
    if (constant_wasps.size() != n_fields) {
        stop("\nERROR: constant_wasps.size() != n_fields\n");
    }

    if (pred_rate.size() != n_fields) {
        stop("\nERROR: pred_rate.size() != n_fields\n");
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
    one_negative_check(sigma_x, "sigma_x");
    one_negative_check(sigma_y, "sigma_y");
    one_negative_check(rho, "rho");
    one_negative_check(extinct_N, "extinct_N");
    one_negative_check(a, "a");
    one_negative_check(k, "k");
    one_negative_check(h, "h");

    // objects containing doubles that must be >= 0
    negative_check<std::vector<double>>(K, "K");
    negative_check<std::vector<double>>(K_y, "K_y");
    negative_check<arma::mat>(attack_surv, "attack_surv");
    for (uint32 i = 0; i < leslie_mat.size(); i++) {
        negative_check<arma::cube>(leslie_mat[i], "leslie_mat");
    }
    for (uint32 i = 0; i < aphid_density_0.size(); i++) {
        negative_check<arma::cube>(aphid_density_0[i], "aphid_density_0");
    }
    negative_check<std::vector<double>>(pred_rate, "pred_rate");
    negative_check<arma::mat>(mum_density_0, "mum_density_0");
    negative_check<arma::vec>(rel_attack, "rel_attack");
    negative_check<std::vector<double>>(wasp_density_0, "wasp_density_0");
    negative_check<std::vector<double>>(wasp_field_attract, "wasp_field_attract");


    // doubles / double vectors that must be >= 0 and <= 1
    one_non_prop_check(mum_smooth, "mum_smooth");
    one_non_prop_check(sex_ratio, "sex_ratio");
    one_non_prop_check(wasp_disp_m0, "wasp_disp_m0");
    // one_non_prop_check(wasp_disp_m1, "wasp_disp_m1");  // << this can be anything
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

    // This vector can have zeros, but can't have all zeros:
    double wfa_sum = std::accumulate(wasp_field_attract.begin(),
                                     wasp_field_attract.end(), 0.0);
    if (wfa_sum <= 0) {
        std::string msg = "\nERROR: wasp_field_attract sums to <= 0.\n";
        stop(msg.c_str());
    }


    return;
}





//[[Rcpp::export]]
List sim_pseudogameofclones_cpp(const uint32& n_reps,
                       const uint32& n_fields,
                       const uint32& max_t,
                       const uint32& save_every,
                       const std::vector<double>& K,
                       const std::vector<double>& K_y,
                       const arma::mat& attack_surv,
                       const bool& aphid_demog_error,
                       const bool& wasp_demog_error,
                       const double& sigma_x,
                       const double& sigma_y,
                       const double& rho,
                       const double& extinct_N,
                       const std::vector<std::string>& aphid_name,
                       const std::vector<arma::cube>& leslie_mat,
                       const std::vector<arma::cube>& aphid_density_0,
                       const std::vector<double>& alate_b0,
                       const std::vector<double>& alate_b1,
                       const double& alate_field_disp_p,
                       const std::vector<uint32>& field_disp_start,
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
                       const double& wasp_disp_m0,
                       const double& wasp_disp_m1,
                       const std::vector<double>& wasp_field_attract,
                       const double& sex_ratio,
                       const std::vector<double>& s_y,
                       const std::vector<bool>& constant_wasps,
                       const std::vector<uint32>& perturb_when,
                       const std::vector<uint32>& perturb_where,
                       const std::vector<uint32>& perturb_who,
                       const std::vector<double>& perturb_how,
                       const bool& sep_adults,
                       uint32 n_threads,
                       const bool& show_progress) {

    uint32 n_lines = aphid_name.size();
    uint32 n_plants = aphid_density_0.size();

    check_args(n_reps, n_lines, n_fields, n_plants,
               max_t, save_every,
               K, K_y,
               attack_surv, aphid_demog_error, wasp_demog_error,
               sigma_x, sigma_y, rho, extinct_N, aphid_name,
               leslie_mat, aphid_density_0, alate_b0, alate_b1,
               field_disp_start,
               living_days,
               pred_rate, mum_density_0, mum_smooth,
               rel_attack, a, k, h,
               wasp_density_0, wasp_delay, wasp_disp_m0, wasp_disp_m1,
               wasp_field_attract,
               sex_ratio, s_y, constant_wasps,
               perturb_when, perturb_where, perturb_who, perturb_how,
               n_threads);

    Progress prog_bar(max_t * n_reps, show_progress);
    std::vector<int> status_codes(n_threads, 0);

    // Generate seeds for random number generators (1 set of seeds per rep)
    const std::vector<std::vector<uint64>> seeds = mt_seeds(n_reps);

    std::vector<RepSummary> summaries(n_reps, RepSummary(sep_adults));

    /*
     Setting info for all AllFields objects.
     I'm doing it this way so that I can store all info (including stages)
     for the end of each rep, which also allows an easier time of repeating
     things after a custom perturbation.
     */
    XPtr<std::vector<AllFields>> all_fields_vec_xptr(
            new std::vector<AllFields>(n_reps), true);
    std::vector<AllFields>& all_fields_vec(*all_fields_vec_xptr);
    // Fill in first AllFields object, then copy that to the rest:
    all_fields_vec[0] = AllFields(
        n_fields, wasp_delay,
        sigma_x, sigma_y, rho, aphid_demog_error, wasp_demog_error,
        K, K_y,
        attack_surv,
        aphid_name, leslie_mat, aphid_density_0,
        alate_b0, alate_b1,
        field_disp_start,
        living_days, pred_rate, extinct_N,
        mum_density_0, mum_smooth,
        rel_attack, a, k, h, wasp_density_0, sex_ratio,
        s_y, constant_wasps,
        alate_field_disp_p, wasp_disp_m0, wasp_disp_m1,
        wasp_field_attract, seeds[0]);
    for (uint32 i = 1; i < n_reps; i++){
        all_fields_vec[i] = all_fields_vec[0];
        all_fields_vec[i].reseed(seeds[i]);
    }

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

    // Parallelize the Loop
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (uint32 i = 0; i < n_reps; i++) {
        if (status_code != 0) continue;
        one_rep__(summaries[i], all_fields_vec[i], i,
                  max_t, save_every,
                  perturb_when, perturb_where, perturb_who, perturb_how,
                  prog_bar, status_code);
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
                                _["wasps"] = summ.wasp_N),
                            _["all_info_xptr"] = all_fields_vec_xptr);

    return out;
}










// Make new `XPtr<std::vector<AllFields>>` object with new parameters
// before restarting simulation.
//[[Rcpp::export]]
SEXP restart_fill_other_pars(SEXP all_fields_in_ptr,
                             const std::vector<double>& K,
                             const std::vector<double>& alate_b0,
                             const std::vector<double>& alate_b1,
                             const double& alate_field_disp_p,
                             const std::vector<double>& K_y,
                             const std::vector<double>& s_y,
                             const double& a,
                             const double& k,
                             const double& h,
                             const double& wasp_disp_m0,
                             const double& wasp_disp_m1,
                             const std::vector<double>& wasp_field_attract,
                             const double& mum_smooth,
                             const std::vector<double>& pred_rate) {

    XPtr<std::vector<AllFields>> all_fields_vec_in_xptr(all_fields_in_ptr);
    std::vector<AllFields>& all_fields_vec_in(*all_fields_vec_in_xptr);

    uint32 n_reps = all_fields_vec_in.size();
    if (n_reps == 0) stop("\nERROR: n_reps == 0\n");

    XPtr<std::vector<AllFields>> all_fields_vec_xptr(
            new std::vector<AllFields>(n_reps), true);
    std::vector<AllFields>& all_fields_vec(*all_fields_vec_xptr);
    for (uint32 i = 0; i < n_reps; i++){
        all_fields_vec[i] = all_fields_vec_in[i];
    }

    // Checking dimensions of vectors
    uint32 n_lines = all_fields_vec[0].n_lines();
    uint32 n_fields = all_fields_vec[0].size();
    uint32 n_plants = all_fields_vec[0].n_plants();
    if (n_lines == 0) stop("\nERROR: n_lines == 0\n");
    if (n_fields == 0) stop("\nERROR: n_fields == 0\n");
    if (n_plants == 0) stop("\nERROR: n_plants == 0\n");
    if (alate_b0.size() != n_lines) {
        stop("\nERROR: alate_b0.size() != n_lines\n");
    }
    if (alate_b1.size() != n_lines) {
        stop("\nERROR: alate_b1.size() != n_lines\n");
    }
    if (K.size() != n_fields) {
        stop("\nERROR: K.size() != n_fields\n");
    }
    if (K_y.size() != n_fields) {
        stop("\nERROR: K_y.size() != n_fields\n");
    }
    if (s_y.size() != n_fields) {
        stop("\nERROR: s_y.size() != n_fields\n");
    }
    if (pred_rate.size() != n_fields) {
        stop("\nERROR: pred_rate.size() != n_fields\n");
    }
    if (wasp_field_attract.size() != n_fields) {
        stop("\nERROR: wasp_field_attract.size() != n_fields\n");
    }
    // This vector can have zeros, but can't have all zeros:
    double wfa_sum = std::accumulate(wasp_field_attract.begin(),
                                     wasp_field_attract.end(), 0.0);
    if (wfa_sum <= 0) {
        std::string msg = "\nERROR: wasp_field_attract sums to <= 0.\n";
        stop(msg.c_str());
    }

    for (AllFields& fields : all_fields_vec) {

        fields.set_new_pars(K, alate_b0, alate_b1, alate_field_disp_p,
                            K_y, s_y, a, k, h, wasp_disp_m0, wasp_disp_m1,
                            wasp_field_attract,
                            mum_smooth, pred_rate);

    }


    return all_fields_vec_xptr;

}





//[[Rcpp::export]]
List restart_experiments_cpp(SEXP all_fields_ptr,
                             const uint32& max_t,
                             const uint32& save_every,
                             const bool& stage_ts_out,
                             const bool& sep_adults,
                             const bool& show_progress,
                             const std::vector<uint32>& perturb_when,
                             const std::vector<uint32>& perturb_where,
                             const std::vector<uint32>& perturb_who,
                             const std::vector<double>& perturb_how) {


    XPtr<std::vector<AllFields>> all_fields_vec_xptr(all_fields_ptr);
    std::vector<AllFields>& all_fields_vec(*all_fields_vec_xptr);

    uint32 n_reps = all_fields_vec.size();

    // No perturbations allowed here:
    Progress prog_bar(max_t * n_reps, show_progress);
    int status_code = 0;

    std::vector<RepSummary> summaries(n_reps, RepSummary(sep_adults));

    std::vector<std::vector<AllStageInfo>> stage_ts(n_reps);

    if (stage_ts_out) {
        for (uint32 i = 0; i < n_reps; i++) {
            if (status_code != 0) break;
            one_rep__(summaries[i], stage_ts[i], all_fields_vec[i], i,
                      max_t, save_every,
                      perturb_when, perturb_where, perturb_who, perturb_how,
                      prog_bar, status_code);
        }
    } else {
        for (uint32 i = 0; i < n_reps; i++) {
            if (status_code != 0) break;
            one_rep__(summaries[i], all_fields_vec[i], i,
                      max_t, save_every,
                      perturb_when, perturb_where, perturb_who, perturb_how,
                      prog_bar, status_code);
        }
    }

    if (status_code != 0) stop("User interruption.");


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
                                _["wasps"] = summ.wasp_N),
                            _["all_info_xptr"] = all_fields_vec_xptr);
    if (stage_ts_out) {
        List stage_ts_list;
        if (n_reps > 1) {
            stage_ts_list = List(n_reps);
            for (uint32 i = 0; i < n_reps; i++) {
                const std::vector<AllStageInfo>& stage_ts_i(stage_ts[i]);
                List stage_ts_list_i = List(stage_ts_i.size());
                for (uint32 j = 0; j < stage_ts_i.size(); j++) {
                    stage_ts_list_i(j) = stage_ts_i[j].to_data_frame();
                }
                stage_ts_list(i) = stage_ts_list_i;
            }
        } else {
            const std::vector<AllStageInfo>& stage_ts_i(stage_ts[0]);
            stage_ts_list = List(stage_ts_i.size());
            for (uint32 j = 0; j < stage_ts_i.size(); j++) {
                stage_ts_list(j) = stage_ts_i[j].to_data_frame();
            }
        }
        out["stage_ts"] = stage_ts_list;
    }

    return out;

}

