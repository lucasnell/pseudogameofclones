
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

// #define RCPPTHREAD_OVERRIDE_COUT 1    // std::cout override
// #define RCPPTHREAD_OVERRIDE_CERR 1    // std::cerr override
// #define RCPPTHREAD_OVERRIDE_THREAD 1  // std::thread override
#include <RcppThread.h>         // multithreading



#include "pseudogameofclones_types.hpp"  // integer types
#include "aphids.hpp"           // aphid classes
#include "patches.hpp"          // field classes
#include "pcg.hpp"              // mt_seeds seed_pcg fxns







//' Check that the number of threads doesn't exceed the number available.
//'
//' @noRd
//'
inline void thread_check(uint32& n_threads) {

    if (n_threads == 0) n_threads = 1;

    uint32 max_threads = std::thread::hardware_concurrency();

    if (n_threads > max_threads) {
        std::string mt_str = std::to_string(max_threads);
        std::string err_msg = "\nThe number of requested threads (" +
            std::to_string(n_threads) +
            ") exceeds the max available on the system (" + mt_str + ").";
        stop(err_msg.c_str());
    }

    return;
}






// Calculate the number of rows per rep.
void calc_rep_rows(uint32& n_rows,
                   uint32& n_rows_wasps,
                   const uint32& max_t,
                   const uint32& save_every,
                   const uint32& n_lines,
                   const uint32& n_fields,
                   const bool& sep_adults) {

    // Rows for wasps is just the number of time points * number of fields:
    n_rows_wasps = (max_t / save_every) + 1;
    if (max_t % save_every > 0) n_rows_wasps++;
    n_rows_wasps *= n_fields;

    // For aphids, you also include lines, and type (alate, apterous,
    // parasitized). If `sep_adults == TRUE`, then you have separate
    // rows for adult vs juvenile apterous and adult vs juvenile alates.
    n_rows = n_lines * n_rows_wasps;
    if (sep_adults) {
        /*
         for (1) juvenile alate vs (2) adult alate vs
         (3) juvenile apterous vs (4) adult apterous vs (5) parasitized
         */
        n_rows *= 5;
    } else n_rows *= 3;  // for (1) alate vs (2) apterous vs (3) parasitized
    n_rows += n_rows_wasps;  // for mummies

    return;
}




struct RepSummary {

    bool sep_adults;

    std::vector<uint32> rep;
    std::vector<uint32> time;
    std::vector<uint32> field;
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
          rep(), time(), field(), line(), type(), N(),
          wasp_rep(), wasp_time(), wasp_field(), wasp_N(), r() {};

    void reserve(const uint32& rep_,
                 const uint32& max_t,
                 const uint32& save_every,
                 const uint32& n_lines,
                 const uint32& n_fields) {
        uint32 n_rows, n_rows_wasps;
        calc_rep_rows(n_rows, n_rows_wasps, max_t, save_every,
                      n_lines, n_fields, sep_adults);
        rep.reserve(n_rows);
        time.reserve(n_rows);
        field.reserve(n_rows);
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

            for (uint32 i = 0; i < field.size(); i++) {
                const AphidPop& aphid(field[i]);
                if (sep_adults) {
                    append_living_aphids__(t, k, aphid.aphid_name,
                                           aphid.total_adult_alates(),
                                           aphid.total_juven_alates(),
                                           aphid.total_adult_apterous(),
                                           aphid.total_juven_apterous(),
                                           aphid.paras.total_aphids());
                } else {
                    append_living_aphids__(t, k, aphid.aphid_name,
                                           aphid.alates.total_aphids(),
                                           aphid.apterous.total_aphids(),
                                           aphid.paras.total_aphids());
                }
            }
            append_mummies__(t, k, field.total_mummies());

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
                                       const std::string& l,
                                       const double& N_ala,
                                       const double& N_apt,
                                       const double& N_par) {

        for (uint32 i = 0; i < 3; i++) {
            rep.push_back(r);
            time.push_back(t);
            field.push_back(c);
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
                                 const double& N_mum) {
        rep.push_back(r);
        time.push_back(t);
        field.push_back(c);
        line.push_back("");
        type.push_back("mummy");
        N.push_back(N_mum);
        return;
    }
};










void one_rep__(RepSummary& summary,
               std::vector<AllStageInfo>& stage_ts,
               AllFields& fields,
               const uint32& rep,
               const uint32& max_t,
               const uint32& save_every,
               const std::deque<PerturbInfo> perturbs_,
               RcppThread::ProgressBar& prog_bar,
               const bool& show_progress,
               const bool& stage_ts_out) {

    uint32 n_lines = fields.n_lines();
    uint32 n_fields = fields.size();

    summary.reserve(rep, max_t, save_every, n_lines, n_fields);

    if (stage_ts_out) stage_ts.reserve(max_t);

    std::deque<PerturbInfo> perturbs(perturbs_);

    summary.push_back(0, fields);

    for (uint32 t = 1; t <= max_t; t++) {

        // Basic updates for perturbations, dispersal, and population dynamics.
        // Returns true if all fields are empty.
        bool all_empty = fields.update(t, perturbs);

        if (t % save_every == 0 || t == max_t) summary.push_back(t, fields);

        // If all fields are empty, then stop this rep.
        if (all_empty) break;

        if (stage_ts_out) stage_ts.push_back(fields.out_all_info());

        if (show_progress) prog_bar++;

        if (t % 10 == 0) RcppThread::checkUserInterrupt();

    }


    return;

}






//[[Rcpp::export]]
List sim_pseudogameofclones_cpp(SEXP fields_ptr,
                                SEXP perturb_ptr,
                                const uint32& n_reps,
                                const uint32& max_t,
                                const uint32& save_every,
                                const bool& sep_adults,
                                uint32 n_threads,
                                const bool& show_progress) {

    /*
     These are the starting conditions for all fields.
     (Note that it's not a vector of AllFields objects like
      `all_fields_vec` below.)
     I don't want to change these in case the `AllFields` object in R
     is used again.
     */
    XPtr<AllFields> fields_xptr(fields_ptr);
    const AllFields& fields0(*fields_xptr);

    XPtr<std::deque<PerturbInfo>> perturb_xptr(perturb_ptr);
    const std::deque<PerturbInfo>& perturbs(*perturb_xptr);

    // Check that # threads isn't too high:
    thread_check(n_threads);

    RcppThread::ProgressBar prog_bar(n_reps * max_t, 1);

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
            new std::vector<AllFields>(n_reps, fields0), true);
    std::vector<AllFields>& all_fields_vec(*all_fields_vec_xptr);
    for (uint32 i = 0; i < n_reps; i++) {
        all_fields_vec[i].reseed(seeds[i]);
    }

    // Have to create this to make one_rep__ compatible with `restart_experiments`
    std::vector<std::vector<AllStageInfo>> stage_ts(n_reps);

    // Parallelized loop
    RcppThread::parallelFor(0, n_reps, [&] (uint32 i) {
        one_rep__(summaries[i], stage_ts[i], all_fields_vec[i], i,
                  max_t, save_every, perturbs,
                  prog_bar, show_progress, false);
    }, n_threads);


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
    if (n_lines == 0) stop("\nERROR: n_lines == 0\n");
    if (n_fields == 0) stop("\nERROR: n_fields == 0\n");
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
                             SEXP perturb_ptr,
                             const uint32& max_t,
                             const uint32& save_every,
                             const bool& stage_ts_out,
                             const bool& sep_adults,
                             const bool& show_progress,
                             uint32 n_threads) {


    XPtr<std::vector<AllFields>> all_fields_vec_xptr(all_fields_ptr);
    std::vector<AllFields>& all_fields_vec(*all_fields_vec_xptr);

    XPtr<std::deque<PerturbInfo>> perturb_xptr(perturb_ptr);
    const std::deque<PerturbInfo>& perturbs(*perturb_xptr);

    uint32 n_reps = all_fields_vec.size();

    // Check that # threads isn't too high:
    thread_check(n_threads);

    // No perturbations allowed here:
    RcppThread::ProgressBar prog_bar(max_t * n_reps, 1);

    std::vector<RepSummary> summaries(n_reps, RepSummary(sep_adults));

    std::vector<std::vector<AllStageInfo>> stage_ts(n_reps);

    // Parallelized loop
    RcppThread::parallelFor(0, n_reps, [&] (uint32 i) {
        one_rep__(summaries[i], stage_ts[i], all_fields_vec[i], i,
                  max_t, save_every, perturbs,
                  prog_bar, show_progress, stage_ts_out);
    }, n_threads);

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

