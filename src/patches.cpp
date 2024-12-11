
#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <functional>           // std::greater
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "pseudogameofclones_types.hpp"  // integer types
#include "patches.hpp"          // field classes



using namespace Rcpp;






/*
 Iterate one time step, not including alate and adult-wasp dispersal.
 */
void OneField::update(pcg32& eng) {

    wasps.x = total_unpar_aphids();
    double old_mums = wasps.mummies.Y.back();

    z = total_aphids();

    S = 1 / (1 + z / K);
    S_y = 1 / (1 + z / K_y);

    empty = true;

    double nm = 0; // newly mummified

    for (uint32 i = 0; i < aphids.size(); i++) {

        /*
         Update population, including dispersal and (optionally) process error.
         Also return # newly mummified from that line
         */
        nm += aphids[i].update(this, &wasps, eng);

        // Adjust for potential extinction or re-colonization:
        extinct_colonize(i);
    }

    wasps.mummies.update(pred_rate, nm);
    double mums = arma::accu(wasps.mummies.Y);
    if (mums < extinct_N) wasps.mummies.Y.fill(0);

    // Lastly update adult wasps (if constant_wasps = false):
    if (!constant_wasps) {
        wasps.update(old_mums, eng);
        if (wasps.Y < extinct_N) wasps.Y = 0;
    }

    return;

}



// Remove aphid dispersers from this field:
arma::mat OneField::remove_field_dispersers(const double& disp_prop) {

    uint32 n_lines = aphids.size();
    uint32 n_stages = aphids[0].alates.X.n_elem;

    arma::mat D = arma::zeros<arma::mat>(n_stages, n_lines);

    arma::vec Di;
    for (uint32 i = 0; i < n_lines; i++) {
        Di = aphids[i].remove_field_dispersers(disp_prop);
        D.col(i) += Di;
    }

    return D;

}

// Add aphid dispersers from another field:
void OneField::add_field_dispersers(const arma::mat& D) {

    uint32 n_lines = aphids.size();

    for (uint32 i = 0; i < n_lines; i++) {
        aphids[i].alates.X += D.col(i);
    }

    return;

}








void AllFields::do_perturb(std::deque<PerturbInfo>& perturbs,
                           const uint32& t) {

    if (perturbs.empty()) return;

    uint32 nl = n_lines();

    uint32 i = 0;
    while (i < perturbs.size()) {
        const PerturbInfo& pert(perturbs[i]);
        if (pert.time != t) break;
        double mult = pert.multiplier;
        OneField& field(fields[pert.field]);
        if (pert.index < nl) { // aphids (non-parasitized)
            AphidPop& aphids(field.aphids[pert.index]);
            aphids.apterous.clear(mult);
            aphids.alates.clear(mult);
            if (aphids.total_aphids() < extinct_N) aphids.clear();
            if ((field.total_aphids() + field.total_mummies()) == 0) {
                field.empty = true;
            }
        } else if (pert.index == nl) { // mummies + parasitized aphids
            field.wasps.mummies.clear(mult);
            if (arma::accu(field.wasps.mummies.Y) < extinct_N) field.wasps.mummies.Y.fill(0);
            for (AphidPop& aphids : field.aphids) {
                aphids.paras.clear(mult);
                if (aphids.total_aphids() < extinct_N) aphids.clear();
            }
            if ((field.total_aphids() + field.total_mummies()) == 0) {
                field.empty = true;
            }
        } else { // adult wasps
            field.wasps.Y *= mult;
            if (field.wasps.Y < extinct_N) field.wasps.Y = 0;
        }
        i++;
    }

    if (i > 0) {
        perturbs.erase(perturbs.begin(), perturbs.begin() + i);
    }

    return;
}







// Update for one time step.
// Returns true if all fields are empty
bool AllFields::update(const uint32& t,
                       std::deque<PerturbInfo>& perturbs) {

    do_perturb(perturbs, t);

    // Alate and wasp dispersal across fields:
    across_field_disp_alates();
    across_field_disp_wasps();

    bool all_empty = true;

    for (OneField& field : fields) {

        field.wasps.add_Y_0_check(t);

        field.update(eng);

        if (all_empty) {
            if (!field.empty) {
                all_empty = false;
                break;
            }
        }

    }

    return all_empty;
}








AllStageInfo AllFields::out_all_info() const {

    AllStageInfo out;
    out.reserve(total_stages);

    for (uint32 k = 0; k < fields.size(); k++) {

        const OneField& field(fields[k]);

        // Adult wasps:
        out.push_back(k+1, "", "wasp", 1, field.wasps.Y);

        // Mummies:
        for (uint32 i = 0; i < field.wasps.mummies.Y.n_elem; i++) {
            out.push_back(k+1, "", "mummy", i+1, field.wasps.mummies.Y[i]);
        }

        // Everything else:
        for (uint32 i = 0; i < field.size(); i++) {

            const AphidPop& aphid(field[i]);

            for (uint32 j = 0; j < aphid.apterous.X.n_elem; j++) {
                out.push_back(k+1, aphid.aphid_name, "apterous",
                              j+1, aphid.apterous.X[j]);
            }
            for (uint32 j = 0; j < aphid.alates.X.n_elem; j++) {
                out.push_back(k+1, aphid.aphid_name, "alate",
                              j+1, aphid.alates.X[j]);
            }
            for (uint32 j = 0; j < aphid.paras.X.n_elem; j++) {
                out.push_back(k+1, aphid.aphid_name, "parasitized",
                              j+1, aphid.paras.X[j]);
            }

        }

    }

    return out;
}



// It's assumed this vector is in the same order as the outputs are above!
// Make sure this happens from the R side.
void AllFields::from_vector(std::vector<double>& N) {

    if (N.size() != total_stages) {
        std::string err_msg("ERROR:\n");
        err_msg += "The vector used to update internal C++ object ";
        err_msg += "`AllFields` is not of the required size (";
        err_msg += std::to_string(total_stages);
        err_msg += "). This is perhaps because you removed row(s) from ";
        err_msg += "the dataframe output from here.";
        Rcpp::stop(err_msg.c_str());
    }


    uint32 i = 0;
    for (OneField& field : fields) {
        field.wasps.Y = N[i];
        i++;
        for (double& y : field.wasps.mummies.Y) {
            y = N[i];
            i++;
        }
        for (AphidPop& aphid : field.aphids) {
            for (double& x : aphid.apterous.X) {
                x = N[i];
                i++;
            }
            for (double& x : aphid.alates.X) {
                x = N[i];
                i++;
            }
            for (double& x : aphid.paras.X) {
                x = N[i];
                i++;
            }
        }
    }

    return;

}













//[[Rcpp::export]]
List fields_to_data_frames(SEXP all_fields_ptr) {

    XPtr<std::vector<AllFields>> all_fields_vec_xptr(all_fields_ptr);
    const std::vector<AllFields>& all_fields_vec(*all_fields_vec_xptr);

    uint32 n_reps = all_fields_vec.size();

    List out(n_reps);

    for (uint32 i = 0; i < n_reps; i++) {
        const AllFields& all_fields(all_fields_vec[i]);
        AllStageInfo all_info = all_fields.out_all_info();
        out[i] = all_info.to_data_frame();
    }

    return out;

}



//[[Rcpp::export]]
void fields_from_vectors(SEXP all_fields_ptr,
                         std::vector<std::vector<double>>& N_vecs) {

    XPtr<std::vector<AllFields>> all_fields_vec_xptr(all_fields_ptr);
    std::vector<AllFields>& all_fields_vec(*all_fields_vec_xptr);

    uint32 n_reps = all_fields_vec.size();

    if (N_vecs.size() != n_reps) stop("`N_vecs` not correct size.");

    for (uint32 i = 0; i < n_reps; i++) {
        AllFields& all_fields(all_fields_vec[i]);
        all_fields.from_vector(N_vecs[i]);
    }

    return;

}





void AllFields::set_new_pars(const std::vector<double>& K_,
                             const std::vector<double>& alate_b0_,
                             const std::vector<double>& alate_b1_,
                             const double& alate_field_disp_p_,
                             const std::vector<double>& K_y_,
                             const std::vector<double>& s_y_,
                             const double& a_,
                             const double& k_,
                             const double& h_,
                             const double& wasp_disp_m0_,
                             const double& wasp_disp_m1_,
                             const std::vector<double>& wasp_field_attract_,
                             const double& mum_smooth_,
                             const std::vector<double>& pred_rate_) {

    this->alate_field_disp_p = alate_field_disp_p_;
    this->wasp_disp_m0 = wasp_disp_m0_;
    this->wasp_disp_m1 = wasp_disp_m1_;
    this->wasp_field_attract = wasp_field_attract_;
    //' Make sure `wasp_field_attract` sums to 1 (negative values
    //' and a sum <= 0 are already checked for in sim_fields_cpp):
    double wfa_sum = std::accumulate(wasp_field_attract.begin(),
                                     wasp_field_attract.end(), 0.0);
    for (double& x : wasp_field_attract) x /= wfa_sum;

    for (uint32 i = 0; i < fields.size(); i++) {

        OneField& field(fields[i]);

        field.K = K_[i];
        field.K_y = K_y_[i];
        field.wasps.s_y = s_y_[i];
        field.wasps.attack.a = a_;
        field.wasps.attack.k = k_;
        field.wasps.attack.h = h_;
        /*
         Since we're restarting this simulation, we don't want wasps to be
         added when we would normally have added them in the original
         simulation:
         */
        field.wasps.Y_0 = 0;

        field.wasps.mummies.smooth = mum_smooth_;
        field.pred_rate = pred_rate_[i];

        for (uint32 k = 0; k < field.aphids.size(); k++) {

            AphidPop& aphid(field.aphids[k]);
            aphid.apterous.alate_b0_ = alate_b0_[k];
            aphid.apterous.alate_b1_ = alate_b1_[k];

        }

    }

    return;

}





//[[Rcpp::export]]
SEXP make_field_ptr(const bool& aphid_demog_error,
                    const arma::mat& aphid_density_0,
                    const bool& wasp_demog_error,
                    const std::vector<double>& wasp_density_0,
                    const std::vector<uint32>& wasp_delay,
                    const std::vector<double>& mummy_density_0,
                    const bool& environ_error,
                    SEXP aphids_ptr,
                    SEXP wasp_ptr,
                    const uint32& n_fields,
                    const std::vector<double>& K,
                    const std::vector<double>& K_y,
                    const std::vector<double>& pred_rate,
                    const double& extinct_N,
                    const std::vector<bool>& constant_wasps,
                    const double& alate_field_disp_p,
                    const double& wasp_disp_m0,
                    const double& wasp_disp_m1,
                    const std::vector<double>& wasp_field_attract,
                    const arma::vec& new_rel_attack) {


    if (aphid_density_0.n_rows != n_fields)
        stop("aphid_density_0.n_rows != n_fields");
    if (wasp_density_0.size() != n_fields)
        stop("wasp_density_0.size() != n_fields");
    if (wasp_delay.size() != n_fields)
        stop("wasp_delay.size() != n_fields");
    if (mummy_density_0.size() != n_fields)
        stop("mummy_density_0.size() != n_fields");
    if (K.size() != n_fields)
        stop("K.size() != n_fields");
    if (K_y.size() != n_fields)
        stop("K_y.size() != n_fields");
    if (pred_rate.size() != n_fields)
        stop("pred_rate.size() != n_fields");
    if (constant_wasps.size() != n_fields)
        stop("constant_wasps.size() != n_fields");
    if (wasp_field_attract.size() != n_fields)
        stop("wasp_field_attract.size() != n_fields");

    XPtr<std::vector<AphidPop>> aphids_xptr(aphids_ptr);
    const std::vector<AphidPop>& aphids(*aphids_xptr);

    if (new_rel_attack.n_elem != aphids[0].apterous.leslie().n_rows) {
        stop("`new_rel_attack` is wrong length");
    }

    uint32 n_lines = aphids.size();

    if (aphid_density_0.n_cols != n_lines)
        stop("aphid_density_0.n_cols != n_lines");

    XPtr<WaspPop> wasp_xptr(wasp_ptr);
    const WaspPop& wasps(*wasp_xptr);

    XPtr<AllFields> fields_xptr(
            new AllFields(n_fields, aphids, wasps, K, K_y, pred_rate,
                          extinct_N, constant_wasps, alate_field_disp_p,
                          wasp_disp_m0, wasp_disp_m1, wasp_field_attract), true);
    AllFields& fields(*fields_xptr);

    for (uint32 i = 0; i < n_fields; i++) {

        OneField& field(fields[i]);
        for (uint32 j = 0; j < n_lines; j++) {
            field[j].adjust_error_X0(aphid_demog_error,
                                     environ_error,
                                     aphid_density_0(i,j));
        }

        field.wasps.adjust_error_Y0(wasp_demog_error, environ_error, wasp_delay[i],
                                    wasp_density_0[i], mummy_density_0[i]);
        field.wasps.attack.rel_attack = new_rel_attack;

    }


    return fields_xptr;

}

/*
 Extract the following elements from an AllFields object (for use in
 `restarted_field_ptr` below):
    - aphid_demog_error
    - wasp_demog_error
    - K
    - K_y
    - pred_rate
    - constant_wasps
    - alate_field_disp_p
    - wasp_disp_m0
    - wasp_disp_m1
    - wasp_field_attract

 It's assumed that all reps are the same.

*/
//[[Rcpp::export]]
List get_field_pars(SEXP all_fields_in_ptr) {

    XPtr<std::vector<AllFields>> all_fields_vec_in_xptr(all_fields_in_ptr);
    if (all_fields_vec_in_xptr->empty()) stop("\nn_reps == 0\n");
    const AllFields& repr_fields((*all_fields_vec_in_xptr)[0]);

    uint32 n_fields = repr_fields.size();
    if (n_fields == 0) stop("\nn_fields == 0\n");

    // parameters that have one value per `AllFields` object (or are defined
    //  at this level):
    bool aphid_demog_error = repr_fields[0].aphids[0].demog_error;
    bool wasp_demog_error = repr_fields[0].wasps.demog_error;
    double alate_field_disp_p = repr_fields.alate_field_disp_p;
    double wasp_disp_m0 = repr_fields.wasp_disp_m0;
    double wasp_disp_m1 = repr_fields.wasp_disp_m1;
    std::vector<double> wasp_field_attract = repr_fields.wasp_field_attract;

    // parameters that can vary by field:
    std::vector<double> K;
    std::vector<double> K_y;
    std::vector<double> pred_rate;
    std::vector<bool> constant_wasps;
    K.reserve(n_fields);
    K_y.reserve(n_fields);
    pred_rate.reserve(n_fields);
    constant_wasps.reserve(n_fields);
    for (const OneField& field : repr_fields) {
        K.push_back(field.K);
        K_y.push_back(field.K_y);
        pred_rate.push_back(field.pred_rate);
        constant_wasps.push_back(field.constant_wasps);
    }



    List out = List::create(_["aphid_demog_error"] = aphid_demog_error,
                            _["wasp_demog_error"] = wasp_demog_error,
                            _["K"] = K,
                            _["K_y"] = K_y,
                            _["pred_rate"] = pred_rate,
                            _["constant_wasps"] = constant_wasps,
                            _["alate_field_disp_p"] = alate_field_disp_p,
                            _["wasp_disp_m0"] = wasp_disp_m0,
                            _["wasp_disp_m1"] = wasp_disp_m1,
                            _["wasp_field_attract"] = wasp_field_attract);

    return out;

}



// note #2: this differs from `make_field_ptr` in that this creates a vector
// of AllFields objects, one per simulation rep
//
//[[Rcpp::export]]
SEXP restarted_field_ptr(SEXP all_fields_in_ptr,
                         const bool& aphid_demog_error,
                         const bool& wasp_demog_error,
                         const std::vector<double>& K,
                         const std::vector<double>& K_y,
                         const std::vector<double>& pred_rate,
                         const std::vector<bool>& constant_wasps,
                         const double& alate_field_disp_p,
                         const double& wasp_disp_m0,
                         const double& wasp_disp_m1,
                         std::vector<double> wasp_field_attract) {

    XPtr<std::vector<AllFields>> all_fields_vec_in_xptr(all_fields_in_ptr);
    const std::vector<AllFields>& all_fields_vec_in(*all_fields_vec_in_xptr);

    uint32 n_reps = all_fields_vec_in.size();
    if (n_reps == 0) stop("\nn_reps == 0\n");

    // copy input XPtr to a new XPtr
    XPtr<std::vector<AllFields>> all_fields_vec_xptr(
            new std::vector<AllFields>(), true);
    std::vector<AllFields>& all_fields_vec(*all_fields_vec_xptr);
    all_fields_vec.reserve(n_reps);
    for (const AllFields& in_fields : all_fields_vec_in){
        all_fields_vec.push_back(in_fields);
    }

    // For checking dimensions of vectors
    uint32 n_lines = all_fields_vec[0].n_lines();
    uint32 n_fields = all_fields_vec[0].size();
    if (n_lines == 0) stop("\nn_lines == 0\n");
    if (n_fields == 0) stop("\nn_fields == 0\n");

    // This vector can't have values < 0:
    double wfa_min = *min_element(wasp_field_attract.begin(),
                                  wasp_field_attract.end());
    if (wfa_min < 0) stop("\nmin(wasp_field_attract) < 0.\n");
    // It can have zeros, but can't have all zeros:
    double wfa_sum = std::accumulate(wasp_field_attract.begin(),
                                     wasp_field_attract.end(), 0.0);
    if (wfa_sum <= 0) stop("\nwasp_field_attract sums to <= 0.\n");
    for (double& x : wasp_field_attract) x /= wfa_sum;


    // ------------*
    // Now fill parameters
    for (AllFields& fields : all_fields_vec) {
        for (uint32 i = 0; i < fields.size(); i++) {
            OneField& field(fields[i]);
            for (AphidPop& a : field) a.demog_error = aphid_demog_error;
            field.wasps.demog_error = wasp_demog_error;
            field.K = K[i];
            field.K_y = K_y[i];
            field.pred_rate = pred_rate[i];
            field.constant_wasps = constant_wasps[i];
        }
        fields.alate_field_disp_p = alate_field_disp_p;
        fields.wasp_disp_m0 = wasp_disp_m0;
        fields.wasp_disp_m1 = wasp_disp_m1;
        fields.wasp_field_attract = wasp_field_attract;
    }


    return all_fields_vec_xptr;

}







//[[Rcpp::export]]
SEXP make_perturb_ptr(const std::vector<uint32>& perturb_when,
                      const std::vector<uint32>& perturb_where,
                      const std::vector<uint32>& perturb_who,
                      const std::vector<double>& perturb_how) {

    uint32 n_perturbs = perturb_when.size();
    if (perturb_where.size() != n_perturbs)
        stop("perturb_where.size() != perturb_when.size()");
    if (perturb_who.size() != n_perturbs)
        stop("perturb_who.size() != perturb_when.size()");
    if (perturb_how.size() != n_perturbs)
        stop("perturb_how.size() != perturb_when.size()");

    XPtr<std::deque<PerturbInfo>> perturbs_xptr(
            new std::deque<PerturbInfo>(n_perturbs), true);
    std::deque<PerturbInfo>& perturbs(*perturbs_xptr);

    for (uint32 i = 0; i < n_perturbs; i++) {
        perturbs[i] = PerturbInfo(perturb_when[i], perturb_where[i],
                                  perturb_who[i], perturb_how[i]);
    }

    return perturbs_xptr;

}
