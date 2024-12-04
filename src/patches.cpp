
#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <functional>           // std::greater
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "pseudogameofclones_types.hpp"  // integer types
#include "patches.hpp"          // field and plant classes



using namespace Rcpp;









/*
 Adjust for potential extinction or re-colonization:
 */
void OnePlant::extinct_colonize(const uint32& i) {

    AphidPop& ap(aphids[i]);

    double N = ap.total_aphids();

    if (N < extinct_N) {
        ap.clear();
    } else {
        empty = false;
        ap.extinct = false;
    }

    return;
}




/*
 Iterate one time step.
 */
void OnePlant::update(const WaspPop* wasps,
                      pcg32& eng) {

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
        nm += aphids[i].update(this, wasps, eng);

        // Adjust for potential extinction or re-colonization:
        extinct_colonize(i);
    }

    mummies.update(pred_rate, nm);
    double mums = arma::accu(mummies.Y);
    if (mums < extinct_N) mummies.Y.fill(0);

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
            for (OnePlant& p : field.plants) {
                AphidPop& aphids(p.aphids[pert.index]);
                aphids.apterous.clear(mult);
                aphids.alates.clear(mult);
                if (aphids.total_aphids() < extinct_N) aphids.clear();
                if ((p.total_aphids() + p.total_mummies()) == 0) {
                    p.empty = true;
                }
            }
        } else if (pert.index == nl) { // mummies + parasitized aphids
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

    if (i > 0) {
        perturbs.erase(perturbs.begin(), perturbs.begin() + i);
    }

    return;
}







// Update for one time step.
// Returns true if all fields/plants are empty
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
            for (const OnePlant& p : field.plants) {
                if (!p.empty) {
                    all_empty = false;
                    break;
                }
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
        out.push_back(k+1, 0, "", "wasp", 1, field.wasps.Y);


        // Everything but adult wasps:
        for (uint32 j = 0; j < field.size(); j++) {

            const OnePlant& plant(field[j]);

            // Mummies:
            for (uint32 ii = 0; ii < plant.mummies.Y.n_elem; ii++) {
                out.push_back(k+1, j+1, "", "mummy", ii+1, plant.mummies.Y[ii]);
            }

            // Everything but mummies:
            for (uint32 i = 0; i < plant.size(); i++) {

                const AphidPop& aphid(plant[i]);

                for (uint32 ii = 0; ii < aphid.apterous.X.n_elem; ii++) {
                    out.push_back(k+1, j+1, aphid.aphid_name, "apterous",
                                  ii+1, aphid.apterous.X[ii]);
                }
                for (uint32 ii = 0; ii < aphid.alates.X.n_elem; ii++) {
                    out.push_back(k+1, j+1, aphid.aphid_name, "alate",
                                  ii+1, aphid.alates.X[ii]);
                }
                for (uint32 ii = 0; ii < aphid.paras.X.n_elem; ii++) {
                    out.push_back(k+1, j+1, aphid.aphid_name, "parasitized",
                                  ii+1, aphid.paras.X[ii]);
                }

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
        for (OnePlant& plant : field.plants) {
            for (double& y : plant.mummies.Y) {
                y = N[i];
                i++;
            }
            for (AphidPop& aphid : plant.aphids) {
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
    //' and a sum <= 0 are already checked for in sim_pseudogameofclones_cpp):
    double wfa_sum = std::accumulate(wasp_field_attract.begin(),
                                     wasp_field_attract.end(), 0.0);
    for (double& x : wasp_field_attract) x /= wfa_sum;

    for (uint32 i = 0; i < fields.size(); i++) {

        OneField& field(fields[i]);

        field.K_ = K_[i];
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

        for (uint32 j = 0; j < field.plants.size(); j++) {

            OnePlant& plant(field.plants[j]);

            plant.K = K_[i];
            plant.K_y = K_y_[i];
            plant.mummies.smooth = mum_smooth_;
            plant.pred_rate = pred_rate_[i];

            for (uint32 k = 0; k < plant.aphids.size(); k++) {

                AphidPop& aphid(plant.aphids[k]);
                aphid.apterous.alate_b0_ = alate_b0_[k];
                aphid.apterous.alate_b1_ = alate_b1_[k];

            }

        }

    }

    return;

}


