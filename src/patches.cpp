
#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <functional>           // std::greater
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "gameofclones_types.hpp"  // integer types
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





void OnePlant::update_z_wilted() {

    z = total_aphids();

    // Once it turns wilted, it stays that way until cleared
    if (wilted_) return;
    // Wilting process is ignored if wilted_N == 0
    if (wilted_N > 0) {
        wilted_ = z > wilted_N;
    }

    return;
}



/*
 Iterate one time step, after calculating dispersal numbers
 */
void OnePlant::update(const arma::cube& emigrants,
                      const arma::cube& immigrants,
                      const WaspPop* wasps,
                      pcg32& eng) {

    update_z_wilted();

    S = 1 / (1 + z / K);
    S_y = 1 / (1 + z / K_y);

    empty = true;

    double nm = 0; // newly mummified

    for (uint32 i = 0; i < aphids.size(); i++) {


        /*
         Update population, including dispersal and (optionally) process error.
         Also return # newly mummified from that line
         */
        nm += aphids[i].update(this, wasps, emigrants.slice(i).col(this_j),
                               immigrants.slice(i).col(this_j), eng);

        if (wilted_) {
            aphids[i].apterous.X *= wilted_mort;
            aphids[i].alates.X *= wilted_mort;
            aphids[i].paras.X *= wilted_mort;
        }

        // Adjust for potential extinction or re-colonization:
        extinct_colonize(i);

    }

    mummies.update(pred_rate, nm);
    double mums = arma::accu(mummies.Y);
    if (mums < extinct_N) mummies.Y.fill(0);
    if (max_mum_density > 0 && mums > max_mum_density) {
        mummies.Y *= (max_mum_density / mums);
    }

    age++;

    return;

}




// Do the actual clearing of plants while avoiding extinction
template <typename T>
inline void OneField::do_clearing(std::vector<PlantClearingInfo<T>>& clear_plants,
                                    double& remaining,
                                    std::vector<bool>& wilted,
                                    const double& clear_surv,
                                    pcg32& eng) {

    int n_plants = plants.size();
    int n_wilted = std::accumulate(wilted.begin(), wilted.end(), 0);

    double K, K_y, wilted_mort;

    if (clear_plants.size() == 0 && n_wilted < n_plants) return;

    /*
     This is to avoid removing no plants when they're all dying.
     Clearing the least-populated plant here means that there's a healthy plant
     for aphids to disperse to, which hopefully reduces the chance of extinction
     before the next time plants get checked.
     */
    if (clear_plants.size() == 0 && n_wilted == n_plants) {
        for (const OnePlant& p : plants) {
            double N = p.total_aphids();
            clear_plants.push_back(PlantClearingInfo<T>(p.this_j, N, true,
                                                         static_cast<T>(N)));
        }
        // Sort ascending by `N`:
        std::sort(clear_plants.begin(), clear_plants.end());

        remaining = 0;

        for (uint32 i = 1; i < n_plants; i++) {
            remaining += clear_plants.back().N;
            clear_plants.pop_back();
        }


        set_K(K, K_y, eng);
        set_wilted_mort(wilted_mort, eng);

        if (clear_surv > 0) {
            plants[clear_plants.front().ind].clear(K, K_y, wilted_mort,
                                        clear_surv);
        } else plants[clear_plants.front().ind].clear(K, K_y, wilted_mort);

        return;

    }

    /*
     If removing all these plants cause total extinction, we need to remove
     some indices from `inds`.
     We'll sort (high to low) the vector by the threshold variable (age or # aphids),
     then remove the last item in the vector until we have aphids again.
     */
    if (remaining == 0) {

        // Sort descending:
        std::sort(clear_plants.begin(), clear_plants.end(),
                  std::greater<PlantClearingInfo<T>>());

        while (remaining == 0) {
            remaining += clear_plants.back().N;
            clear_plants.pop_back();
        }

    }

    for (uint32 i = 0; i < clear_plants.size(); i++) {

        set_K(K, K_y, eng);
        set_wilted_mort(wilted_mort, eng);

        if (clear_surv > 0) {
            plants[clear_plants.front().ind].clear(K, K_y, wilted_mort,
                                        clear_surv);
        } else plants[clear_plants.front().ind].clear(K, K_y, wilted_mort);

    }


    return;
}




// Clear plants by either a maximum age or total abundance
void OneField::clear_plants(const uint32& max_age,
                               const double& clear_surv,
                               pcg32& eng) {

    std::vector<PlantClearingInfo<uint32>> clear_plants;
    clear_plants.reserve(plants.size());

    // # remaining aphids for non-cleared plants
    double remaining = 0;

    // keeping track of dying plants:
    std::vector<bool> wilted;
    wilted.reserve(plants.size());

    for (const OnePlant& p : plants) {
        double N = p.total_aphids();
        wilted.push_back(p.wilted());
        if (p.age > max_age || p.wilted()) {
            double surv_remaining = N * clear_surv;
            if (surv_remaining < extinct_N) surv_remaining = 0;
            remaining += surv_remaining;
            clear_plants.push_back(PlantClearingInfo<uint32>(p.this_j, N, wilted.back(),
                                                              p.age));
        } else remaining += N;
    }

    do_clearing<uint32>(clear_plants, remaining, wilted, clear_surv, eng);

    return;
}


void OneField::clear_plants(const double& max_N,
                               const double& clear_surv,
                               pcg32& eng) {

    std::vector<PlantClearingInfo<double>> clear_plants;
    clear_plants.reserve(plants.size());

    double remaining = 0;

    std::vector<bool> wilted;
    wilted.reserve(plants.size());

    for (const OnePlant& p : plants) {
        double N = p.total_aphids();
        wilted.push_back(p.wilted());
        if (N > max_N || p.wilted()) {
            double surv_remaining = N * clear_surv;
            if (surv_remaining < extinct_N) surv_remaining = 0;
            remaining += surv_remaining;
            clear_plants.push_back(PlantClearingInfo<double>(p.this_j, N, wilted.back(),
                                                              N));
        } else remaining += N;
    }

    do_clearing<double>(clear_plants, remaining, wilted, clear_surv, eng);

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
                       std::deque<PerturbInfo>& perturbs,
                       std::deque<uint32>& check_for_clear) {

    do_perturb(perturbs, t);

    // Alate and wasp dispersal across fields:
    if (!check_for_clear.empty() && t == check_for_clear.front()) {
        across_field_disp_alates();
        across_field_disp_wasps();
    }

    bool all_empty = true;

    for (OneField& field : fields) {

        field.wasps.add_Y_0_check(t);

        if (disp_error) {
            field.calc_dispersal(eng);
        } else field.calc_dispersal();

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


// Clear plants by age or total N thresholds (or neither)
void AllFields::clear_plants(const uint32& t,
                             std::deque<uint32>& check_for_clear) {
    if (check_for_clear.empty()) return;
    if (t != check_for_clear.front()) return;
    check_for_clear.pop_front();
    if (max_N > 0) {
        for (OneField& field : fields) {
            field.clear_plants(max_N, clear_surv, eng);
        }
    } else if (max_age > 0) {
        for (OneField& field : fields) {
            field.clear_plants(max_age, clear_surv, eng);
        }
    }
    return;
}






List AllFields::to_list() const {

    std::vector<uint32> field_out;
    std::vector<uint32> plant_out;
    std::vector<std::string> line_out;
    std::vector<std::string> type_out;
    std::vector<uint32> stage_out;
    std::vector<double> N_out;

    field_out.reserve(total_stages);
    plant_out.reserve(total_stages);
    line_out.reserve(total_stages);
    type_out.reserve(total_stages);
    stage_out.reserve(total_stages);
    N_out.reserve(total_stages);


    for (uint32 k = 0; k < fields.size(); k++) {

        const OneField& field(fields[k]);

        // Adult wasps:
        field_out.push_back(k+1);
        plant_out.push_back(0);
        line_out.push_back("");
        type_out.push_back("wasp");
        stage_out.push_back(1);
        N_out.push_back(field.wasps.Y);

        // Everything but adult wasps:
        for (uint32 j = 0; j < field.size(); j++) {

            const OnePlant& plant(field[j]);

            // Mummies:
            for (uint32 ii = 0; ii < plant.mummies.Y.n_elem; ii++) {
                field_out.push_back(k+1);
                plant_out.push_back(j+1);
                line_out.push_back("");
                type_out.push_back("mummy");
                stage_out.push_back(ii+1);
                N_out.push_back(plant.mummies.Y[ii]);
            }

            // Everything but mummies:
            for (uint32 i = 0; i < plant.size(); i++) {

                const AphidPop& aphid(plant[i]);

                uint32 n_stages = aphid.apterous.X.n_elem +
                    aphid.alates.X.n_elem +
                    aphid.paras.X.n_elem;

                for (uint32 ii = 0; ii < n_stages; ii++) {
                    field_out.push_back(k+1);
                    plant_out.push_back(j+1);
                    line_out.push_back(aphid.aphid_name);
                }

                for (uint32 ii = 0; ii < aphid.apterous.X.n_elem; ii++) {
                    type_out.push_back("apterous");
                    stage_out.push_back(ii+1);
                    N_out.push_back(aphid.apterous.X[ii]);
                }
                for (uint32 ii = 0; ii < aphid.alates.X.n_elem; ii++) {
                    type_out.push_back("alate");
                    stage_out.push_back(ii+1);
                    N_out.push_back(aphid.alates.X[ii]);
                }
                for (uint32 ii = 0; ii < aphid.paras.X.n_elem; ii++) {
                    type_out.push_back("parasitized");
                    stage_out.push_back(ii+1);
                    N_out.push_back(aphid.paras.X[ii]);
                }

            }

        }

    }

    List out = List::create(
        _["field"] = field_out,
        _["plant"] = plant_out,
        _["line"] = line_out,
        _["type"] = type_out,
        _["stage"] = stage_out,
        _["N"] = N_out);

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
List fields_to_list(SEXP all_fields_ptr) {

    XPtr<std::vector<AllFields>> all_fields_vec_xptr(all_fields_ptr);
    const std::vector<AllFields>& all_fields_vec(*all_fields_vec_xptr);

    uint32 n_reps = all_fields_vec.size();

    List out(n_reps);

    for (uint32 i = 0; i < n_reps; i++) {
        const AllFields& all_fields(all_fields_vec[i]);
        out[i] = all_fields.to_list();
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





void AllFields::set_new_pars(const double& K_,
                             const std::vector<double>& alate_b0_,
                             const std::vector<double>& alate_b1_,
                             const double& alate_field_disp_p_,
                             const std::vector<double>& K_y_mult_,
                             const std::vector<double>& s_y_,
                             const double& a_,
                             const double& k_,
                             const double& h_,
                             const double& wasp_disp_m0_,
                             const double& wasp_disp_m1_,
                             const std::vector<double>& wasp_field_attract_,
                             const double& mum_smooth_,
                             const std::vector<double>& pred_rate_,
                             const uint32& max_plant_age_,
                             const double& clear_surv_) {

    this->alate_field_disp_p = alate_field_disp_p_;
    this->wasp_disp_m0 = wasp_disp_m0_;
    this->wasp_disp_m1 = wasp_disp_m1_;
    this->wasp_field_attract = wasp_field_attract_;
    //' Make sure `wasp_field_attract` sums to 1 (negative values
    //' and a sum <= 0 are already checked for in sim_gameofclones_cpp):
    double wfa_sum = std::accumulate(wasp_field_attract.begin(),
                                     wasp_field_attract.end(), 0.0);
    for (double& x : wasp_field_attract) x /= wfa_sum;
    this->max_age = max_plant_age_;
    this->clear_surv = clear_surv_;

    for (uint32 i = 0; i < fields.size(); i++) {

        OneField& field(fields[i]);

        field.mean_K_ = K_;
        field.K_y_mult = K_y_mult_[i];
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

            plant.K = K_;
            plant.K_y = K_ * K_y_mult_[i];
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


