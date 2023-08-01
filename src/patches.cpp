
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
void OnePlant::update(const AllAphidPlantDisps& emigrants,
                      const AllAphidPlantDisps& immigrants,
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
        nm += aphids[i].update(this, wasps, emigrants[i], immigrants[i], eng);

        if (wilted_) {
            aphids[i].apterous.X *= (1 - wilted_mort);
            aphids[i].alates.X *= (1 - wilted_mort);
            aphids[i].paras.X *= (1 - wilted_mort);
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







void OneField::clear_plants(const double& clear_surv,
                            const uint32& min_rm,
                            pcg32& eng) {

    std::vector<WiltedPlantsInfo> wilted_plants;

    // We can remove at max half the plants:
    uint32 max_rm = std::max(1U, static_cast<uint32>(plants.size() / 2));
    if (min_rm > max_rm) max_rm = min_rm;

    wilted_plants.reserve(max_rm);

    for (const OnePlant& p : plants) {
        if (p.wilted()) {
            WiltedPlantsInfo tmp_wpi(p.this_j, p.total_aphids());
            wilted_plants.push_back(tmp_wpi);
        }
    }

    if (wilted_plants.size() == 0 && min_rm == 0) return;

    /*
     If we have greater wilted plants than the max, we need to remove
     some from the `wilted_plants` vector.
     */
    if (wilted_plants.size() > max_rm) {
        /*
         Sort by ascending `N` so that we clear the plants with the
         fewest aphids:
         */
        std::sort(wilted_plants.begin(), wilted_plants.end());
        while (wilted_plants.size() > max_rm) wilted_plants.pop_back();
    }
    /*
     If requested to remove a minimum number of plants from each field on
     this day but not enough are wilted, add non-wilted plants that have
     the fewest aphids.
     */
    if (wilted_plants.size() < min_rm) {
        // Pretend non-wilted plants are wilted to sort by `N` as I do when
        // `wilted_plants.size() > max_rm` above:
        std::deque<WiltedPlantsInfo> non_wilted;  // use deque for .pop_front()
        for (const OnePlant& p : plants) {
            if (!p.wilted()) {
                WiltedPlantsInfo tmp_wpi(p.this_j, p.total_aphids());
                non_wilted.push_back(tmp_wpi);
            }
        }
        // Now sort by ascending `N`
        std::sort(non_wilted.begin(), non_wilted.end());
        // Add non-wilted to the end of `wilted_plants` until `min_rm`
        // requirement is met
        while (wilted_plants.size() < min_rm && ! non_wilted.empty()) {
            wilted_plants.push_back(non_wilted.front());
            non_wilted.pop_front();
        }
    }

    double K, K_y;

    for (uint32 i = 0; i < wilted_plants.size(); i++) {

        OnePlant& plant_i(plants[wilted_plants[i].ind]);

        set_K(K, K_y, eng);

        if (clear_surv > 0) {
            plant_i.clear(K, K_y, clear_surv);
        } else plant_i.clear(K, K_y);

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

        field.calc_plant_dispersal();

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


// Clear plants
void AllFields::clear_plants(const uint32& t,
                             std::deque<uint32>& check_for_clear,
                             std::deque<std::pair<uint32, uint32>>& extra_plant_removals) {
    if (check_for_clear.empty()) return;
    if (t != check_for_clear.front()) return;
    check_for_clear.pop_front();
    uint32 min_rm = 0;
    if (!extra_plant_removals.empty() && t >= extra_plant_removals.front().first) {
        min_rm = extra_plant_removals.front().second;
        extra_plant_removals.pop_front();
    }
    for (OneField& field : fields) {
        field.clear_plants(clear_surv, min_rm, eng);
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
                             const std::vector<double>& pred_rate_) {

    this->alate_field_disp_p = alate_field_disp_p_;
    this->wasp_disp_m0 = wasp_disp_m0_;
    this->wasp_disp_m1 = wasp_disp_m1_;
    this->wasp_field_attract = wasp_field_attract_;
    //' Make sure `wasp_field_attract` sums to 1 (negative values
    //' and a sum <= 0 are already checked for in sim_gameofclones_cpp):
    double wfa_sum = std::accumulate(wasp_field_attract.begin(),
                                     wasp_field_attract.end(), 0.0);
    for (double& x : wasp_field_attract) x /= wfa_sum;

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


