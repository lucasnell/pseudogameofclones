
#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <functional>           // std::greater
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "clonewars_types.hpp"  // integer types
#include "patches.hpp"          // field and plant classes
#include "math.hpp"             // combine_leslies fxn



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
 Carrying capacity for plant.
 It depends on the Leslie matrix for each line's apterous and alates.
 (I'm not including parasitized aphids because they shouldn't be too numerous.)
 The overall carrying capacity is weighted based on each line's abundance.
 */
double OnePlant::carrying_capacity() const {

    arma::vec cc(aphids.size());
    arma::vec Ns(aphids.size());
    double total_N = 0;

    arma::mat L;
    arma::cx_vec eigval;
    double ev;

    for (uint32 i = 0; i < aphids.size(); i++) {

        Ns[i] = aphids[i].total_aphids();

        total_N += Ns[i];

        combine_leslies(L,
                        aphids[i].apterous.leslie(),
                        aphids[i].alates.leslie(),
                        aphids[i].apterous.alate_prop(this),
                        aphids[i].alates.alate_plant_disp_p(),
                        aphids[i].alates.disp_mort(),
                        aphids[i].alates.disp_start());

        eigval = arma::eig_gen( L );
        ev = eigval.max().real();
        cc[i] = (ev - 1) * K;
    }

    double avg_cc;

    if (total_N > 0) {
        avg_cc = arma::accu(cc % Ns / total_N);
    } else avg_cc = arma::mean(cc);

    return avg_cc;
}





void OnePlant::update_z_wilted() {

    z = total_aphids();

    // Once it turns wilted, it stays that way until cleared
    if (wilted_) return;
    // Wilting process is ignored if wilted_prop > 1
    if (wilted_prop <= 1) {
        wilted_ = carrying_capacity() >= (K * wilted_prop);
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


        // Update population, including process error and dispersal.
        // Also return # newly mummified from that line
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
// Same but minus stochasticity
void OnePlant::update(const arma::cube& emigrants,
                      const arma::cube& immigrants,
                      const WaspPop* wasps) {

    update_z_wilted();

    S = 1 / (1 + z / K);
    S_y = 1 / (1 + z / K_y);

    empty = true;

    double nm = 0; // newly mummified

    for (uint32 i = 0; i < aphids.size(); i++) {

        nm += aphids[i].update(this, wasps, emigrants.slice(i).col(this_j),
                               immigrants.slice(i).col(this_j));

        if (wilted_) {
            aphids[i].apterous.X *= wilted_mort;
            aphids[i].alates.X *= wilted_mort;
            aphids[i].paras.X *= wilted_mort;
        }

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
