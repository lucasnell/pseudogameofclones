
#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <functional>           // std::greater
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "clonewars_types.hpp"  // integer types
#include "patches.hpp"          // patch classes
#include "math.hpp"             // combine_leslies fxn
#include "pcg.hpp"              // runif_ fxns



using namespace Rcpp;









/*
 Adjust for potential extinction or re-colonization:
 */
void OnePatch::extinct_colonize(const uint32& i) {

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
 Carrying capacity for patch.
 It depends on the Leslie matrix for each line's apterous and alates.
 (I'm not including parasitized aphids because they shouldn't be too numerous.)
 The overall carrying capacity is weighted based on each line's abundance.
 */
double OnePatch::carrying_capacity() const {

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
                        aphids[i].alates.disp_rate(),
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





void OnePatch::update_z_wilted() {

    z = total_aphids();

    // Once it turns wilted, it stays that way until cleared
    if (wilted_) return;
    wilted_ = carrying_capacity() >= (K * death_prop);

    return;
}



/*
 Iterate one time step, after calculating dispersal numbers
 */
void OnePatch::update(const arma::cube& emigrants,
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
            aphids[i].apterous.X *= death_mort;
            aphids[i].alates.X *= death_mort;
            aphids[i].paras.X *= death_mort;
        }

        // Adjust for potential extinction or re-colonization:
        extinct_colonize(i);

    }

    mummies.update(pred_rate, nm);
    if (arma::accu(mummies.Y) < extinct_N) mummies.Y.fill(0);

    age++;

    return;

}
// Same but minus stochasticity
void OnePatch::update(const arma::cube& emigrants,
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
            aphids[i].apterous.X *= death_mort;
            aphids[i].alates.X *= death_mort;
            aphids[i].paras.X *= death_mort;
        }

        extinct_colonize(i);

    }

    mummies.update(pred_rate, nm);
    if (arma::accu(mummies.Y) < extinct_N) mummies.Y.fill(0);

    age++;

    return;

}





// Do the actual clearing of patches while avoiding extinction
template <typename T>
inline void AllPatches::do_clearing(std::vector<PatchClearingInfo<T>>& clear_patches,
                                    double& remaining,
                                    std::vector<bool>& wilted,
                                    const double& clear_surv,
                                    pcg32& eng) {

    int n_patches = patches.size();
    int n_wilted = std::accumulate(wilted.begin(), wilted.end(), 0);

    double K, K_y, death_mort;

    if (clear_patches.size() == 0 && n_wilted < n_patches) return;

    /*
     This is to avoid removing no plants when they're all dying.
     Clearing the least-populated plant here means that there's a healthy plant
     for aphids to disperse to, which hopefully reduces the chance of extinction
     before the next time plants get checked.
     */
    if (clear_patches.size() == 0 && n_wilted == n_patches) {
        for (const OnePatch& p : patches) {
            double N = p.total_aphids();
            clear_patches.push_back(PatchClearingInfo<T>(p.this_j, N, true,
                                                         static_cast<T>(N)));
        }
        // Sort ascending by `N`:
        std::sort(clear_patches.begin(), clear_patches.end());

        remaining = 0;

        for (uint32 i = 1; i < n_patches; i++) {
            remaining += clear_patches.back().N;
            clear_patches.pop_back();
        }


        set_K(K, K_y, eng);
        set_death_mort(death_mort, eng);

        if (clear_surv > 0) {
            patches[clear_patches.front().ind].clear(K, K_y, death_mort,
                                        clear_surv);
        } else patches[clear_patches.front().ind].clear(K, K_y, death_mort);

        return;

    }

    /*
     If removing all these patches cause total extinction, we need to remove
     some indices from `inds`.
     We'll sort (high to low) the vector by the threshold variable (age or # aphids),
     then remove the last item in the vector until we have aphids again.
     */
    if (remaining == 0) {

        // Sort descending:
        std::sort(clear_patches.begin(), clear_patches.end(),
                  std::greater<PatchClearingInfo<T>>());

        while (remaining == 0) {
            remaining += clear_patches.back().N;
            clear_patches.pop_back();
        }

    }

    for (uint32 i = 0; i < clear_patches.size(); i++) {

        set_K(K, K_y, eng);
        set_death_mort(death_mort, eng);

        if (clear_surv > 0) {
            patches[clear_patches.front().ind].clear(K, K_y, death_mort,
                                        clear_surv);
        } else patches[clear_patches.front().ind].clear(K, K_y, death_mort);

    }


    return;
}




// Clear patches by either a maximum age or total abundance
void AllPatches::clear_patches(const uint32& max_age,
                               const double& clear_surv,
                               pcg32& eng) {

    std::vector<PatchClearingInfo<uint32>> clear_patches;
    clear_patches.reserve(patches.size());

    // # remaining aphids for non-cleared patches
    double remaining = 0;

    // keeping track of dying plants:
    std::vector<bool> wilted;
    wilted.reserve(patches.size());

    for (const OnePatch& p : patches) {
        double N = p.total_aphids();
        wilted.push_back(p.wilted());
        if (p.age > max_age || p.wilted()) {
            clear_patches.push_back(PatchClearingInfo<uint32>(p.this_j, N, wilted.back(),
                                                              p.age));
        } else remaining += N;
    }

    do_clearing<uint32>(clear_patches, remaining, wilted, clear_surv, eng);

    return;
}


void AllPatches::clear_patches(const double& max_N,
                               const double& clear_surv,
                               pcg32& eng) {

    std::vector<PatchClearingInfo<double>> clear_patches;
    clear_patches.reserve(patches.size());

    double remaining = 0;

    std::vector<bool> wilted;
    wilted.reserve(patches.size());

    for (const OnePatch& p : patches) {
        double N = p.total_aphids();
        wilted.push_back(p.wilted());
        if (N > max_N || p.wilted()) {
            clear_patches.push_back(PatchClearingInfo<double>(p.this_j, N, wilted.back(),
                                                              N));
        } else remaining += N;
    }

    do_clearing<double>(clear_patches, remaining, wilted, clear_surv, eng);

    return;
}
