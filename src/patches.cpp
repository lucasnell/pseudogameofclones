
#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "clonewars_types.hpp"  // integer types
#include "patches.hpp"          // patch classes
#include "pcg.hpp"              // runif_ fxns



using namespace Rcpp;



// Do the actual clearing of patches while avoiding extinction
template <typename T>
inline void AllPatches::do_clearing(std::vector<PatchClearingInfo<T>>& clear_patches,
                                    double& remaining,
                                    pcg32& eng) {

    if (clear_patches.size() == 0) return;

    /*
     If removing all these patches cause total extinction, we need to remove
     some indices from `inds`.
     We'll sort (high to low) the vector by the threshold variable (age or # aphids),
     then remove the last item in the vector until we have aphids again.
     */
    if (remaining == 0) {

        std::sort(clear_patches.begin(), clear_patches.end());

        while (remaining == 0) {
            remaining += clear_patches.back().N;
            clear_patches.pop_back();
        }

    }

    if (K_error) {
        for (uint32 i = 0; i < clear_patches.size(); i++) {
            double K_ = mu_K_ + sigma_K_ * norm_distr(eng);
            patches[clear_patches[i].ind].clear(K_);
        }
    } else {
        for (uint32 i = 0; i < clear_patches.size(); i++) {
            patches[clear_patches[i].ind].clear();
        }
    }

    return;
}




// Clear patches by either a maximum age or total abundance
void AllPatches::clear_patches(const uint32& max_age,
                               pcg32& eng) {

    std::vector<PatchClearingInfo<uint32>> clear_patches;
    clear_patches.reserve(patches.size());

    // # remaining aphids for non-cleared patches
    double remaining = 0;

    for (const OnePatch& p : patches) {
        double N = p.total_aphids();
        if (p.age > max_age) {
            clear_patches.push_back(PatchClearingInfo<uint32>(p.this_j, N, p.age));
        } else remaining += N;
    }

    do_clearing<uint32>(clear_patches, remaining, eng);

    return;
}


void AllPatches::clear_patches(const double& max_N,
                               pcg32& eng) {

    std::vector<PatchClearingInfo<double>> clear_patches;
    clear_patches.reserve(patches.size());

    double remaining = 0;

    for (const OnePatch& p : patches) {
        double N = p.total_aphids();
        if (N > max_N) {
            clear_patches.push_back(PatchClearingInfo<double>(p.this_j, N, N));
        } else remaining += N;
    }

    do_clearing<double>(clear_patches, remaining, eng);

    return;
}
