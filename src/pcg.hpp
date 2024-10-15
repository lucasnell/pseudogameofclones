#ifndef __GAMEOFCLONES_PCG_H
#define __GAMEOFCLONES_PCG_H

#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <pcg/pcg_extras.hpp>  // pcg 128-bit integer type
#include <pcg/pcg_random.hpp> // pcg prng

#include "pseudogameofclones_types.hpp"

using namespace Rcpp;


typedef pcg_extras::pcg128_t uint128;


namespace pcg {
    const double max = static_cast<double>(pcg32::max());
    const long double max64 = static_cast<long double>(pcg64::max());
}



/*
 ========================

 Seeding

 ========================
 */

/*
 For consistent results from multi-thread operations, you should call `mt_seeds`
 when outside multi-thread mode, then create a PRNG (`pcg32` object) once
 inside multi-thread mode, one PRNG per thread.
 For each simulation rep, use a `std::vector<uint64>` vector (from inside
 the object output from `mt_seeds`) to seed the PRNG using `seed_pcg`.
 Using seeds for each rep (not thread) means that if you run `set.seed(999)`
 in R before running a process that can be run across multiple threads,
 you'll get the same output no matter how many threads you use.
 */

// To sample for seeds before multi-thread operations
inline std::vector<std::vector<uint64>> mt_seeds(const uint32& n_reps) {

    std::vector<std::vector<uint64>> sub_seeds(n_reps, std::vector<uint64>(4));

    for (uint32 i = 0; i < n_reps; i++) {
        // These are 32-bit integers cast as 64-bit for downstream compatibility
        sub_seeds[i] = as<std::vector<uint64>>(Rcpp::runif(4,0,4294967296));
    }

    return sub_seeds;
}


// Change seed of existing rng from four 32-bit seeds (casted to 64-bit)
inline void seed_pcg(pcg32& eng, const std::vector<uint64>& sub_seeds) {

    uint64 seed1;
    uint64 seed2;
    // Converting to two 64-bit seeds for pcg32
    seed1 = (sub_seeds[0]<<32) + sub_seeds[1];
    seed2 = (sub_seeds[2]<<32) + sub_seeds[3];

    eng.seed(seed1, seed2);

    return;
}


/*
 ========================
 Number generation
 ========================
 */

// uniform in range (0,1)
inline double runif_01(pcg32& eng) {
    return (static_cast<double>(eng()) + 1) / (pcg::max + 2);
}

// uniform in range (a,b)
inline double runif_ab(pcg32& eng, const double& a, const double& b) {
    return a + ((static_cast<double>(eng()) + 1) / (pcg::max + 2)) * (b - a);
}


#endif
