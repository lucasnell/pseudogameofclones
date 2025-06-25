
// [[Rcpp::depends(dqrng, BH, RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <pcg_random.hpp>
#include <vector>
#include <algorithm>

using namespace Rcpp;

namespace pcg {
    const double max = static_cast<double>(pcg32::max());
    const long double max64 = static_cast<long double>(pcg64::max());
}

// uniform in range (0,1)
inline double runif_01(pcg32& eng) {
    return (static_cast<double>(eng()) + 1) / (pcg::max + 2);
}

typedef uint_fast32_t uint32;

/*
 Convert from 2D (x and y) to 1D.
 It assumes that x and y coordinates are sorted first by y coordinates, then
 by x coordinates.
 It also assumes 0-based indices.
 */
uint32 to_1d(const int& x, const int& y, const int& x_size) {
    uint32 k = (y * x_size + x);
    return k;
}


//'
//' Proportion of alates and plants that are inoculated, for a given landscape.
//'
//' ASSUMPTIONS OF THIS FUNCTION:
//'
//' 1) Vectors `land`, `alate`, `to`, `x`, and `y` are all sorted first by
//'    `land`, then by `alate`, then by time (not an argument).
//'
//' 2) For `to`, 1 = just virus, 2 = just Pseudomonas, 3 = none, and 4 = both.
//'    Thus, 1 and 4 are virus-infected, and 2 and 3 are not.
//'
//' 3) On the timeline of these simulations, even if an alate inoculates
//'    a plant with a virus, the virus doesn't have enough time to replicate
//'    inside the plant to cause a subsequent alate to get infected.
//'
//' 4) All landscapes have the same number of alates and landscape dimensions
//'    (`x_size` and `y_size`).
//'
//'
//[[Rcpp::export]]
arma::mat inoc_sims(const int& x_size,
                    const int& y_size,
                    const double& delta_a,
                    const double& delta_p,
                    const std::vector<int>& land,
                    const std::vector<int>& alate,
                    const std::vector<int>& to,
                    const std::vector<int>& x,
                    const std::vector<int>& y,
                    const std::vector<double>& rnds) {
    //                           const uint32& seed) {

    if (land.size() == 0) stop("land.size() == 0");
    if (land.size() != alate.size()) stop("land.size() != alate.size()");
    if (land.size() != to.size()) stop("land.size() != to.size()");
    if (land.size() != x.size()) stop("land.size() != x.size()");
    if (land.size() != y.size()) stop("land.size() != y.size()");

    uint32 n = land.size();

    std::vector<double>::const_iterator u = rnds.begin();
    // pcg32 rng(seed);

    // Focal landscape and alate:
    int this_land = land[0];
    int this_alate = alate[0];
    if (this_land != 1) stop("land should be sorted from 1 to max(land)");
    if (this_alate != 1) stop("alate should be sorted from 1 to max(alates)");
    // Go through once to get size of output and test that land and alate
    // range from 1 to their max and all integers in between.
    uint32 n_lands = 1;
    uint32 n_alates = 1;
    bool new_land;
    for (uint32 i = 0; i < n; i++) {
        new_land = land[i] != this_land;
        if (new_land) {
            if ((land[i] - this_land) != 1)
                stop("land should be sorted from 1 to max(land) and include all values from 1 to max(land)");
            if (this_alate != static_cast<int>(n_alates))
                stop("all landscapes should  have the same number of alates");
            if (alate[i] == this_alate && n_alates > 1)
                stop("alate should switch when land does if max(alates) > 1");
            this_land = land[i];
            n_lands++;
        }
        if (alate[i] != this_alate) {
            if ((alate[i] - this_alate) != 1 && !new_land)
                stop("alate should be sorted from 1 to max(alates) and include all values from 1 to max(alates)");
            this_alate = alate[i];
            if (land[i] == 1) {
                n_alates++;
            } else if (this_alate > n_alates)
                stop("alate should be sorted from 1 to n_alates");
        }
        // Also checking `to`, `x`, and `y` for valid values:
        if (to[i] < 1 || to[i] > 4) stop("to contains item(s) outside 1:4");
        if (x[i] < 1 || x[i] > x_size) stop("x contains item(s) outside 1:x_size");
        if (y[i] < 1 || y[i] > y_size) stop("y contains item(s) outside 1:y_size");
    }


    uint32 n_plants = x_size * y_size;

    /*
     Vector of which plants have been inoculated. It's important to keep
     track of individual plants to avoid counting the same plants twice if
     they've been visited by more than one inoculated alate.
     */
    std::vector<bool> plant_inocs(n_plants, false);

    // Total alates and plants inoculated for each landscape:
    arma::mat inocs(n_lands, 2U, arma::fill::zeros);

    // Whether current alate is infected:
    bool alate_inf = false;
    // Whether current plant is infected:
    bool plant_inf = false;
    // index for which plant the x and y refers to:
    uint32 k;
    // Reset focal landscape and alate:
    this_land = land[0];
    this_alate = alate[0];



    for (uint32 i = 0; i < n; i++) {

        if (land[i] != this_land) {
            this_land = land[i];
            std::fill(plant_inocs.begin(), plant_inocs.end(), 0);
        }

        if (alate[i] != this_alate) {
            alate_inf = false;
            this_alate = alate[i];
        }

        /*
         For `plant_inf` on the next line, I'm not considering `plant_inocs[k]`
         because that would cause alates within one time step to cause
         subsequent alates to get infected.
         We assume that the virus needs more time than this to replicate
         inside the plant to infect alates.
         */
        plant_inf = (to[i] == 1 || to[i] == 4);

        if (!alate_inf && plant_inf) {
            // Uninoculated alate gets inoculated:
            if (*u < delta_a) {
                alate_inf = true;
                inocs(this_land-1, 0) += 1;
            }
            u++;
            if (u == rnds.end()) stop("u == rnds.end()");
        } else if (alate_inf && !plant_inf) {
            k = to_1d(x[i]-1, y[i]-1, x_size);
            if (!plant_inocs[k]) {
                // Uninoculated plant gets inoculated:
                if (*u < delta_p) {
                    plant_inf = true;
                    plant_inocs[k] = true;
                    inocs(this_land-1, 1) += 1;
                }
                u++;
                if (u == rnds.end()) stop("u == rnds.end()");
            }
        }
    }

    inocs.col(0) /= static_cast<double>(n_alates);
    inocs.col(1) /= static_cast<double>(n_plants);

    return inocs;

}
