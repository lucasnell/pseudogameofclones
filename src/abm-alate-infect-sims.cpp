#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <pcg/pcg_random.hpp>   // pcg prng

#include "pseudogameofclones_types.hpp"  // integer types
#include "pcg.hpp"                       // runif_01, seed_rng functions

using namespace Rcpp;


/*
 Convert from 2D (x and y) to 1D.
 It assumes that x and y coordinates are sorted first by y coordinates, then
 by x coordinates.
 (^ This doesn't really matter, as long as I'm consistent, which I am here
    since I'm only using this function.)
 It also assumes 0-based indices.
 */
uint32 to_1d(const int& x, const int& y, const int& x_size) {
    uint32 k = (y * x_size + x);
    return k;
}


void check_ais_args(uint32& n_lands,
                    uint32& n_alates,
                    const int& x_size,
                    const int& y_size,
                    const double& delta_a,
                    const double& delta_p,
                    const std::vector<int>& land,
                    const std::vector<int>& alate,
                    const std::vector<bool>& infected,
                    const std::vector<int>& x,
                    const std::vector<int>& y) {

    if (x_size < 2 || x_size > 1e9) stop("x_size < 2 || x_size > 1e9");
    if (y_size < 2 || y_size > 1e9) stop("y_size < 2 || y_size > 1e9");
    if (delta_a < 0 || delta_a > 1) stop("delta_a < 0 || delta_a > 1");
    if (delta_p < 0 || delta_p > 1) stop("delta_p < 0 || delta_p > 1");

    if (land.size() == 0) stop("land.size() == 0");
    if (land.size() != alate.size()) stop("land.size() != alate.size()");
    if (land.size() != infected.size()) stop("land.size() != infected.size()");
    if (land.size() != x.size()) stop("land.size() != x.size()");
    if (land.size() != y.size()) stop("land.size() != y.size()");

    uint32 n = land.size();

    // Focal landscape and alate:
    int this_land = land[0];
    int this_alate = alate[0];
    if (this_land != 1) stop("land should start with 1");
    if (this_alate != 1) stop("alate should start with 1");
    // Go through once to get size of output and test that land and alate
    // range from 1 to their max and all integers in between.
    n_lands = 1;
    n_alates = 1;
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
        // Also checking `x` and `y` for valid values:
        if (x[i] < 1 || x[i] > x_size) stop("x contains item(s) outside 1:x_size");
        if (y[i] < 1 || y[i] > y_size) stop("y contains item(s) outside 1:y_size");
    }


    return;


}


//' Proportion of alates and plants that are inoculated, for a given landscape.
//'
//' ASSUMPTIONS OF THIS FUNCTION:
//'
//' 1) Vectors `land`, `alate`, `infected`, `x`, and `y` are all sorted first by
//'    `land`, then by `alate`, then by time (not an argument).
//'
//' 2) On the timeline of these simulations, even if an alate inoculates
//'    a plant with a virus, the virus doesn't have enough time to replicate
//'    inside the plant to cause a subsequent alate to get infected.
//'
//' 3) All landscapes have the same number of alates and landscape dimensions
//'    (`x_size` and `y_size`).
//'
//'
//' @param x_size Single integer indicating x dimension of landscape.
//'     Locations must be between `1` and `x_size`.
//' @param y_size Single integer indicating y dimension of landscape.
//'     Locations must be between `1` and `y_size`.
//' @param delta_a Single numeric indicating the probability that an
//'     uninoculated alate is loaded with a virus if it interacts with an
//'     inoculated plant.
//' @param delta_p Single numeric indicating the probability that an
//'     uninoculated plant is loaded with a virus if it interacts with an
//'     inoculated alate.
//' @param land Integer vector indicating the focal landscape.
//' @param alate Integer vector indicating the focal alate.
//' @param infected Logical vector indicating whether focal plant is
//'     infected with virus.
//' @param x Integer vector indicating location of the focal plant.
//'     All values must be between `1` and `x_size`.
//' @param y Integer vector indicating location of the focal plant.
//'     All values must be between `1` and `y_size`.
//'
//' @export
//'
//[[Rcpp::export]]
arma::mat alate_infect_sims(const int& x_size,
                            const int& y_size,
                            const double& delta_a,
                            const double& delta_p,
                            const std::vector<int>& land,
                            const std::vector<int>& alate,
                            const std::vector<bool>& infected,
                            const std::vector<int>& x,
                            const std::vector<int>& y) {

    uint32 n_lands, n_alates;

    // Check args and set `n_lands` and `n_alates`:
    check_ais_args(n_lands, n_alates, x_size, y_size, delta_a, delta_p,
                   land, alate, infected, x, y);

    uint32 n = land.size();

    pcg32 rng;
    seed_pcg(rng);

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
    // Focal landscape and alate:
    int this_land = land[0];
    int this_alate = alate[0];
    // ~ U(0,1) for inoculation sampling:
    double u;

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
        plant_inf = infected[i];

        if (!alate_inf && plant_inf) {
            // Uninoculated alate gets inoculated:
            u = runif_01(rng);
            if (u < delta_a) {
                alate_inf = true;
                inocs(this_land-1, 0) += 1;
            }
        } else if (alate_inf && !plant_inf) {
            k = to_1d(x[i]-1, y[i]-1, x_size);
            if (!plant_inocs[k]) {
                // Uninoculated plant gets inoculated:
                u = runif_01(rng);
                if (u < delta_p) {
                    plant_inf = true;
                    plant_inocs[k] = true;
                    inocs(this_land-1, 1) += 1;
                }
            }
        }
    }

    inocs.col(0) /= static_cast<double>(n_alates);
    inocs.col(1) /= static_cast<double>(n_plants);

    return inocs;

}
