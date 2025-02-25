
#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <pcg/pcg_random.hpp>   // pcg prng

#include "pseudogameofclones_types.hpp"  // integer types
#include "pcg.hpp"              // runif_01, seed_rng functions
#include "abm-targets.hpp"      // DimensionConverter and LocationSampler classes



using namespace Rcpp;



//' @export
//' @noRd
//[[Rcpp::export]]
DataFrame target_type_sims(const int& x_size,
                           const int& y_size,
                           const arma::mat& corr,
                           const arma::ivec& n_samples,
                           const bool& fill_all = false) {

    if (x_size <= 0) stop("x_size must be > 0");
    if (y_size <= 0) stop("y_size must be > 0");
    if (corr.n_elem == 0) stop("corr cannot be empty");
    if (!corr.is_symmetric()) stop("corr must be symmetric");
    if (arma::any(arma::vectorise(corr) < 0)) stop("corr cannot contain values < 0");
    if (n_samples.n_elem != corr.n_rows)
        stop("length(n_samples) == nrow(corr) must be true");
    if (arma::any(n_samples <= 0)) stop("n_samples cannot contain values <= 0");
    if (arma::any(n_samples > x_size * y_size))
        stop("n_samples cannot contain values > x_size * y_size");

    uint32 n_types = corr.n_rows;
    uint32 n_points = x_size * y_size;

    // One location sampler for each type:
    std::vector<LocationSampler> samplers(n_types, LocationSampler(n_points));

    // Convert coordinates between 1D and 2D:
    DimensionConverter dim_conv(x_size, y_size);

    // Object collecting samples:
    std::vector<std::vector<uint32>> samps(n_points);
    // Note: not reserving storage here. I do it below when an item in `samps`
    // is assigned its first item. Doing it this way saves a bit of memory.

    // Number of points assigned to 1 or more types (used to define output below)
    uint32 n_used_pts = 0;

    int32_t total_samps = arma::accu(n_samples);
    // # sims done for each type (also doubles as indices for output):
    arma::ivec sims_done(n_types, arma::fill::zeros);
    uint32 x, y, k;
    std::vector<uint32> neighbors;
    neighbors.reserve(9); // highest number of neighbors possible

    pcg32 eng;
    seed_pcg(eng);

    while (total_samps > 0) {
        for (uint32 i = 0; i < n_types; i++) {

            if (sims_done(i) >= n_samples(i)) continue;

            k = samplers[i].sample(eng);
            dim_conv.to_2d(x, y, k); // assign new x and y based on k

            // reserve memory so that it doesn't have to be moved after this:
            if (samps[k].empty()) {
                n_used_pts++;
                samps[k].reserve(n_types);
            }
            // add to output:
            samps[k].push_back(i+1); //+1 to convert to R's 1-based indexing

            // Adjust sampling probabilities:
            dim_conv.get_neighbors(neighbors, k); // fill neighbors vector
            for (uint32 j = 0; j < n_types; j++) {
                samplers[j].update_weights(neighbors, corr(i,j));
            }
            /*
             Note: You don't have to update `k`th prob to zero after the call
             to `update_weights` on `neighbors` even though the latter also
             updates `k` because `update_weights` never updates weights that
             are already set to zero.
             */
            samplers[i].update_weights(k, 0.0);

            // Iterate sample numbers:
            sims_done(i)++;
            total_samps--;
        }
    }


    if (fill_all) n_used_pts = n_points;

    // Create output dataframe:
    List out_type(n_used_pts); // this has to be created separately & added later
    DataFrame out = DataFrame::create(
        _["x"] = IntegerVector(n_used_pts),
        _["y"] = IntegerVector(n_used_pts));
    // References to columns:
    IntegerVector out_x = out[0];
    IntegerVector out_y = out[1];

    bool samps_empty;
    if (fill_all) {
        std::vector<uint32> filler(1, n_types+1);
        for (uint32 k = 0; k < samps.size(); k++) {
            const uint32& i(k); // for consistently when !fill_all below
            samps_empty = samps[k].empty();
            if (samps_empty) {
                out_type(i) = filler;
            } else {
                std::sort(samps[k].begin(), samps[k].end());
                out_type(i) = samps[k];
            }
            dim_conv.to_2d(out_x(i), out_y(i), k);
        }
    } else {
        uint32 i = 0;
        for (uint32 k = 0; k < samps.size(); k++) {
            samps_empty = samps[k].empty();
            if (!samps_empty) {
                std::sort(samps[k].begin(), samps[k].end());
                out_type(i) = samps[k];
                dim_conv.to_2d(out_x(i), out_y(i), k);
                i++;
            }
        }
    }
    out["type"] = out_type;

    out.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});
    // `row.names` is required for playing nice with list column!
    out.attr("row.names") = Rcpp::seq(1, n_used_pts);

    return out;

}
