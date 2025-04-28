
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



//' Simulate targets of different types
//'
//' Simulate locations of targets of different types along an evenly spaced
//' grid of integers, where the placement of one target can affect subsequent
//' placement of other targets.
//' Target locations are drawn from all combinations of `1` to
//' `x_size-1` and `1` to `y_size-1`.
//' I'm subtracting 1 because the landscape in `searcher_sims` is defined by
//' the bounds `c(0, x_size)` and `c(0, y_size)`, so this subtraction
//' keeps a gap of 1 between the most extreme locations and the bounds
//' of the landscape.
//' This prevents targets from being located on the bounds.
//'
//'
//' @param x_size Single integer indicating x dimension of search area.
//'     Locations will be drawn from `1` to `x_size-1`.
//'     See description above for why.
//' @param y_size Single integer indicating y dimension of search area.
//'     Locations will be drawn from `1` to `y_size-1`.
//'     See description above for why.
//' @param wt_mat Square numeric matrix indicating how sample weighting on
//'     neighboring locations is affected by a target of each type being
//'     placed in a particular spot.
//'     Item `wt_mat[i,j]` indicates the effect of target type `i` on
//'     subsequent samplings of target type `j`.
//'     Values above 1 cause neighboring locations to be more likely to be
//'     sampled later, while values below 1 cause them to be less likely
//'     sampled.
//'     For locations that have been adjusted using `wt_mat` multiple times,
//'     the weights are multiplied by each other
//'     (e.g., `w *= wt_mat[1,3]` at time `t`, then
//'     `w *= wt_mat[1,2]` at time `t+1`).
//'     All weights start with values of `1`.
//'     The matrix should have the same number of rows and columns as the
//'     number of target types.
//' @param n_samples Integer vector indicating the number of samples to
//'     produce per target type. The length should equal the number of
//'     target types, and all values should be `>=1`.
//' @param allow_overlap Single logical indicating whether to allow
//'     targets to be more than one type.
//'     An example of a target being multiple types would be a plant that
//'     hosts an epiphytic bacteria and is also infected with a virus.
//'     Defaults to `TRUE`.
//' @param fill_all Single logical indicating whether to have a target at all
//'     points in the landscape. If `TRUE`, then locations where no targets
//'     were simulated will be assigned to a target type `length(n_samples)+1`.
//'     Defaults to `TRUE`.
//'
//'
//' @return A [`tibble`][tibble::tbl_df] with columns `x`, `y`, and `type`.
//'
//'
//' @export
//'
//[[Rcpp::export]]
DataFrame target_type_sims(int x_size,
                           int y_size,
                           const arma::mat& wt_mat,
                           const arma::ivec& n_samples,
                           const bool& allow_overlap = true,
                           const bool& fill_all = true) {

    if (x_size < 2) stop("x_size must be >= 2");
    if (y_size < 2) stop("y_size must be >= 2");
    if (wt_mat.n_elem == 0) stop("wt_mat cannot be empty");
    if (!wt_mat.is_square()) stop("wt_mat must be square");
    if (arma::any(arma::vectorise(wt_mat) < 0)) stop("wt_mat cannot contain values < 0");
    if (n_samples.n_elem != wt_mat.n_rows)
        stop("length(n_samples) == nrow(wt_mat) must be true");
    if (arma::any(n_samples <= 0)) stop("n_samples cannot contain values <= 0");

    // I'm reducing these by 1 bc I want to sample from 1 to floor(x_size)
    // and floor(y_size). See description in docs above for why.
    x_size--;
    y_size--;

    if (arma::any(n_samples > x_size * y_size))
        stop("n_samples cannot contain values > (x_size-1) * (y_size-1)");

    uint32 n_types = wt_mat.n_rows;
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
    if (!allow_overlap && total_samps > n_points) {
        stop("sum(n_samples) cannot be > (x_size-1) * (y_size-1) when allow_overlap = FALSE");
    }

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
                samplers[j].update_weights(neighbors, wt_mat(i,j));
                if (!allow_overlap) samplers[j].update_weights(k, 0.0);
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
            // Convert from 0- to 1-based indexing:
            out_x(i)++;
            out_y(i)++;
        }
    } else {
        uint32 i = 0;
        for (uint32 k = 0; k < samps.size(); k++) {
            samps_empty = samps[k].empty();
            if (!samps_empty) {
                std::sort(samps[k].begin(), samps[k].end());
                out_type(i) = samps[k];
                dim_conv.to_2d(out_x(i), out_y(i), k);
                out_x(i)++;
                out_y(i)++;
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
