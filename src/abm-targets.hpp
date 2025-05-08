#ifndef __PSEUDOGAMEOFCLONES_ABM_TARGETS_H
#define __PSEUDOGAMEOFCLONES_ABM_TARGETS_H

#include <RcppArmadillo.h>
#include <vector>
#include <deque>
#include <string>
#include <limits>
#include <pcg/pcg_random.hpp>   // pcg prng


#include "pseudogameofclones_types.hpp"
#include "pcg.hpp"              // runif_ab fxn


using namespace Rcpp;





/*
 Class to convert back and forth between 2D to 1D indices.
 It assumes that x and y coordinates are sorted first by y coordinates, then
 by x coordinates.
 It also assumes 0-based indices.
 An example matrix of x and y coordinates might start like this:

 #>      x y
 #> [0,] 0 0
 #> [1,] 1 0
 #> [2,] 2 0
 #> [3,] 3 0
 #> [4,] 4 0
 #> [5,] 5 0

 */
class DimensionConverter {

    std::vector<uint32> neigh_x;
    std::vector<uint32> neigh_y;

    uint32 x_size;
    uint32 y_size;

public:

    DimensionConverter(const uint32& x_size_, const uint32& y_size_)
        : neigh_x(), neigh_y(), x_size(x_size_), y_size(y_size_) {
        neigh_x.reserve(3);
        neigh_y.reserve(3);
    }

    // Convert from 1D to 2D:
    void to_2d(uint32& x, uint32& y, const uint32& k) const {
        x = k - y_size * (k / y_size);
        y = k / y_size;
        return;
    }
    // Overloaded for signed ints (for use with Rcpp::IntegerVector)
    void to_2d(int& x, int& y, const uint32& k) const {
        x = k - y_size * (k / y_size);
        y = k / y_size;
        return;
    }
    // Convert from 2D to 1D:
    void to_1d(uint32& k, const uint32& x, const uint32& y) const {
        k = (y * x_size + x);
        return;
    }
    /*
     Return indices (in 1D) for all neighbors based on a 1D input coordinate.
     NOTE:
        - It clears `indices` before adding to it
        - It also returns an index for the focal point
     */
    void get_neighbors(std::vector<uint32>& indices,
                       const uint32& k) {
        uint32 x0 = k - y_size * (k / y_size);
        uint32 y0 = k / y_size;
        indices.clear();
        neigh_x.clear();
        neigh_y.clear();
        if (x0 > 0) neigh_x.push_back(x0-1);
        neigh_x.push_back(x0);
        if (x0 < x_size-1) neigh_x.push_back(x0+1);
        if (y0 > 0) neigh_y.push_back(y0-1);
        neigh_y.push_back(y0);
        if (y0 < x_size-1) neigh_y.push_back(y0+1);
        uint32 k_out;
        for (const uint32& x : neigh_x) {
            for (const uint32& y : neigh_y) {
                to_1d(k_out, x, y);
                indices.push_back(k_out);
            }
        }
        return;
    }


};




/*
 Sample locations from a vector of probabilities, where probabilities
 get updated as samples are chosen.
 */
class LocationSampler {

    // This is to avoid infinite weights.
    // It's still a massive number (1.797693e+302 on my machine):
    const double max_wt = std::numeric_limits<double>::max() / 1e6;

    arma::vec weights;
    arma::vec cs_probs;
    uint32 n;

    bool needs_recalc;

    // Make `cs_probs` into a vector that's the cumulative sum of
    // weights / sum(weights). The last value in `cs_probs` should always be 1.
    void calc_cumsum() {
        double p_sum = arma::accu(weights);
        if (p_sum <= 0) {
            weights.ones();
            p_sum = weights.n_elem;
        }
        cs_probs(0) = weights(0) / p_sum;
        for (uint32 i = 1; i < n; i++) {
            cs_probs(i) = cs_probs(i-1) + weights(i) / p_sum;
        }
        needs_recalc = false;
        return;
    }

public:

    LocationSampler(const uint32& n_)
        : weights(n_, arma::fill::ones),
          cs_probs(n_, arma::fill::none),
          n(n_),
          needs_recalc(true) {
        calc_cumsum();
    }

    /*
     Update weights by multiplying by a given value, but only if the weight
     wasn't already changed to zero (which happens if its location was
     already sampled).
     These functions do NOT re-calculate cumulative sums because
     that is done later inside `sample`.
     They both also change `needs_recalc` to true if they actually change
     one or more probabilities.
     */
    void update_weights(const std::vector<uint32>& indices,
                        const double& wt_val) {
        if (wt_val == 1) return;
        uint32 n_changed = 0;
        for (const uint32& k : indices) {
            if (weights(k) > 0 && weights(k) < max_wt) {
                weights(k) *= wt_val;
                if (weights(k) > max_wt) weights(k) = max_wt;
                n_changed++;
            }
        }
        if (n_changed > 0) needs_recalc = true;
        return;
    }
    void update_weights(const uint32& k,
                        const double& wt_val) {
        if (weights(k) > 0 && weights(k) < max_wt && wt_val != 1) {
            weights(k) *= wt_val;
            if (weights(k) > max_wt) weights(k) = max_wt;
            needs_recalc  = true;
        }
        return;
    }

    // Check to see if probabilities needs re-calculated, then
    // do weighted sampling for an index from 0 to (n-1):
    uint32 sample(pcg32& eng) {
        if (needs_recalc) calc_cumsum();
        double u = runif_01(eng);
        uint32 k = 0;
        while (k < n && cs_probs(k) < u) k++;
        return k;
    }



};


class LandSimmer {

    arma::mat wt_mat;
    arma::ivec n_samples;
    uint32 n_types;
    uint32 n_points;
    int x_size;
    int y_size;
    bool allow_overlap;

    // One location sampler for each type:
    std::vector<LocationSampler> samplers;

    // Convert coordinates between 1D and 2D:
    DimensionConverter dim_conv;

    // Object collecting samples:
    std::vector<std::vector<uint32>> samps;

    // Random number generator
    pcg32 eng;

    // Number of points assigned to 1 or more types (used to define output)
    uint32 n_used_pts;


public:

    // Note: not reserving storage for each item in `samps` upon initialization.
    // I do it below when an item in `samps` is assigned its first item.
    // Doing it this way saves a bit of memory.

    LandSimmer(const arma::mat& wt_mat_,
               const arma::ivec& n_samples_,
               int x_size_,
               int y_size_,
               const bool& allow_overlap_)
    : wt_mat(wt_mat_),
      n_samples(n_samples_),
      n_types(wt_mat_.n_rows),
      n_points(x_size_ * y_size_),
      x_size(x_size_),
      y_size(y_size_),
      allow_overlap(allow_overlap_),
      samplers(n_types, LocationSampler(n_points)),
      dim_conv(x_size, y_size),
      samps(n_points),
      eng(),
      n_used_pts(0U) {
        seed_pcg(eng);
      }


    void run(RcppThread::ProgressBar& prog_bar, const bool& show_progress) {

        int32_t total_samps = arma::accu(n_samples);
        uint32 n = 0;// for iterating

        // # sims done for each type (also doubles as indices for output):
        arma::ivec sims_done(n_types, arma::fill::zeros);
        uint32 x, y, k;
        std::vector<uint32> neighbors;
        neighbors.reserve(9); // highest number of neighbors possible

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
                n++;

                if (show_progress) prog_bar++;
                if (n % 10 == 0) RcppThread::checkUserInterrupt();
            }
        }

        return;
    }

    DataFrame create_output(const bool& fill_all) {

        if (fill_all) n_used_pts = n_points;

        // Create output dataframe:
        DataFrame out_df = DataFrame::create(
            _["x"] = IntegerVector(n_used_pts),
            _["y"] = IntegerVector(n_used_pts),
            _["type"] = CharacterVector(n_used_pts));
        // References to columns:
        IntegerVector out_x = out_df[0];
        IntegerVector out_y = out_df[1];
        CharacterVector out_type = out_df[2];

        if (fill_all) {
            for (uint32 k = 0; k < samps.size(); k++) {
                out_type(k) = "";
                if (!samps[k].empty()) {
                    std::sort(samps[k].begin(), samps[k].end());
                    for (uint32 s = 0; s < samps[k].size(); s++) {
                        out_type(k) += std::to_string(samps[k][s]);
                        if (s < samps[k].size() - 1) out_type(k) += "_";
                    }
                }
                dim_conv.to_2d(out_x(k), out_y(k), k);
                // Convert from 0- to 1-based indexing:
                out_x(k)++;
                out_y(k)++;
            }
        } else {
            uint32 i = 0;
            for (uint32 k = 0; k < samps.size(); k++) {
                if (!samps[k].empty()) {
                    std::sort(samps[k].begin(), samps[k].end());
                    out_type(i) = "";
                    for (uint32 s = 0; s < samps[k].size(); s++) {
                        out_type(k) += std::to_string(samps[k][s]);
                        if (s < samps[k].size() - 1) out_type(k) += "_";
                    }
                    dim_conv.to_2d(out_x(i), out_y(i), k);
                    out_x(i)++;
                    out_y(i)++;
                    i++;
                }
            }
        }

        out_df.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});

        return out_df;


    }



};



#endif
