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


// class AliasSampler {
// public:
//     AliasSampler() : Prob(), Alias(), n(0) {};
//     AliasSampler(const std::vector<double>& probs)
//         : Prob(probs.size()), Alias(probs.size()), n(probs.size()) {
//         arma::vec p(probs);
//         construct(p);
//     }
//     AliasSampler(arma::vec probs)
//         : Prob(probs.n_elem), Alias(probs.n_elem), n(probs.n_elem) {
//         construct(probs);
//     }
//     // Copy constructor
//     AliasSampler(const AliasSampler& other)
//         : Prob(other.Prob), Alias(other.Alias), n(other.n) {}
//
//     // Actual alias sampling
//     inline uint32 sample(pcg64& eng) const {
//         // Fair dice roll from n-sided die
//         uint32 i = runif_01(eng) * n;
//         // uniform in range (0,1)
//         double u = runif_01(eng);
//         if (u < Prob[i]) return(i);
//         return Alias[i];
//     };
//
// private:
//     std::vector<double> Prob;
//     std::vector<uint32> Alias;
//     uint32 n;
//
//
//     void construct(arma::vec& p) {
//
//         p /= arma::accu(p);  // make sure they sum to 1
//         p *= n;
//
//         std::deque<uint32> Small;
//         std::deque<uint32> Large;
//         for (uint32 i = 0; i < n; i++) {
//             if (p(i) < 1) {
//                 Small.push_back(i);
//             } else Large.push_back(i);
//         }
//
//         uint32 l, g;
//         while (!Small.empty() && !Large.empty()) {
//             l = Small.front();
//             Small.pop_front();
//             g = Large.front();
//             Large.pop_front();
//             Prob[l] = p(l);
//             Alias[l] = g;
//             p(g) = (p(g) + p(l)) - 1;
//             if (p(g) < 1) {
//                 Small.push_back(g);
//             } else Large.push_back(g);
//         }
//         while (!Large.empty()) {
//             g = Large.front();
//             Large.pop_front();
//             Prob[g] = 1;
//         }
//         while (!Small.empty()) {
//             l = Small.front();
//             Small.pop_front();
//             Prob[l] = 1;
//         }
//
//         return;
//     }
// };



#endif
