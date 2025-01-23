#ifndef __PSEUDOGAMEOFCLONES_ABM_H
#define __PSEUDOGAMEOFCLONES_ABM_H

#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <memory>
#include <pcg/pcg_random.hpp>   // pcg prng

#ifndef RCPPTHREAD_OVERRIDE_COUT
#define RCPPTHREAD_OVERRIDE_COUT 1    // std::cout override
#endif
#ifndef RCPPTHREAD_OVERRIDE_CERR
#define RCPPTHREAD_OVERRIDE_CERR 1    // std::cerr override
#endif
// #ifndef RCPPTHREAD_OVERRIDE_THREAD
// #define RCPPTHREAD_OVERRIDE_THREAD 1  // std::thread override
// #endif
#include <RcppThread.h>         // multithreading

#include "pseudogameofclones_types.hpp"
#include "pcg.hpp"              // runif_ab fxn


using namespace Rcpp;



//
// New dx or dy that includes effects of reflecting off bound(s).
//
inline double reflect(const double& start,
                      double d,
                      const arma::vec& bounds) {
    if (bounds.n_elem != 2) stop("bounds.n_elem != 2");
    if (bounds(1) <= bounds(0)) stop("bounds(1) <= bounds(0)");
    arma::uvec past_bounds = {(start + d) < bounds(0), (start + d) > bounds(1)};
    double bound;
    while (arma::any(past_bounds)) {
        bound = past_bounds(0) ? bounds(0) : bounds(1);
        d = 2 * (bound - start) - d;
        past_bounds(0) = (start + d) < bounds(0);
        past_bounds(1) = (start + d) > bounds(1);
    }
    return d;
}

//
// Calculate Euclidean distance between two points
//
inline double distance(const double& x1, const double& y1,
                       const double& x2, const double& y2) {
    return std::sqrt(std::pow((x1 - x2), 2U) + std::pow((y1 - y2), 2U));
}


class TargetInfo {

    arma::mat target_xy;
    std::vector<uint32> type_map;
    std::vector<double> l_star_;
    std::vector<double> l_i_;
    std::vector<double> bias_;
    std::vector<uint32> n_stay_;
    std::vector<uint32> n_ignore_;

public:

    TargetInfo(const arma::mat& target_xy_,
                    const std::vector<uint32>& type_map_,
                    const std::vector<double>& l_star__,
                    const std::vector<double>& l_i__,
                    const std::vector<double>& bias__,
                    const std::vector<uint32>& n_stay__,
                    const std::vector<uint32>& n_ignore__)
        : target_xy(target_xy_), type_map(type_map_),
          l_star_(l_star__), l_i_(l_i__),
          bias_(bias__), n_stay_(n_stay__), n_ignore_(n_ignore__) {

        if (target_xy.n_cols != 2) stop("target_xy.n_cols != 2");
        if (target_xy.n_rows != type_map.size())
            stop("target_xy.n_rows != type_map.size()");

        uint32 max_idx = *std::max_element(type_map.begin(), type_map.end());
        if (l_star_.size() != max_idx+1U) stop("wrong length for l_star");
        if (l_i_.size() != max_idx+1U) stop("wrong length for l_i");
        if (bias_.size() != max_idx+1U) stop("wrong length for bias");
        if (n_stay_.size() != max_idx+1U) stop("wrong length for n_stay");
        if (n_ignore_.size() != max_idx+1U) stop("wrong length for n_ignore");

    };

    uint32 n_targets() const {
        return type_map.size();
    }

    /*
     For all below, `i` should be the index for the target_xy matrix row.
     This class maps that index to the index for the target type via
     the `type_map` vector.
     */

    double x(const uint32& i) const {
        return target_xy(i, 0);
    }
    double y(const uint32& i) const {
        return target_xy(i, 1);
    }
    double l_star(const uint32& i) const {
        return l_star_[type_map[i]];
    }
    double l_i(const uint32& i) const {
        return l_i_[type_map[i]];
    }
    double bias(const uint32& i) const {
        return bias_[type_map[i]];
    }
    uint32 n_stay(const uint32& i) const {
        return n_stay_[type_map[i]];
    }
    uint32 n_ignore(const uint32& i) const {
        return n_ignore_[type_map[i]];
    }
    // Overloaded to make arma vectors:
    arma::vec l_star(const arma::uvec& i) const {
        arma::vec out(i.n_elem);
        for (uint32 j = 0; j < i.n_elem; j++) {
            out(j) = l_star_[type_map[i(j)]];
        }
        return out;
    }
    arma::vec l_i(const arma::uvec& i) const {
        arma::vec out(i.n_elem);
        for (uint32 j = 0; j < i.n_elem; j++) {
            out(j) = l_i_[type_map[i(j)]];
        }
        return out;
    }
    arma::vec bias(const arma::uvec& i) const {
        arma::vec out(i.n_elem);
        for (uint32 j = 0; j < i.n_elem; j++) {
            out(j) = bias_[type_map[i(j)]];
        }
        return out;
    }
    arma::uvec n_stay(const arma::uvec& i) const {
        arma::uvec out(i.n_elem);
        for (uint32 j = 0; j < i.n_elem; j++) {
            out(j) = n_stay_[type_map[i(j)]];
        }
        return out;
    }
    arma::uvec n_ignore(const arma::uvec& i) const {
        arma::uvec out(i.n_elem);
        for (uint32 j = 0; j < i.n_elem; j++) {
            out(j) = n_ignore_[type_map[i(j)]];
        }
        return out;
    }


};







class ABMsimulator {

    pcg32 eng;
    const TargetInfo* target_info;
    double delta;
    arma::vec x_bounds;
    arma::vec y_bounds;
    double max_size;
    uint32 n_targets;
    // Vector of distances from searcher to target(s)
    arma::vec l_vec;
    // which target(s) are within l_star (and thus influence trajectory):
    arma::uvec wi_lstar;
    // which target (if any) is within l_i; if none, this is set to 1+n_targets:
    uint32 wi_li;
    // distance to nearest target:
    double l_min;
    // number of time points within l_i of most recent target:
    uint32 stayed;
    // max bias and n_ignore for any type of target:
    double max_bias;
    uint32 max_n_ignore;
    // max time steps to simulate:
    uint32 max_t;
    // `visited` below indicates number of time points since the most
    // recent visit for each target (only starts after leaving within
    // l_i of the target):
    arma::uvec visited;
    // These are objects used temporarily in the simulations themselves:
    double rw_theta;
    double rw_dx;
    double rw_dy;
    double dr_theta;
    double dr_dx;
    double dr_dy;
    bool within_li;
    bool within_lstar;
    bool keep_staying;
    arma::vec wts;

    /*

     Update `wi_lstar` (index / indices for targets within l_star),
     `wi_li` (index for target within l_i),
     `l_min` (distance to nearest target),
     and optionally `l_vec` (distances to all targets),
     including ignoring recently visited target(s).

     `t` is the index for time, indicating which values from
     `x` and `y` to use.

    */
    void update_wi_lvec(const uint32& t,
                        const bool& update_l_vec);

    void stay_at_target(const uint32& t);

    void biased_movement(const uint32& t);

    void iterate(const uint32& t);


public:

    // *==================== outputs ====================*
    arma::vec x;
    arma::vec y;
    // whether on target at given time:
    std::vector<bool> on_target;
    // whether hit a *new* target at given time:
    std::vector<bool> new_target;


    ABMsimulator(const TargetInfo& target_info_,
                 const double& delta_,
                 const double& x_size,
                 const double& y_size,
                 const std::vector<double>& bias,
                 const std::vector<uint32>& n_ignore,
                 const uint32& max_t_,
                 const double& x0,
                 const double& y0,
                 const bool& randomize_xy0)
        : eng(),
          target_info(&target_info_),
          delta(delta_),
          x_bounds({-x_size / 2, x_size / 2}),
          y_bounds({-y_size / 2, y_size / 2}),
          max_size(std::sqrt(std::pow(x_size, 2U) + std::pow(y_size, 2U))),
          n_targets(target_info_.n_targets()),
          l_vec(target_info_.n_targets(), arma::fill::none),
          wi_lstar(),
          wi_li(),
          l_min(),
          stayed(0),
          max_bias(*std::max_element(bias.begin(), bias.end())),
          max_n_ignore(*std::max_element(n_ignore.begin(), n_ignore.end())),
          max_t(max_t_),
          visited(target_info_.n_targets(),
                  arma::fill::value(max_n_ignore * 2U)),
          rw_theta(),
          rw_dx(),
          rw_dy(),
          dr_theta(),
          dr_dx(),
          dr_dy(),
          within_li(),
          within_lstar(),
          keep_staying(),
          wts(),
          x(max_t_ + 1U, arma::fill::none),
          y(max_t_ + 1U, arma::fill::none),
          on_target(max_t_ + 1U, false),
          new_target(max_t_ + 1U, false) {

        std::vector<uint64> seeds = mt_seeds();

        seed_pcg(eng, seeds);
        if (randomize_xy0) {
            x(0) = runif_ab(eng, x_bounds(0), x_bounds(1));
            y(0) = runif_ab(eng, y_bounds(0), y_bounds(1));
        } else {
            x(0) = x0;
            y(0) = y0;
        }

        // Update l_vec, wi_lstar, and l_min:
        update_wi_lvec(0, true);

        // Update on_target and new_target if searcher starts near a target:
        if (wi_li < n_targets && l_min <= target_info->l_i(wi_li)) {
            on_target[0] = true;
            new_target[0] = true;
            visited(wi_li) = 0;
        }

    };


    void run(RcppThread::ProgressBar& prog_bar, const bool& show_progress) {
        for (uint32 t = 0; t < max_t; t++) {
            this->iterate(t);
            if (show_progress) prog_bar++;
            if (t % 10 == 0) RcppThread::checkUserInterrupt();
        }
        return;
    }


};




#endif
