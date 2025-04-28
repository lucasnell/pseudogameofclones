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


// Simple sign function:
template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//
// Iterate x or y and include effects of reflecting off bound(s).
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
    return start + d;
}

//
// Calculate Euclidean distance between two points
//
inline double distance(const double& x1, const double& y1,
                       const double& x2, const double& y2) {
    return std::sqrt(std::pow((x1 - x2), 2U) + std::pow((y1 - y2), 2U));
}




class OneTargetType {

    std::vector<double> l_star_;
    std::vector<double> bias_;
    std::vector<double> abs_bias_;
    double l_int_;
    uint32 n_stay_;
    uint32 n_ignore_;

public:

    OneTargetType(const std::vector<double>& l_star__,
                  const std::vector<double>& bias__,
                  const double& l_int__,
                  const uint32& n_stay__,
                  const uint32& n_ignore__)
        : l_star_(),
          bias_(),
          abs_bias_(),
          l_int_(l_int__),
          n_stay_(n_stay__),
          n_ignore_(n_ignore__) {

        uint32 n = l_star__.size();
        if (bias__.size() != n) stop("bias and l_star sizes don't match");

        l_star_.reserve(n);
        bias_.reserve(n);
        abs_bias_.reserve(n);

        // Create vector of indices that sorts by increasing `l_star`:
        std::vector<uint32> indices(n);
        std::iota(indices.begin(), indices.end(), 0U);
        std::sort(indices.begin(), indices.end(),
                  [&](uint32 A, uint32 B) -> bool {
                      return l_star__[A] < l_star__[B];
                  });
        // Now add `l_star`, `bias`, and `abs_bias` based on this sorting:
        for (const uint32& i : indices) {
            l_star_.push_back(l_star__[i]);
            bias_.push_back(bias__[i]);
            abs_bias_.push_back(std::abs(bias__[i]));
        }

        // Check that increasing `l_star` always coincides with decreasing
        // `abs_bias` (so `abs_bias` must be sorted in decreasing order):
        for (uint32 i = 1; i < n; i++) {
            if (abs_bias_[i-1] < abs_bias_[i]) {
                stop("higher `l_star` must always coincide with lower `abs(bias)`");
            }
        }


    };

    double l_star() const {
        return l_star_.back(); // bc they're sorted, this is the greatest value
    }
    double bias(const double& l_i) const {
        uint32 i = 0;
        while (l_star_[i] < l_i) i++;
        if (i >= bias_.size()) return -100;
        return bias_[i];
    }
    double abs_bias(const double& l_i) const {
        uint32 i = 0;
        while (l_star_[i] < l_i) i++;
        if (i >= abs_bias_.size()) return -100;
        return abs_bias_[i];
    }
    double l_int() const {
        return l_int_;
    }
    uint32 n_stay() const {
        return n_stay_;
    }
    uint32 n_ignore() const {
        return n_ignore_;
    }

};




class TargetInfo {

    arma::mat target_xy;
    std::vector<uint32> type_map_;
    std::vector<OneTargetType> types_;

public:

    TargetInfo(const arma::mat& target_xy_,
                    const std::vector<uint32>& type_map__,
                    const std::vector<std::vector<double>>& l_star_,
                    const std::vector<std::vector<double>>& bias_,
                    const std::vector<double>& l_int_,
                    const std::vector<uint32>& n_stay_,
                    const std::vector<uint32>& n_ignore_)
        : target_xy(target_xy_), type_map_(type_map__),
          types_() {

        if (target_xy.n_cols != 2) stop("target_xy.n_cols != 2");
        if (target_xy.n_rows != type_map_.size())
            stop("target_xy.n_rows != type_map_.size()");

        /*
         The next two lines are done bc this allows me to have target types
         that aren't present in the landscape.
         This is useful for simulations where you remove or replace
         target type(s).
         */
        uint32 n_types = 1U + *std::max_element(type_map_.begin(), type_map_.end());
        if (l_star_.size() > n_types) n_types = l_star_.size();

        if (l_star_.size() < n_types) stop("l_star too short");
        if (bias_.size() < n_types) stop("bias too short");
        if (l_int_.size() < n_types) stop("l_int too short");
        if (n_stay_.size() < n_types) stop("n_stay too short");
        if (n_ignore_.size() < n_types) stop("n_ignore too short");

        types_.reserve(n_types);
        for (uint32 i = 0; i < n_types; i++) {
            types_.push_back(OneTargetType(l_star_[i], bias_[i], l_int_[i],
                                           n_stay_[i], n_ignore_[i]));
        }

    };

    // Total number of targets:
    uint32 n_targets() const {
        return type_map_.size();
    }

    // Number of target types:
    uint32 n_types() const {
        return types_.size();
    }

    // Vector of target types:
    std::vector<uint32> types() const {
        return type_map_;
    }

    /*
     For all below, `i` should be the index for the target_xy matrix row.
     This class maps that index to the index for the target type via
     the `type_map_` vector.
     */

    double x(const uint32& i) const {
        return target_xy(i, 0);
    }
    double y(const uint32& i) const {
        return target_xy(i, 1);
    }
    double l_star(const uint32& i) const {
        return types_[type_map_[i]].l_star();
    }
    double bias(const uint32& i, const double& l_i) const {
        return types_[type_map_[i]].bias(l_i);
    }
    double abs_bias(const uint32& i, const double& l_i) const {
        return types_[type_map_[i]].abs_bias(l_i);
    }
    double l_int(const uint32& i) const {
        return types_[type_map_[i]].l_int();
    }
    uint32 n_stay(const uint32& i) const {
        return types_[type_map_[i]].n_stay();
    }
    uint32 n_ignore(const uint32& i) const {
        return types_[type_map_[i]].n_ignore();
    }


};







class ABMsimulator {

    pcg32 eng;
    const TargetInfo* target_info;
    double d;
    arma::vec x_bounds;
    arma::vec y_bounds;
    double max_size;
    uint32 n_targets;
    // Vector of distances from searcher to target(s)
    arma::vec l_vec;
    // which target (if any) are within l_star (and thus influence trajectory);
    // if none, this is set to 1+n_targets:
    uint32 wi_lstar;
    // which target (if any) is within l_int; if none, this is set to 1+n_targets:
    uint32 wi_li;
    // distance to nearest target:
    double l_min;
    // number of time points within l_int of most recent target:
    uint32 stayed;
    // max n_ignore for any type of target:
    uint32 max_n_ignore;
    // max time steps to simulate:
    uint32 max_t;
    // `visited` below indicates number of time points since the most
    // recent visit for each target (only starts after leaving within
    // l_int of the target):
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

    /*

     Update `wi_lstar` (index for target influencing movement),
     `wi_li` (index for target within l_int),
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
    // which target searcher is on at given time (-1 for none):
    std::vector<sint32> on_target;
    // whether hit a *new* target at given time:
    std::vector<bool> new_target;


    ABMsimulator(const TargetInfo& target_info_,
                 const double& d_,
                 const double& x_size,
                 const double& y_size,
                 const std::vector<std::vector<double>>& bias,
                 const std::vector<uint32>& n_ignore,
                 const uint32& max_t_,
                 const double& x0,
                 const double& y0,
                 const bool& randomize_xy0)
        : eng(),
          target_info(&target_info_),
          d(d_),
          x_bounds({0, x_size}),
          y_bounds({0, y_size}),
          max_size(std::sqrt(std::pow(x_size, 2U) + std::pow(y_size, 2U))),
          n_targets(target_info_.n_targets()),
          l_vec(target_info_.n_targets(), arma::fill::none),
          wi_lstar(),
          wi_li(),
          l_min(),
          stayed(0),
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
          x(max_t_ + 1U, arma::fill::none),
          y(max_t_ + 1U, arma::fill::none),
          on_target(max_t_ + 1U, -1),
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
        if (wi_li < n_targets && l_min <= target_info->l_int(wi_li)) {
            on_target[0] = wi_li;
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
