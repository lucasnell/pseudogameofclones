
#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <pcg/pcg_random.hpp>   // pcg prng
#include <RcppThread.h>         // multithreading

#include "pseudogameofclones_types.hpp"  // integer types
#include "pcg.hpp"              // runif_01 fxn
#include "util.hpp"              // thread_check fxn
#include "abm.hpp"              // classes for ABM



using namespace Rcpp;








/*

 Update `wi_lstar` (index for targets affecting movement),
 `wi_li` (index for target within l_int),
 `l_min` (distance to nearest target),
 and optionally `l_vec` (distances to all targets),
 including ignoring recently visited target(s).

 `t` is the index for time, indicating which values from `x` and `y` to use.

 */
void ABMsimulator::update_wi_lvec(const uint32& t,
                                  const bool& update_l_vec) {

    wi_li = n_targets + 1U;
    wi_lstar = n_targets + 1U;
    l_min = max_size * 2;

    // values used to choose wi_lstar:
    double max_ab = 0;              // abs(bias)
    double min_l = max_size * 2;    // l_i
    double ab_j;
    bool new_wi_lstar;
    std::vector<uint32> ties;
    ties.reserve(std::max(10U, n_targets / 2U));

    const double& xt(x(t));
    const double& yt(y(t));

    for (uint32 j = 0; j < n_targets; j++) {
        if (update_l_vec)
            l_vec(j) = distance(xt, yt, target_info->x(j), target_info->y(j));
        if (visited(j) < target_info->n_ignore(j)) continue;
        if (l_vec(j) <= target_info->l_star(j)) {
            ab_j = target_info->abs_bias(j, l_vec(j));
            if (ab_j >= max_ab) {
                // Update to new target if abs(bias) is greater OR if it's
                // the same and closer:
                new_wi_lstar = ab_j > max_ab || (ab_j == max_ab &&
                    l_vec(j) < min_l);
                // If both of those are a tie, then choose randomly (have to
                // do this later to be unbiased):
                if (!new_wi_lstar && ab_j == max_ab && l_vec(j) == min_l) {
                    ties.push_back(j);
                }
                if (new_wi_lstar) {
                    wi_lstar = j;
                    max_ab = ab_j;
                    min_l = l_vec(j);
                }
            }
        }
        if (l_vec(j) < l_min) {
            l_min = l_vec(j);
            if (l_min <= target_info->l_int(j)) wi_li = j;
        }
    }

    // When choosing wi_lstar, if there are ties for bias and distance,
    // choose among them randomly:
    if (ties.size() > 0) {
        // Filter for indices that are tied with the final minimum values:
        std::vector<uint32> real_ties;
        real_ties.reserve(ties.size() + 1);
        for (const uint32& j : ties) {
            if (l_vec(j) == min_l && target_info->abs_bias(j, l_vec(j)) == max_ab) {
                real_ties.push_back(j);
            }
        }
        if (real_ties.size() > 0) {
            // Also include the item that was assigned to wi_lstar:
            real_ties.push_back(wi_lstar);
            uint32 rnd = runif_01(eng) * real_ties.size();
            wi_lstar = real_ties[rnd];
        }
    }

    return;

}



void ABMsimulator::stay_at_target(const uint32& t) {

    if (l_min > 0) {
        dr_theta = std::atan2((target_info->y(wi_li) - y(t)),
                              (target_info->x(wi_li) - x(t)));
        dr_dx = std::min(d, l_min) * std::cos(dr_theta);
        dr_dy = std::min(d, l_min) * std::sin(dr_theta);
        x(t+1) = x(t) + dr_dx;
        y(t+1) = y(t) + dr_dy;
        l_min = distance(x(t+1), y(t+1),
                         target_info->x(wi_li), target_info->y(wi_li));
        l_vec[wi_li] = l_min;
    } else {
        x(t+1) = x(t);
        y(t+1) = y(t);
    }
    on_target[t+1] = wi_li;
    stayed++;

    return;

}


void ABMsimulator::biased_movement(const uint32& t) {

    if (l_min == 0) {
        // If it's directly on the target, then the directed portion
        // would be to stay in place (this if-else avoids NaNs):
        x(t+1) = x(t);
        y(t+1) = y(t);
    } else {

        double dx, dy;

        // If within l_star (but not directly on it), then include
        // target-directed motion (note: all biases must be != 0).
        if (wi_lstar < n_targets) {

            const uint32& i(wi_lstar);
            const double& l(l_vec(i));
            double b = target_info->bias(i, l_min);
            double ab = target_info->abs_bias(i, l_min);
            if (b < -1 || ab < 0) stop("l_min too big for chosen target");

            // this is used to reverse direction if bias is negative:
            double bias_sign = sgn<double>(b);
            double delta = (b > 0 && l < d) ? l : d;
            dr_theta = std::atan2(bias_sign * (target_info->y(i) - y(t)),
                                  bias_sign * (target_info->x(i) - x(t)));
            dr_dx = delta * std::cos(dr_theta);
            dr_dy = delta * std::sin(dr_theta);

            // combine:
            dx = (1 - ab) * rw_dx + ab * dr_dx;
            dy = (1 - ab) * rw_dy + ab * dr_dy;

        } else {

            // just do random walk if bias is zero:
            dx = rw_dx;
            dy = rw_dy;

        }

        // Negative biases and random walks can result in exceeding bounds,
        // so pass to `reflect`:
        x(t+1) = reflect(x(t), dx, x_bounds);
        y(t+1) = reflect(y(t), dy, y_bounds);
    }

}


void ABMsimulator::iterate(const uint32& t) {

    within_li = wi_li < n_targets && l_min <= target_info->l_int(wi_li);
    keep_staying = wi_li < n_targets && stayed < target_info->n_stay(wi_li);
    // If within l_int and have NOT stayed long enough, deterministically
    // move towards target and go to next iteration
    if (within_li && keep_staying) {
        stay_at_target(t);
        return;
    }
    // If within l_int and have stayed long enough, adjust some things
    // for potentially ignoring target(s), then proceed as normal:
    if (within_li && ! keep_staying) {
        stayed = 0U;
        visited(wi_li) = 0U;
        update_wi_lvec(t, false); // don't update l_vec
    }

    // -----------------*
    // random walk portion of step:
    rw_theta = runif_ab(eng, 0, 2 * M_PI);
    rw_dx = d * std::cos(rw_theta);
    rw_dy = d * std::sin(rw_theta);

    // -----------------*
    within_lstar = wi_lstar < n_targets;
    if (within_lstar) {
        biased_movement(t);
    } else {
        // random walk otherwise:
        x(t+1) = reflect(x(t), rw_dx, x_bounds);
        y(t+1) = reflect(y(t), rw_dy, y_bounds);
    }
    // did it hit a new target? (also update distances for next iteration)
    update_wi_lvec(t+1, true);
    if (wi_li < n_targets && l_min <= target_info->l_int(wi_li)) {
        on_target[t+1] = wi_li;
        new_target[t+1] = true;
        visited(wi_li) = 0;
        stayed = 0;
    }

    visited += 1U;

    return;

}




void check_args(const double& d,
                const uint32& max_t,
                const double& x_size,
                const double& y_size,
                const arma::mat& target_xy,
                std::vector<uint32>& target_types,
                const std::vector<std::vector<double>>& l_star,
                const std::vector<std::vector<double>>& bias,
                const std::vector<double>& l_int,
                const std::vector<uint32>& n_stay,
                const std::vector<uint32>& n_ignore,
                Nullable<NumericVector> xy0,
                const bool& randomize_xy0,
                double& x0,
                double& y0,
                const uint32& n_searchers,
                uint32& n_threads) {

    // Check that # threads isn't too high:
    thread_check(n_threads);

    if (d <= 0) stop("d <= 0");
    if (max_t == 0) stop("max_t == 0");
    if (x_size < 2) stop("x_size < 2");
    if (y_size < 2) stop("y_size < 2");
    if (n_searchers == 0) stop("n_searchers == 0");

    if (target_types.size() == 0) stop("target_types.size() == 0");
    if (target_types.size() != target_xy.n_rows)
        stop("target_types.size() != target_xy.n_rows");

    uint32 n_types = l_star.size();

    if (n_types == 0) stop("length(l_star) == 0");
    if (l_int.size() != n_types) stop("length(l_int) != length(l_star)");
    if (bias.size() != n_types) stop("length(bias) != length(l_star)");
    if (n_stay.size() != n_types) stop("length(n_stay) != length(l_star)");
    if (n_ignore.size() != n_types) stop("length(n_ignore) != length(l_star)");

    auto stop_i = [](std::string message, const uint32& i) {
        message += " on item ";
        message += std::to_string(i);
        stop(message.c_str());
    };

    for (uint32 i = 0; i < n_types; i++) {
        if (l_star[i].size() != bias[i].size()) {
            stop_i("length(l_star[[i]]) != length(bias[[i]])", i);
        }
        for (uint32 j = 0; j < l_star[i].size(); j++) {
            if (l_star[i][j] <= 0) stop_i("l_star <= 0", i);
            if (bias[i][j] < -1 || bias[i][j] > 1)
                stop_i("bias < -1 || bias > 1", i);
            if (bias[i][j] == 0) stop_i("bias == 0", i);
        }
        if (l_int[i] <= 0) stop_i("l_int <= 0", i);
    }

    // Check values and convert from R 1-based to c++ 0-based indices:
    for (uint32& tt : target_types) {
        if (tt > n_types) stop("target_types contains values > length(l_star)");
        tt--;
    }

    arma::vec x_bounds = {0, x_size};
    arma::vec y_bounds = {0, y_size};

    if (target_xy.n_cols != 2) stop("target_xy.n_cols != 2");
    if (target_xy.n_rows < 1) stop("target_xy.n_rows < 1");
    if (! arma::all(target_xy.col(0) > x_bounds[0]))
        stop("! arma::all(target_xy.col(0) > x_bounds[0])");
    if (! arma::all(target_xy.col(0) < x_bounds[1]))
        stop("! arma::all(target_xy.col(0) < x_bounds[1])");
    if (! arma::all(target_xy.col(1) > y_bounds[0]))
        stop("! arma::all(target_xy.col(1) > y_bounds[0])");
    if (! arma::all(target_xy.col(1) < y_bounds[1]))
        stop("! arma::all(target_xy.col(1) < y_bounds[1])");


    // Lastly check for overlapping bounds.
    // I'll check that no distance between targets is less than the sum
    // of their l_int values
    double dist, sum_l_int;
    for (uint32 i = 1; i < target_xy.n_rows; i++) {
        for (uint32 j = 0; j < i; j++) {
            dist = distance(target_xy(i,0), target_xy(i,1),
                            target_xy(j,0), target_xy(j,1));
            sum_l_int = l_int[target_types[i]] + l_int[target_types[j]];
            if (dist <= sum_l_int) {
                std::string err("targets are too close together which would ");
                err += "result in searchers interacting with >1 at a time.";
                stop(err.c_str());
            }
        }
    }

    if (xy0.isNotNull()) {
        if (randomize_xy0) {
            std::string msg("If `randomize_xy0 = TRUE`, this will negate ");
            msg += "your provided starting values for `xy0`";
            stop(msg.c_str());
        }
        NumericVector xy0_vec(xy0);
        if (xy0_vec.size() != 2) {
            stop("xy0_vec must be NULL or a length-2 numeric vector");
        }
        if (xy0_vec(0) <= x_bounds(0) || xy0_vec(0) >= x_bounds(1) ||
            xy0_vec(1) <= y_bounds(0) || xy0_vec(1) >= y_bounds(1)) {
            stop("xy0_vec cannot contain items that reach or exceed bounds");
        }
        x0 = xy0_vec(0);
        y0 = xy0_vec(1);
    }

    return;

}





DataFrame create_output(const uint32& max_t,
                        const uint32& n_searchers,
                        const std::vector<ABMsimulator>& simmers,
                        const TargetInfo& target_info,
                        const bool& summarize) {

    DataFrame out;


    std::vector<uint32> tt = target_info.types();
    std::vector<int> searcher;
    std::vector<int> type;

    if (summarize) {

        uint32 n_types = target_info.n_types();
        std::vector<int> on;
        std::vector<int> hit;

        searcher.reserve(n_types * n_searchers);
        type.reserve(n_types * n_searchers);
        on.reserve(n_types * n_searchers);
        hit.reserve(n_types * n_searchers);

        // Fill with correct `searcher` and `type` values, then with zeros for
        // `on` and `hit` (which will be updated below):
        for (uint32 i = 0; i < n_searchers; i++) {
            for (uint32 j = 0; j < n_types; j++) {
                searcher.push_back(i+1);
                type.push_back(j+1);
                on.push_back(0);
                hit.push_back(0);
            }
        }

        // Now fill for `on` and `hit` values
        uint32 k, j;
        for (uint32 i = 0; i < n_searchers; i++) {
            for (uint32 t = 0; t <= max_t; t++) {
                if (simmers[i].on_target[t] >= 0) {
                    j = tt[simmers[i].on_target[t]];
                    k = j + i * n_types;
                    on[k]++;
                    if (simmers[i].new_target[t]) hit[k]++;
                }
            }
        }

        out = DataFrame::create(_["searcher"] = searcher,
                                _["type"] = type,
                                _["on"] = on,
                                _["hit"] = hit);

    } else {

        std::vector<int> time;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<int> tar;
        std::vector<bool> hit;
        searcher.reserve((max_t+1U) * n_searchers);
        time.reserve((max_t+1U) * n_searchers);
        x.reserve((max_t+1U) * n_searchers);
        y.reserve((max_t+1U) * n_searchers);
        tar.reserve((max_t+1U) * n_searchers);
        type.reserve((max_t+1U) * n_searchers);
        hit.reserve((max_t+1U) * n_searchers);

        for (uint32 i = 0; i < n_searchers; i++) {
            for (uint32 t = 0; t <= max_t; t++) {
                searcher.push_back(i+1);
                time.push_back(t);
                x.push_back(simmers[i].x[t]);
                y.push_back(simmers[i].y[t]);
                tar.push_back(simmers[i].on_target[t] + 1);
                if (simmers[i].on_target[t] < 0) {
                    type.push_back(0);
                } else type.push_back(tt[simmers[i].on_target[t]] + 1);
                hit.push_back(simmers[i].new_target[t]);
            }
        }

        out = DataFrame::create(
            _["searcher"] = searcher,
            _["time"] = time,
            _["x"] = x,
            _["y"] = y,
            _["tar"] = tar,
            _["type"] = type,
            _["hit"] = hit);

    }


    out.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});

    return out;

}







//' Simulate searchers seeking targets.
//'
//' @param d Single numeric indicating search rate of searcher.
//' @param max_t Single integer indicating the number of time steps to simulate.
//' @param x_size Single numeric indicating x dimension of search area.
//'     Bounds will be `c(0, x_size)`.
//' @param y_size Single numeric indicating y dimension of search area.
//'     Bounds will be `c(0, y_size)`.
//' @param target_xy Two-column numeric matrix containing x and y coordinates
//'     for targets.
//' @param target_types Integer vector with the same number of items as the
//'     number of rows in `target_xy` indicating the type of target each
//'     target is. Targets can vary in their `l_star`, `l_int`, `bias`,
//'     `n_stay`, and `n_ignore` parameters (see descriptions below).
//'     This vector should consist of integers from 1 to the number of items
//'     in the arguments `l_star`, `l_int`, `bias`, `n_stay`, and `n_ignore`.
//'     It's allowed that some target types don't show up in this vector
//'     since this can be useful for simulations where you remove or replace
//'     target type(s).
//' @param l_star List where each element is a numeric vector indicating,
//'     for each target type, the distance(s) from searcher to target that
//'     causes targets to bias searcher movement.
//'     Its length should be equal to the number of unique type(s) of targets.
//'     Target types can have multiple `l_star` values if they have, for
//'     example, both a virus that attracts searchers and an epiphytic
//'     bacteria that repels them. If the cues from the virus and bacteria
//'     happen at different spatial scales, then this would result in
//'     multiple `l_star` values. It would also result in multiple `bias`
//'     values for targets of this type.
//'     For multiple values of `l_star` and `bias`, it's assumed that they
//'     are ordered the same (i.e., `l_star[[i]][[j]]` coincides with
//'     `bias[[i]][[j]]`).
//'     It's also required that within the same target type, higher values
//'     of `l_star` coincide with lower values of `abs(bias)`.
//' @param l_int Numeric vector indicating, for each target type,
//'     the distance from searcher to target that causes an interaction.
//'     Its length should be equal to the number of unique type(s) of targets.
//' @param bias List where each item is a numeric vector indicating, for
//'     each target type, the bias(es) the target causes searcher movement
//'     once they're `<= l_star` away from it. Values range from -1 to 1,
//'     with negative values resulting in searchers being repelled.
//'     Values of +1 (-1) cause searchers to move directly towards (away from)
//'     targets once within `l_star`.
//'     Values nearer to zero cause movement that is more similar to
//'     a random walk.
//'     Its length should be equal to the number of unique type(s) of targets.
//'     See the description for parameter `l_star` for an example of why you
//'     might need multiple `bias` values for one target type.
//' @param n_stay Numeric vector indicating, for each target type,
//'     the number of time steps searchers stay at targets when they interact
//'     with them.
//'     Its length should be equal to the number of unique type(s) of targets.
//' @param n_ignore Numeric vector indicating, for each target type,
//'     the number of time steps searchers ignore targets after
//'     interacting with them.
//'     Its length should be equal to the number of unique type(s) of targets.
//' @param xy0 Numeric vector of length 2 indicating starting x and y
//'     coordinates for searcher(s). When `xy0 = NULL` and
//'     `randomize_xy0 = FALSE`, searcher(s) start at the middle of the
//'     search space (x = 0, y = 0).
//'     This function throws an error when `xy0` is provided and
//'     `randomize_xy0 = TRUE` because these are conflicting.
//'     Defaults to `NULL`.
//' @param randomize_xy0 Single logical for whether to randomize starting
//'     coordinates for searcher(s). If `TRUE`, x and y coordinates are
//'     generated from a random uniform distribution ranging from the lower
//'     to upper bound for each dimension.
//'     Defaults to `TRUE`.
//' @param n_searchers Single integer indicating the number of independent
//'     searchers to simulate.
//'     Defaults to `1L`.
//' @param summarize Single logical for whether to summarize output by
//'     searcher. See below for details on how this changes the output.
//'     Defaults to `FALSE`.
//' @param show_progress Single logical for whether to show progress bar.
//'     Defaults to `FALSE`.
//' @param n_threads Single integer for the number of threads to use.
//'     Ignored if `n_searchers == 1`.
//'     Defaults to `1L`.
//'
//' @returns If `summarize = FALSE`, then it outputs a tibble with the columns
//'     `searcher` (searcher number),
//'     `time` (time), `x` (x coordinate), `y` (y coordinate),
//'     `tar` (which target is searcher interacting with (within `l_int`)?),
//'     `type` (which target type is searcher interacting with?),
//'      and
//'     `hit` (logical - is searcher interacting with a new target?).
//'     Columns `tar` and `type` are `0` if the searcher is not on any targets.
//'     If `summarize = TRUE`, then it outputs a tibble with the columns
//'     `searcher` (searcher number),
//'     `type` (which target type is searcher interacting with?),
//'     `on` (integer - how many time steps did the searcher spend on this
//'     type of target?).
//'     and
//'     `hit` (integer - how many times did this searcher interact with
//'     a new target of this type?).
//'
//'
//' @export
//'
//[[Rcpp::export]]
DataFrame searcher_sims(const double& d,
                        const uint32& max_t,
                        const double& x_size,
                        const double& y_size,
                        const arma::mat& target_xy,
                        std::vector<uint32> target_types,
                        const std::vector<std::vector<double>>& l_star,
                        const std::vector<std::vector<double>>& bias,
                        const std::vector<double>& l_int,
                        const std::vector<uint32>& n_stay,
                        const std::vector<uint32>& n_ignore,
                        Nullable<NumericVector> xy0 = R_NilValue,
                        const bool& randomize_xy0 = true,
                        const uint32& n_searchers = 1,
                        const bool& summarize = false,
                        const bool& show_progress = false,
                        uint32 n_threads = 1) {

    double x0 = x_size / 2;
    double y0 = y_size / 2;

    // Check arguments and optionally set x0 and y0:
    check_args(d, max_t, x_size, y_size, target_xy, target_types, l_star, bias,
               l_int, n_stay, n_ignore, xy0, randomize_xy0, x0, y0,
               n_searchers, n_threads);

    const TargetInfo target_info(target_xy, target_types, l_star, bias,
                                 l_int, n_stay, n_ignore);

    // Create one simulator object per searcher.
    // This also generates a seeded RNG for each, using R's runif(...) for
    // reproducibility.
    std::vector<ABMsimulator> simmers;
    simmers.reserve(n_searchers);
    for (uint32 i = 0; i < n_searchers; i++) {
        simmers.push_back(ABMsimulator(target_info, d, x_size, y_size, bias,
                                       n_ignore, max_t, x0, y0, randomize_xy0));
    }

    RcppThread::ProgressBar prog_bar(n_searchers * max_t, 1);

    if (n_threads > 1U && n_searchers > 1U) {
        RcppThread::parallelFor(0, n_searchers, [&] (uint32 i) {
            simmers[i].run(prog_bar, show_progress);
        }, n_threads);
    } else {
        for (uint32 i = 0; i < n_searchers; i++) {
            simmers[i].run(prog_bar, show_progress);
        }
    }


    DataFrame out = create_output(max_t, n_searchers, simmers, target_info,
                                  summarize);

    return out;
}
