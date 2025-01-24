
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

 Update `wi_lstar` (index / indices for targets within l_star),
 `wi_li` (index for target within l_i),
 `l_min` (distance to nearest target),
 and optionally `l_vec` (distances to all targets),
 including ignoring recently visited target(s).

 `t` is the index for time, indicating which values from `x` and `y` to use.

 */
void ABMsimulator::update_wi_lvec(const uint32& t,
                                  const bool& update_l_vec) {

    wi_li = n_targets + 1U;
    l_min = max_size * 2;
    // it's faster than `arma::find(...)` to get size of `wi_lstar` in first
    // loop then create it in another loop:
    uint32 wil_size = 0U;

    const double& xt(x(t));
    const double& yt(y(t));

    for (uint32 j = 0; j < n_targets; j++) {
        if (update_l_vec)
            l_vec(j) = distance(xt, yt, target_info->x(j), target_info->y(j));
        if (visited(j) < target_info->n_ignore(j)) continue;
        if (l_vec(j) <= target_info->l_star(j)) wil_size++;
        if (l_vec(j) < l_min) {
            l_min = l_vec(j);
            if (l_min <= target_info->l_i(j)) wi_li = j;
        }
    }

    wi_lstar.set_size(wil_size);
    if (wil_size > 0) {
        uint32 k = 0;
        for (uint32 j = 0; j < n_targets; j++) {
            if (visited(j) >= target_info->n_ignore(j) &&
                l_vec(j) <= target_info->l_star(j)) {
                wi_lstar(k) = j;
                k++;
            }
        }
    }

    return;

}



void ABMsimulator::stay_at_target(const uint32& t) {

    if (l_min > 0) {
        dr_theta = std::atan2((target_info->y(wi_li) - y(t)),
                              (target_info->x(wi_li) - x(t)));
        dr_dx = std::min(delta, l_min) * std::cos(dr_theta);
        dr_dy = std::min(delta, l_min) * std::sin(dr_theta);
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

        // If within l_star (but not directly on it) and there's
        // target bias, then include target-directed motion
        // (potentially including multiple targets).
        //
        // because we're weighting each target by its l_star / distance:
        wts = target_info->l_star(wi_lstar) / l_vec(wi_lstar);
        wts /= arma::accu(wts);
        // Now calculate get weighted mean for all targets:
        double mean_bias = 0;
        dr_dx = 0;
        dr_dy = 0;
        for (uint32 j = 0; j < wi_lstar.n_elem; j++) {
            const arma::uword& k(wi_lstar(j));
            const double& wt(wts(j));
            const double& l(l_vec(k));
            dr_theta = std::atan2((target_info->y(k) - y(t)),
                                  (target_info->x(k) - x(t)));
            dr_dx += (std::min(delta, l) * std::cos(dr_theta) * wt);
            dr_dy += (std::min(delta, l) * std::sin(dr_theta) * wt);
            mean_bias += target_info->bias(k) * wt;
        }
        // combine:
        x(t+1) = x(t) + (1 - mean_bias) * rw_dx + mean_bias * dr_dx;
        y(t+1) = y(t) + (1 - mean_bias) * rw_dy + mean_bias * dr_dy;
    }

}


void ABMsimulator::iterate(const uint32& t) {

    within_li = wi_li < n_targets && l_min <= target_info->l_i(wi_li);
    keep_staying = wi_li < n_targets && stayed < target_info->n_stay(wi_li);
    // If within l_i and have NOT stayed long enough, deterministically
    // move towards target and go to next iteration
    if (within_li && keep_staying) {
        stay_at_target(t);
        return;
    }
    // If within l_i and have stayed long enough, adjust some things
    // for potentially ignoring target(s), then proceed as normal:
    if (within_li && ! keep_staying) {
        stayed = 0U;
        visited(wi_li) = 0U;
        update_wi_lvec(t, false); // don't update l_vec
    }

    // -----------------*
    // random walk portion of step:
    rw_theta = runif_ab(eng, 0, 2 * M_PI);
    rw_dx = reflect(x(t), delta * std::cos(rw_theta), x_bounds);
    rw_dy = reflect(y(t), delta * std::sin(rw_theta), y_bounds);

    // -----------------*
    within_lstar = wi_lstar.n_elem > 0;
    if (within_lstar && max_bias > 0) {
        biased_movement(t);
    } else {
        // random walk otherwise:
        x(t+1) = x(t) + rw_dx;
        y(t+1) = y(t) + rw_dy;
    }
    // did it hit a new target? (also update distances for next iteration)
    update_wi_lvec(t+1, true);
    if (wi_li < n_targets && l_min <= target_info->l_i(wi_li)) {
        on_target[t+1] = wi_li;
        new_target[t+1] = true;
        visited(wi_li) = 0;
        stayed = 0;
    }

    visited += 1U;

    return;

}




void check_args(const double& delta,
                const uint32& max_t,
                const double& x_size,
                const double& y_size,
                const arma::mat& target_xy,
                std::vector<uint32>& target_types,
                const std::vector<double>& l_star,
                const std::vector<double>& l_i,
                const std::vector<double>& bias,
                const std::vector<uint32>& n_stay,
                const std::vector<uint32>& n_ignore,
                Nullable<NumericVector> xy0,
                const bool& randomize_xy0,
                double& x0,
                double& y0,
                const uint32& n_reps,
                uint32& n_threads) {

    // Check that # threads isn't too high:
    thread_check(n_threads);

    if (delta <= 0) stop("delta <= 0");
    if (max_t == 0) stop("max_t == 0");
    if (x_size <= 0) stop("x_size <= 0");
    if (y_size <= 0) stop("y_size <= 0");
    if (n_reps == 0) stop("n_reps == 0");

    if (target_types.size() == 0) stop("target_types.size() == 0");
    if (target_types.size() != target_xy.n_rows)
        stop("target_types.size() != target_xy.n_rows");

    // Create vector of unique, sorted values from `target_types`:
    std::vector<uint32> unq_targets = target_types;
    std::sort(unq_targets.begin(), unq_targets.end());
    std::vector<uint32>::iterator it;
    it = std::unique(unq_targets.begin(), unq_targets.end());
    unq_targets.resize(std::distance(unq_targets.begin(), it));

    uint32 n_types = unq_targets.size();

    if (l_star.size() != n_types)
        stop("l_star.size() != length(unique(target_types))");
    if (l_i.size() != n_types)
        stop("l_i.size() != length(unique(target_types))");
    if (bias.size() != n_types)
        stop("bias.size() != length(unique(target_types))");
    if (n_stay.size() != n_types)
        stop("n_stay.size() != length(unique(target_types))");
    if (n_ignore.size() != n_types)
        stop("n_ignore.size() != length(unique(target_types))");

    auto stop_i = [](std::string message, const uint32& i) {
        message += " on item ";
        message += std::to_string(i);
        stop(message.c_str());
    };

    for (uint32 i = 0; i < n_types; i++) {
        if (l_star[i] <= 0) stop_i("l_star <= 0", i);
        if (l_i[i] <= 0) stop_i("l_i <= 0", i);
        if (bias[i] < -1 || bias[i] > 1) stop_i("bias < -1 || bias > 1", i);
        if (unq_targets[i] != i+1U) {
            std::string msg("target_types should contain unique values that, ");
            msg += "when sorted, are identical to a vector from 1 to ";
            msg += "the number of unique values. Yours is c(";
            for (uint32& ut : unq_targets) msg += std::to_string(ut) + ", ";
            msg += ").";
            stop(msg);
        }
    }

    // convert from R 1-based to c++ 0-based indices:
    for (uint32& tt : target_types) tt--;

    arma::vec x_bounds = {-1, 1};
    x_bounds *= (x_size / 2);
    arma::vec y_bounds = {-1, 1};
    y_bounds *= (y_size / 2);

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
    // of their l_i values
    double dist, sum_l_i;
    for (uint32 i = 1; i < target_xy.n_rows; i++) {
        for (uint32 j = 0; j < i; j++) {
            dist = distance(target_xy(i,0), target_xy(i,1),
                            target_xy(j,0), target_xy(j,1));
            sum_l_i = l_i[target_types[i]] + l_i[target_types[j]];
            if (dist <= sum_l_i) {
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
                        const uint32& n_reps,
                        const std::vector<ABMsimulator>& simmers,
                        const TargetInfo& target_info) {

    std::vector<uint32> tt = target_info.types();

    std::vector<int> rep;
    std::vector<int> time;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<int> tar;
    std::vector<int> type;
    std::vector<bool> hit;
    rep.reserve((max_t+1U) * n_reps);
    time.reserve((max_t+1U) * n_reps);
    x.reserve((max_t+1U) * n_reps);
    y.reserve((max_t+1U) * n_reps);
    tar.reserve((max_t+1U) * n_reps);
    type.reserve((max_t+1U) * n_reps);
    hit.reserve((max_t+1U) * n_reps);

    for (uint32 i = 0; i < n_reps; i++) {
        for (uint32 j = 0; j <= max_t; j++) {
            rep.push_back(i);
            time.push_back(j);
            x.push_back(simmers[i].x[j]);
            y.push_back(simmers[i].y[j]);
            tar.push_back(simmers[i].on_target[j] + 1);
            if (simmers[i].on_target[j] < 0) {
                type.push_back(0);
            } else type.push_back(tt[simmers[i].on_target[j]] + 1);
            hit.push_back(simmers[i].new_target[j]);
        }
    }

    DataFrame out = DataFrame::create(
        _["rep"] = rep,
        _["time"] = time,
        _["x"] = x,
        _["y"] = y,
        _["tar"] = tar,
        _["type"] = type,
        _["hit"] = hit);

    out.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});

    return out;

}







//' Simulate searchers seeking targets.
//'
//' @param delta Single numeric indicating search rate of searcher.
//' @param max_t Single integer indicating the number of time steps to simulate.
//' @param x_size Single numeric indicating x dimension of search area.
//'     Bounds will be `c(-x_size/2, x_size/2)`.
//' @param y_size Single numeric indicating y dimension of search area.
//'     Bounds will be `c(-y_size/2, y_size/2)`.
//' @param target_xy Two-column numeric matrix containing x and y coordinates
//'     for targets.
//' @param target_types Integer vector with the same number of items as the
//'     number of rows in `target_xy` indicating the type of target each
//'     target is. Targets can vary in their `l_star`, `l_i`, `bias`,
//'     `n_stay`, and `n_ignore` parameters (see descriptions below).
//'     This vector should consist of integers from 1 to the number of items
//'     in the arguments `l_star`, `l_i`, `bias`, `n_stay`, and `n_ignore`.
//' @param l_star Numeric vector indicating, for each target type,
//'     the distance from searcher to target that causes targets to bias
//'     searcher movement.
//'     Its length should be equal to the number of unique type(s) of targets.
//' @param l_i Numeric vector indicating, for each target type,
//'     the distance from searcher to target that causes an interaction.
//'     Its length should be equal to the number of unique type(s) of targets.
//' @param bias Numeric vector indicating, for each target type,
//'     the bias the target causes searcher movement once they're `<= l_star`
//'     away from it. Values range from -1 to 1,
//'     with negative values resulting in searchers being repelled.
//'     Values of +1 (-1) cause searchers to move directly towards (away from)
//'     targets once within `l_star`.
//'     Values nearer to zero cause movement that is more similar to
//'     a random walk.
//'     Its length should be equal to the number of unique type(s) of targets.
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
//' @param n_reps Single integer indicating the number of independent
//'     repetitions to run. This can be thought of as the number of searchers
//'     on the landscape if searchers are independent of one another.
//'     Defaults to `1L`.
//' @param show_progress Single logical for whether to show progress bar.
//'     Defaults to `FALSE`.
//' @param n_threads Single integer for the number of threads to use.
//'     For now, this is ignored because the multithreading is causing a
//'     strange warning I cannot explain.
//'     Defaults to `1L`.
//'
//' @returns A tibble with the columns `rep` (repetition number),
//'     `time` (time), `x` (x coordinate), `y` (y coordinate),
//'     `tar` (which target is searcher interacting with (within `l_i`)?),
//'     `type` (which target type is searcher interacting with?),
//'      and
//'     `hit` (logical - is searcher interacting with a new target?).
//'     Columns `tar` and `type` are `0` if the searcher is not on any targets.
//'
//'
//' @export
//'
//[[Rcpp::export]]
DataFrame searcher_sims(const double& delta,
                        const uint32& max_t,
                        const double& x_size,
                        const double& y_size,
                        const arma::mat& target_xy,
                        std::vector<uint32> target_types,
                        const std::vector<double>& l_star,
                        const std::vector<double>& l_i,
                        const std::vector<double>& bias,
                        const std::vector<uint32>& n_stay,
                        const std::vector<uint32>& n_ignore,
                        Nullable<NumericVector> xy0 = R_NilValue,
                        const bool& randomize_xy0 = true,
                        const uint32& n_reps = 1,
                        const bool& show_progress = false,
                        uint32 n_threads = 1) {

    double x0 = 0;
    double y0 = 0;


    // Check arguments and optionally set x0 and y0:
    check_args(delta, max_t, x_size, y_size, target_xy, target_types, l_star, l_i,
               bias, n_stay, n_ignore, xy0, randomize_xy0, x0, y0,
               n_reps, n_threads);

    const TargetInfo target_info(target_xy, target_types, l_star,
                                 l_i, bias, n_stay, n_ignore);

    // Create one simulator object per rep.
    // This also generates a seeded RNG for each, using R's runif(...) for
    // reproducibility.
    std::vector<ABMsimulator> simmers;
    simmers.reserve(n_reps);
    for (uint32 i = 0; i < n_reps; i++) {
        simmers.push_back(ABMsimulator(target_info, delta, x_size, y_size, bias,
                                       n_ignore, max_t, x0, y0, randomize_xy0));
    }

    RcppThread::ProgressBar prog_bar(n_reps * max_t, 1);

    // Parallelized loop
    RcppThread::parallelFor(0, n_reps, [&] (uint32 i) {
        simmers[i].run(prog_bar, show_progress);
    }, n_threads);


    DataFrame out = create_output(max_t, n_reps, simmers, target_info);

    return out;
}
