

#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iterator>

using namespace Rcpp;



class TargetTypesInfo {


    std::vector<uint32_t> type_map;
    std::vector<double> l_star_;
    std::vector<double> l_i_;
    std::vector<double> bias_;
    std::vector<uint32_t> n_stay_;
    std::vector<uint32_t> n_ignore_;

public:

    TargetTypesInfo(const std::vector<uint32_t>& type_map_,
                    const std::vector<double>& l_star__,
                    const std::vector<double>& l_i__,
                    const std::vector<double>& bias__,
                    const std::vector<uint32_t>& n_stay__,
                    const std::vector<uint32_t>& n_ignore__)
        : type_map(type_map_), l_star_(l_star__), l_i_(l_i__),
          bias_(bias__), n_stay_(n_stay__), n_ignore_(n_ignore__) {

        uint32_t max_idx = *std::max_element(type_map.begin(), type_map.end());
        if (l_star_.size() != max_idx+1U) stop("wrong length for l_star");
        if (l_i_.size() != max_idx+1U) stop("wrong length for l_i");
        if (bias_.size() != max_idx+1U) stop("wrong length for bias");
        if (n_stay_.size() != max_idx+1U) stop("wrong length for n_stay");
        if (n_ignore_.size() != max_idx+1U) stop("wrong length for n_ignore");

    };

    /*
     For all below, `i` should be the index for the target_xy matrix row.
     This class maps that index to the index for the target type via
     the `type_map` vector.
     */
    double l_star(const uint32_t& i) const {
        return l_star_[type_map[i]];
    }
    double l_i(const uint32_t& i) const {
        return l_i_[type_map[i]];
    }
    double bias(const uint32_t& i) const {
        return bias_[type_map[i]];
    }
    uint32_t n_stay(const uint32_t& i) const {
        return n_stay_[type_map[i]];
    }
    uint32_t n_ignore(const uint32_t& i) const {
        return n_ignore_[type_map[i]];
    }
    // Overloaded to make arma vectors:
    arma::vec l_star(const arma::uvec& i) const {
        arma::vec out(i.n_elem);
        for (uint32_t j = 0; j < i.n_elem; j++) {
            out(j) = l_star_[type_map[i(j)]];
        }
        return out;
    }
    arma::vec l_i(const arma::uvec& i) const {
        arma::vec out(i.n_elem);
        for (uint32_t j = 0; j < i.n_elem; j++) {
            out(j) = l_i_[type_map[i(j)]];
        }
        return out;
    }
    arma::vec bias(const arma::uvec& i) const {
        arma::vec out(i.n_elem);
        for (uint32_t j = 0; j < i.n_elem; j++) {
            out(j) = bias_[type_map[i(j)]];
        }
        return out;
    }
    arma::uvec n_stay(const arma::uvec& i) const {
        arma::uvec out(i.n_elem);
        for (uint32_t j = 0; j < i.n_elem; j++) {
            out(j) = n_stay_[type_map[i(j)]];
        }
        return out;
    }
    arma::uvec n_ignore(const arma::uvec& i) const {
        arma::uvec out(i.n_elem);
        for (uint32_t j = 0; j < i.n_elem; j++) {
            out(j) = n_ignore_[type_map[i(j)]];
        }
        return out;
    }


};



//
// New dx or dy that includes effects of reflecting off bound(s).
//
double reflect(const double& start,
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
double distance(const double& x1, const double& y1,
                const double& x2, const double& y2) {
    return std::sqrt(std::pow((x1 - x2), 2U) + std::pow((y1 - y2), 2U));
}



//
// Update `wi_lstar` (index / indices for targets within l_star),
// `wi_li` (index for target within l_i),
// `l_min` (distance to nearest target),
// and `l_vec` (distances to all targets),
// including ignoring recently visited target(s).
//
void update_wi_lvec(arma::uvec& wi_lstar,
                    uint32_t& wi_li,
                    double& l_min,
                    arma::vec& l_vec,
                    const double& xt,
                    const double& yt,
                    const arma::mat& target_xy,
                    const arma::uvec& visited,
                    const TargetTypesInfo& target_info,
                    const double& max_size) {

    wi_li = l_vec.n_elem + 1U;
    l_min = max_size * 2;
    // it's faster than `arma::find(...)` to get size of `wi_lstar` in first
    // loop then create it in another loop:
    uint32_t wil_size = 0U;

    for (uint32_t j = 0; j < l_vec.n_elem; j++) {
        l_vec(j) = distance(xt, yt, target_xy(j,0), target_xy(j,1));
        if (visited(j) < target_info.n_ignore(j)) continue;
        if (l_vec(j) <= target_info.l_star(j)) wil_size++;
        if (l_vec(j) < l_min) {
            l_min = l_vec(j);
            if (l_vec(j) <= target_info.l_i(j)) wi_li = j;
        }
    }

    wi_lstar.set_size(wil_size);
    if (wil_size > 0) {
        uint32_t k = 0;
        for (uint32_t j = 0; j < l_vec.n_elem; j++) {
            if (visited(j) >= target_info.n_ignore(j) &&
                l_vec(j) <= target_info.l_star(j)) {
                wi_lstar(k) = j;
                k++;
            }
        }
    }

    return;

}


//
// Same as above, but without updating `l_vec`
//
void update_wi(arma::uvec& wi_lstar,
               uint32_t& wi_li,
                double& l_min,
                const arma::vec& l_vec,
                const arma::uvec& visited,
                const TargetTypesInfo& target_info,
                const double& max_size) {

    wi_li = l_vec.n_elem + 1U;
    l_min = max_size * 2;
    // it's faster than `arma::find(...)` to get size of `wi_lstar` in first
    // loop then create it in another loop:
    uint32_t wil_size = 0U;

    for (uint32_t j = 0; j < l_vec.n_elem; j++) {
        if (visited(j) < target_info.n_ignore(j)) continue;
        if (l_vec(j) <= target_info.l_star(j)) wil_size++;
        if (l_vec(j) < l_min) {
            l_min = l_vec(j);
            if (l_vec(j) <= target_info.l_i(j)) wi_li = j;
        }
    }

    wi_lstar.set_size(wil_size);
    if (wil_size > 0) {
        uint32_t k = 0;
        for (uint32_t j = 0; j < l_vec.n_elem; j++) {
            if (visited(j) >= target_info.n_ignore(j) &&
                l_vec(j) <= target_info.l_star(j)) {
                wi_lstar(k) = j;
                k++;
            }
        }
    }

    return;

}



void check_args(const double& delta,
                const uint32_t& maxt,
                const double& x_size,
                const double& y_size,
                const arma::mat& target_xy,
                const std::vector<uint32_t>& target_types,
                const std::vector<double>& l_star,
                const std::vector<double>& l_i,
                const std::vector<double>& bias,
                const std::vector<uint32_t>& n_stay,
                const std::vector<uint32_t>& n_ignore,
                Nullable<NumericVector> xy0,
                double& x0,
                double& y0) {

    if (delta <= 0) stop("delta <= 0");
    if (x_size <= 0) stop("x_size <= 0");
    if (y_size <= 0) stop("y_size <= 0");

    if (target_types.size() == 0) stop("target_types.size() == 0");
    if (target_types.size() != target_xy.n_rows)
        stop("target_types.size() != target_xy.n_rows");

    // Create vector of unique, sorted values from `target_types`:
    std::vector<uint32_t> unq_targets = target_types;
    std::sort(unq_targets.begin(), unq_targets.end());
    std::vector<uint32_t>::iterator it;
    it = std::unique(unq_targets.begin(), unq_targets.end());
    unq_targets.resize(std::distance(unq_targets.begin(), it));

    uint32_t n_types = unq_targets.size();

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

    auto stop_i = [](std::string message, const uint32_t& i) {
        message += " on item ";
        message += std::to_string(i);
        stop(message.c_str());
    };

    for (uint32_t i = 0; i < n_types; i++) {
        if (l_star[i] <= 0) stop_i("l_star <= 0", i);
        if (l_i[i] <= 0) stop_i("l_i <= 0", i);
        if (bias[i] < -1 || bias[i] > 1) stop_i("bias < -1 || bias > 1", i);
        if (unq_targets[i] != i+1U) {
            std::string msg("target_types should contain unique values that, ");
            msg += "when sorted, are identical to a vector from 1 to ";
            msg += "the number of unique values. Yours is c(";
            for (uint32_t& ut : unq_targets) msg += std::to_string(ut) + ", ";
            msg += ").";
            stop(msg);
        }
    }

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
    // The lower triangle of this matrix contains distances between targets
    // minus the sum of their l_i values
    arma::mat target_dists(target_xy.n_rows, target_xy.n_rows, arma::fill::none);
    for (uint32_t i = 1; i < target_xy.n_rows; i++) {
        for (uint32_t j = 0; j < i; j++) {
            target_dists(i,j) = distance(target_xy(i,0), target_xy(i,1),
                                         target_xy(j,0), target_xy(j,1));
            target_dists(i,j) -= l_i[target_types[i]-1U];
            target_dists(i,j) -= l_i[target_types[j]-1U];
        }
    }
    arma::uvec lower_tri = arma::trimatl_ind(arma::size(target_dists), -1);
    arma::vec dists = target_dists(lower_tri);
    if (arma::min(dists) <= 0) {
        std::string err("targets are too close together which would ");
        err += "result in searchers interacting with >1 at a time.";
        stop(err.c_str());
    }

    if (xy0.isNotNull()) {
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




// [[Rcpp::export]]
DataFrame bias_bound_rw_cpp(const double& delta,
                            const uint32_t& maxt,
                            const double& x_size,
                            const double& y_size,
                            const arma::mat& target_xy,
                            std::vector<uint32_t> target_types,
                            const std::vector<double>& l_star,
                            const std::vector<double>& l_i,
                            const std::vector<double>& bias,
                            const std::vector<uint32_t>& n_stay,
                            const std::vector<uint32_t>& n_ignore,
                            Nullable<NumericVector> xy0 = R_NilValue) {

    double x0 = 0;
    double y0 = 0;

    // Check arguments and optionally set x0 and y0:
    check_args(delta, maxt, x_size, y_size, target_xy, target_types, l_star, l_i,
               bias, n_stay, n_ignore, xy0, x0, y0);

    // convert from R 1-based to c++ 0-based indices:
    for (uint32_t& tt : target_types) tt--;

    const TargetTypesInfo target_info(target_types, l_star, l_i, bias, n_stay,
                                      n_ignore);

    arma::vec x_bounds = {-x_size / 2, x_size / 2};
    arma::vec y_bounds = {-y_size / 2, y_size / 2};
    double max_size = std::sqrt(std::pow(x_size, 2U) + std::pow(y_size, 2U));

    // --------------------------------------------*
    // CREATE OUTPUTS
    // X and Y coordinates of searcher:
    arma::vec x(maxt + 1U);
    arma::vec y(maxt + 1U);
    x(0) = x0;
    y(0) = y0;

    // whether on target at given time:
    std::vector<bool> on_target(maxt + 1U, false);
    // whether hit a *new* target at given time:
    std::vector<bool> new_target(maxt + 1U, false);

    // --------------------------------------------*
    // OTHER USEFUL OBJECTS:
    uint32_t n_targets = target_xy.n_rows;
    // Vector of distances from searcher to target(s) (will be updated below)
    arma::vec l_vec(n_targets, arma::fill::none);
    // which target(s) are within l_star (and thus influence trajectory):
    arma::uvec wi_lstar; // (will be updated below)
    // which target (if any) is within l_i; if none, this is set to 1+n_targets:
    uint32_t wi_li; // (will be updated below)
    // distance to nearest target:
    double l_min; // (will be updated below)
    // number of time points within l_i of most recent target:
    uint32_t stayed = 0;
    // max bias for any type of target:
    double max_bias = *std::max_element(bias.begin(), bias.end());
    // `visited` below indicates number of time points since the most
    // recent visit for each target (only starts after leaving within
    // l_i of the target):
    arma::uvec visited(l_vec.n_elem, arma::fill::value(
            (*std::max_element(n_ignore.begin(), n_ignore.end())) * 2U));
    // Update l_vec, wi_lstar, and l_min:
    update_wi_lvec(wi_lstar, wi_li, l_min, l_vec, x(0), y(0), target_xy, visited,
                   target_info, max_size);
    // Update on_target and new_target if searcher starts near a target:
    if (wi_li < n_targets && l_min <= target_info.l_i(wi_li)) {
        on_target[0] = true;
        new_target[0] = true;
        visited(wi_li) = 0;
    }

    double rw_theta, rw_dx, rw_dy, dr_theta, dr_dx, dr_dy;
    bool within_li, within_lstar, keep_staying;
    arma::vec wts;


    for (uint32_t t = 0; t < maxt; t++) {
        within_li = wi_li < n_targets && l_min <= target_info.l_i(wi_li);
        keep_staying = wi_li < n_targets && stayed < target_info.n_stay(wi_li);
        // If within l_i and have NOT stayed long enough, deterministically
        // move towards target and go to next iteration
        if (within_li && keep_staying) {
            if (l_min > 0) {
                dr_theta = std::atan2((target_xy(wi_li,1) - y(t)),
                                      (target_xy(wi_li,0) - x(t)));
                dr_dx = std::min(delta, l_min) * std::cos(dr_theta);
                dr_dy = std::min(delta, l_min) * std::sin(dr_theta);
                x(t+1) = x(t) + dr_dx;
                y(t+1) = y(t) + dr_dy;
                l_min = distance(x(t+1), y(t+1),
                                 target_xy(wi_li,0), target_xy(wi_li,1));
                l_vec[wi_li] = l_min;
            } else {
                x(t+1) = x(t);
                y(t+1) = y(t);
            }
            on_target[t+1] = true;
            stayed++;
            continue;
        }
        // If within l_i and have stayed long enough, adjust some things
        // for potentially ignoring target(s), then proceed as normal:
        if (within_li && ! keep_staying) {
            stayed = 0U;
            visited(wi_li) = 0U;
            update_wi(wi_lstar, wi_li, l_min, l_vec, visited, target_info,
                      max_size);
        }

        // -----------------*
        // random walk portion of step:
        rw_theta = runif(1, 0, 2 * M_PI)(0);
        rw_dx = reflect(x[t], delta * std::cos(rw_theta), x_bounds);
        rw_dy = reflect(y[t], delta * std::sin(rw_theta), y_bounds);

        // -----------------*
        within_lstar = wi_lstar.n_elem > 0;
        if (within_lstar && max_bias > 0) {
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
                wts = target_info.l_star(wi_lstar) / l_vec(wi_lstar);
                wts /= arma::accu(wts);
                // Now calculate get weighted mean for all targets:
                double mean_bias = 0;
                dr_dx = 0;
                dr_dy = 0;
                for (uint32_t j = 0; j < wi_lstar.n_elem; j++) {
                    const arma::uword& k(wi_lstar(j));
                    const double& wt(wts(j));
                    const double& l(l_vec(k));
                    dr_theta = std::atan2((target_xy(k,1) - y(t)),
                        (target_xy(k,0) - x(t)));
                    dr_dx += (std::min(delta, l) * std::cos(dr_theta) * wt);
                    dr_dy += (std::min(delta, l) * std::sin(dr_theta) * wt);
                    mean_bias += target_info.bias(k) * wt;
                }
                // combine:
                x(t+1) = x(t) + (1 - mean_bias) * rw_dx + mean_bias * dr_dx;
                y(t+1) = y(t) + (1 - mean_bias) * rw_dy + mean_bias * dr_dy;
            }
        } else {
            // random walk otherwise:
            x(t+1) = x(t) + rw_dx;
            y(t+1) = y(t) + rw_dy;
        }
        // did it hit a new target? (also update distances for next iteration)
        update_wi_lvec(wi_lstar, wi_li, l_min, l_vec, x(t+1), y(t+1),
                       target_xy, visited, target_info, max_size);
        if (wi_li < n_targets && l_min <= target_info.l_i(wi_li)) {
            on_target[t+1] = true;
            new_target[t+1] = true;
            visited(wi_li) = 0;
            stayed = 0;
        }

        visited += 1U;

    }


    DataFrame out = DataFrame::create(
        _["t"] = Rcpp::seq(0U, maxt),
        _["x"] = x,
        _["y"] = y,
        _["on"] = on_target,
        _["hit"] = new_target);

    out.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});

    return out;
}
