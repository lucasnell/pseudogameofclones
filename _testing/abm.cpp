
//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <vector>
#include <math.h>

using namespace Rcpp;



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
// Update `i` (index for closest target), `l` (distance to nearest target),
// and `l_vec` (distances to all targets),
// including ignoring recently visited target(s).
//
void update_i_l_lvec(uint32_t& i,
                     double& l,
                     arma::vec& l_vec,
                     const double& xt,
                     const double& yt,
                     const arma::mat& target_xy,
                     const arma::uvec& visited,
                     const uint32_t& n_ignore,
                     const double& max_size) {

    i = 0;
    l = 2000000; // max_size * 2;

    for (uint32_t j = 0; j < l_vec.n_elem; j++) {
        l_vec(j) = distance(xt, yt, target_xy(j,0), target_xy(j,1));
        if (visited(j) < n_ignore) continue;
        if (l_vec(j) < l) {
            l = l_vec(j);
            i = j;
        }
    }

    return;

}


//
// Same as above, but without updating `l_vec`
//
void update_i_l(uint32_t& i,
                double& l,
                const arma::vec& l_vec,
                const arma::uvec& visited,
                const uint32_t& n_ignore,
                const double& max_size) {

    i = 0;
    l = 2000000; // max_size * 2;

    for (uint32_t j = 0; j < l_vec.n_elem; j++) {
        if (visited(j) < n_ignore) continue;
        if (l_vec(j) < l) {
            l = l_vec(j);
            i = j;
        }
    }

    return;

}




void check_args(const double& delta,
                const uint32_t& maxt,
                const double& x_size,
                const double& y_size,
                const arma::mat& target_xy,
                const double& l_star,
                const double& l_i,
                const double& bias,
                const uint32_t& n_stay,
                const uint32_t& n_ignore) {

    if (delta <= 0) stop("delta <= 0");
    if (x_size <= 0) stop("x_size <= 0");
    if (y_size <= 0) stop("y_size <= 0");
    if (l_star <= 0) stop("l_star <= 0");
    if (l_i <= 0) stop("l_i <= 0");
    if (bias < 0 || bias > 1) stop("bias < 0 || bias > 1");

    arma::vec x_bounds = {-1, 1};
    x_bounds *= (x_size / 2);
    arma::vec y_bounds = {-1, 1};
    y_bounds *= (y_size / 2);

    if (target_xy.n_cols != 2) stop("target_xy.n_cols != 2");
    if (target_xy.n_rows < 1) stop("target_xy.n_rows < 1");
    if (! arma::all(target_xy.col(0) > x_bounds[0])) {
        stop("! arma::all(target_xy.col(0) > x_bounds[0])");
    }
    if (! arma::all(target_xy.col(0) < x_bounds[1])) {
        stop("! arma::all(target_xy.col(0) < x_bounds[1])");
    }
    if (! arma::all(target_xy.col(1) > y_bounds[0])) {
        stop("! arma::all(target_xy.col(1) > y_bounds[0])");
    }
    if (! arma::all(target_xy.col(1) < y_bounds[1])) {
        stop("! arma::all(target_xy.col(1) < y_bounds[1])");
    }
    uint32_t n_targets = target_xy.n_rows;
    // Check for overlapping bounds:
    arma::mat target_dists(n_targets, n_targets, arma::fill::zeros);
    for (uint32_t i = 1; i < n_targets; i++) {
        for (uint32_t j = 0; j < i; j++) {
            target_dists(i,j) = distance(target_xy(i,0), target_xy(i,1),
                         target_xy(j,0), target_xy(j,1));
        }
    }
    arma::vec dists = target_dists(arma::trimatl_ind(arma::size(target_dists), -1));
    if (arma::min(dists) <= l_i) {
        std::string err("targets are too close together which would ");
        err += "result in searchers interacting with >1 at a time.";
        stop(err.c_str());
    }

    return;

}




// [[Rcpp::export]]
DataFrame bias_bound_rw_cpp(const double& delta,
                            const uint32_t& maxt,
                            const double& x_size,
                            const double& y_size,
                            const arma::mat& target_xy,
                            const double& l_star,
                            const double& l_i,
                            const double& bias,
                            const uint32_t& n_stay,
                            Nullable<IntegerVector> n_ignore = R_NilValue,
                            const bool& random_xy0 = true) {

    uint32_t n_ignore_ = static_cast<uint32_t>(std::ceil(2.0 * l_star / delta));
    if (n_ignore.isNotNull()) {
        IntegerVector niv(n_ignore);
        if (niv.size() != 1 || niv(0) < 0) {
            stop("n_ignore must be NULL or a single non-negative integer");
        }
        n_ignore_ = niv(0);
    }

    check_args(delta, maxt, x_size, y_size, target_xy, l_star, l_i,
               bias, n_stay, n_ignore_);

    arma::vec x_bounds = {-1, 1};
    x_bounds *= (x_size / 2);
    arma::vec y_bounds = {-1, 1};
    y_bounds *= (y_size / 2);
    double max_size = std::sqrt(std::pow(x_size, 2U) + std::pow(y_size, 2U));

    // --------------------------------------------*
    // CREATE OUTPUTS
    // X and Y coordinates of searcher:
    arma::vec x(maxt + 1U);
    arma::vec y(maxt + 1U);
    if (random_xy0) {
        x(0) = runif(1, x_bounds(0), x_bounds(1))(0);
        y(0) = runif(1, y_bounds(0), y_bounds(1))(0);
    } else{
        x(0) = 0.0;
        y(0) = 0.0;
    }
    // whether on target at given time:
    std::vector<bool> on_target(maxt + 1U, false);
    // whether hit a *new* target at given time:
    std::vector<bool> new_target(maxt + 1U, false);

    // --------------------------------------------*
    // OTHER USEFUL OBJECTS:
    uint32_t n_targets = target_xy.n_rows;
    // Vector of distances from searcher to target(s) (will be updated below)
    arma::vec l_vec(n_targets, arma::fill::none);
    // which target is closest (choose first one if there are ties):
    uint32_t i; // (will be updated below)
    // distance to nearest target:
    double l; // (will be updated below)
    // number of time points within l_i of most recent target:
    uint32_t stayed = 0;
    // `visited` below indicates number of time points since the most
    // recent visit for each target (only starts after leaving within
    // l_i of the target):
    arma::uvec visited(l_vec.n_elem, arma::fill::value(n_ignore_ * 2U));
    // Update l_vec, i, and l:
    update_i_l_lvec(i, l, l_vec, x(0), y(0), target_xy, visited,
                    n_ignore_, max_size);
    // Update on_target and new_target if searcher starts near a target:
    if (l <= l_i) {
        on_target[0] = true;
        new_target[0] = true;
        visited(i) = 0;
    }

    double rw_theta, rw_dx, rw_dy, dr_theta, dr_dx, dr_dy;


    for (uint32_t t = 0; t < maxt; t++) {
        // If within l_i and have NOT stayed long enough, stay in the spot
        // and go to next iteration
        if (l <= l_i && stayed < n_stay) {
            x(t+1) = x(t);
            y(t+1) = y(t);
            on_target[t+1] = true;
            stayed++;
            continue;
        }
        // If within l_i and have stayed long enough, adjust some things
        // for potentially ignoring target(s), then proceed as normal:
        if (l <= l_i && stayed >= n_stay) {
            stayed = 0U;
            visited(i) = 0U;
            update_i_l(i, l, l_vec, visited, n_ignore_, max_size);
        }

        // -----------------*
        // random walk portion of step:
        rw_theta = runif(1, 0, 2 * M_PI)(0);
        rw_dx = reflect(x[t], delta * std::cos(rw_theta), x_bounds);
        rw_dy = reflect(y[t], delta * std::sin(rw_theta), y_bounds);

        // -----------------*
        if (l <= l_star && bias > 0) {
            if (l == 0) {
                // If it's directly on the target, then the directed portion
                // would be to stay in place (this if-else avoids NaNs):
                dr_dx = 0;
                dr_dy = 0;
            } else {
                // If within l_star (but not directly on it) and there's
                // target bias, then include target-directed motion:
                dr_theta = std::atan2((target_xy(i,1) - y(t)),
                    (target_xy(i,0) - x(t)));
                dr_dx = std::min(delta, l) * std::cos(dr_theta);
                dr_dy = std::min(delta, l) * std::sin(dr_theta);
            }
            // combine:
            x(t+1) = x(t) + (1 - bias) * rw_dx + bias * dr_dx;
            y(t+1) = y(t) + (1 - bias) * rw_dy + bias * dr_dy;
        } else {
            // random walk otherwise:
            x(t+1) = x(t) + rw_dx;
            y(t+1) = y(t) + rw_dy;
        }
        // did it hit a new target? (also update distances for next iteration)
        update_i_l_lvec(i, l, l_vec, x(t+1), y(t+1), target_xy, visited,
                        n_ignore_, max_size);
        if (l <= l_i) {
            on_target[t+1] = true;
            new_target[t+1] = true;
            visited(i) = 0;
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
