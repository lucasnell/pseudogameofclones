
#'
#' Agent-based model of searchers (aphids) and targets (plants)
#' where some targets attract searchers while others repel them.
#' I use this to describe the functional form of interaction rates in response
#' to both attracting and repelling.
#'

library(tidyverse)
library(gganimate)

# ggplot2 theme:
theme_set(theme_classic() %+replace%
              theme(strip.background = element_blank(),
                    strip.text = element_text(size = 11),
                    legend.background = element_blank(),
                    plot.title = element_text(size = 14, hjust = 0)))

if (interactive()) {
    setHook(packageEvent("grDevices", "onLoad"),
            function(...) grDevices::quartz.options(width = 4, height = 4,
                                                    pointsize = 10))
    options("device" = "quartz")
    grDevices::graphics.off()
}



# unbounded random walk:
random_walk <- function(delta, maxt) {
    x <- numeric(maxt + 1L)
    y <- numeric(maxt + 1L)
    for (t in 1:maxt) {
        theta <- runif(1, 0, 2 * pi)
        x[t+1] <- x[t] + delta * cos(theta)
        y[t+1] <- y[t] + delta * sin(theta)
    }
    return(tibble(t = 0:maxt, x = x, y = y))
}



# random_walk(1, 1000) |>
#     ggplot(aes(x,y)) +
#     geom_hline(yintercept = 0, color = "gray70") +
#     geom_vline(xintercept = 0, color = "gray70") +
#     geom_path(linewidth = 1) +
#     geom_point(size = 2) +
#     coord_equal()




#
# New dx or dy that includes effects of reflecting off bound(s).
#
reflect <- function(start, d, bounds) {
    stopifnot(length(bounds) == 2)
    stopifnot(bounds[2] > bounds[1])
    past_bounds <- c((start + d) < bounds[1], (start + d) > bounds[2])
    while (any(past_bounds)) {
        bound <- bounds[past_bounds]
        d <- 2 * (bound - start) - d
        past_bounds <- c((start + d) < bounds[1], (start + d) > bounds[2])
    }
    return(d)
}

#
# Calculate Euclidean distance between two points
#
distance <- function(x1, y1, x2, y2) sqrt((x1 - x2)^2 + (y1 - y2)^2)


# random walk with boundaries:
bound_rw <- function(delta, maxt, x_size, y_size) {

    x_bounds <- c(-1, 1) * x_size / 2
    y_bounds <- c(-1, 1) * y_size / 2

    x <- numeric(maxt + 1L)
    y <- numeric(maxt + 1L)

    for (t in 1:maxt) {
        theta <- runif(1, 0,  2 * pi)
        dx <- reflect(x[t], delta * cos(theta), x_bounds)
        dy <- reflect(y[t], delta * sin(theta), y_bounds)
        x[t+1] <- x[t] + dx
        y[t+1] <- y[t] + dy
    }

    tibble(t = 0:maxt, x = x, y = y)
}



xs <- 10; ys <- 10

bound_rw(1, 100, xs, ys) |>
    ggplot(aes(x,y)) +
    geom_hline(yintercept = c(-1, 1) * xs / 2, color = "gray70") +
    geom_vline(xintercept = c(-1, 1) * ys / 2, color = "gray70") +
    geom_path(linewidth = 1) +
    geom_point(size = 2) +
    coord_equal()





#' Random walk with boundaries and bias towards targets.
#'
#' Searchers can only interact (when distance is <= l_i) with 1 target, but
#' they can be influenced (when distance is <= l_star) with multiple targets.
#' When multiple targets influence trajectories, the effect each has on the
#' searcher is proportional to the distance from the searcher divided by
#' the l_star for that target.
#'
bias_bound_rw <- function(delta, maxt, x_size, y_size, target_xy,
                          l_star, l_i, bias,
                          n_stay = 5L,
                          n_ignore = NULL,
                          xy0 = NULL) {

    x_bounds <- c(-1, 1) * x_size / 2
    y_bounds <- c(-1, 1) * y_size / 2

    stopifnot(bias >= 0 && bias <= 1)
    stopifnot(ncol(target_xy) == 2)
    stopifnot(nrow(target_xy) >= 1)
    stopifnot(all(target_xy[,1] > x_bounds[1]))
    stopifnot(all(target_xy[,1] < x_bounds[2]))
    stopifnot(all(target_xy[,2] > y_bounds[1]))
    stopifnot(all(target_xy[,2] < y_bounds[2]))

    if (is.null(n_ignore)) n_ignore <- ceiling(2 * l_star / delta)

    # Check for overlapping bounds:
    target_dists <- matrix(0, nrow(target_xy), nrow(target_xy))
    for (i in 2:nrow(target_xy)) {
        for (j in 1:(i-1L)) {
            target_dists[i,j] <- distance(target_xy[i,1], target_xy[i,2],
                                          target_xy[j,1], target_xy[j,2])
        }
    }
    if (min(target_dists[lower.tri(target_dists)]) <= l_i) {
        stop("targets are too close together which would result in searchers ",
             "interacting with >1 at a time.")
    }

    # These are outputs:
    x <- numeric(maxt + 1L)
    y <- numeric(maxt + 1L)
    if (!is.null(xy0)) {
        stopifnot(is.numeric(xy0) && length(xy0) == 2)
        stopifnot(all(xy0[1] > x_bounds[1] & xy0[1] < x_bounds[2]))
        stopifnot(all(xy0[2] > y_bounds[1] & xy0[2] < y_bounds[2]))
        x[1] <- xy0[1]
        y[1] <- xy0[2]
    }
    on_target <- logical(maxt + 1L) # whether on target at given time
    new_target <- logical(maxt + 1L) # whether hit a *new* target at given time
    # Vector of distances from searcher to target(s)
    l_vec <- apply(target_xy, 1, \(xy) distance(x[1], y[1], xy[1], xy[2]))
    # which target(s) are within l_star (and thus influence trajectory):
    wi_lstar <- which(l_vec <= l_star)
    # which target is within l_i (there can be a max of 1 at a time):
    wi_li <- which(l_vec <= l_i)
    # minimum distance to a target:
    l_min <- ifelse(length(wi_lstar) == 0, l_star * 2, min(l_vec[wi_lstar]))
    # number of time points within l_i of most recent target:
    stayed <- 0L
    # `visited` below indicates number of time points since the most
    # recent visit for each target (only starts after leaving within
    # l_i of the target):
    visited <- rep(as.integer(n_ignore * 2), length(l_vec))
    # If starting out on target, adjust these:
    if (l_min <= l_i) {
        on_target[1] <- TRUE
        new_target[1] <- TRUE
        visited[wi_li] <- TRUE
    }

    for (t in 1:maxt) {
        # If within l_i and have NOT stayed long enough, deterministically
        # move towards target and go to next iteration
        if (l_min <= l_i && stayed < n_stay) {
            if (l_min > 0) {
                dr_theta <- atan2((target_xy[wi_li,2] - y[t]),
                                  (target_xy[wi_li,1] - x[t]))
                dr_dx <- unname(min(delta, l_min) * cos(dr_theta))
                dr_dy <- unname(min(delta, l_min) * sin(dr_theta))
                x[t+1] <- x[t] + dr_dx
                y[t+1] <- y[t] + dr_dy
                l_min <- distance(x[t+1], y[t+1], target_xy[wi_li,1],
                                  target_xy[wi_li,2])
                l_vec[wi_li] <- l_min
            } else {
                x[t+1] <- x[t]
                y[t+1] <- y[t]
            }
            on_target[t+1] <- TRUE
            stayed <- stayed + 1L
            next
        }
        # If within l_i and have stayed long enough, adjust some things
        # for potentially ignoring target(s), then proceed as normal:
        if (l_min <= l_i && stayed >= n_stay) {
            stayed <- 0L
            visited[wi_li] <- 0L
            wi_lstar <- which(l_vec <= l_star & visited >= n_ignore)
            wi_li <- which(l_vec <= l_i & visited >= n_ignore)
            l_min <- ifelse(length(wi_lstar) == 0, l_star * 2,
                            min(l_vec[wi_lstar]))
        }

        # -----------------*
        # random walk portion of step:
        rw_theta <- runif(1, 0, 2 * pi)
        rw_dx <- reflect(x[t], delta * cos(rw_theta), x_bounds)
        rw_dy <- reflect(y[t], delta * sin(rw_theta), y_bounds)

        # -----------------*
        if (l_min <= l_star && bias > 0) {

            if (l_min == 0) {

                # If it's directly on the target, then the directed portion
                # would be to stay in place (this if-else avoids NaNs):
                dr_dx <- 0
                dr_dy <- 0

            } else {

                # If within l_star (but not directly on it) and there's
                # target bias, then include target-directed motion:
                dr_dxy <- sapply(wi_lstar, \(k) {
                    l <- l_vec[[k]]
                    dr_theta <- atan2((target_xy[k,2] - y[t]),
                                         (target_xy[k,1] - x[t]))
                    dr_dx <- unname(min(delta, l) * cos(dr_theta))
                    dr_dy <- unname(min(delta, l) * sin(dr_theta))
                    return(c(x = dr_dx, y = dr_dy))
                })
                wts <- (l_star / l_vec[wi_lstar]) / sum(l_star / l_vec[wi_lstar])
                dr_dx <- sum(dr_dxy[1,] * wts)
                dr_dy <- sum(dr_dxy[2,] * wts)
            }
            # combine:
            x[t+1] <- x[t] + (1 - bias) * rw_dx + bias * dr_dx
            y[t+1] <- y[t] + (1 - bias) * rw_dy + bias * dr_dy

        } else {

            # random walk otherwise:
            x[t+1] <- x[t] + rw_dx
            y[t+1] <- y[t] + rw_dy

        }

        # did it hit a new target? (also update distances for next iteration)
        l_vec <- apply(target_xy, 1, \(xy) distance(x[t+1], y[t+1], xy[1], xy[2]))
        wi_lstar <- which(l_vec <= l_star & visited >= n_ignore)
        wi_li <- which(l_vec <= l_i & visited >= n_ignore)
        l_min <- ifelse(length(wi_lstar) == 0, l_star * 2,
                        min(l_vec[wi_lstar]))
        if (l_min <= l_i) {
            on_target[t+1] <- TRUE
            new_target[t+1] <- TRUE
            visited[i] <- 0L
            stayed <- 0L
        }

        visited <- visited + 1L

    }

    tibble(t = 0:maxt, x = x, y = y, on = on_target, hit = new_target)
}




Rcpp::sourceCpp("_testing/abm.cpp")



xs <- 11; ys <- 11
target_xy <- crossing(x = -5:5, y = -5:5) |>
    as.matrix()

crossing(x = -5:5, y = -5:5) |>
    mutate(t = factor((1:n()) %% 2 == 0L)) |>
    ggplot(aes(x,y)) +
    geom_hline(yintercept = c(-1, 1) * xs / 2, color = "gray70") +
    geom_vline(xintercept = c(-1, 1) * ys / 2, color = "gray70") +
    geom_point(aes(color = t), size = 3) +
    coord_equal()



target_xy2 <- rbind(c(-3, 0), c(3, 0))

.seed <- sample.int(2^31-1, 1); set.seed(.seed); print(.seed)
x <- bias_bound_rw(delta = 0.5, maxt = 100, x_size = xs, y_size = ys,
                   target_xy = target_xy2,
                   l_star = 5, l_i = 1, bias = 1,
                   n_stay = 5L,
                   n_ignore = 1e6,
                   xy0 = c(-1.9, 0))
                   # xy0 = runif(2, c(-xs/2, -ys/2), c(xs/2, ys/2)))

p <- x |>
    filter(t < 10) |>
    ggplot(aes(x,y)) +
    geom_hline(yintercept = c(-1, 1) * xs / 2, color = "gray70") +
    geom_vline(xintercept = c(-1, 1) * ys / 2, color = "gray70") +
    geom_point(data = as_tibble(target_xy2) |> set_names(c("x", "y")),
               color = "red", size = 3) +
    geom_path(linewidth = 1) +
    geom_point(size = 2) +
    coord_equal()

anim <- p +
    labs(title = "Time: {frame_along}") +
    transition_reveal(t) +
    ease_aes("linear")

anim






# anim_save("~/Desktop/anim.gif")





Rcpp::cppFunction(depends = "RcppArmadillo",
"arma::uvec arma_find1(const arma::vec& x, const std::vector<double>& y) {
    arma::uvec z;
    for (size_t i = 0; i < y.size(); i++) z = arma::find(x > y[i]);
    return z;
}
")
Rcpp::cppFunction(depends = "RcppArmadillo",
"arma::uvec arma_find2(const arma::vec& x, const std::vector<double>& y) {
    arma::uvec z;
    size_t n;
    for (size_t i = 0; i < y.size(); i++) {
        n = 0;
        for (const double& xx : x) if (xx > y[i]) n++;
        z.set_size(n);
        n = 0;
        for (size_t j = 0; j < x.n_elem; j++) {
            if (x(j) > y[i]) {
                z(n) = j;
                n++;
            }
        }
    }
    return z;
}
")
Rcpp::cppFunction(depends = "RcppArmadillo",
"arma::uvec arma_find3(const arma::vec& x, const std::vector<double>& y) {
    arma::uvec z;
    std::vector<arma::uword> zv;
    zv.reserve(x.n_elem);
    for (size_t i = 0; i < y.size(); i++) {
        zv.clear();
        for (arma::uword j = 0; j < x.n_elem; j++) {
            if (x(j) > y[i]) zv.push_back(j);
        }
        z.set_size(zv.size());
        for (size_t j = 0; j < zv.size(); j++) z(j)  = zv[j];
    }
    return z;
}
")


x <- runif(100)
y <- runif(1000)
real <- which(x > y[1]) - 1L

identical(as.integer(arma_find1(x, y[1])), real) &&
    identical(as.integer(arma_find2(x, y[1])), real) &&
    identical(as.integer(arma_find3(x, y[1])), real)

identical(arma_find1(x, y), arma_find2(x, y)) &&
    identical(arma_find1(x, y), arma_find3(x, y))

microbenchmark::microbenchmark(a = arma_find1(x, y),
                               b = arma_find2(x, y),
                               c = arma_find3(x, y))


