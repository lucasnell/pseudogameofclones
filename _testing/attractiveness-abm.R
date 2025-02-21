
#'
#' Agent-based model of searchers (aphids) and targets (plants)
#' where some targets attract searchers while others repel them.
#' I use this to describe the functional form of interaction rates in response
#' to both attracting and repelling.
#'

library(tidyverse)
library(gganimate)
library(pseudogameofclones)

# number of threads:
.n_threads <- max(1L, parallel::detectCores() - 2L)


if (interactive()) {
    setHook(packageEvent("grDevices", "onLoad"),
            function(...) grDevices::quartz.options(width = 4, height = 4,
                                                    pointsize = 10))
    options("device" = "quartz")
    grDevices::graphics.off()
}




xs <- 101L; ys <- 101L

target_xy <- crossing(x = 1:(xs-1L), y = 1:(ys-1L)) |>
    as.matrix()
# target_types <- as.integer((1:nrow(target_xy)) %% 2 == 0L) + 1L
target_types <- lapply(1:(xs-1L), \(i) if (i %% 2L == 0) rep(1:2, (ys-1)/2)
                       else rep(2:1, (ys-1)/2)) |>
    do.call(what = c)

# target_xy <- cbind(x = c(-3, -3), y = c(0, 0))
# target_types <- 1:2

target_p <- as_tibble(target_xy) |>
    mutate(type = factor(target_types,
                         labels = c("virus", "bacteria"))) |>
    ggplot(aes(x,y)) +
    geom_hline(yintercept = c(0, xs), color = "gray70") +
    geom_vline(xintercept = c(0, ys), color = "gray70") +
    geom_point(aes(color = type), shape = 19, size = 1) +
    scale_color_manual(NULL, values = c("dodgerblue", "firebrick")) +
    coord_equal()
# target_p

# Info for target types: virus- and Pseudomonas-infected, respectively:
target_info <- list(l_star = c(5, 1),
                    l_i = c(0.2, 0.2),
                    bias = c(0.3, -0.6),
                    n_stay = c(5L, 5L),
                    n_ignore = c(1e6, 1e6))
target_info$n_ignore  <- ceiling(2 * target_info$l_star / 0.5)

# Change these if only two targets (used for testing)
if (nrow(target_xy) == 2) {
    target_info$l_star <- c(5, 5)
    target_info$l_i <- c(1, 1)
}


# Simulations with bacteria:
b_sims <- searcher_sims(d = 2, max_t = 100, x_size = xs, y_size = ys,
                        target_xy = target_xy, target_types = target_types,
                        l_star = target_info$l_star,
                        l_i = target_info$l_i,
                        bias = target_info$bias,
                        n_stay = target_info$n_stay,
                        n_ignore = target_info$n_ignore,
                        n_reps = 100L,
                        n_threads = .n_threads)
# no bacteria sims:
nb_sims <- searcher_sims(d = 2, max_t = 100, x_size = xs, y_size = ys,
                         target_xy = target_xy[target_types == 1,],
                         target_types = target_types[target_types == 1],
                         l_star = target_info$l_star[1],
                         l_i = target_info$l_i[1],
                         bias = target_info$bias[1],
                         n_stay = target_info$n_stay[1],
                         n_ignore = target_info$n_ignore[1],
                         n_reps = 100L,
                         n_threads = .n_threads)


p <- b_sims |>
    mutate(rep = factor(rep)) |>
    # filter(t < 100) |>
    ggplot(aes(x,y)) +
    geom_hline(yintercept = c(0, xs), color = "gray70") +
    geom_vline(xintercept = c(0, ys), color = "gray70") +
    geom_point(data = as_tibble(target_xy) |>
                   mutate(type = factor(target_types,
                                        labels = c("virus", "bacteria"))),
               aes(fill = type), shape = 21, size = 1, stroke = 0) +
    geom_path(aes(color = rep), linewidth = 1, alpha = 0.5) +
    geom_point(aes(color = rep), size = 2, alpha = 0.5) +
    scale_color_viridis_d(NULL, guide = "none") +
    scale_fill_manual(NULL, values = c("dodgerblue", "goldenrod"),
                      guide = "none") +
    coord_equal()

anim <- p +
    labs(title = "Time: {frame_along}") +
    transition_reveal(time) +
    ease_aes("linear")

# animate(anim, nframes = length(unique(b_sims$time)))


# anim_save("~/Desktop/searchers.gif", anim, nframes = length(unique(b_sims$time)),
#           height = 5, width = 5, units = "in", res = 300)

# b_sims |>
#     filter(type > 0) |>
#     mutate(type = factor(type, labels = c("virus", "bacteria"))) |>
#     group_by(type, rep) |>
#     summarize(on = n(), hit = sum(hit), .groups = "drop") |>
#     pivot_longer(on:hit) |>
#     ggplot(aes(type, value)) +
#     geom_hline(yintercept = 0, linewidth = 0.5, color = "gray70") +
#     geom_jitter(aes(color = type), width = 0.25, height = 0, shape = 1) +
#     facet_wrap(~ name, scales = "free_y") +
#     scale_color_manual(NULL, values = c("dodgerblue", "goldenrod"),
#                        guide = "none")



bind_rows(b_sims |> mutate(bact = "w"),
          nb_sims |> mutate(bact = "wo")) |>
    mutate(bact = factor(bact, levels = c("w", "wo"))) |>
    filter(type == 1) |>
    group_by(rep, bact) |>
    summarize(on = n(), hit = sum(hit), .groups = "drop") |>
    pivot_longer(on:hit) |>
    ggplot(aes(bact, value)) +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "gray70") +
    geom_jitter(width = 0.25, height = 0, shape = 1) +
    ylab("Number of events") +
    stat_summary(geom = "point", fun = mean, size = 3, color = "firebrick") +
    facet_wrap(~ name, scales = "free_y")



library(tidyverse)
library(pseudogameofclones)




# ========================================================*
# Mean difference between two random variables generated by sampling
# integer locations without replacement
# ========================================================*

# ---------------------*
# in one dimension:
sapply(1:100e3, \(i) {
    x <- sample.int(50, 2, replace = FALSE)
    return(abs(diff(x)))
}) |> mean()
avg_diff_1d <- function(n) {
    d <- 0:(n-1)
    return((2 * sum(d * (n - d))) / (n * (n - 1)))
}
avg_diff_1d(50)


# ---------------------*
# in two dimensions:
sapply(1:100e3, \(i) {
    y_size <- x_size <- 50L
    i <- sample.int(x_size * y_size, 2, replace = FALSE)
    x <- i - y_size * ((i - 1L) %/% y_size)
    y <- (i - 1L) %/% y_size + 1L
    return(sqrt(diff(x)^2 + diff(y)^2))
}) |> mean()
avg_diff_2d <- function(nx, ny) {
    dx = 1:(nx-1)
    dy = 1:(ny-1)
    d_sum <- (4 * sum(sqrt(rep(dx, ny-1)^2 + rep(dy, each = nx-1)^2) *
                          (nx - rep(dx, ny-1)) * (ny - rep(dy, each = nx-1)))) +
        (2 * nx * sum(dy * (ny - dy))) + (2 * ny * sum(dx * (nx - dx)))
    return(d_sum / (nx * ny * (nx * ny - 1)))
}
avg_diff_2d(nx = 50L, ny = 50L)




target_type_sims_R <- function(x_size,
                               y_size,
                               corr,
                               n_samples) {

    stopifnot(is.numeric(corr) && length(corr) > 0)
    stopifnot(is.matrix(corr) && isSymmetric(corr) && all(corr >= 0))
    stopifnot(is.numeric(n_samples) && length(n_samples) == nrow(corr) &&
                  all(n_samples <= x_size * y_size) && all(n_samples > 0))

    n_types <- nrow(corr)

    nz <- x_size * y_size

    # Vector of locations:
    Z <- mapply(\(x, y) c(x, y), rep(1:x_size, y_size),
                rep(1:y_size, each = x_size), SIMPLIFY = FALSE)

    # Distance matrix:
    M <- matrix(0, nz, nz)
    for (i in 1:(nz-1)) for (j in (i+1):nz) {
        M[i, j] <- M[j, i] <- sqrt(sum((Z[[i]] - Z[[j]])^2))
    }

    # Vector of sampling probabilities for each type:
    P <- lapply(1:n_types, \(i) rep(1.0, nz))

    # Output data frame of xy coordinates by type:
    samp_df <- tibble(x = do.call(c, lapply(Z,\(z) z[1])),
                      y = do.call(c, lapply(Z,\(z) z[2])),
                      type = rep(list(NULL), nz))

    # Now sample for points by type:
    n_samps <- rep(0L, n_types)
    while (sum(n_samps) < sum(n_samples)) {
        for (j in 1:n_types) {
            if (n_samps[j] >= n_samples[j]) next
            k <- sample.int(nz, 1, prob = P[[j]])
            # assign type to output:
            if (length(samp_df$type[[k]]) >= n_types)
                stop("length(samp_df$type[[k]]) >= n_types")
            samp_df$type[[k]] <- c(samp_df$type[[k]], j)
            # adjust sampling probabilities:
            nearby <- M[k,] <= 1.5
            P[[j]][nearby & P[[j]] > 0] <- P[[j]][nearby & P[[j]] > 0] * corr[j,j]
            for (i in (1:n_types)[-j]) {
                P[[i]][nearby & P[[i]] > 0] <- P[[i]][nearby & P[[i]] > 0] * corr[i,j]
            }
            P[[j]][[k]] <- 0
            n_samps[j] <- n_samps[j] + 1L
        }
    }

    samp_df <- samp_df |>
        filter(map_lgl(type, \(x) !is.null(x))) |>
        mutate(type = map_chr(type, \(x) paste(sort(x), collapse = "_")))

    return(samp_df)

}



x_size = 50L
y_size = 50L
corr = rbind(c(1, 10),
             c(10, 10))
# corr <- c(0,0,0)
# corr_d = 2 * c(1, 1, 1)
# n_samples <- NULL
n_samples <- c(1000, 1000)

samp_df <- target_type_sims(x_size,
                 y_size,
                 corr,
                 n_samples) |>
    group_by(x, y) |>
    summarize(type = paste(type, collapse = "_"), .groups = "drop")


samp_df2 <- target_type_sims_R(x_size,
                               y_size,
                               corr,
                               n_samples)



samp_df |>
    group_by(type) |>
    summarize(n = n())
samp_df2 |>
    group_by(type) |>
    summarize(n = n())


samp_df |>
    filter(map_lgl(type, \(x) !is.null(x))) |>
    mutate(type = map_chr(type, \(x) paste(sort(x), collapse = "_"))) |>
    ggplot(aes(x, y, color = type)) +
    geom_point() +
    scale_color_viridis_d(begin = 0.1, end = 0.9) +
    coord_equal(xlim = c(0, x_size), ylim = c(0, y_size))


samp_diff_df <- samp_df |>
    split(~ type) |>
    map_dfr(\(x, n = 1000) {
        # Choose n random pairs:
        cm <- t(combn(nrow(x), 2))
        if (nrow(cm) > n) cm <- cm[sample.int(nrow(cm), n),]
        diffs <- sqrt((x$x[cm[,1]] - x$x[cm[,2]])^2 + (x$y[cm[,1]] - x$y[cm[,2]])^2)
        return(tibble(type = x$type[1], diff = diffs))
    })

samp_diff_df |>
    ggplot(aes(type, diff)) +
    geom_jitter(alpha = 0.25, height = 0, width = 0.2, shape = 1) +
    geom_hline(yintercept = avg_diff_2d(nx = x_size, ny = y_size),
               linetype = "22", color = "red", linewidth = 1) +
    stat_summary(fun = mean, geom = "point", size = 3, color = "red")






samp_listw <- cbind(samp_df$x, samp_df$y) |>
    knearneigh(k = 8) |>
    knn2nb() |>
    nb2listw(style="B")

joincount.multi(factor(samp_df$type), samp_listw)








