
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
target_info <- list(l_star = c(5, 1) |> as.list(),
                    l_i = c(0.2, 0.2),
                    bias = c(0.3, -0.6) |> as.list(),
                    n_stay = c(5L, 5L),
                    n_ignore = c(1e6, 1e6))
target_info$n_ignore  <- ceiling(2 * unlist(target_info$l_star) / 0.5)

# Change these if only two targets (used for testing)
if (nrow(target_xy) == 2) {
    target_info$l_star <- c(5, 5) |> as.list()
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
                        n_searchers = 100L,
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
                         n_searchers = 100L,
                         n_threads = .n_threads)


p <- b_sims |>
    mutate(searcher = factor(searcher)) |>
    # filter(t < 100) |>
    ggplot(aes(x,y)) +
    geom_hline(yintercept = c(0, xs), color = "gray70") +
    geom_vline(xintercept = c(0, ys), color = "gray70") +
    geom_point(data = as_tibble(target_xy) |>
                   mutate(type = factor(target_types,
                                        labels = c("virus", "bacteria"))),
               aes(fill = type), shape = 21, size = 1, stroke = 0) +
    geom_path(aes(color = searcher), linewidth = 1, alpha = 0.5) +
    geom_point(aes(color = searcher), size = 2, alpha = 0.5) +
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
#     group_by(type, searcher) |>
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
    group_by(searcher, bact) |>
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
# Simulating targets ----
# ========================================================*



x_size = 50L
y_size = 50L
corr = rbind(c(1, 1.2),
             c(1.2, 2))
n_samples <- c(1000, 1000)

samp_df <- target_type_sims(x_size,
                 y_size,
                 corr,
                 n_samples) |>
    mutate(type = map_chr(type, \(x) paste(x, collapse = "_")))


samp_df |>
    group_by(type) |>
    summarize(n = n())



samp_df |>
    filter(map_lgl(type, \(x) !is.null(x))) |>
    mutate(type = map_chr(type, \(x) paste(sort(x), collapse = "_"))) |>
    ggplot(aes(x, y, color = type)) +
    geom_point() +
    scale_color_viridis_d(begin = 0.1, end = 0.9) +
    coord_equal(xlim = c(0, x_size), ylim = c(0, y_size))



library(spdep)

samp_listw <- cbind(samp_df$x, samp_df$y) |>
    knearneigh(k = 8) |>
    knn2nb() |>
    nb2listw(style="B")

joincount.multi(factor(samp_df$type), samp_listw)








