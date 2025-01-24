
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




xs <- 101; ys <- 101

target_xy <- crossing(x = -(floor(xs / 2)):(floor(xs / 2)),
                      y = -(floor(ys / 2)):(floor(ys / 2))) |>
    as.matrix()
target_types <- as.integer((1:nrow(target_xy)) %% 2 == 0L) + 1L

# target_xy <- cbind(x = c(-3, -3), y = c(0, 0))
# target_types <- 1:2

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
b_sims <- searcher_sims(delta = 2, max_t = 100, x_size = xs, y_size = ys,
                        target_xy = target_xy, target_types = target_types,
                        l_star = target_info$l_star,
                        l_i = target_info$l_i,
                        bias = target_info$bias,
                        n_stay = target_info$n_stay,
                        n_ignore = target_info$n_ignore,
                        n_reps = 100L,
                        n_threads = .n_threads)
# no bacteria sims:
nb_sims <- searcher_sims(delta = 2, max_t = 100, x_size = xs, y_size = ys,
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
    geom_hline(yintercept = c(-1, 1) * xs / 2, color = "gray70") +
    geom_vline(xintercept = c(-1, 1) * ys / 2, color = "gray70") +
    geom_point(data = as_tibble(target_xy) |>
                   mutate(type = factor(target_types,
                                        labels = c("virus", "bacteria"))),
               aes(fill = type), shape = 21, size = 3) +
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
    stat_summary(geom = "point", fun = mean, size = 3, color = "firebrick") +
    facet_wrap(~ name, scales = "free_y")


