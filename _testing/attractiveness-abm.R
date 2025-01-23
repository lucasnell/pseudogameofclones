
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




xs <- 11; ys <- 11

target_xy <- crossing(x = -5:5, y = -5:5) |>
    as.matrix()
target_types <- as.integer((1:nrow(target_xy)) %% 2 == 0L) + 1L

# target_xy <- cbind(x = c(-3, -3), y = c(0, 0))
# target_types <- 1:2

if (nrow(target_xy) == 2) {
    l_star <- c(5, 5)
    l_i <- c(1, 1)
} else {
    l_star <- c(2, 0.5)
    l_i <- c(0.1, 0.1)
}

x <- searcher_sims(delta = 0.5, max_t = 100, x_size = xs, y_size = ys,
                   target_xy = target_xy, target_types = target_types,
                   l_star = l_star,
                   l_i = l_i,
                   bias = c(0.5, -0.1),
                   n_stay = c(5L, 5L),
                   # n_ignore = ceiling(2 * l_star / 0.5),
                   n_ignore = c(1e6, 1e6),
                   xy0 = NULL,
                   randomize_xy0 = TRUE,
                   n_reps = 100L,
                   show_progress = FALSE,
                   n_threads = .n_threads)

p <- x |>
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

animate(anim, nframes = length(unique(x$time)))


anim_save("~/Desktop/searchers.gif", anim, nframes = length(unique(x$time)),
          height = 5, width = 5, units = "in", res = 300)

x |>
    filter(type > 0) |>
    mutate(type = factor(type, labels = c("virus", "bacteria"))) |>
    group_by(type, rep) |>
    summarize(on = n(), hit = sum(hit), .groups = "drop") |>
    pivot_longer(on:hit) |>
    ggplot(aes(type, value)) +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "gray70") +
    geom_jitter(aes(color = type), width = 0.25, height = 0, shape = 1) +
    facet_wrap(~ name, scales = "free_y") +
    scale_color_manual(NULL, values = c("dodgerblue", "goldenrod"),
                       guide = "none")


