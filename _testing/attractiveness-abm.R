
#'
#' Agent-based model of searchers (aphids) and targets (plants)
#' where some targets attract searchers while others repel them.
#' I use this to describe the functional form of interaction rates in response
#' to both attracting and repelling.
#'


source("_testing/_preamble.R")

library(gganimate)



# color palette for target types:
type_pal <- c(viridis(3, begin = 0.1, end = 0.9)[c(1,3,2)], "gray80") |>
    set_names(c("bacteria", "virus", "bacteria + virus", "none"))



xs <- 101L; ys <- 101L

target_sims <- target_type_sims(xs, ys,
                                wt_mat = rbind(c(1, 1), c(1, 1)),
                                n_samples = c(5000, 5000))
#' Convert to form usable in `searcher_sims`, where it's required that
#' there's only one type per location.
#' This means we'll encode a new 4th type referring to both 1 & 2 (the 3rd type
#' refers to neither virus nor bacteria so cannot be involved with multiples).
target_types <- target_sims[["type"]] |>
    map_int(\(x) {
        if (length(x) == 1L) return(x)
        else return(4L)
    })
# Same but with no bacteria
# This means that type 2 (bacteria) is changed to type 3 (none) and
# type 4 (both) to type 1 (virus)
nb_target_types <- case_when(target_types == 4L ~ 1L,
                             target_types == 2L ~ 3L,
                             .default = target_types)

#' Info for target types (in this order):
#' 1. virus-infected
#' 2. Pseudomonas-infected
#' 3. none
#' 4. both
l_star <- list(v = 5,   b = 2,    n = 0.5)
bias <-   list(v = 0.2, b = -0.8, n = 0.1)
target_info <- list(l_star = list(l_star$v, l_star$b, l_star$n, c(l_star$v, l_star$b)),
                    bias =   list(bias$v, bias$b, bias$n, c(bias$v, bias$b)),
                    l_i = rep(0.2, 4L),
                    n_stay = rep(5L, 4L),
                    # n_stay = rep(1e6L, 4L),
                    n_ignore = rep(1e6, 4L))
target_info$n_ignore  <- rep(10L, 4L) # ceiling(2 * map_dbl(target_info$l_star, max) / 2)

.d <- 0.1


# Simulations with bacteria:
b_sims <- searcher_sims(d = .d, max_t = 1000, x_size = xs, y_size = ys,
                        target_xy = target_sims[,c("x","y")] |> as.matrix(),
                        target_types = target_types,
                        l_star = target_info$l_star,
                        l_i = target_info$l_i,
                        bias = target_info$bias,
                        n_stay = target_info$n_stay,
                        n_ignore = target_info$n_ignore,
                        n_searchers = 100L,
                        n_threads = .n_threads)
# no bacteria sims:
nb_sims <- searcher_sims(d = .d, max_t = 1000, x_size = xs, y_size = ys,
                         target_xy = target_sims[,c("x","y")] |> as.matrix(),
                         target_types = nb_target_types,
                         l_star = target_info$l_star,
                         l_i = target_info$l_i,
                         bias = target_info$bias,
                         n_stay = target_info$n_stay,
                         n_ignore = target_info$n_ignore,
                         n_searchers = 100L,
                         n_threads = .n_threads)



bind_rows(b_sims |> mutate(bact = "w"),
          nb_sims |> mutate(bact = "wo")) |>
    filter(type != 0L) |> # these are rows where searchers aren't on targets
    mutate(bact = factor(bact, levels = c("w", "wo"),
                         labels = c("*Pseudomonas*", "no *Pseudomonas*")),
           infect = factor(type == 1L | type == 4L, levels = c(TRUE, FALSE),
                           labels = c("infected", "uninfected"))) |>
    group_by(searcher, bact, infect) |>
    summarize(hit = sum(hit), .groups = "drop") |>
    ggplot(aes(bact, hit, color = bact)) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "gray70") +
    geom_point(shape = 1, alpha = 0.25,
               position = position_jitterdodge(jitter.width = 0.25,
                                               dodge.width = 0.75)) +
    stat_summary(geom = "point", fun = mean, size = 4, shape = 18,
                 position = position_dodge(width = 0.75)) +
    ylab("Number of hits on target") +
    facet_wrap( ~ infect, ncol = 1) +
    scale_color_viridis_d(NULL, guide = "none", option = "plasma", begin = 0.2, end = 0.8) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_markdown(color = "black", size = 10),
          legend.text = element_markdown())




# p <- b_sims |>
#     mutate(searcher = factor(searcher)) |>
#     # filter(t < 100) |>
#     ggplot(aes(x,y)) +
#     geom_hline(yintercept = c(0, xs), color = "gray70") +
#     geom_vline(xintercept = c(0, ys), color = "gray70") +
#     geom_point(data = target_sims |>
#                               mutate(type = type |>
#                                          map(\(x) c("bacteria", "virus", "none")[x]) |>
#                                          map_chr(\(x) paste(x, collapse = " + ")) |>
#                                          factor(levels = c("bacteria", "virus",
#                                                            "bacteria + virus", "none"))),
#                aes(fill = type), shape = 21, size = 1, stroke = 0) +
#     geom_path(aes(color = searcher), linewidth = 1, alpha = 0.5) +
#     geom_point(aes(color = searcher), size = 2, alpha = 0.5) +
#     scale_color_viridis_d(NULL, guide = "none", option = "plasma") +
#     scale_fill_manual(NULL, values = type_pal) +
#     coord_equal() +
#     guides(fill = guide_legend(override.aes = list(size = 3)))
#
#
# # p
#
# anim <- p +
#     labs(title = "Time: {frame_along}") +
#     transition_reveal(time) +
#     ease_aes("linear")
#
# # animate(anim, nframes = length(unique(b_sims$time)))
#
#
# anim_save("~/Desktop/searchers.gif", anim, nframes = length(unique(b_sims$time)),
#           height = 5, width = 5, units = "in", res = 300)
#
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







# ========================================================*
# Test searcher sims ----
# ========================================================*


test_sims <- searcher_sims(d = 2, max_t = 100, x_size = 11L, y_size = 11L,
                           target_xy = rbind(c(2L, 5L), c(8L, 5L)),
                           target_types = c(1, 2),
                           l_star = list(c(5, 2), c(5, 2)),
                           bias = list(c(0.3, -0.7), c(-0.3, 0.7)),
                           l_i = c(0.5, 0.5),
                           n_stay = c(5, 5),
                           n_ignore = c(10, 10),
                           n_searchers = 1L)
p <- test_sims |>
    mutate(searcher = factor(searcher)) |>
    # filter(t < 100) |>
    ggplot(aes(x,y)) +
    geom_hline(yintercept = c(0, 11), color = "gray70") +
    geom_vline(xintercept = c(0, 11), color = "gray70") +
    geom_path(aes(color = searcher), linewidth = 1, alpha = 0.5) +
    geom_point(aes(color = searcher), size = 2, alpha = 0.5) +
    geom_point(data = tibble(type = factor(c(1, 2)),
                             x = c(2, 8),
                             y = c(5, 5)),
               aes(fill = type), shape = 21, size = 5, stroke = 0) +
    scale_color_viridis_d(NULL, guide = "none", option = "plasma", begin = 0.4) +
    scale_fill_viridis_d(NULL, guide = "none", end = 0.5) +
    # guides(fill = guide_legend(override.aes = list(size = 3))) +
    coord_equal()

# p


anim <- p +
    labs(title = "Time: {frame_along}") +
    transition_reveal(time) +
    ease_aes("linear")


animate(anim, nframes = length(unique(test_sims$time)), fps = 5)





# ========================================================*
# Simulating targets ----
# ========================================================*


library(tidyverse)
library(pseudogameofclones)


x_size = 50L
y_size = 50L
n_samples <- c(1000, 1000)

samp_df <- target_type_sims(x_size,
                 y_size,
                 wt_mat,
                 n_samples) |>
    mutate(type = type |>
               map(\(x) c("bacteria", "virus", "none")[x]) |>
               map_chr(\(x) paste(x, collapse = " + ")) |>
               factor(levels = c("bacteria", "virus", "bacteria + virus", "none")))


samp_df |>
    group_by(type) |>
    summarize(n = n())



samp_df |>
    filter(map_lgl(type, \(x) !is.null(x))) |>
    ggplot(aes(x, y, color = type)) +
    geom_point() +
    scale_color_manual(NULL, values = type_pal) +
    coord_equal(xlim = c(0, x_size), ylim = c(0, y_size)) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())



library(spdep)

samp_list <- cbind(samp_df$x, samp_df$y) |>
    knearneigh(k = 8) |>
    knn2nb() |>
    nb2listw(style="B")

joincount.multi(factor(samp_df$type), samp_listw)

