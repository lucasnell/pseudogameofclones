library(tidyverse)
library(pseudogameofclones)
library(viridisLite)
library(gganimate)


if (interactive()) {
    setHook(packageEvent("grDevices", "onLoad"),
            function(...) grDevices::quartz.options(width = 4, height = 4,
                                                    pointsize = 10))
    options("device" = "quartz")
    grDevices::graphics.off()
}



# color palette for target types:
type_pal <- c(viridis(3, begin = 0.1, end = 0.9), "gray80")[c(2,1,4,3)] |>
    set_names(c("virus", "pseudo", "none", "both"))

xs <- ys <- 134L  # approximates a hectare with 0.75 m plant spacing

#' Below converts the `type` column to a form usable in `alate_search_sims`,
#' where it's required that there's only one type per location.
#' Additionally, there should be a max of 4 types: `"virus"`, `"pseudo"`,
#' `"none"`, or `"both"` (in that order).
target_sims <- target_type_sims(xs+1L, ys+1L, #<< bc it sims locations 1:(*_size-1)
                                wt_mat = rbind(c(1, 1), c(1, 1)),
                                n_samples = rep(round(xs*ys * 0.5), 2)) |>
    getElement(1) |>
    mutate(type = factor(type, levels = c("1", "2", "", "1_2")) |>
               as.integer(),
           #' Same but with no bacteria
           #' This means that `"pseudo"` is changed to `"none"` and
           #' `"both"` to `"virus"`
           nb_type = case_when(type == 2L ~ 3L,
                               type == 4L ~ 1L,
                               .default = type))


{
    t0 <- Sys.time()
    sim_df <- alate_search_sims(max_t = 100,
                                plant_xy = target_sims[,c("x","y")] |> as.matrix(),
                                plant_types = target_sims$type,
                                alpha = 0.5,
                                beta = -1,
                                epsilon = 0.5,
                                n_alates = 100e3L,
                                summarize = "types",
                                n_threads = 6)
    t1 <- Sys.time()
    print(t1 - t0); rm(t0, t1)
}


#







# =============================================================================*
# =============================================================================*
# old code ----
# =============================================================================*
# =============================================================================*

#' From "The Role of Aphid Behaviour in the Epidemiology of Potato Virus Y:
#' a Simulation Study" by Thomas Nemecek (1993; p. 72), dispersal distances
#' follow a Weibull distribution with shape = 0.6569 and scale = 9.613.
#'
#' We'll use the median of this distribution for our dispersal radius,
#' where alates have to move to a patch somewere within it.
#'
#' I'm dividing by 0.75 to convert from meters to plant locations that are
#' 0.75 meters apart (typical spacing for pea).
#'
radius <- qweibull(0.5, 0.6569, 9.613) / 0.75

#' For a given distance or below above the current plant, how many plants
#' to the left or right are within the radius?
sapply(1:floor(radius), \(d) sqrt(radius^2 - d^2))

# color palette for target types:
type_pal <- c(viridis(3, begin = 0.1, end = 0.9), "gray80")[c(2,1,4,3)] |>
    set_names(c("virus", "pseudo", "none", "both"))

xs <- ys <- 134L  # approximates a hectare with 0.75 m plant spacing

#' Below converts the `type` column to a form usable in `alate_search_sims`,
#' where it's required that there's only one type per location.
#' Additionally, there should be a max of 4 types: `"virus"`, `"pseudo"`,
#' `"none"`, or `"both"` (in that order).
target_sims <- target_type_sims(xs+1L, ys+1L, #<< bc it sims locations 1:(*_size-1)
                                wt_mat = rbind(c(1, 1), c(1, 1)),
                                n_samples = rep(round(xs*ys * 0.5), 2)) |>
    mutate(type = type |>
               map_int(\(x) {
                   if (length(x) == 1L) return(x)
                   else return(4L)
                   }),
           #' Same but with no bacteria
           #' This means that `"pseudo"` is changed to `"none"` and
           #' `"both"` to `"virus"`
           nb_type = case_when(type == 2L ~ 3L,
                               type == 4L ~ 1L,
                               .default = type))


#' Making matrices of out types:
tt_mat <- nb_tt_mat <- matrix(NA_integer_, xs, ys)
for (i in 1:nrow(target_sims)) {
    x <- target_sims$x[[i]]
    y <- target_sims$y[[i]]
    tt_mat[x,y] <- target_sims$type[[i]]
    nb_tt_mat[x,y] <- target_sims$nb_type[[i]]
}; rm(i, x, y)



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ----




max_t = 100
target_types = tt_mat
radius = 7.336451
alpha = 0.5
beta = -1
w = 0.2
epsilon = 0.5

# xy0 = NULL
# randomize_xy0 = TRUE
# n_searchers = 1

stopifnot(ncol(target_types) > 0 && nrow(target_types) > 0)

n_x <- ncol(target_types)
n_y <- nrow(target_types)

stopifnot(all(target_types %in% 1:4))

# Binary variable for whether a plant has virus:
v <- apply(target_types, c(1,2),
           \(tt) as.integer(tt == 1L | tt == 4L))
# Binary variable for whether a plant has *Pseudomonas*:
b <- apply(target_types, c(1,2),
           \(tt) as.integer(tt == 2L | tt == 4L))

# change in x and y for neighbors (those within radius)
neigh_dxdy <- lapply(-floor(radius):floor(radius), \(dx) {
    ddy <- floor(sqrt(radius^2 - dx^2))
    dy <- -ddy:ddy
    return(lapply(dy, \(.dy) c(dx, .dy)))
}) |>
    do.call(what = c) |>
    do.call(what = rbind)






# Location of aphid:
x <- rep(NA_real_, max_t+1)
y <- rep(NA_real_, max_t+1)
# Types of plants the aphid traveled to and from:
to <- rep(NA_integer_, max_t+1)
from <- rep(NA_integer_, max_t+1)

x[1] <- sample.int(n_x, 1)
y[1] <- sample.int(n_y, 1)
to[1] <- target_types[x[1], y[1]]
from[1] <- target_types[x[1], y[1]]  # to avoid NAs

for (t in 1:max_t) {

    # feed_p <- w
    # if (v[x[t],y[t]]) feed_p <- feed_p * epsilon
    #
    # if (runif(1) < feed_p) {
    #     # x[(t+1):length(x)] <- x[t]
    #     # y[(t+1):length(y)] <- y[t]
    #     # to[(t+1):length(x)] <- to[t]
    #     # from[(t+1):length(y)] <- from[t]
    #     break
    # }

    # Sampling weights:
    s_wts <- apply(neigh_dxdy, 1,
                   \(dxdy) {
                       .x <- x[t] + dxdy[1]
                       .y <- y[t] + dxdy[2]
                       if (.x <= 0 | .x > n_x) return(0)
                       if (.y <= 0 | .y > n_y) return(0)
                       return(exp(alpha * v[.x,.y] + beta * b[.x,.y]))
                   })

    neigh_i <- sample.int(length(s_wts), 1, prob = s_wts / sum(s_wts))

    x[t+1] <- x[t] + neigh_dxdy[neigh_i,1]
    y[t+1] <- y[t] + neigh_dxdy[neigh_i,2]
    to[t+1] <- target_types[x[t+1], y[t+1]]
    from[t+1] <- to[t]

}



tibble(x = x, y = y) |>
    mutate(type = map2_chr(x, y, \(i, j) target_types[i,j]),
           type = factor(type, levels = 1:4,
                         labels = c("virus", "pseudo", "none", "both"))) |>
    group_by(type) |>
    summarize(hits = n())

tibble(from = from, to = to) |>
    mutate(across(everything(), \(x) {
        factor(x, levels = 1:4,
               labels = c("virus", "pseudo", "none", "both"))
    })) |>
    group_by(from, to, .drop = FALSE) |>
    summarize(n = n(), .groups = "drop")


tibble(from = from, to = to) |>
    mutate(across(everything(), \(x) {
        factor(x, levels = 1:4,
               labels = c("virus", "pseudo", "none", "both"))
    })) |>
    ggplot(aes(from, to)) +
    geom_bin_2d() +
    scale_fill_viridis_c()





p <- tibble(time = 0:max_t, x = x, y = y) |>
    mutate(searcher = 1L) |>
    mutate(searcher = factor(searcher)) |>
    # filter(t < 100) |>
    ggplot(aes(x,y)) +
    geom_hline(yintercept = c(0, n_x+1L), color = "gray70") +
    geom_vline(xintercept = c(0, n_x+1L), color = "gray70") +
    geom_point(data = target_sims |>
                              mutate(type = type |>
                                         factor(levels = c("virus", "pseudo", "none", "both"))),
               aes(fill = type), shape = 21, size = 1, stroke = 0) +
    geom_path(aes(color = searcher), linewidth = 1, alpha = 0.5) +
    geom_point(aes(color = searcher), size = 2, alpha = 0.5) +
    scale_color_viridis_d(NULL, guide = "none", option = "plasma") +
    scale_fill_manual(NULL, values = type_pal) +
    coord_equal() +
    guides(fill = guide_legend(override.aes = list(size = 3)))


# p

anim <- p +
    labs(title = "Time: {frame_along}") +
    transition_reveal(time) +
    ease_aes("linear")

animate(anim, nframes = max_t+1, fps = 5)


# anim_save("~/Desktop/searchers.gif", anim, nframes = length(unique(b_sims$time)),
#           height = 5, width = 5, units = "in", res = 300)









