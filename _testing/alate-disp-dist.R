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



{
    t0 <- Sys.time()
    target_sims <- target_type_sims(134L+1L, 134L+1L, #<< bc it sims locations 1:(*_size-1)
                                    n_lands = 1,
                                    n_threads = 1,
                                    wt_mat = rbind(c(1, 1), c(1, 1)),
                                    n_samples = rep(round(134L*134L * 0.5), 2))
    t1 <- Sys.time()
    print(t1 - t0); rm(t0, t1)
}
# 5.513737 secs



n_samps <- c(100, 200)
wt_mat <- rbind(c(2, 1),
                c(1, 3))
n_x = 100L
n_y = 50L

stopifnot(ncol(wt_mat) == nrow(wt_mat))
stopifnot(length(n_samps) == nrow(wt_mat))
n_types <- length(n_samps)
type_wts <- rep(1, n_types)
total_pts <- n_x * n_y
out_list <- rep(list(rep(list(numeric()), n_y)), n_x)
neigh_diffs <- cbind(x = rep(-1:1, 3), y = rep(-1:1, each = 3))
coll_samps <- rep(0L, n_types)
wts <- rep(0, n_types)
probs <- rep(0, n_types)
#' `m` samples to take
#' `n` total number of sampling options
#' `i` signifies the (0-based) index of the current element
#' `k` counts the number of elements already selected for the sample.
p_fun <- function(m, n, i, k) {
    return((m - k) / (n - i))
}


# Just for testing:
out_list[[1]][[1]] <- 1
out_list[[1]][[2]] <- 2
out_list[[1]][[3]] <- 1:2
out_list[[2]][[3]] <- 2


# for (i in 1:n_x) {
#     for (j in 1:n_y) {
i = 2; j = 2



for (k in sample.int(n_types)) {
    probs[k] <- p_fun(n_samps[k], total_pts, (i-1L) * n_y + j - 1L, coll_samps[k])
    wts[k] <- 1
    for (id in -1:1) for (jd in -1:1) {
        ii <- i + id
        if (ii < 1 || ii > n_x) next
        jj <- j + jd
        if (jj < 1 || jj > n_y) next
        nk_vec <- out_list[[ii]][[jj]]
        if (length(nk_vec) == 0) next
        for (nk in nk_vec) wts[k] <- wts[k] * wt_mat[k, nk]
    }
}

probs
# wts <- wts / mean(wts)

probs |> sum()
(probs * wts / mean(wts)) |> sum()
(probs * wts / sum(wts)) |> sum()


probs |> prod()
(probs * wts / mean(wts)) |> prod()
(probs * wts / sum(wts)) |> prod()






n_x <- 50
n_y <- 20
xy_mat <- expand.grid(x = 1:n_x, y = 1:n_y)
n <- n_x * n_y
dmat <- matrix(0, n, n)
# Takes ~ 9 sec with n=1000
for (i in 1:(n-1)) for (j in (i+1):n) {
    d <- sqrt((xy_mat[i,1] - xy_mat[j,1])^2 + (xy_mat[i,2] - xy_mat[j,2])^2)
    dmat[i,j] <- dmat[j,i] <- d
}; rm(i, j)

Sigma <- 1 / dmat^3
diag(Sigma) <- 3
ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
all(ev >= -1e-6 * abs(ev[1L])) # needs to be TRUE

# Takes ~1.5 sec
rn <- MASS::mvrnorm(1, rep(0, n), Sigma)


Rcpp::cppFunction(
"double pnorm0(const double& x, const double& s) {
    double y = 0.5 * erfc(-1 * x / (s * M_SQRT2));
    return y;
}"
)


for (i in 1:1000) {
    xx <- rnorm(1)
    ss <- runif(1, 0, 2)
    if (!all.equal(pnorm0(xx, ss), pnorm(xx, sd = ss))) stop("ugh")
}; cat("YAY\n"); rm(i, xx, ss)


ru <- pnorm(rn, rep(0, n), sqrt(diag(Sigma)))

# hist(rn)
# hist(ru)

p <- c(0.4, 0.4, 0.1)
cs_p <- cumsum(p)
rt <- rep(length(p)+1, n)
rt[ru <= cs_p[1]] <- 1
rt[ru > cs_p[1] & ru <= cs_p[2]] <- 2
rt[ru > cs_p[2] & ru <= cs_p[3]] <- 3
table(rt) / n


xy_mat |>
    as_tibble() |>
    mutate(type = factor(rt)) |>
    ggplot(aes(x, y, fill = type)) +
    geom_raster() +
    scale_fill_viridis_d()

xy_mat |>
    as_tibble() |>
    mutate(type = ru) |>
    ggplot(aes(x, y, fill = type)) +
    geom_raster() +
    scale_fill_viridis_c()

wts <- dmat
diag(wts) <- 0

ape::Moran.I(rn, wts)


library(spdep)

samp_list <- xy_mat |>
    as.matrix() |>
    knearneigh(k = 8) |>
    knn2nb() |>
    nb2listw(style="B")

joincount.multi(factor(rt), samp_list)






# coll_samps[k] <- coll_samps[k] - 1L


neigh_mat <- t(matrix(c(i, j), 2, 8)) + neigh_diffs
neigh_mat <- neigh_mat[neigh_mat[,1] > 0 & neigh_mat[,2] > 0, ]
neigh_mat <- neigh_mat[neigh_mat[,1] <= n_x & neigh_mat[,2] <= n_y, ]
neigh_mat

#     }
# }; rm(i, j)





# color palette for target types:
type_pal <- c(viridis(3, begin = 0.1, end = 0.9), "gray80")[c(2,1,4,3)] |>
    set_names(c("virus", "pseudo", "none", "both"))

xs <- ys <- 134L  # approximates a hectare with 0.75 m plant spacing

#' Below converts the `type` column to a form usable in `alate_search_sims`,
#' where it's required that there's only one type per location.
#' Additionally, there should be a max of 4 types: `"virus"`, `"pseudo"`,
#' `"none"`, or `"both"` (in that order).
target_sims <- target_type_sims(xs+1L, ys+1L, #<< bc it sims locations 1:(*_size-1)
                                # n_lands = 100,
                                # n_threads = 1,
                                wt_mat = rbind(c(1, 1), c(1, 1)),
                                n_samples = rep(round(xs*ys * 0.5), 2)) |>
    getElement(1) |>
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









