
source("_testing/_preamble.R")

library(gganimate)



xs <- ys <- 134L  # approximates a hectare with 0.75 m plant spacing
n_plants <- xs * ys  # total number of plants; used below

# folder containing files created here for time-consuming sims:
alate_rds_folder <- "_testing/alate-sims-rds"


#' #' Below converts the `type` column to a form usable in `alate_search_sims`,
#' #' where it's required that there's only one type per location.
#' #' Additionally, there should be a max of 4 types: `"virus"`, `"pseudo"`,
#' #' `"none"`, or `"both"` (in that order).
#' target_sims <- target_type_sims(xs+1L, ys+1L, #<< bc it sims locations 1:(*_size-1)
#'                                 wt_mat = rbind(c(1, 1), c(1, 1)),
#'                                 n_samples = rep(round(xs*ys * 0.5), 2)) |>
#'     getElement(1) |>
#'     mutate(type = factor(type, levels = c("1", "2", "", "1_2")) |>
#'                as.integer(),
#'            #' Same but with no bacteria
#'            #' This means that `"pseudo"` is changed to `"none"` and
#'            #' `"both"` to `"virus"`
#'            nb_type = case_when(type == 2L ~ 3L,
#'                                type == 4L ~ 1L,
#'                                .default = type))
#'
#'
#'
#' sim_df <- alate_search_sims(max_t = 100,
#'                             plant_xy = target_sims[,c("x","y")] |> as.matrix(),
#'                             plant_types = target_sims$type,
#'                             alpha = 0.5,
#'                             beta = -1,
#'                             epsilon = 0.5,
#'                             n_alates = 100,
#'                             # summarize = "types",
#'                             n_threads = 6)
#'
#' plot_search_sims <- function(alate_sims,
#'                              target_sims,
#'                              anim_return = TRUE,
#'                              anim_start = FALSE,
#'                              anim_fps = 5) {
#'
#'     x_max <- max(alate_sims$x)
#'     y_max <- max(alate_sims$y)
#'
#'     p <- alate_sims |>
#'         mutate(across(c(alate, to, from), factor)) |>
#'         ggplot(aes(x,y)) +
#'         geom_rect(aes(xmin = 0, xmax = x_max+1, ymin = 0, ymax = y_max+1),
#'                   fill = NA, color = "black", linewidth = 0.75) +
#'         geom_point(data = target_sims |>
#'                        mutate(type = factor(type,
#'                                             labels = c("virus", "pseudo",
#'                                                        "none", "both"))),
#'                    aes(fill = type), shape = 21, size = 1, stroke = 0) +
#'         geom_path(aes(color = alate), linewidth = 1, alpha = 0.5) +
#'         geom_point(aes(color = alate), size = 2) +
#'         scale_color_viridis_d(NULL, guide = "none", option = "plasma") +
#'         scale_fill_manual(NULL, values = type_pal) +
#'         guides(fill = guide_legend(override.aes = list(size = 3))) +
#'         coord_equal(clip = "off", expand = FALSE) +
#'         theme(axis.line = element_blank())
#'
#'     if (!anim_return) return(p)
#'
#'     anim <- p +
#'         labs(title = "Time: {frame_along}") +
#'         transition_reveal(time) +
#'         ease_aes("linear") +
#'         theme(plot.title = element_text(margin = margin(0,0,0,b=12)))
#'
#'     if (anim_start) {
#'         animate(anim, nframes = max(alate_sims$time)+1, fps = anim_fps)
#'     }
#'
#'     return(anim)
#' }
#'
#'
#' dt <- 10
#'
#' list("different types repel" = c(same = 1, diff = 1/dt),
#'      "different types attract" = c(same = 1, diff = dt),
#'      "same types attract" = c(same = dt, diff = 1),
#'      "same types repel" = c(same = 1/dt, diff = 1),  ## same as neutral
#'      "same repel, diff. attract" = c(same = 1/dt, diff = dt),
#'      "same attract, diff. repel" = c(same = dt, diff = 1/dt),
#'      "neutral" = c(same = 1, diff = 1)) |>
#'     imap(\(x, n) {
#'         wt_mat <- rbind(c(x[["same"]], x[["diff"]]),
#'                         c(x[["diff"]], x[["same"]]))
#'         t_sim_df <- target_type_sims(xs+1L, ys+1L,
#'                                     wt_mat,
#'                                     n_samples = rep(round(xs*ys * 0.4), 2)) |>
#'             getElement(1) |>
#'             mutate(type = factor(type, levels = c("1", "2", "", "1_2"),
#'                                  labels = c("virus", "pseudo", "none", "both")))
#'         t_sim_df |>
#'             ggplot(aes(x, y, color = type)) +
#'             geom_rect(aes(xmin = 0, xmax = xs+1L, ymin = 0, ymax = ys+1),
#'                       fill = NA, color = "black", linewidth = 0.75) +
#'             geom_point(size = 0.5) +
#'             scale_color_manual(NULL, values = type_pal) +
#'             coord_equal(clip = "off", expand = FALSE) +
#'             guides(color = guide_legend(override.aes = list(size = 3))) +
#'             ggtitle(n) +
#'             theme(axis.line = element_blank(),
#'                   axis.title = element_blank(),
#'                   axis.text = element_blank(),
#'                   axis.ticks = element_blank(),
#'                   plot.title = element_text(margin = margin(0,0,0,b=12)))
#'     }) |>
#'     c(list(nrow = 2, guides = "collect")) |>
#'     do.call(what = wrap_plots)







# =============================================================================*
# =============================================================================*
# Simulations ----
# =============================================================================*
# =============================================================================*


#' Variables to vary:
#'
#' 1. density of each plant type
#' 2. plant type landscape, one of the following:
#'     a. different types repel: same = 1, diff < 1
#'     b. different types attract: same = 1, diff > 1
#'     c. same types attract: same > 1 diff = 1
#'     d. same attract, diff. repel: same > 1, diff < 1
#'     e. neutral: same = 1, diff = 1
#' 3. `alpha`: effect of virus infection on alate alighting
#' 4. `beta`: effect of *Pseudomonas* infection on alate alighting
#' 5. `epsilon`: effect of virus infection on alate acceptance
#'
#'
#' Descriptions of other parameters below:
#'
#' `vd`: starting density of virus
#' `pd`: density of *Pseudomonas*
#' `land_type`: type of landscape (letters follow from above)
#' `n_lands`: separate landscapes per iteration
#' `n_alates`: number of alates per landscape
#'


one_alate_sim <- function(vd, pd, land_type, alpha, beta, epsilon, n_alates, seed,
                          n_lands = 12L) {

    stopifnot(length(vd) == 1 && inherits(vd, "numeric") && vd > 0 && vd < 1)
    stopifnot(length(pd) == 1 && inherits(pd, "numeric") && pd >= 0 && pd <= 1)
    stopifnot(is.character(land_type) && length(land_type) == 1)
    stopifnot(land_type %in% letters[1:5])
    stopifnot(length(alpha) == 1 && inherits(alpha, "numeric"))
    stopifnot(length(beta) == 1 && inherits(beta, "numeric"))
    stopifnot(length(epsilon) == 1 && inherits(epsilon, "numeric") &&
                  epsilon > 0 && epsilon * formals(alate_search_sims)[["w"]] < 1)
    stopifnot(length(seed) == 1 && is.numeric(seed))

    dt <- 2

    same_diff <- switch(land_type,
                        a = c(1,  1/dt),  # different types repel
                        b = c(1,  dt),  # different types attrac
                        c = c(dt, 1),  # same types attract
                        d = c(dt, 1/dt),  # same attract, diff. repel
                        e = c(1,  1))  # neutral

    wt_mat <- rbind(same_diff, rev(same_diff))
    densities <- round(xs*ys * c(vd, pd))
    # Then happens when Pseudomonas is not present:
    if (densities[2] == 0) {
        densities <- densities[1]
        wt_mat <- wt_mat[1,1,drop=FALSE]
    }

    set.seed(seed)
    target_sims <- target_type_sims(xs, ys,
                                    wt_mat = wt_mat,
                                    n_samples = densities,
                                    n_lands = n_lands,
                                    n_threads = .n_threads) |>
        map(mutate, type = as.integer(factor(type, levels = c("1", "2", "", "1_2"))))


    alate_df <- map_dfr(1:length(target_sims), \(i) {
        alate_search_sims(max_t = 100,
                          plant_xy = target_sims[[i]][,c("x","y")] |> as.matrix(),
                          plant_types = target_sims[[i]][["type"]],
                          alpha = alpha,
                          beta = beta,
                          epsilon = epsilon,
                          n_alates = n_alates,
                          n_threads = .n_threads) |>
            mutate(land = i)
    })

    return(alate_df)
}


set.seed(1372060266)
alate_sim_combos <- crossing(.vd = c(0.1, 0.3, 0.5),
                                  .pd = c(0, 0.1, 0.5, 0.9),
                                  .land_type = letters[1:5],
                                  .alpha = c(-1, 0, 1),
                                  .beta = c(-1, 0, 1),
                                  .epsilon = c(0.25, 1, 2),
                                  .n_alates = c(10, 100, 1000, 10e3)) |>
    # So landscapes from individual runs can be replicated:
    mutate(.seed = sample.int(2^31-1, n()),
           i = 1:n())

if (!dir.exists(alate_rds_folder) ||
    length(list.files(alate_rds_folder, "*.rds")) < nrow(alate_sim_combos)) {

    if (!dir.exists(alate_rds_folder)) dir.create(alate_rds_folder)

    # Takes ~4.4 hrs
    x <- alate_sim_combos |>
        pmap(\(.vd, .pd, .land_type, .alpha, .beta, .epsilon, .n_alates, .seed, i) {
            file_name <- sprintf("%s/sim-%04i.rds", alate_rds_folder, i)
            if (file.exists(file_name)) return(invisible(NULL))
            out <- one_alate_sim(.vd, .pd, .land_type, .alpha, .beta, .epsilon,
                                 .n_alates, .seed) |>
                mutate(vd = .vd, pd = .pd, land_type = .land_type, alpha = .alpha,
                       beta = .beta, epsilon = .epsilon, n_alates = .n_alates) |>
                select(vd, pd, land_type, alpha, beta, epsilon, n_alates,
                       land, alate, everything())
            attr(out, "seed") <- .seed
            write_rds(out, file_name, compress = "xz")
            return(invisible(NULL))
        }, .progress = TRUE)
    rm(x)
}


# probability of virus loading onto alate and plant:
d_a = 0.5
d_p = 0.5


sim_df <- list.files(alate_rds_folder, "*.rds", full.names = TRUE)[[1]] |>
    read_rds()

# 1. Final locations of alates (1 for each `land` and `alate`)
sim_df |>
    group_by(land, alate) |>
    filter(time == max(time)) |>
    ungroup()



# 2. Proportion of alates and plants that are inoculated (1 for each `land`)
n_lands <- sim_df$land |> unique() |> length()
n_alates <- sim_df$n_alates[[1]]
# n_plants



Rcpp::sourceCpp("_testing/alate-sims.cpp")
set.seed(1452903413); U <- runif(nrow(sim_df)*2)

ui <- 1

z <- sim_df |>
    split(~ land) |>
    map(\(d) {
        d |>
            split(~ alate) |>
            map(\(di) {
                to <- di$to
                # outputs: whether alate inoculated, list of plants inoculated:
                a_inoc <- FALSE
                p_inocs <- c()

                # When did it probe onto an infected plant(s)?
                on_inf <- which(to == 1L | to == 4L)
                if (length(on_inf) == 0) {
                    return(tibble(a_inoc = FALSE, p_inocs = list(c())))
                }

                # When, if ever, did it get inoculated with virus?
                inoc_a_ind <- length(to) + 1
                for (i in on_inf) {
                    # u <- runif(1)
                    u <- U[ui]
                    ui <<- ui + 1
                    if (ui > length(U)) stop("ui > length(U)")
                    if (u < d_a) {
                        inoc_a_ind <- i
                        a_inoc <- TRUE
                        break
                    }
                }
                if (!a_inoc) return(tibble(a_inoc = FALSE, p_inocs = list(c())))

                # When, if ever, did it probe uninfected plants?
                on_uninf <- which(to == 2L | to == 3L)
                if (length(on_uninf) == 0 || all(on_uninf < inoc_a_ind)) {
                    return(tibble(a_inoc = a_inoc, p_inocs = list(c())))
                }
                on_uninf <- on_uninf[on_uninf > inoc_a_ind]
                for (i in on_uninf) {
                    # u <- runif(1)
                    u <- U[ui]
                    ui <<- ui + 1
                    if (ui > length(U)) stop("ui > length(U)")
                    k <- ((di$y[i]-1) * xs + (di$x[i]-1))
                    if (u < d_p) p_inocs <- c(p_inocs, k)
                }

                return(tibble(a_inoc = a_inoc, p_inocs = list(p_inocs)))

            }) |>
            list_rbind() |>
            summarize(a_inoc = mean(a_inoc),
                      p_inoc = length(unique(do.call(c, p_inocs))) / n_plants) |>
            mutate(land = d$land[[1]]) |>
            select(land, everything())
    }) |>
    list_rbind()









# Alate inoculations ----
#'
#' Count the number of times each alate probes a virus-infected plant.
#'
alate_inf_encs <- alate_sims |>
    group_by(vd, pd, land_type, alpha, beta, epsilon, land, alate) |>
    summarize(n_encs = (\(z) sum(z == 1L))(to), .groups = "drop")


#'
#' Count the number of times each alate probes a virus-infected plant before
#' each uninfected plant it probes.
#'

# d <- alate_sims[1:4,]
d <- alate_sims |>
    filter(vd == 0.5, pd == 0, land_type == "a", alpha == 0, beta == 0,
           epsilon == 1, land == 1, alate == 4)

to <- d$to
from <- d$from

# When did it probe onto an infected and uninfected plant(s)?
on_inf <- which(to == 1L)
on_uninf <- which(to == 2L | to == 3L)

if (length(on_inf) == 0 || all(on_uninf < on_inf[1])) {
    e_jk <- 0L
} else {

    map_int(on_uninf, \(x) sum(on_inf < x))

    # on_uninf <- on_uninf[on_uninf > on_inf[1]]
    # out[["n_pass"]] <-
}

microbenchmark::microbenchmark(a = !any(d$to == 1),
                               b = sum(d$to == 1) == 0)

on_inf
on_uninf



# Spatial clustering ----










#' # =============================================================================*
#' # =============================================================================*
#' # old code ----
#' # =============================================================================*
#' # =============================================================================*
#'
#' #' From "The Role of Aphid Behaviour in the Epidemiology of Potato Virus Y:
#' #' a Simulation Study" by Thomas Nemecek (1993; p. 72), dispersal distances
#' #' follow a Weibull distribution with shape = 0.6569 and scale = 9.613.
#' #'
#' #' We'll use the median of this distribution for our dispersal radius,
#' #' where alates have to move to a patch somewere within it.
#' #'
#' #' I'm dividing by 0.75 to convert from meters to plant locations that are
#' #' 0.75 meters apart (typical spacing for pea).
#' #'
#' radius <- qweibull(0.5, 0.6569, 9.613) / 0.75
#'
#' #' For a given distance or below above the current plant, how many plants
#' #' to the left or right are within the radius?
#' sapply(1:floor(radius), \(d) sqrt(radius^2 - d^2))
#'
#' # color palette for target types:
#' type_pal <- c(viridis(3, begin = 0.1, end = 0.9), "gray80")[c(2,1,4,3)] |>
#'     set_names(c("virus", "pseudo", "none", "both"))
#'
#' xs <- ys <- 134L  # approximates a hectare with 0.75 m plant spacing
#'
#' #' Below converts the `type` column to a form usable in `alate_search_sims`,
#' #' where it's required that there's only one type per location.
#' #' Additionally, there should be a max of 4 types: `"virus"`, `"pseudo"`,
#' #' `"none"`, or `"both"` (in that order).
#' target_sims <- target_type_sims(xs+1L, ys+1L, #<< bc it sims locations 1:(*_size-1)
#'                                 wt_mat = rbind(c(1, 1), c(1, 1)),
#'                                 n_samples = rep(round(xs*ys * 0.5), 2)) |>
#'     mutate(type = type |>
#'                map_int(\(x) {
#'                    if (length(x) == 1L) return(x)
#'                    else return(4L)
#'                    }),
#'            #' Same but with no bacteria
#'            #' This means that `"pseudo"` is changed to `"none"` and
#'            #' `"both"` to `"virus"`
#'            nb_type = case_when(type == 2L ~ 3L,
#'                                type == 4L ~ 1L,
#'                                .default = type))
#'
#'
#' #' Making matrices of out types:
#' tt_mat <- nb_tt_mat <- matrix(NA_integer_, xs, ys)
#' for (i in 1:nrow(target_sims)) {
#'     x <- target_sims$x[[i]]
#'     y <- target_sims$y[[i]]
#'     tt_mat[x,y] <- target_sims$type[[i]]
#'     nb_tt_mat[x,y] <- target_sims$nb_type[[i]]
#' }; rm(i, x, y)
#'
#'
#'
#' # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ----
#' # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ----
#'
#'
#'
#'
#' max_t = 100
#' target_types = tt_mat
#' radius = 7.336451
#' alpha = 0.5
#' beta = -1
#' w = 0.2
#' epsilon = 0.5
#'
#' # xy0 = NULL
#' # randomize_xy0 = TRUE
#' # n_searchers = 1
#'
#' stopifnot(ncol(target_types) > 0 && nrow(target_types) > 0)
#'
#' n_x <- ncol(target_types)
#' n_y <- nrow(target_types)
#'
#' stopifnot(all(target_types %in% 1:4))
#'
#' # Binary variable for whether a plant has virus:
#' v <- apply(target_types, c(1,2),
#'            \(tt) as.integer(tt == 1L | tt == 4L))
#' # Binary variable for whether a plant has *Pseudomonas*:
#' b <- apply(target_types, c(1,2),
#'            \(tt) as.integer(tt == 2L | tt == 4L))
#'
#' # change in x and y for neighbors (those within radius)
#' neigh_dxdy <- lapply(-floor(radius):floor(radius), \(dx) {
#'     ddy <- floor(sqrt(radius^2 - dx^2))
#'     dy <- -ddy:ddy
#'     return(lapply(dy, \(.dy) c(dx, .dy)))
#' }) |>
#'     do.call(what = c) |>
#'     do.call(what = rbind)
#'
#'
#'
#'
#'
#'
#' # Location of aphid:
#' x <- rep(NA_real_, max_t+1)
#' y <- rep(NA_real_, max_t+1)
#' # Types of plants the aphid traveled to and from:
#' to <- rep(NA_integer_, max_t+1)
#' from <- rep(NA_integer_, max_t+1)
#'
#' x[1] <- sample.int(n_x, 1)
#' y[1] <- sample.int(n_y, 1)
#' to[1] <- target_types[x[1], y[1]]
#' from[1] <- target_types[x[1], y[1]]  # to avoid NAs
#'
#' for (t in 1:max_t) {
#'
#'     # feed_p <- w
#'     # if (v[x[t],y[t]]) feed_p <- feed_p * epsilon
#'     #
#'     # if (runif(1) < feed_p) {
#'     #     # x[(t+1):length(x)] <- x[t]
#'     #     # y[(t+1):length(y)] <- y[t]
#'     #     # to[(t+1):length(x)] <- to[t]
#'     #     # from[(t+1):length(y)] <- from[t]
#'     #     break
#'     # }
#'
#'     # Sampling weights:
#'     s_wts <- apply(neigh_dxdy, 1,
#'                    \(dxdy) {
#'                        .x <- x[t] + dxdy[1]
#'                        .y <- y[t] + dxdy[2]
#'                        if (.x <= 0 | .x > n_x) return(0)
#'                        if (.y <= 0 | .y > n_y) return(0)
#'                        return(exp(alpha * v[.x,.y] + beta * b[.x,.y]))
#'                    })
#'
#'     neigh_i <- sample.int(length(s_wts), 1, prob = s_wts / sum(s_wts))
#'
#'     x[t+1] <- x[t] + neigh_dxdy[neigh_i,1]
#'     y[t+1] <- y[t] + neigh_dxdy[neigh_i,2]
#'     to[t+1] <- target_types[x[t+1], y[t+1]]
#'     from[t+1] <- to[t]
#'
#' }
#'
#'
#'
#' tibble(x = x, y = y) |>
#'     mutate(type = map2_chr(x, y, \(i, j) target_types[i,j]),
#'            type = factor(type, levels = 1:4,
#'                          labels = c("virus", "pseudo", "none", "both"))) |>
#'     group_by(type) |>
#'     summarize(hits = n())
#'
#' tibble(from = from, to = to) |>
#'     mutate(across(everything(), \(x) {
#'         factor(x, levels = 1:4,
#'                labels = c("virus", "pseudo", "none", "both"))
#'     })) |>
#'     group_by(from, to, .drop = FALSE) |>
#'     summarize(n = n(), .groups = "drop")
#'
#'
#' tibble(from = from, to = to) |>
#'     mutate(across(everything(), \(x) {
#'         factor(x, levels = 1:4,
#'                labels = c("virus", "pseudo", "none", "both"))
#'     })) |>
#'     ggplot(aes(from, to)) +
#'     geom_bin_2d() +
#'     scale_fill_viridis_c()
#'
#'
#'
#'
#'
#' p <- tibble(time = 0:max_t, x = x, y = y) |>
#'     mutate(searcher = 1L) |>
#'     mutate(searcher = factor(searcher)) |>
#'     # filter(t < 100) |>
#'     ggplot(aes(x,y)) +
#'     geom_hline(yintercept = c(0, n_x+1L), color = "gray70") +
#'     geom_vline(xintercept = c(0, n_x+1L), color = "gray70") +
#'     geom_point(data = target_sims |>
#'                               mutate(type = type |>
#'                                          factor(levels = c("virus", "pseudo", "none", "both"))),
#'                aes(fill = type), shape = 21, size = 1, stroke = 0) +
#'     geom_path(aes(color = searcher), linewidth = 1, alpha = 0.5) +
#'     geom_point(aes(color = searcher), size = 2, alpha = 0.5) +
#'     scale_color_viridis_d(NULL, guide = "none", option = "plasma") +
#'     scale_fill_manual(NULL, values = type_pal) +
#'     coord_equal() +
#'     guides(fill = guide_legend(override.aes = list(size = 3)))
#'
#'
#' # p
#'
#' anim <- p +
#'     labs(title = "Time: {frame_along}") +
#'     transition_reveal(time) +
#'     ease_aes("linear")
#'
#' animate(anim, nframes = max_t+1, fps = 5)
#'
#'
#' # anim_save("~/Desktop/searchers.gif", anim, nframes = length(unique(b_sims$time)),
#' #           height = 5, width = 5, units = "in", res = 300)
#'
