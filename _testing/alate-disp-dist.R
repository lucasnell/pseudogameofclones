library(tidyverse)
library(pseudogameofclones)

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


xs <- ys <- 135L  # approximates a hectare with 0.75 m plant spacing

#' Below converts the `type` column to a form usable in `alate_searching`,
#' where it's required that there's only one type per location.
#' Additionally, there should be a max of 4 types: `"virus"`, `"pseudo"`,
#' `"none"`, or `"both"`.
target_sims <- target_type_sims(xs, ys,
                                wt_mat = rbind(c(1, 1), c(1, 1)),
                                n_samples = rep(round(xs*ys * 0.5), 2)) |>
    mutate(type = type |>
               map_chr(\(x) {
                   opts <- c("virus", "pseudo", "none", "both")
                   if (length(x) == 1L) return(opts[x])
                   else return(opts[4L])
                   }),
           #' Same but with no bacteria
           #' This means that `"pseudo"` is changed to `"none"` and
           #' `"both"` to `"virus"`
           nb_types = case_when(type == "pseudo" ~ "none",
                                type == "both" ~ "virus",
                                .default = type))






# # choose_plant <- function(x, b0, b1) {
# # choose_plant <- function(x, p) {
# choose_plant <- function(x, zeta) {
#     stopifnot(all(x >= 0))
#     # # # convert from units of m to row spacing (30" or ~0.75 m)
#     # # x <- x / 0.75
#     # # Allows passing multiple `b1` values for multiple powers of `x`:
#     # X <- matrix(x, length(x), length(b1))
#     # if (length(b1) > 1) for (k in 2:length(b1)) X[,k] <- X[,k]^k
#     # Z <- X %*% b1
#     # z <- inv_logit(b0 + Z)
#     # # z <- numeric(length(x))
#     # # z[1] <- p[1]
#     # # for (i in 2:length(x)) {
#     # #     z[i] <- z[i-1] * p[i]
#     # # }
#     # z <- exp(zeta * x^2)
#     z <- 1 / x^(zeta)
#     z <- z / sum(z)
#     return(z)
# }
#
#
# dists <- 1:20
#
# # From "The Role of Aphid Behaviour in the Epidemiology of Potato Virus Y:
# # a Simulation Study" by Thomas Nemecek (1993; p. 72)
# # shape = 0.6569
# # scale = 9.613
# obs <- map_dbl(dists, \(x) diff(pweibull(c(x-1, x), 0.6569, 9.613)))
#
#
#
#
#
# weib_fit_fun <- function(pars) {
#     # b0 <- - exp(pars[1])
#     # b1 <- pars[-1]
#     # fit <- choose_plant(dists, b0, b1)
#     # # p <- inv_logit(pars)
#     # # fit <- choose_plant(dists, p)
#     obs <- obs / sum(obs)
#     fit <- choose_plant(dists, pars)
#     return(sum(abs(obs - fit)))
# }
# # (weib_op <- optim(c(1, rep(1, 1)), weib_fit_fun))
#
# # weib_op <- aphidsync::winnowing_optim(weib_fit_fun,
# #                                       rep(-10, 3),
# #                                       rep(10, 3),
# #                                       polished_control = list(maxit = 10e3L,
# #                                                               reltol = 1e-10),
# #                                       n_bevals = 1000L)
# # (weib_op <- weib_op[[1]])
#
# # (weib_op <- optim(runif(45), weib_fit_fun, control = list(maxit = 100e3L, abstol = 1e-10)))
# (weib_op <- optimize(weib_fit_fun, c(0, 3)))
#
# # curve(sapply(x, \(xx) weib_fit_fun(xx)), -1, 0)
# plot(dists, choose_plant(dists, 0.9250706), ylim = c(0, 0.25), type = "l"); lines(dists, obs / sum(obs), lty = "22", col = "red")
#
#
#
# # (weib_pars <- c(- exp(weib_op$par[1]), weib_op$par[-1]))
# # (weib_pars <- inv_logit(weib_op$par))
# # plot(weib_pars, type = "l")
#
#
# weib_cdf <- function(x, L = 9.613, k = 0.6569) {
#     z0 <- 1 - exp(-((x-1)/L)^k)
#     z <- 1 - exp(-(x/L)^k)
#     return(z - z0)
# }
#
# {
#     plot(dists, obs, type = "l")
#     # lines(dists, choose_plant(dists, b0 = weib_pars[1], b1 = weib_pars[-1]),
#     lines(dists, choose_plant(dists, weib_op$minimum),
#           # lines(dists, choose_plant(dists, p = weib_pars),
#           # curve(weib_cdf(x), add = TRUE,
#           col = "red", lty = "22")
# }
