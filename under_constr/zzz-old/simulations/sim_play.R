

# ======================================================================================`
# ======================================================================================`

# This file plays around with simulations.

# ======================================================================================`
# ======================================================================================`



suppressPackageStartupMessages({
    library(clonewars)
    library(tidyverse)
})


if (.Platform[['GUI']] != "X11") source(".Rprofile")




# parameters ----
#'
#' Fecundities for alates (from `mike_alate_fecundity.R`):
#'

fecunds <- list(UT3 = c(0.00000000, 0.21232513, 1.20283952, 1.79183374, 2.08507755,
                        2.20126995, 2.19932939, 2.14685547, 2.05808361, 1.87638612,
                        1.66680035, 1.26676500, 0.77241969, 0.57007240, 0.44494360,
                        0.34887375, 0.23972804, 0.12908435, 0.09602138, 0.07894355,
                        0.06138690, 0.04335143, 0.00000000, 0.00000000, 0.00000000,
                        0.00000000, 0.00000000),
                `WI-L4` = c(0.00000000, 2.33605889, 3.19374758, 3.31982305, 3.40349764,
                            3.43087162, 3.37017379, 3.26930565, 3.12693780, 2.84847149,
                            2.52952167, 1.92221193, 1.17203425, 0.86498836, 0.67512275,
                            0.52935274, 0.36374361, 0.19586192, 0.14569489, 0.11978240,
                            0.09314339, 0.06577786, 0.00000000, 0.00000000, 0.00000000,
                            0.00000000, 0.00000000))


#'
#' Parameters to generate the probabilities of alate production
#' (from `alate_production_assays.R`):
#'

alate_b0_ <- c(UT3 = -5.480921, `WI-L4` = -1.075751)
alate_b1_ <- c(UT3 = 0.0002434254, `WI-L4` = 0.0002434254)


#'
#' The rest of the dispersal parameters are the same:
#'

disp_rate_ = rep(0.8, 2)
disp_mort_ = rep(0.4, 2)
disp_start_ = rep(sum(head(dev_times$instar_days$lowT, -1)), 2)




#'
#' Leslie matrices for UT3 alates, WI-L4 alates, and apterous.
#' I assume that both lines' apterous aphids are the same and have values in between
#' the `high` and `low` growth rate lines that were used in Meisner et al. (2014).
#'

UT3_alates <- leslie_matrix(dev_times$instar_days$lowT, populations$surv_juv$low,
                        populations$surv_adult$low, fecunds[["UT3"]])
WIL4_alates <- leslie_matrix(dev_times$instar_days$lowT, populations$surv_juv$low,
                             populations$surv_adult$low, fecunds[["WI-L4"]])
apterous <- leslie_matrix(dev_times$instar_days$lowT,
                          mean(do.call(c, populations$surv_juv)),
                          colMeans(do.call(rbind, populations$surv_adult)),
                          colMeans(do.call(rbind, populations$repro)))

#'
#' Construction arrays that contain Leslie matrices for both apterous and alates:
#'
L_UT3 <- array(c(apterous, UT3_alates), c(dim(apterous), 2))
L_WIL4 <- array(c(apterous, WIL4_alates), c(dim(apterous), 2))
leslie_mat_ <- list(L_UT3, L_WIL4)

#'
#' Average and SD for carrying capacity is from the population-growth assays.
#' See `under_constr/simulations/carrying_capacity.R` for more info on why I'm using
#' the part starting with `Re(eigen(apterous, ...`.
#'
load_data(filter_pars = NULL) %>%
    group_by(line, rep) %>%
    summarize(mn = max(N)) %>%
    .[["mn"]] %>%
    {mean_K_ <<- mean(.); sd_K_ <<- sd(.)}
mean_K_ <- mean_K_ / (Re(eigen(apterous, FALSE, TRUE)$values[[1]]) - 1)
sd_K_ <- sd_K_ / (Re(eigen(apterous, FALSE, TRUE)$values[[1]]) - 1)



#'
#' Plant death parameters.
#' First two are for log-normal distribution of ages at which plants die.
#' Second two are for beta distribution of plant-mortality-induced effects on
#' growth rates.
#' See `under_constr/_assays/plant_death.R` for more info.
#'
death_prop_ <- 0.8
shape1_death_mort_ <- 3.736386
shape2_death_mort_ <- 5.777129

#'
#' Starting densities: rows are aphid stages, columns are types (alate vs apterous),
#' and slices are aphid lines.
#' (Just four 4th instars to start.)
#' `aphid_density_0_` is replicated for all `n_patches` patches.
#'
n_patches <- 4
density_0 <- array(0, c(nrow(apterous), 2, 2))
density_0[sum(head(dev_times$instar_days$highT, -1)), 1, ] <- 4
aphid_density_0_ <- replicate(n_patches, density_0, simplify = FALSE)




# one-rep function ----

sim_for_max_N <- function(.max_N, no_error = FALSE, ...) {

    args <- list(n_reps = 100,
                 max_plant_age = 0,
                 max_N = .max_N,
                 check_for_clear = seq(3, 100, 3),
                 max_t = 100,
                 save_every = 1,
                 mean_K = mean_K_,
                 sd_K = sd_K_,
                 death_prop = death_prop_,
                 shape1_death_mort = shape1_death_mort_,
                 shape2_death_mort = shape2_death_mort_,
                 disp_error = TRUE,
                 demog_error = TRUE,
                 sigma_x = environ$sigma_x,
                 rho = environ$rho,
                 extinct_N = 1,
                 aphid_name = c("competitor (UT3)",
                                "disperser (WI_L4)"),
                 leslie_mat = as.name("leslie_mat_"),
                 aphid_density_0 = as.name("aphid_density_0_"),
                 alate_b0 = alate_b0_,
                 alate_b1 = alate_b1_,
                 disp_rate = disp_rate_,
                 disp_mort = disp_mort_,
                 disp_start = disp_start_,
                 pred_rate = rep(0, n_patches),
                 n_threads = 1,
                 show_progress = FALSE)

    if (no_error) {
        args$sd_K <- 0
        args$shape1_death_mort = shape1_death_mort_ /
            (shape1_death_mort_ + shape2_death_mort_)
        args$shape2_death_mort = 0
        args$disp_error = FALSE
        args$demog_error = FALSE
        args$sigma_x = 0
    }

    other_args <- list(...)

    if (length(other_args) > 0) args[names(other_args)] <- other_args[names(other_args)]

    sim_df <- do.call(clonewars:::sim_clonewars_cpp, args)

    sim_df <- sim_df %>%
        as_tibble() %>%
        mutate(max_N = .max_N)
    return(sim_df)
}


# actual simulations ----

# set.seed(7890345)
# sim_df <- map_dfr(seq(200, 2000, 100), sim_for_max_N) %>%
#     select(max_N, everything()) %>%
#     mutate(max_N = factor(max_N))




# # To check whether the last patch has different # of aphids on it,
# # which shouldn't happen when everything's deterministic, starting densities
# # are the same among patches, plants don't die, and there is no clearing of plants.
# sim_for_max_N(.max_N = 1000, n_reps = 1,
#               no_error = TRUE,
#               max_t = 100,
#               check_for_clear = integer(0),
#               death_prop = 2,
#               mean_K = mean_K_ * 5) %>%
#     group_by(time, line) %>%
#     summarize(N0 = sum(N[patch == 0]),
#               N3 = sum(N[patch == 3]),
#               dN = N3 - N0) %>%
#     filter(dN != 0) %>%
#     print()


Z <- sim_for_max_N(.max_N = 2500, n_reps = 1,
                   no_error = TRUE,
                   max_t = 100,
                   mean_K = mean_K_ * 5)

Z %>%
    # group_by(rep, time, patch, line) %>%
    # summarize(N = sum(N)) %>%
    # ungroup() %>%
    # ggplot(aes(time, N, color = line)) +
    ggplot(aes(time, N, linetype = type, color = line)) +
    geom_line() +
    facet_wrap(~ patch) +
    scale_color_manual(values = c("firebrick", "dodgerblue")) +
    scale_linetype_manual(values = 2:1)



Z %>%
    group_by(rep, time, patch) %>%
    summarize(N = sum(N)) %>%
    ungroup() %>%
    ggplot(aes(time, N)) +
    geom_line() +
    facet_wrap(~ patch)











# # sim_df %>%
# #     group_by(max_N, rep, time, line) %>%
# #     summarize(N = sum(N)) %>%
# #     ungroup() %>%
# #     filter(rep < 9) %>%
# #     ggplot(aes(time, (N))) +
# #     geom_line(aes(color = line)) +
# #     facet_wrap(~ rep, nrow = 3) +
# #     scale_color_manual(values = c("dodgerblue", "firebrick"))
#
#
#
# sim_summ_df <- sim_df %>%
#     group_by(max_N, rep) %>%
#     filter(time == max(time)) %>%
#     group_by(max_N, rep, line) %>%
#     summarize(N = sum(N)) %>%
#     group_by(max_N, rep) %>%
#     summarize(comp = N[line == "competitor (UT3)"],
#               disp = N[line == "disperser (WI_L4)"]) %>%
#     mutate(outcome = case_when(
#         comp == 0  & disp == 0 ~ "extinct",
#         comp > 0   & disp == 0 ~ "UT3",
#         comp == 0  & disp > 0 ~ "WI-L4",
#         comp > 0   & disp > 0 ~ "coexist",
#         TRUE ~ "huh"
#     )) %>%
#     select(max_N, rep, outcome) %>%
#     group_by(max_N, outcome) %>%
#     summarize(prop = n() / 100) %>%
#     ungroup() %>%
#     mutate(outcome = factor(outcome, levels = c("extinct", "UT3", "WI-L4", "coexist"))) %>%
#     arrange(max_N, outcome)
#
#
#
# sim_summ_df %>%
#     mutate(max_N = as.numeric(paste(max_N))) %>%
#     ggplot(aes(max_N, prop)) +
#     geom_point(aes(color = outcome)) +
#     geom_line(aes(color = outcome)) +
#     scale_color_brewer(palette = "Dark2")
#
#
#
#
#
#
#
#
#
#
#
# # # =======================================================================================`
# # # =======================================================================================`
# # # =======================================================================================`
# # # =======================================================================================`
# # # =======================================================================================`
# # # =======================================================================================`
# # # =======================================================================================`
# # # =======================================================================================`
# #
# #
# #
# # source("under_constr/simulations/_preamble.R")
# #
# #
# # sim_env <- new.env()
# #
# #
# # with(sim_env, {
# #
# #     stan_fit <- readr::read_rds("inst/extdata/stan_fit.rds")
# #
# #     line_names <-
# #         clonewars::load_data() %>%
# #         .[["line"]] %>%
# #         levels()
# #     n_plants <- 8
# #     n_lines <- 8
# #     N_0 <- matrix(rep(6, n_lines * n_plants), n_plants, n_lines)
# #     max_t <- 1000
# #     # R <- rep(mean(apply(rstan::extract(stan_fit, "R", permuted = FALSE), 3, mean)), n_lines)
# #     R <- apply(rstan::extract(stan_fit, "R", permuted = FALSE), 3, mean)
# #     # A <- rep(mean(apply(rstan::extract(stan_fit, "A", permuted = FALSE), 3, mean)), n_lines)
# #     A <- apply(rstan::extract(stan_fit, "A", permuted = FALSE), 3, mean)
# #     D_vec <- clonewars::disp_estimates$b0
# #     # D_vec <- rep(mean(clonewars::disp_estimates$b0), n_lines)
# #     # Remove process and dispersal error to see patterns more easily
# #     process_error <- 0
# #     disp_error <- FALSE
# #
# #     plant_mort_0 <- rep(clonewars::plant_death$after_max_mort_coefs$b0[[1]], n_lines)
# #     plant_mort_1 <- rep(clonewars::plant_death$after_max_mort_coefs$b1[[1]], n_lines)
# #
# #     extinct_N <- 0
# #
# #     # Because we're not including any stochasticity, we can run it only once
# #     n_cages <- 1
# #     n_cores <- 1
# #
# #     plant_death_age_mean <- 1e6
# #     plant_death_age_sd <- 0.0
# #     repl_age <- 100
# #     repl_times <- seq(4, max_t, 4) - 1
# #
# #     # Modifying dispersal rates:
# #     D_mods <- c(0.010866458, 0.014979071, 0.007795972, 0.001171823,
# #                 0.011079937, 0.013445825, 0.027356186, 0.006618402)
# #     D_mods <- D_mods / mean(D_mods)
# #     D_vec <- mean(D_vec) + log(D_mods)
# #
# #     sim <- function(repl_threshold_, max_t_ = 1000) {
# #         # if (missing(repl_age_)) repl_age_ <- repl_age
# #         # lines_ <- pools[[i]]
# #         lines_ <- 1:n_lines
# #         simi <- cwsims:::sim_cages(n_cages = n_cages,
# #                                     N_0 = N_0[,lines_, drop = FALSE],
# #                                     max_t = max_t_,
# #                                     R = R[lines_],
# #                                     A = A[lines_],
# #                                     D_vec = D_vec[lines_],
# #                                     process_error = process_error,
# #                                     disp_error = disp_error,
# #                                     plant_mort_0 = plant_mort_0[lines_],
# #                                     plant_mort_1 = plant_mort_1[lines_],
# #                                     plant_death_age_mean = plant_death_age_mean,
# #                                     plant_death_age_sd = plant_death_age_sd,
# #                                     repl_times = repl_times,
# #                                     repl_age = repl_age,
# #                                     extinct_N = extinct_N,
# #                                     repl_threshold = repl_threshold_,
# #                                     n_cores = n_cores,
# #                                     line_names = lines_)
# #
# #         return(simi)
# #     }
# #
# #     rm(stan_fit)
# #
# # })
# #
# #
# # sims <- sim_env$sim(repl_threshold_ = 500) %>%
# #     .[["N"]] %>%
# #     dplyr::select(-pool) %>%
# #     arrange(date, line) %>%
# #     group_by(date) %>%
# #     mutate(prop = N / sum(N),
# #            prop_end = cumsum(prop),
# #            prop_start = lag(prop_end, default = 0)) %>%
# #     ungroup()
# #
# #
# # sims %>%
# #     mutate(line = factor(line, levels = 1:sim_env$n_plants,
# #                          labels = sprintf("%.4f | %.1f", exp(sim_env$D_vec),
# #                                           1 / sim_env$A))) %>%
# #     # filter(date < 500) %>%
# #     ggplot(aes(date, N, color = line)) +
# #     geom_line(size = 0.75) +
# #     scale_color_brewer(palette = "Dark2")
# #
# #
# #
# #
# #
# #
# #
# # sim_env$sim(repl_threshold_ = 8000, max_t_ = 1000) %>%
# #     .[["N"]] %>%
# #     dplyr::select(-pool) %>%
# #     arrange(date, line) %>%
# #     # group_by(date) %>%
# #     # mutate(prop = N / sum(N),
# #     #        prop_end = cumsum(prop),
# #     #        prop_start = lag(prop_end, default = 0)) %>%
# #     # ungroup() %>%
# #     mutate(line = factor(line, levels = 1:sim_env$n_plants,
# #                          labels = sprintf("%.4f | %.1f", exp(sim_env$D_vec),
# #                                           1 / sim_env$A))) %>%
# #     # filter(date < 500) %>%
# #     ggplot(aes(date, N, color = line)) +
# #     geom_line(size = 0.75) +
# #     scale_color_brewer(palette = "Dark2")
