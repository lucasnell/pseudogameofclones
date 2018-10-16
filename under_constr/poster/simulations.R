
# ======================================================================================
# ======================================================================================

# This file does only the simulations and no analyses of those simulations.

# ======================================================================================
# ======================================================================================


suppressPackageStartupMessages({
    library(clonewars)
    library(cwsims)
})





sim_env <- new.env()


with(sim_env, {

    stan_fit <- readr::read_rds("data-raw/stan_fit.rds")

    line_names <-
        clonewars::load_data() %>%
        .[["line"]] %>%
        levels()
    n_plants <- 8
    n_lines <- 8
    N_0 <- matrix(rep(6, n_lines * n_plants), n_plants, n_lines)
    max_t <- 180
    R <- apply(rstan::extract(stan_fit, "R", permuted = FALSE), 3, mean) %>%
        as.numeric()
    A <- apply(rstan::extract(stan_fit, "A", permuted = FALSE), 3, mean) %>%
        as.numeric()
    D_vec <- clonewars::disp_estimates$b0
    process_error <- apply(rstan::extract(stan_fit, "s_epsilon", permuted = FALSE),
                           3, mean) %>%
        as.numeric()
    plant_mort_0 <- clonewars::plant_death$after_max_mort_coefs$b0
    plant_mort_1 <- clonewars::plant_death$after_max_mort_coefs$b1
    plant_death_age_mean <- clonewars::plant_death$until_max_summ$max_mean
    plant_death_age_sd <- clonewars::plant_death$until_max_summ$max_sd
    repl_times <- seq(4, max_t, 4) - 1
    repl_age <- 3
    extinct_N <- 1
    n_cages <- 100
    n_cores <- parallel::detectCores() - 1

    # # All combinations of pools from 1 to 8
    # pools <- map(1:sim_env$n_lines, ~ combn(1:n_lines, .x) %>%
    #                  t() %>%
    #                  split(1:nrow(.)) %>%
    #                  set_names(NULL)) %>%
    #     flatten()

    # Running longer to try to find patterns
    max_t <- 500
    # n_cages <- 100
    # repl_times <- seq(1, max_t, 1) - 1
    plant_death_age_mean <- 1e6
    plant_death_age_sd <- 0.0
    repl_age <- 0

    # Also want to remove process error to see patterns more easily
    process_error <- 0
    # process_error <- process_error / 10
    disp_error <- FALSE
    # Because we're not including any stochasticity, we can run it only once
    n_cages <- 1
    n_cores <- 1

    sim <- function(repl_threshold_) {
        # if (missing(repl_age_)) repl_age_ <- repl_age
        # lines_ <- pools[[i]]
        lines_ <- 1:n_lines
        simi <- cwsims:::sim_cages(n_cages = n_cages,
                                    N_0 = N_0[,lines_, drop = FALSE],
                                    max_t = max_t,
                                    R = R[lines_],
                                    A = A[lines_],
                                    D_vec = D_vec[lines_],
                                    process_error = process_error,
                                    disp_error = disp_error,
                                    plant_mort_0 = plant_mort_0[lines_],
                                    plant_mort_1 = plant_mort_1[lines_],
                                    plant_death_age_mean = plant_death_age_mean,
                                    plant_death_age_sd = plant_death_age_sd,
                                    repl_times = repl_times,
                                    repl_age = repl_age,
                                    extinct_N = extinct_N,
                                    repl_threshold = repl_threshold_,
                                    n_cores = n_cores,
                                    line_names = lines_)

        return(simi)
    }

    rm(stan_fit)

})


# readr::write_rds(sim_env, "data-raw/sim_env.rds")


# surv_by_thresh <- function(thresh) {
#
#     sims_ <- sim_env$sim(repl_threshold_ = thresh)
#     pool_sims_N <- sims_$N %>%
#         dplyr::select(-pool) %>%
#         mutate(thresh = thresh1)
#
#     pool_sims_N %>%
#         filter(date == max(date)) %>%
#         group_by(rep, line) %>%
#         summarize(s = ifelse(N == 0, 0L, 1L)) %>%
#         ungroup() %>%
#         spread(line, s) %>%
#         dplyr::select(-rep) %>%
#         as.matrix() %>%
#         colMeans()
# }
#
# threshes <- seq(350, 575, 10)
# t0 <- Sys.time()
# set.seed(549489+2)
# sbt <- map(threshes, surv_by_thresh)
# Sys.time() - t0; rm(t0)
#
#
# do.call(rbind, sbt) %>%
#     as_data_frame() %>%
#     mutate(thresh = as.integer(threshes)) %>%
#     filter_at(vars(-thresh), all_vars(. != 1.0))
#
# do.call(rbind, sbt) %>%
#     as_data_frame() %>%
#     mutate(thresh = as.integer(threshes)) %>%
#     gather("line", "p", -thresh, factor_key = TRUE) %>%
#     ggplot(aes(thresh, p, color = line)) +
#     geom_line() +
#     geom_point() +
#     scale_color_brewer(palette = "Dark2") +
#     facet_wrap(~ line) +
#     NULL


test_threshes <- function(threshes) {
    sim_list <- lapply(threshes,
           function(thresh_) {
               sims_ <- sim_env$sim(repl_threshold_ = thresh_)
               sims_$N <- sims_$N %>%
                   dplyr::select(-pool) %>%
                   mutate(thresh = thresh_)
               sims_$X <- sims_$X %>%
                   dplyr::select(-pool) %>%
                   mutate(thresh = thresh_)
               sims_$Z <- sims_$Z %>%
                   dplyr::select(-pool) %>%
                   mutate(thresh = thresh_)
               return(sims_)
           })
    pool_sims <- list()
    pool_sims$N <- bind_rows(lapply(sim_list, function(x) x$N)) %>%
        dplyr::select(-rep) %>%
        dplyr::mutate_at(vars(line, thresh), factor) %>%
        dplyr::select(thresh, everything())
    pool_sims$X <- bind_rows(lapply(sim_list, function(x) x$X)) %>%
        dplyr::select(-rep) %>%
        dplyr::mutate_at(vars(plant, thresh), factor) %>%
        dplyr::select(thresh, everything())
    pool_sims$Z <- bind_rows(lapply(sim_list, function(x) x$Z)) %>%
        dplyr::select(-rep) %>%
        dplyr::mutate_at(vars(plant, thresh), factor) %>%
        dplyr::select(thresh, everything())
    return(pool_sims)
}




