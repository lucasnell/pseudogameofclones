
# ======================================================================================
# ======================================================================================

# This file plays around with simulations.

# ======================================================================================
# ======================================================================================


suppressPackageStartupMessages({
    library(clonewars)
    library(cwsims)
})


source("under_constr/simulations/_preamble.R")


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
    max_t <- 1000
    # R <- rep(mean(apply(rstan::extract(stan_fit, "R", permuted = FALSE), 3, mean)), n_lines)
    R <- apply(rstan::extract(stan_fit, "R", permuted = FALSE), 3, mean)
    # A <- rep(mean(apply(rstan::extract(stan_fit, "A", permuted = FALSE), 3, mean)), n_lines)
    A <- apply(rstan::extract(stan_fit, "A", permuted = FALSE), 3, mean)
    D_vec <- clonewars::disp_estimates$b0
    # D_vec <- rep(mean(clonewars::disp_estimates$b0), n_lines)
    # Remove process and dispersal error to see patterns more easily
    process_error <- 0
    disp_error <- FALSE

    plant_mort_0 <- rep(clonewars::plant_death$after_max_mort_coefs$b0[[1]], n_lines)
    plant_mort_1 <- rep(clonewars::plant_death$after_max_mort_coefs$b1[[1]], n_lines)

    extinct_N <- 0

    # Because we're not including any stochasticity, we can run it only once
    n_cages <- 1
    n_cores <- 1

    plant_death_age_mean <- 1e6
    plant_death_age_sd <- 0.0
    repl_age <- 100
    repl_times <- seq(4, max_t, 4) - 1

    # Modifying dispersal rates:
    D_mods <- c(0.010866458, 0.014979071, 0.007795972, 0.001171823,
                0.011079937, 0.013445825, 0.027356186, 0.006618402)
    D_mods <- D_mods / mean(D_mods)
    D_vec <- mean(D_vec) + log(D_mods)

    sim <- function(repl_threshold_, max_t_ = 1000) {
        # if (missing(repl_age_)) repl_age_ <- repl_age
        # lines_ <- pools[[i]]
        lines_ <- 1:n_lines
        simi <- cwsims:::sim_cages(n_cages = n_cages,
                                    N_0 = N_0[,lines_, drop = FALSE],
                                    max_t = max_t_,
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


sims <- sim_env$sim(repl_threshold_ = 500) %>%
    .[["N"]] %>%
    dplyr::select(-pool) %>%
    arrange(date, line) %>%
    group_by(date) %>%
    mutate(prop = N / sum(N),
           prop_end = cumsum(prop),
           prop_start = lag(prop_end, default = 0)) %>%
    ungroup()


sims %>%
    mutate(line = factor(line, levels = 1:sim_env$n_plants,
                         labels = sprintf("%.4f | %.1f", exp(sim_env$D_vec),
                                          1 / sim_env$A))) %>%
    # filter(date < 500) %>%
    ggplot(aes(date, N, color = line)) +
    geom_line(size = 0.75) +
    scale_color_brewer(palette = "Dark2")







sim_env$sim(repl_threshold_ = 8000, max_t_ = 1000) %>%
    .[["N"]] %>%
    dplyr::select(-pool) %>%
    arrange(date, line) %>%
    # group_by(date) %>%
    # mutate(prop = N / sum(N),
    #        prop_end = cumsum(prop),
    #        prop_start = lag(prop_end, default = 0)) %>%
    # ungroup() %>%
    mutate(line = factor(line, levels = 1:sim_env$n_plants,
                         labels = sprintf("%.4f | %.1f", exp(sim_env$D_vec),
                                          1 / sim_env$A))) %>%
    # filter(date < 500) %>%
    ggplot(aes(date, N, color = line)) +
    geom_line(size = 0.75) +
    scale_color_brewer(palette = "Dark2")
