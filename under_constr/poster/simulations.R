
# ======================================================================================
# ======================================================================================

# This file does only the simulations and no analyses of those simulations.

# ======================================================================================
# ======================================================================================


suppressPackageStartupMessages({
    library(clonewars)
})





sim_env <- new.env()


with(sim_env, {

    stan_fit <- read_rds("data-raw/stan_fit.rds")

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
    extinct_N <- 6
    n_cages <- 1000
    n_cores <- parallel::detectCores() - 1

    # All combinations of pools from 1 to 8
    pools <- map(1:sim_env$n_lines, ~ combn(1:n_lines, .x) %>%
                     t() %>%
                     split(1:nrow(.)) %>%
                     set_names(NULL)) %>%
        flatten()

    # Running longer to try to find patterns
    max_t <- 2000
    n_cages <- 100
    repl_times <- seq(4, max_t, 4) - 1

    # Also want to remove process error to see patterns more easily
    process_error <- 0

    sim <- function(i) {
        lines_ <- pools[[i]]
        simi <- cwsims:::sim_cages(n_cages = n_cages,
                                    N_0 = N_0[,lines_, drop = FALSE],
                                    max_t = max_t,
                                    R = R[lines_],
                                    A = A[lines_],
                                    D_vec = D_vec[lines_],
                                    process_error = process_error,
                                    plant_mort_0 = plant_mort_0[lines_],
                                    plant_mort_1 = plant_mort_1[lines_],
                                    plant_death_age_mean = plant_death_age_mean,
                                    plant_death_age_sd = plant_death_age_sd,
                                    repl_times = repl_times,
                                    repl_age = repl_age,
                                    extinct_N = extinct_N,
                                    n_cores = n_cores,
                                    line_names = lines_)

        return(simi)
    }

    rm(stan_fit)

readr::write_rds(sim_env, "data-raw/sim_env.rds")

})

library(progress)

pb <- progress_bar$new(
    format = "  simulating [:bar] :percent in :elapsed",
    total = with(sim_env, sum(sapply(pools, length))),
    clear = FALSE, width = options("width")$width)


# Takes ~11 min
pool_sims_N <- rep(list(NA), length(sim_env$pools))
pool_sims_X <- rep(list(NA), length(sim_env$pools))
pool_sims_Z <- rep(list(NA), length(sim_env$pools))
set.seed(549489)
for (i in 1:length(sim_env$pools)) {

    n_ <- length(sim_env$pools[[i]])
    sims_ <- sim_env$sim(i)

    pool_sims_N[[i]] <- sims_$N
    pool_sims_X[[i]] <- sims_$X
    pool_sims_Z[[i]] <- sims_$Z

    pb$tick(n_)

}; rm(i, n_, sims_)


# Takes ~30 sec
pool_sims_N <- bind_rows(pool_sims_N)
invisible(gc())
pool_sims_X <- bind_rows(pool_sims_X)
invisible(gc())
pool_sims_Z <- bind_rows(pool_sims_Z)
invisible(gc())

pool_sims_N <- pool_sims_N %>%
    mutate(pool = as.integer(pool))
pool_sims_Z <- pool_sims_Z %>%
    mutate(pool = as.integer(pool))

# readr::write_rds(pool_sims_N, "data-raw/pool_sims_N.rds")
# readr::write_rds(pool_sims_X, "data-raw/pool_sims_X.rds")
# readr::write_rds(pool_sims_Z, "data-raw/pool_sims_Z.rds")


# Takes ~30 sec
readr::write_rds(pool_sims_N, "data-raw/pool_sims_N.rds")
readr::write_rds(pool_sims_X, "data-raw/pool_sims_X.rds")
readr::write_rds(pool_sims_Z, "data-raw/pool_sims_Z.rds")

# Takes ~20 sec
pool_sims_N <- pool_sims_N %>%
    group_by(pool, rep, line) %>%
    summarize(N = mean(N)) %>%
    ungroup()
readr::write_rds(pool_sims_N, "data-raw/pool_sims_N_by_linerep.rds")




