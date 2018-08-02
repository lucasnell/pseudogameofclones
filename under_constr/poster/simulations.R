
# ======================================================================================
# ======================================================================================

# This file does only the simulations and no analyses of those simulations.

# ======================================================================================
# ======================================================================================


suppressPackageStartupMessages({
    library(clonewars)
    library(bigmemory)
})




stan_fit <- read_rds("data-raw/stan_fit.rds")

sim_env <- new.env()


with(sim_env, {

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
    D_binom <- clonewars::disp_estimates$binom
    D_nb <- clonewars::disp_estimates$nb
    process_error <- apply(rstan::extract(stan_fit, "s_epsilon", permuted = FALSE),
                           3, mean) %>%
        as.numeric()
    plant_mort_0 <- clonewars::plant_death$after_max_mort_coefs$inter
    plant_mort_1 <- clonewars::plant_death$after_max_mort_coefs$date
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

    # Creating objects to use C++ function directly
    D_mat <- as.matrix(cbind(D_binom[,c("b0", "b1", "b2")],
                             D_nb[,c("b0", "theta")]))
    colnames(D_mat) <- NULL

    sim <- function(i) {
        lines_ <- pools[[i]]
        simi <- cwsims:::sim_cages_(n_cages = n_cages,
                                    N_0 = N_0[,lines_, drop = FALSE],
                                    max_t = max_t,
                                    R = R[lines_],
                                    A = A[lines_],
                                    D_mat = D_mat[lines_,, drop = FALSE],
                                    process_error = process_error,
                                    plant_mort_0 = plant_mort_0[lines_],
                                    plant_mort_1 = plant_mort_1[lines_],
                                    plant_death_age_mean = plant_death_age_mean,
                                    plant_death_age_sd = plant_death_age_sd,
                                    repl_times = repl_times,
                                    repl_age = repl_age,
                                    extinct_N = extinct_N,
                                    n_cores = n_cores,
                                    by_cage = FALSE,
                                    show_progress = FALSE)
        for (j in 1:length(simi)) {
            dims_ <- c(n_plants, length(lines_), max_t + 1)
            date_ = rep(1:dims_[3] - 1L, each = prod(dims_[1:2]))
            plant_ = rep(1:dims_[1], prod(dims_[2:3]))
            line_ = rep(1:dims_[2], each = dims_[1]) %>% rep(dims_[3])
            simi[[j]] <- data_frame(plant = plant_,
                                 line = line_,
                                 date = date_,
                                 N = round(simi[[j]], 4),
                                 pool = paste(lines_, collapse = ""),
                                 rep = j) %>%
                arrange(rep, plant, line, date) %>%
                mutate(line = map_int(line, ~ lines_[.x]))
        }
        return(simi)
    }
})



# pb <- progress::progress_bar$new(
#     format = "  writing [:bar] :percent in :elapsed",
#     total = with(sim_env, length(pools) * n_cages),
#     clear = FALSE, width = options("width")$width)
#
# # Takes ~2 hrs
#
# set.seed(549489)
# for (i in 1:length(sim_env$pools)) {
#
#     simi <- sim_env$sim(i)
#
#     for (j in 1:length(simi)) {
#         readr::write_csv(simi[[j]], path = "/Volumes/750gb/__clonewars/pool_sims.csv",
#                          append = !(i == 1 & j == 1),
#                          col_names = (i == 1 & j == 1))
#         pb$tick()
#     }
#
# }; # rm(i, j, simi)


# ----------
# Old version:
# ----------
# # Takes ~6.5 min
# set.seed(549489)
# pool_sims <- map(1:length(sim_env$pools), sim_env$sim)
# readr::write_rds(pool_sims,
#                  path = "data-raw/pool_sims.rds",
#                  compress = "gz")
# # Takes ~1.75 min
# pool_sims <- readr::read_rds("data-raw/pool_sims.rds")
#
#
# ----------
# New version:
# ----------
# pb <- progress::progress_bar$new(
#     format = "  writing [:bar] :percent in :elapsed",
#     total = sum(map_int(pool_sims, length)),
#     clear = FALSE, width = options("width")$width)
#
# # n_plants, n_lines, n_times
# dims_ <- with(sim_env, c(n_plants, n_lines, max_t+1))
# for (i in 1:length(pool_sims)) {
#
#     # Should be same dims for all in this pool:
#     pool_ <- paste(sim_env$pools[[i]], collapse = "")
#     dims_[2] <- nchar(pool_)
#
#     for (j in 1:length(pool_sims[[i]])) {
#         date_ = rep(1:dims_[3] - 1L, each = prod(dims_[1:2]))
#         plant_ = rep(1:dims_[1], prod(dims_[2:3]))
#         line_ = rep(1:dims_[2], each = dims_[1]) %>% rep(dims_[3])
#         out_df <- data_frame(plant = plant_,
#                              line = line_,
#                              date = date_,
#                              N = round(as.numeric(pool_sims[[i]][[j]]), 4),
#                              pool = pool_,
#                              rep = j) %>%
#             arrange(rep, plant, line, date) %>%
#             mutate(line = map_int(line, ~ sim_env$pools[[i]][.x]))
#
#         readr::write_csv(out_df, path = "data-raw/pool_sims.csv", append = TRUE,
#                          col_names = (i == 1 & j == 1))
#
#         pb$tick()
#     }
# }; # rm(dims_, pool_, out_df, date_, plant_, line_, i, j)

# plant,line,date,N,pool,rep
# Takes ~ 2 hrs
# t0 <- Sys.time()
# pool_sims <- read.big.matrix("/Volumes/750gb/__clonewars/pool_sims.csv",
#     header = TRUE, type = "double", shared = TRUE,
#     backingfile = "sims.bin",
#     backingpath = "/Volumes/750gb/__clonewars",
#     descriptor = "sims.desc")
# Sys.time() - t0; rm(t0)

# Takes ~0.008757114 sec
desc <- dget("/Volumes/750gb/__clonewars/sims.desc")
pool_sims <- attach.big.matrix(desc)
head(pool_sims)
rm(pool_sims)



# pool_sims <- pool_sims %>%
#     group_by(pool, rep, line) %>%
#     summarize(N = mean(N)) %>%
#     ungroup()
# readr::write_rds(pool_sims, path = "data-raw/pool_sim_mean_no_error.rds")


# pool_sims <- readr::read_rds("data-raw/pool_sim_mean.rds")
# pool_sims <- readr::read_rds("data-raw/pool_sim_mean_no_error.rds")


# pool_sims <- pool_sims %>%
#     mutate(line = factor(line,
#                          levels = seq_along(levels(growth$line)),
#                          labels = levels(growth$line)),
#            rep = factor(rep),
#            pool_size = nchar(pool),
#            pool = factor(pool, levels = map_chr(pools, ~ paste(.x, collapse = "")))) %>%
#     identity()

# {
#     # This should have zero rows bc having total extinction makes no sense:
#     pool_sims %>%
#         filter(date == pool_max_t) %>%
#         group_by(pool, rep) %>%
#         summarize(N = sum(N)) %>%
#         ungroup() %>%
#         filter(N == 0) %>%
#         nrow() %>%
#         `==`(0) %>%
#         print()
#     # This should also have zero rows bc it means aphids spontaneously appeared
#     pool_sims %>%
#         group_by(pool, rep, line) %>%
#         filter(N == 0, dplyr::lead(N) > 0) %>%
#         ungroup() %>%
#         nrow() %>%
#         `==`(0) %>%
#         print()
# }


