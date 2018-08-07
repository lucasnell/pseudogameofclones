
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
    D_mat[,4] <- exp(D_mat[,4])

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
        dims_ <- c(n_plants, length(lines_), max_t + 1)
        for (j in 1:length(simi)) {
            date_ = rep(1:dims_[3] - 1L, each = prod(dims_[1:2]))
            plant_ = rep(1:dims_[1], prod(dims_[2:3]))
            line_ = rep(1:dims_[2], each = dims_[1]) %>% rep(dims_[3])
            simi[[j]] <- data_frame(plant = plant_,
                                 line = line_,
                                 date = date_,
                                 N = round(simi[[j]], 4),
                                 pool = paste(lines_, collapse = ""),
                                 rep = j) %>%
                arrange(rep, line, plant, date) %>%
                mutate(line = map_int(line, ~ lines_[.x]))
        }
        return(simi)
    }
}); rm(stan_fit)


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
#         readr::write_csv(simi[[j]], path = "data-raw/big/pool_sims.csv",
#                          append = !(i == 1 & j == 1),
#                          col_names = (i == 1 & j == 1))
#         pb$tick()
#     }
#
# }; rm(i, j, simi)
#
# gc()
#
# # plant,line,date,N,pool,rep
# # Takes ~ 2 hrs
# cat(sprintf("Started at %s\n", Sys.time()))
# pool_sims <- read.big.matrix("data-raw/big/pool_sims.csv",
#     header = TRUE, type = "double", shared = TRUE,
#     backingfile = "sims.bin",
#     backingpath = "data-raw/big",
#     descriptor = "sims.desc")




# # Takes ~0.008757114 sec
# desc <- dget("data-raw/big/sims.desc")
# pool_sims <- attach.big.matrix(desc)
# head(pool_sims, 20)
# nrow(pool_sims)
#
#
# # Takes ~ 10 min
# group_tree <- cwsims:::make_group_tree(pool_sims@address,
#                 pool_sizes = sapply(sim_env$pools, length),
#                 n_reps = sim_env$n_cages,
#                 n_plants = sim_env$n_plants,
#                 max_t = sim_env$max_t)
#
# # Takes ~0.25 sec
# cwsims:::save_group_tree(group_tree, filename = "data-raw/big/sims_group_tree.cereal")
#
# # Takes ~0.5 sec
# group_tree <- cwsims:::load_group_tree(filename = "data-raw/big/sims_group_tree.cereal")
#
#
# # Each below takes ~3.5 min
# by_plant_mean <- cwsims:::grouped_mean(pool_sims@address,
#                                    group_tree,
#                                    by_plant = TRUE,
#                                    by_date = FALSE)
# readr::write_rds(by_plant_mean, "data-raw/big/by_plant_mean.rds")
#
# by_plant_zeros <- cwsims:::grouped_mean(pool_sims@address,
#                                     group_tree,
#                                     by_plant = TRUE,
#                                     by_date = FALSE,
#                                     zeros = TRUE)
# readr::write_rds(by_plant_zeros, "data-raw/big/by_plant_zeros.rds")
#
# by_cage_mean <- cwsims:::grouped_mean(pool_sims@address,
#                                   group_tree,
#                                   by_plant = FALSE,
#                                   by_date = TRUE)
# readr::write_rds(by_cage_mean, "data-raw/big/by_cage_mean.rds")
#
# by_cage_zeros <- cwsims:::grouped_mean(pool_sims@address,
#                                    group_tree,
#                                    by_plant = FALSE,
#                                    by_date = TRUE,
#                                    zeros = TRUE)
# readr::write_rds(by_cage_zeros, "data-raw/big/by_cage_zeros.rds")


#
# sims_bydate <- readr::read_rds("data-raw/big/by_cage_mean.rds")
# sims_bydate <- as_data_frame(sims_bydate)
# colnames(sims_bydate) <- c("pool", "rep", "line", "date", "N")
#
# for (i in c("pool", "rep", "line", "date")) {
#     sims_bydate[[i]] <- as.integer(sims_bydate[[i]])
#     print(i)
# }; rm(i)
#
# zeros <- readr::read_rds("data-raw/big/by_cage_zeros.rds")
# zeros <- zeros[,5]
#
# sims_bydate$Z <- zeros
# rm(zeros); invisible(gc())
#
# sims_bydate <- mutate(sims_bydate, line = factor(line, levels = 1:8,
#                                                  labels = clonewars::load_data() %>%
#                                                      .[["line"]] %>%
#                                                      levels()))
#
# readr::write_rds(sims_bydate, path = "data-raw/big/sims_bydate.rds")
# rm(sims_bydate); invisible(gc())
#
# sims_byplant <- readr::read_rds("data-raw/big/by_plant_mean.rds") %>%
#     as_data_frame() %>%
#     set_names(c("pool", "rep", "line", "plant", "N")) %>%
#     mutate(Z = {readr::read_rds("data-raw/big/by_plant_zeros.rds")}[,5],
#            line = factor(line, levels = 1:8,
#                          labels = clonewars::load_data() %>%
#                              .[["line"]] %>%
#                              levels())) %>%
#     mutate_at(vars(pool, rep, plant), funs(as.integer))
#
# readr::write_rds(sims_byplant, path = "data-raw/big/sims_byplant.rds")

sims_bydate <- readr::read_rds("data-raw/big/sims_bydate.rds")
sims_byplant <- readr::read_rds("data-raw/big/sims_byplant.rds")


# mean total (i.e., among all lines) N on plants
# mean time a plant is empty

