
# ======================================================================================
# ======================================================================================

# This file does only the simulations and no analyses of those simulations.

# ======================================================================================
# ======================================================================================


suppressPackageStartupMessages({
    library(clonewars)
})
source(".Rprofile")





sim_env <- new.env()


with(sim_env, {

    stan_fit <- read_rds("data-raw/stan_fit.rds")

    n_plants <- 8
    n_lines <- 8
    N_0 <- matrix(rep(6, n_lines * n_plants), n_plants, n_lines)
    max_t <- 180
    R <- apply(rstan::extract(stan_fit, "R", permuted = FALSE), 3, mean) %>%
        as.numeric() %>%
        .[1] %>%
        rep(n_lines)
    A <- apply(rstan::extract(stan_fit, "A", permuted = FALSE), 3, mean) %>%
        as.numeric() %>%
        .[1] %>%
        rep(n_lines)
    D_binom <- clonewars::disp_estimates$binom %>%
        .[1,] %>%
        list() %>%
        rep(n_lines) %>%
        bind_rows()
    D_nb <- clonewars::disp_estimates$nb %>%
        .[1,] %>%
        list() %>%
        rep(n_lines) %>%
        bind_rows()
    process_error <- apply(rstan::extract(stan_fit, "s_epsilon", permuted = FALSE),
                           3, mean) %>%
        as.numeric() %>%
        .[1] %>%
        rep(n_lines)
    plant_mort_0 <- clonewars::plant_death$after_max_mort_coefs$inter %>%
        .[1] %>%
        rep(n_lines)
    plant_mort_1 <- clonewars::plant_death$after_max_mort_coefs$date %>%
        .[1] %>%
        rep(n_lines)
    plant_death_age_mean <- clonewars::plant_death$until_max_summ$max_mean
    plant_death_age_sd <- clonewars::plant_death$until_max_summ$max_sd
    repl_times <- seq(4, max_t, 4) - 1
    repl_age <- 3
    extinct_N <- 1
    n_cages <- 1000
    n_cores <- parallel::detectCores() - 1

    # All combinations of pools from 1 to 8
    pools <- map(1:n_lines, ~ combn(1:n_lines, .x) %>%
                     t() %>%
                     split(1:nrow(.)) %>%
                     set_names(NULL)) %>%
        flatten()

    # Running longer to try to find patterns
    max_t <- 2000
    n_cages <- 10  # only 10 bc this is just a test
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
        simi <- cwsims:::sim_cages(n_cages = n_cages,
                                    N_0 = N_0[,lines_, drop = FALSE],
                                    max_t = max_t,
                                    R = R[lines_],
                                    A = A[lines_],
                                    D_binom = D_binom[lines_,, drop = FALSE],
                                    D_nb = D_nb[lines_,, drop = FALSE],
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

})



library(progress)

pb <- progress_bar$new(
    format = "  simulating [:bar] :percent in :elapsed",
    total = with(sim_env, sum(sapply(pools, length))),
    clear = FALSE, width = options("width")$width)


# Takes ~1 min
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



pool_sims_N <- bind_rows(pool_sims_N) %>%
    mutate(pool_size = as.integer(floor(log10(pool)) + 1L))
invisible(gc())
pool_sims_X <- bind_rows(pool_sims_X) %>%
    mutate(pool_size = as.integer(floor(log10(pool)) + 1L))
invisible(gc())
pool_sims_Z <- bind_rows(pool_sims_Z) %>%
    mutate(pool_size = as.integer(floor(log10(pool)) + 1L))
invisible(gc())

pool_sims <- pool_sims_N %>%
    group_by(pool_size, pool, rep, line) %>%
    summarize(N = mean(N)) %>%
    ungroup()


pool_sims_Z %>%
    group_by(pool_size, pool, rep) %>%
    summarize(Z = mean(Z)) %>%
    ungroup() %>%
    ggplot(aes(pool_size, Z)) +
    geom_point(alpha = 0.25, shape = 1, aes(color = factor(pool_size)),
               position = position_jitter(width = 0.25, height = 0)) +
    stat_summary(fun.y = mean, geom = "point", size = 5) +
    scale_color_brewer(palette = "RdYlBu", guide = FALSE) +
    ylab("Mean proportion of empty time") +
    xlab("Pool size")


pool_sims_N %>%
    filter(pool == 1234, date > 1500) %>%
    ggplot(aes(date, N)) +
    geom_line(aes(color = factor(line), group = factor(rep)), alpha = 0.5) +
    facet_wrap(~ factor(line))



library(Perc)
library(stringr)
# https://cran.r-project.org/web/packages/Perc/vignettes/Perc.html



# All pairwise combinations:
pw_combs <- combn(1:sim_env$n_lines, 2) %>%
    t() %>%
    split(1:nrow(.)) %>%
    set_names(NULL)

conf_list <- map(2:sim_env$n_lines, function(ss) {

    df_ <- pool_sims %>%
        filter(pool_size == ss) %>%
        dplyr::select(-pool_size) %>%
        mutate(line = as.integer(line) %>% paste(),
               pool = pool %>% paste())

    map_dfr(pw_combs, function(cc) {
        df_ %>%
            filter(str_detect(pool, sprintf("(?=.*%s)(?=.*%s)", cc[1], cc[2])),
                   line %in% cc) %>%
            group_by(pool, rep) %>%
            summarize(win = line[N == max(N)],
                      lose = line[N == min(N)]) %>%
            ungroup() %>%
            summarize(winner = list(c(win[1], lose[1]) %>% sort()),
                      loser = list(c(lose[1], win[1]) %>% sort() %>% rev()),
                      freq = list(c(sum(win == winner[[1]][1]),
                                    sum(win == winner[[1]][2])))) %>%
            unnest()
    }) %>%
        mutate_all(list(as.numeric)) %>%
        as.matrix()

})

conf_mats <- map(conf_list, as.conflictmat, weighted = TRUE)

# conf_trans <- transitivity(conf_mat)

# maxLength > 2 doesn't change much for any of them:
dom_probs <- map(conf_mats, conductance, maxLength = 2)


# Takes ~3 min
set.seed(6328321)
s_ranks <- map(dom_probs, ~ simRankOrder(.x$p.hat))



rank_df <- map_dfr(1:length(s_ranks),
                   function(j) {
                       data_frame(line = paste(s_ranks[[j]]$BestSimulatedRankOrder$ID) %>%
                                      as.integer() %>%
                                      map_chr(~ letters[.x]),
                                  ranking = 1:length(line),
                                  pool_size = j+1)
                   }) %>%
    mutate(line = factor(line),
           ranking = factor(ranking, levels = sim_env$n_lines:1))

rank_df %>%
    # ggplot(aes(pool_size, ranking)) +
    ggplot(aes(pool_size, as.integer(ranking))) +
    # geom_raster(aes(fill = line)) +
    geom_line(aes(color = line), size = 1) +
    geom_point(aes(color = line), size = 3) +
    scale_fill_brewer("aphid\nline:", palette = "RdYlBu") +
    scale_color_brewer("aphid\nline:", palette = "RdYlBu") +
    scale_x_continuous("Pool size", breaks = 2:sim_env$n_lines) +
    scale_y_continuous("Ranking", breaks = 1:sim_env$n_lines, labels = sim_env$n_lines:1) +
    NULL


pool_sims_Z %>%
    group_by(pool_size, pool, rep) %>%
    summarize(Z = mean(Z)) %>%
    ungroup() %>%
    ggplot(aes(pool_size, Z)) +
    geom_point(alpha = 0.25, shape = 1, aes(color = factor(pool_size)),
               position = position_jitter(width = 0.25, height = 0)) +
    stat_summary(fun.y = mean, geom = "point", size = 5) +
    scale_color_brewer(palette = "RdYlBu", guide = FALSE) +
    ylab("Mean proportion of empty time") +
    xlab("Pool size")


pool_sims_X %>%
    group_by(pool_size, pool, rep) %>%
    summarize(X = mean(X)) %>%
    ungroup() %>%
    ggplot(aes(pool_size, X)) +
    geom_point(alpha = 0.25, shape = 1, aes(color = factor(pool_size)),
               position = position_jitter(width = 0.25, height = 0)) +
    stat_summary(fun.y = mean, geom = "point", size = 5) +
    scale_color_brewer(palette = "RdYlBu", guide = FALSE) +
    ylab("Mean log1p(# aphids)") +
    xlab("Pool size")




#


