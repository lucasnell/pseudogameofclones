# !diagnostics off


source("under_constr/poster/_preamble.R")


library(ape)
library(stringr)


ggplot2::theme_set(
    ggplot2::theme_get() +
        theme(axis.ticks = element_line(),
              axis.title = element_text(),
              axis.text = element_text(),
              strip.text = element_text())
)



growth <-
    load_data(impute_fxn = impute, filter_pars = NULL) %>%
    mutate(line = paste(line)) %>%
    bind_rows(clonewars:::load_pz_data(impute_fxn = impute, filter_pars = NULL)) %>%
    mutate_at(vars(line, rep), funs(factor)) %>%
    # Now filter out early part of each time series, before N > 6
    # N <= 6 is when the stochasticity associated with only starting with 2 adults
    # appears to be strongest
    group_by(line, rep) %>%
    filter(1:n() >= which(N > 6)[1]) %>%
    mutate(date = date - min(date),
           disp = ifelse(is.na(disp), 0, disp),
           N = round(N),
           X = log(N),
           pN = N - disp,  # numbers of aphids on plant
           dD = disp - lag(disp, default = 0), # number of new dispersed aphids
           dD = ifelse(dD < 0, 0, dD)) %>%
    ungroup()





stan_fit <- read_rds("data-raw/stan_fit.rds")

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

# lines_ <- 1:8

# All pairwise combinations:
pw_combs <- combn(1:n_lines, 2) %>%
    t() %>%
    split(1:nrow(.)) %>%
    set_names(NULL)
# Function to do simulations with pairwise combinations:
pw_sims <- function(line1_, line2_) {
    lines_ <- c(line1_, line2_)
    sim_df <- cwsims::sim_cages(n_cages, N_0[,lines_], max_t, R[lines_], A[lines_],
                                D_binom[lines_,], D_nb[lines_,], process_error,
                                plant_mort_0[lines_], plant_mort_1[lines_],
                                plant_death_age_mean, plant_death_age_sd,
                                repl_times, repl_age, extinct_N, n_cores,
                                by_cage = TRUE) %>%
        mutate(pair = sprintf("%i_%i", line1_, line2_), X = log(N),
               line = map_int(line, ~ lines_[.x]))
    return(sim_df)
}

# # Takes ~3 min
# set.seed(654651985L)
# pair_sims <- map_dfr(pw_combs, ~ pw_sims(.x[1], .x[2]))
# pair_sims <- pair_sims %>%
#     mutate(line = factor(line,
#                          levels = seq_along(levels(growth$line)),
#                          labels = levels(growth$line)),
#            rep = factor(rep),
#            pair = factor(pair)) %>%
#     arrange(pair, rep, date, line) %>%
#     # group_by(pair, rep, date) %>%
#     # mutate(prop = N / sum(N),
#     #        prop_end = cumsum(prop),
#     #        prop_start = lag(prop_end, default = 0)) %>%
#     # ungroup() %>%
#     identity()
# readr::write_rds(pair_sims, path = "data-raw/pair_sims.rds")

pair_sims <- readr::read_rds("data-raw/pair_sims.rds")


# {
#     # This should have zero rows bc having total extinction makes no sense:
#     pair_sims %>%
#         filter(date == max_t) %>%
#         group_by(pair, rep) %>%
#         summarize(N = sum(N)) %>%
#         ungroup() %>%
#         filter(N == 0) %>%
#         nrow() %>%
#         `==`(0) %>%
#         print()
#     # This should also have zero rows bc it means aphids spontaneously appeared
#     pair_sims %>%
#         group_by(pair, rep, line) %>%
#         filter(N == 0, dplyr::lead(N) > 0) %>%
#         ungroup() %>%
#         nrow() %>%
#         `==`(0) %>%
#         print()
# }


pair_ranks <- pair_sims %>%
    group_by(pair, line) %>%
    summarize(N = mean(N)) %>%
    group_by(pair) %>%
    mutate(opponent = rev(line), diff = (N - rev(N)) / N) %>%
    ungroup() %>%
    select(opponent, line, N, diff) %>%
    add_row(opponent = pair_sims$line %>% unique() %>% sort(),
            line = pair_sims$line %>% unique() %>% sort(),
            N = NA, diff = 0) %>%
    arrange(opponent, line) %>%
    group_by(opponent) %>%
    # Descending rank:
    mutate(rank = rank(-diff)) %>%
    ungroup()




# ================================================================================
# ================================================================================

# Pooled simulations

# ================================================================================
# ================================================================================

# All combinations of pools from 1 to 8
pools <- map(1:n_lines, ~ combn(1:n_lines, .x) %>%
        t() %>%
        split(1:nrow(.)) %>%
        set_names(NULL)) %>%
    flatten()

# Running longer to try to find patterns
pool_max_t <- 2000
pool_n_cages <- 100
pool_repl_times <- seq(4, pool_max_t, 4) - 1

# Function to do simulations with pools of varying combinations:
pool_sim_fun <- function(lines_) {
    sim_df <- cwsims::sim_cages(pool_n_cages, N_0[,lines_, drop = FALSE], pool_max_t,
                                R[lines_], A[lines_],
                                D_binom[lines_,, drop = FALSE],
                                D_nb[lines_,, drop = FALSE],
                                process_error = 0,
                                plant_mort_0[lines_], plant_mort_1[lines_],
                                plant_death_age_mean, plant_death_age_sd,
                                pool_repl_times, repl_age, extinct_N, n_cores,
                                by_cage = TRUE) %>%
        mutate(pool = paste(lines_, collapse = ""),
               X = log(N),
               line = map_int(line, ~ lines_[.x]))
    return(sim_df)
}

# # Takes ~10 min
# set.seed(549489)
# pool_sims <- map_dfr(pools, pool_sim_fun)
# readr::write_rds(pool_sims, path = "data-raw/pool_sims.rds")
# pool_sims <- pool_sims %>%
#     group_by(pool, rep, line) %>%
#     summarize(N = mean(N)) %>%
#     ungroup()
# readr::write_rds(pool_sims, path = "data-raw/pool_sim_mean.rds")


# set.seed(549489+1)
# pool_sims <- map_dfr(pools, pool_sim_fun)
# readr::write_rds(pool_sims, path = "data-raw/pool_sims_no_error.rds")
# pool_sims <- readr::read_rds("data-raw/pool_sims_no_error.rds")
# pool_sims <- pool_sims %>%
#     group_by(pool, rep, line) %>%
#     summarize(N = mean(N)) %>%
#     ungroup()
# readr::write_rds(pool_sims, path = "data-raw/pool_sim_mean_no_error.rds")


# pool_sims <- readr::read_rds("data-raw/pool_sim_mean.rds")
pool_sims <- readr::read_rds("data-raw/pool_sim_mean_no_error.rds")


pool_sims <- pool_sims %>%
    mutate(line = factor(line,
                         levels = seq_along(levels(growth$line)),
                         labels = levels(growth$line)),
           rep = factor(rep),
           pool_size = nchar(pool),
           pool = factor(pool, levels = map_chr(pools, ~ paste(.x, collapse = "")))) %>%
    identity()

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




pool_ranks <- pool_sims %>%
    filter(pool_size == 8) %>%
    group_by(pool, rep) %>%
    mutate(N = N / sum(N),
           rank = rank(-N)) %>%
    ungroup() %>%
    dplyr::select(-pool, -pool_size) %>%
    # group_by(line) %>%
    # summarize(N = mean(N)) %>%
    # ungroup() %>%
    # # Descending rank by mean N:
    # mutate(rank = rank(-N)) %>%
    identity()


pool_ranks %>%
    ggplot(aes(line, rank, color = line)) +
    geom_point(alpha = 0.5, shape = 1,
               position = position_jitter(width = 0.25, height = 0.1)) +
    scale_color_brewer(palette = "Dark2", guide = FALSE)




pair_ranks <- pool_sims %>%
    filter(pool_size == 2) %>%
    dplyr::select(-pool_size) %>%
    mutate(line = as.integer(line),
           opponent = str_remove_all(pool, paste(line)) %>% as.integer()) %>%
    group_by(pool, rep) %>%
    mutate(N = N / sum(N)) %>%
    ungroup() %>%
    group_by(opponent, line) %>%
    summarize(N = mean(N)) %>%
    ungroup() %>%
    select(opponent, line, N) %>%
    add_row(opponent = 1:n_lines, line = 1:n_lines, N = NA) %>%
    arrange(opponent, line) %>%
    # Now filling in values for opponents:
    group_by(opponent) %>%
    mutate(N = ifelse(is.na(N), 1 - mean(N, na.rm = TRUE), N)) %>%
    # Descending rank:
    mutate(rank = rank(-N)) %>%
    ungroup() %>%
    mutate(line = factor(line, levels = 1:n_lines, labels = levels(pool_ranks$line))) %>%
    identity()


# rank_plot <-
pair_ranks %>%
    ggplot(aes(rank, as.integer(opponent))) +
    geom_vline(xintercept = c(1, 8), linetype = 1, size = 0.5) +
    geom_hline(yintercept = 1:8, linetype = 3, size = 0.25) +
    geom_vline(data = pool_ranks %>%
                   group_by(line) %>%
                   summarize(N = mean(N)) %>%
                   ungroup() %>%
                   # Descending rank by mean N:
                   mutate(rank = rank(-N)),
                   aes(xintercept = rank),
               size = 2, linetype = 1, color = palette$default_primary) +
    geom_path(aes(color = line), size = 1, linetype = 1,
               color = palette$accent) +
    geom_point(aes(color = line), size = 3, shape = 16,
               color = palette$accent) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    scale_linetype_manual(values = c(1, 3), guide = FALSE) +
    facet_wrap(~ line, nrow = 2, scales = "free_x") +
    # reversed to make line #1 up top:
    scale_y_reverse("Opponent",
                    breaks = 1:n_lines,
                    labels = levels(pair_ranks$line)) +
    scale_x_continuous("Rank", breaks = 1:n_lines) +
    theme(axis.line.y.left = element_blank(), axis.ticks.y = element_blank(),
          # axis.title = element_blank(),
          # axis.text = element_blank(),
          # strip.text = element_blank(),
          panel.spacing.y = unit(2, "cm")) +
    NULL



# ggsave(filename = "figs/line_ranks.pdf", plot = rank_plot,
#        width = 20, height = 10, units = "cm", bg = "white", useDingbats = FALSE)



# meanN_plot <-
pool_sims %>%
    mutate(pool_size = factor(pool_size)) %>%
    filter(pool_size != 1) %>%
    group_by(pool_size, pool, line) %>%
    summarize(n = mean(N)) %>%
    ungroup() %>%
    ggplot(aes(as.integer(line), n)) +
    geom_point(data = pool_sims %>%
                   mutate(pool_size = factor(pool_size)) %>%
                   filter(pool_size == 1) %>%
                   group_by(pool_size, pool, line) %>%
                   summarize(n = mean(N)) %>%
                   ungroup(),
               aes(color = pool_size, group = pool), size = 3, alpha = 0.25) +
    geom_line(aes(color = pool_size, group = pool), size = 1, alpha = 0.25) +
    # geom_point(aes(color = factor(pool_size)), size = 3, shape = 1,
    #            position = position_jitter(height = 0, width = 0.25)) +
    # stat_summary(aes(color = factor(pool_size)), geom = "line",
    #              # size = 6, shape = 16, alpha = 0.75,
    #              size = 2,
    #              fun.y = mean) +
    scale_y_continuous("mean(N)", breaks = 1:n_lines, labels = levels(pool_sims$line)) +
    xlab("Aphid line") +
    scale_color_brewer(palette = "Dark2") +
    facet_wrap(~ pool_size, scales = "free_y", nrow = 2) +
    theme(# axis.ticks.y = element_blank(),
          # axis.title = element_blank(),
          # axis.text = element_blank(),
          # strip.text = element_blank(),
          NULL
          ) +
    NULL





# ggsave(filename = "figs/line_survs.pdf", plot = surv_plot,
#        width = 20, height = 8, units = "cm", bg = "white", useDingbats = FALSE)





# ================================================================================
# ================================================================================

# Cluster analysis

# ================================================================================
# ================================================================================




pool_mats <- pool_sims %>%
    arrange(pool_size) %>%
    split(.$pool_size) %>%
    map(function(df_) {

        pool_size_ <- df_$pool_size[1]

        df__ <- df_ %>%
            dplyr::select(-pool_size) %>%
            mutate(line = as.integer(line) %>% paste(),
                   pool = paste(pool),
                   others = str_remove_all(pool, line)) %>%
            arrange(pool, rep, line)

        others_ <- combn(n_lines, pool_size_ - 1) %>%
            t() %>%
            split(1:nrow(.)) %>%
            set_names(NULL) %>%
            map(paste)

        if (pool_size_ > 1) {
            na_rows__ <- expand.grid(N = NA_real_,
                                     rep = factor(1, levels = levels(df__$rep)),
                                     line = df__$line %>% unique() %>% sort(),
                                     pool = sort(unique(df__$pool)),
                                     stringsAsFactors = FALSE) %>%
                tbl_df() %>%
                mutate(others = map2(pool, line,
                                     function(pp, ll) {
                                         if (str_detect(pp, ll)) {
                                             lgl <- map_lgl(others_,
                                                            ~ any(str_detect(ll, .x)))
                                         } else {
                                             lgl <- map_lgl(others_,
                                                            ~ all(str_detect(pp, .x)))
                                         }
                                         strs <- others_[lgl]
                                         map_chr(strs, ~ str_c(.x, collapse = ""))
                                     })) %>%
                unnest()

            df__ <- bind_rows(df__, na_rows__)
        }

        M <- df__ %>%
            arrange(line, others, rep) %>%
            group_by(line, others) %>%
            summarize(N = mean(N, na.rm = TRUE)) %>%
            ungroup() %>%
            mutate(N = ifelse(is.nan(N), NA, N)) %>%
            spread(others, N) %>%
            .[, -1] %>%
            as.matrix()

        if (pool_size_ == df_$line %>% levels() %>% length()) M <- cbind(M[!is.na(M)])

        rownames(M) <- letters[as.integer(unique(df_$line))]

        return(M)
    })


# euclidean distances:
pool_ds <- map(pool_mats, dist)

# complete linkage method for hierarchical clustering:
pool_hc <- map(pool_ds, hclust)

# Converting to a multiPhylo object:
pool_phy <- map(pool_hc, ape::as.phylo) %>%
    do.call(what = c)

par(mfrow = c(2, 4))
for (i in 1:n_lines) {
    # pdf(sprintf("figs/cluster_%02i.pdf", i), width = 16 / 2.54, height = 16 / 2.54,
    #     useDingbats = FALSE)
    pool_phy[[i]] %>%
        plot(
            type = "unrooted",
            main = sprintf("pool size: %i", i),
            cex = 1.5,
            rotate.tree = 0
        )
    # dev.off()
}


library(Perc)
# https://cran.r-project.org/web/packages/Perc/vignettes/Perc.html


conf_list <- map(2:n_lines, function(ss) {

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
        mutate_all(funs(as.numeric)) %>%
        as.matrix()

})


conf_mats <- map(conf_list, as.conflictmat, weighted = TRUE)

# conf_trans <- transitivity(conf_mat)

# maxLength > 2 doesn't change much for any of them:
dom_probs <- map(conf_mats, conductance, maxLength = 2)


# Takes ~3 min
s_ranks <- map(dom_probs, ~ simRankOrder(.x$p.hat))

map_dfr(2:n_lines,
        ~ data_frame(id = paste(s_ranks[[.x-1]]$BestSimulatedRankOrder$ID),
                     ranking = 1:length(id),
                     pool_size = .x)) %>%
    # spread(key = pool_size, value = id) %>%
    mutate_at(vars(pool_size, id), factor) %>%
    mutate(ranking = factor(ranking, levels = n_lines:1)) %>%
    # ggplot(aes(pool_size, ranking)) +
    ggplot(aes(as.integer(pool_size), as.integer(ranking))) +
    # geom_raster(aes(fill = id)) +
    geom_line(aes(color = id), size = 1) +
    geom_point(aes(color = id), size = 3) +
    scale_fill_brewer(palette = "RdYlBu") +
    scale_color_brewer(palette = "RdYlBu") +
    NULL


i <- 5

# s_ranks[[i]]$BestSimulatedRankOrder
# dom_probs[[i]]$p.hat
# valueConverter(dom_probs[[i]]$p.hat)

dyadicLongConverter(dom_probs[[i]]$p.hat) %>%
    mutate_at(vars(ID1, ID2), paste) %>%
    bind_rows(dyadicLongConverter(dom_probs[[i]]$p.hat) %>%
                  mutate_at(vars(ID1, ID2), paste) %>%
                  mutate(ID1_ = ID2, ID2 = ID1, ID1 = ID1_) %>%
                  dplyr::select(-ID1_)) %>%
    distinct(ID1, ID2, .keep_all = TRUE) %>%
    dplyr::rename(certainty = RankingCertainty) %>%
    .[, c("ID1", "ID2", "certainty")] %>%
    tbl_df() %>%
    mutate(ID1 = factor(ID1, levels = paste(s_ranks[[i]]$BestSimulatedRankOrder$ID)),
           ID2 = factor(ID2, levels = rev(paste(s_ranks[[i]]$BestSimulatedRankOrder$ID)))) %>%
    ggplot(aes(ID1, ID2), drop = FALSE) +
    geom_raster(aes(fill = certainty)) +
    scale_fill_gradient(low = palette$dark_primary, high = palette$light_primary,
                        breaks = seq(0.5, 0.9, 0.1), limits = c(0.5, 1.0)) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    NULL
