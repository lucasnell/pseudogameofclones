# !diagnostics off


source("under_constr/poster/_preamble.R")
source("under_constr/poster/simulations.R")

suppressPackageStartupMessages({
    library(ape)
    library(stringr)
})


ggplot2::theme_set(
    ggplot2::theme_get() +
        theme(axis.ticks = element_line(),
              axis.title = element_text(),
              axis.text = element_text(),
              strip.text = element_text())
)

lines_ <- sim_env$pools[[100]]

set.seed(9)
simi1 <- with(sim_env, {
    cwsims::sim_cages(n_cages = 10,
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
                                  by_cage = TRUE,
                                  show_progress = FALSE) %>%
            mutate(pool = paste(lines_, collapse = "") %>% as.integer()) %>%
            arrange(rep, line, date) %>%
            mutate(line = map_int(line, ~ lines_[.x])) %>%
            dplyr::select(pool, rep, everything())
})

set.seed(9)
simi4 <- with(sim_env, {
    cwsims::sim_cages(n_cages = 10,
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
                                  by_cage = FALSE,
                                  show_progress = FALSE) %>%
            mutate(pool = paste(lines_, collapse = "") %>% as.integer()) %>%
            arrange(rep, line, date) %>%
            mutate(line = map_int(line, ~ lines_[.x])) %>%
            dplyr::select(pool, rep, everything())
})

#     sim_small2 <- function() {
#     map_dfr(pools, function(lines_) {

set.seed(9)
simi <- with(sim_env, {
        cwsims:::sim_cages_(n_cages = 10,
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
})
set.seed(9)
simi3 <- with(sim_env, {
        cwsims:::sim_cages_(n_cages = 10,
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
                                    by_cage = TRUE,
                                    show_progress = FALSE)
})

dims_ <- with(sim_env, c(n_plants, length(lines_), max_t + 1))

for (j in 1:length(simi)) {
    date_ = rep(1:dims_[3] - 1L, each = prod(dims_[1:2]))
    plant_ = rep(1:dims_[1], prod(dims_[2:3]))
    line_ = rep(1:dims_[2], each = dims_[1]) %>% rep(dims_[3])
    simi[[j]] <- data_frame(plant = plant_,
                            line = line_,
                            date = date_,
                            N = simi[[j]],
                            pool = paste(lines_, collapse = "") %>% as.integer(),
                            rep = j) %>%
        arrange(rep, line, plant, date) %>%
        mutate(line = map_int(line, ~ lines_[.x])) %>%
        group_by(pool, rep, line, date) %>%
        summarize(N = sum(N)) %>%
        ungroup()
}

for (j in 1:length(simi3)) {
    simi3[[j]] <- cwsims:::simplify_cage(simi3[[j]], j, with(sim_env, N_0[,lines_, drop = FALSE]),
                  sim_env$max_t, by_cage_ = TRUE) %>%
        mutate(pool = paste(lines_, collapse = "") %>% as.integer(),
               line = map_int(line, ~ lines_[.x])) %>%
        arrange(pool, rep, line, date) %>%
        dplyr::select(pool, rep, line, date, N)
}

simi <- bind_rows(simi)
simi3 <- bind_rows(simi3)

simi1; simi
all.equal(simi1, simi)

simi4 %>%
    group_by(rep, line) %>%
    summarize(N = mean(N)) %>%
    group_by(rep) %>%
    mutate(rank = rank(-N)) %>%
    ungroup() %>%
    ggplot(aes(line, rank)) +
    geom_point(alpha = 0.5, shape = 1, position = position_jitter(width = 0.25))




#         return(bind_rows(simi))
#     })}}
# )



set.seed(549489)
test_df <- sim_env$sim_small()
test_df <- test_df %>%
    mutate(pool = as.integer(pool), pool_size = as.integer(floor(log10(pool)) + 1L)) %>%
    dplyr::select(pool_size, pool, rep, line, date, N) %>%
    arrange(pool_size, pool, rep, line, date)
set.seed(549489)
test_df2 <- sim_env$sim_small2()
test_df2 <- test_df2 %>%
    mutate(pool = as.integer(pool), pool_size = as.integer(floor(log10(pool)) + 1L)) %>%
    dplyr::select(pool_size, pool, rep, line, date, N) %>%
    arrange(pool_size, pool, rep, line, date)

# all.equal(test_df$N, test_df2$N)



set.seed(549489)
for (i in 1:length(sim_env$pools)) {

    simi <- sim_env$sim(i) %>%
        map(~ )

    for (j in 1:length(simi)) {
        readr::write_csv(simi[[j]], path = "data-raw/big/pool_sims.csv",
                         append = !(i == 1 & j == 1),
                         col_names = (i == 1 & j == 1))
        pb$tick()
    }

}; rm(i, j, simi)

sims_bydate <- sims_bydate %>%
    mutate(pool_size = as.integer(floor(log10(pool)) + 1L))
sims_byplant <- sims_byplant %>%
    mutate(pool_size = as.integer(floor(log10(pool)) + 1L))

n_lines <- sims_byplant$line %>% unique() %>% length()
n_plants <- sims_byplant$plant %>% unique() %>% length()

sims_byline <- sims_bydate %>%
    group_by(pool, rep, line, pool_size) %>%
    summarize(N = mean(N), Z = mean(Z)) %>%
    ungroup()
# sims_byline <- sims_byplant %>%
#     group_by(pool, rep, line, pool_size) %>%
#     summarize(N = mean(N), Z = mean(Z)) %>%
#     ungroup()

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





pair_ranks <- sims_byline %>%
    filter(pool_size == 2) %>%
    group_by(pool, line) %>%
    summarize(N = mean(N)) %>%
    group_by(pool) %>%
    mutate(opponent = rev(line), diff = (N - rev(N)) / N) %>%
    ungroup() %>%
    select(opponent, line, N, diff) %>%
    add_row(opponent = growth$line %>% unique() %>% sort(),
            line = growth$line %>% unique() %>% sort(),
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



pool_ranks <- sims_byline %>%
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




pair_ranks <- sims_byline %>%
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
sims_byline %>%
    mutate(pool_size = factor(pool_size)) %>%
    filter(pool_size != 1) %>%
    group_by(pool_size, pool, line) %>%
    summarize(n = mean(N)) %>%
    ungroup() %>%
    ggplot(aes(as.integer(line), n)) +
    geom_point(data = sims_byline %>%
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
    scale_y_continuous("mean(N)", breaks = 1:n_lines, labels = levels(sims_byline$line)) +
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




pool_mats <- sims_byline %>%
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







# ================================================================================
# ================================================================================

# Dominance analysis

# ================================================================================
# ================================================================================


library(Perc)
# https://cran.r-project.org/web/packages/Perc/vignettes/Perc.html



# All pairwise combinations:
pw_combs <- combn(1:n_lines, 2) %>%
    t() %>%
    split(1:nrow(.)) %>%
    set_names(NULL)

sims_byline <- test_df2 %>%
    group_by(pool_size, pool, rep, line) %>%
    summarize(N = mean(N)) %>%
    ungroup() %>%
    mutate(line = factor(line, levels = 1:n_lines,
                         labels = clonewars::load_data() %>% .[["line"]] %>% levels()))

# ss = 3
# cc = pw_combs[[1]]
# rm(ss, cc, df_)

conf_list <- map(2:n_lines, function(ss) {

    df_ <- sims_byline %>%
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
set.seed(6328321)
s_ranks <- map(dom_probs, ~ simRankOrder(.x$p.hat))
# readr::write_rds(s_ranks, "data-raw/s_ranks.rds")
# s_ranks <- readr::read_rds("data-raw/s_ranks.rds")



rank_df <- map_dfr(1:length(s_ranks),
                   function(j) {
                       data_frame(line = paste(s_ranks[[j]]$BestSimulatedRankOrder$ID) %>%
                                      as.integer() %>%
                                      map_chr(~ letters[.x]),
                                  ranking = 1:length(line),
                                  pool_size = j+1)
                   }) %>%
    mutate(line = factor(line),
           ranking = factor(ranking, levels = n_lines:1),
           surv = map_dbl(line,
                          function(i) {
                              i_ <- levels(sims_byline$line)[as.integer(i)]
                              ind <- disp_estimates$binom$line == i_
                              stopifnot(sum(ind) == 1)
                              pd <- 2
                              gtools::inv.logit(
                                  plant_death$after_max_mort_coefs$inter[ind] +
                                      plant_death$after_max_mort_coefs$date[ind] * pd)
                          }),
           disp = map_dbl(line,
                          function(i) {
                              i_ <- levels(sims_byline$line)[as.integer(i)]
                              ind <- disp_estimates$binom$line == i_
                              stopifnot(sum(ind) == 1)
                              xmat <- as.matrix(disp_estimates$binom[ind, -1])
                              # Days past plant death:
                              pd <- 0
                              gtools::inv.logit(xmat %*% rbind(1, pd, pd^2)) *
                                  exp(disp_estimates$nb$b0[ind]) * 500
                          }),
           R = map_dbl(line, ~ sim_env$R[as.integer(.x)]),
           A = map_dbl(line, ~ sim_env$A[as.integer(.x)]))


rank_df %>%
    # ggplot(aes(pool_size, ranking)) +
    ggplot(aes(pool_size, as.integer(ranking))) +
    # geom_raster(aes(fill = line)) +
    geom_line(aes(color = line), size = 1) +
    geom_point(aes(color = line), size = 3) +
    scale_fill_brewer("aphid\nline:", palette = "RdYlBu") +
    scale_color_brewer("aphid\nline:", palette = "RdYlBu") +
    scale_x_continuous("Pool size", breaks = 2:n_lines) +
    scale_y_continuous("Ranking", breaks = 1:n_lines, labels = n_lines:1) +
    NULL




rank_df %>%
    mutate(ranking = ranking %>% paste() %>% as.integer()) %>%
    group_by(pool_size) %>%
    summarize(corr = cor(ranking, surv))
rank_df %>%
    mutate(pool_size = factor(pool_size)) %>%
    ggplot(aes(surv, ranking, color = line)) +
    geom_point(size = 3) +
    scale_color_brewer(palette = "RdYlBu") +
    facet_wrap(~pool_size, nrow = 3, scales = "free") +
    ylab("Ranking") +
    xlab("Expected survival at 2 days past plant death") +
    NULL


rank_df %>%
    mutate(ranking = ranking %>% paste() %>% as.integer()) %>%
    group_by(pool_size) %>%
    summarize(corr = cor(ranking, disp))
rank_df %>%
    mutate(pool_size = factor(pool_size)) %>%
    ggplot(aes(disp, ranking, color = line)) +
    geom_point(size = 3) +
    scale_color_brewer(palette = "RdYlBu") +
    facet_wrap(~pool_size, nrow = 3, scales = "free") +
    ylab("Ranking") +
    xlab("Expected dispersal at 0 days past plant death and N = 500") +
    NULL




# rank_df %>%
#     mutate(ranking = ranking %>% paste() %>% as.integer()) %>%
#     group_by(pool_size) %>%
#     summarize(corr = cor(ranking, R))
# rank_df %>%
#     mutate(pool_size = factor(pool_size)) %>%
#     ggplot(aes(R, ranking, color = line)) +
#     geom_point(size = 3) +
#     scale_color_brewer(palette = "RdYlBu") +
#     facet_wrap(~pool_size, nrow = 3, scales = "free") +
#     ylab("Ranking") +
#     xlab("Growth rate") +
#     NULL


rank_df %>%
    mutate(ranking = ranking %>% paste() %>% as.integer()) %>%
    group_by(pool_size) %>%
    summarize(corr = cor(ranking, A))
rank_df %>%
    mutate(pool_size = factor(pool_size)) %>%
    ggplot(aes(A, ranking, color = line)) +
    geom_point(size = 3) +
    scale_color_brewer(palette = "RdYlBu") +
    facet_wrap(~pool_size, nrow = 3, scales = "free") +
    ylab("Ranking") +
    scale_x_continuous("Density dependence", breaks = seq(0.0019, 0.0021, 0.0001)) +
    NULL




with(rank_df %>% distinct(line, A, disp), cor(A, disp))
rank_df %>%
    distinct(line, A, disp) %>%
    ggplot(aes(A, disp, color = line)) +
    geom_point(size = 3) +
    scale_color_brewer(palette = "RdYlBu") +
    ylab("Expected dispersal at 0 days past plant death and N = 500") +
    scale_x_continuous("Density dependence", breaks = seq(0.0019, 0.0021, 0.0001)) +
    NULL








# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================


# s_ranks[[i]]$BestSimulatedRankOrder
# dom_probs[[i]]$p.hat
# valueConverter(dom_probs[[i]]$p.hat)

dom_dfs <- map(1:length(dom_probs),
              function(i) {
                  order_ <- s_ranks[[i]]$BestSimulatedRankOrder$ID %>%
                      paste() %>%
                      as.integer() %>%
                      map_chr(~ letters[.x])
                  df_ <- dyadicLongConverter(dom_probs[[i]]$p.hat) %>%
                      mutate_at(vars(ID1, ID2), paste) %>%
                      mutate_at(vars(ID1, ID2), as.integer) %>%
                      mutate_at(vars(ID1, ID2),
                                function(y) map_chr(y, ~ letters[.x])) %>%
                      dplyr::rename(certainty = RankingCertainty) %>%
                      .[, c("ID1", "ID2", "certainty")] %>%
                      tbl_df()
                  df_rev <- df_ %>%
                      mutate(ID1_ = ID2, ID2 = ID1, ID1 = ID1_) %>%
                      dplyr::select(-ID1_)

                  bind_rows(df_, df_rev) %>%
                      mutate(ID1 = factor(ID1, levels = order_),
                             ID2 = factor(ID2, levels = order_ %>% rev()))
              })

dom_plots <- map(1:length(dom_dfs),
                 function(i) {
                     dom_dfs %>%
                         .[[i]] %>%
                         ggplot(aes(ID1, ID2)) +
                         geom_raster(aes(fill = certainty)) +
                         scale_fill_gradient(low = palette$dark_primary,
                                             high = palette$light_primary,
                                             breaks = seq(0.5, 0.9, 0.1),
                                             limits = c(0.5, 1.0),
                                             guide = FALSE) +
                         scale_x_discrete(drop = FALSE) +
                         scale_y_discrete(drop = FALSE) +
                         NULL
                 })



library(grid)
library(gridExtra)
library(ggpubr)


# Extract the legend. Returns a gtable
leg <- dom_dfs %>%
    .[[1]] %>%
    {
        ggplot(., aes(ID1, ID2)) +
            geom_raster(aes(fill = certainty)) +
            scale_fill_gradient(low = palette$dark_primary,
                                high = palette$light_primary,
                                breaks = seq(0.5, 0.9, 0.1),
                                limits = c(0.5, 1.0)) +
            scale_x_discrete(drop = FALSE) +
            scale_y_discrete(drop = FALSE)
    } %>%
    get_legend() %>%
    as_ggplot() %>%
    {. + theme(plot.margin = margin(0, 0, 0, 0))} %>%
    ggplotGrob()

grob_list <- map(dom_plots, ggplotGrob)
grob_list <- c(grob_list, list(leg))


grid.newpage()
gridExtra::grid.arrange(grobs = grob_list,
                        layout_matrix = rbind(c(1,2,3),
                                              c(4,5,6),
                                              c(7,8,NA)))
