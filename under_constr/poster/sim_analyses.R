# !diagnostics off


source("under_constr/poster/_preamble.R")
source("under_constr/poster/simulations.R")

ggplot2::theme_set(
    ggplot2::theme_get() +
        theme(axis.ticks = element_line(),
              axis.title = element_text(),
              axis.text = element_text(),
              strip.text = element_text())
)

sims <- test_threshes(c(100, 1000))

pool_sims <- sims$N %>%
    arrange(thresh, date, line) %>%
    group_by(thresh, date) %>%
    mutate(prop = N / sum(N),
           prop_end = cumsum(prop),
           prop_start = lag(prop_end, default = 0)) %>%
    ungroup()


pool_sims %>%
    # filter(date < 500) %>%
    ggplot(aes(date, N, color = line)) +
    geom_line() +
    facet_grid(thresh ~ ., scales = "free_y") +
    scale_color_brewer(palette = "Dark2", guide = FALSE)


# survs <- pool_sims %>%
#     # filter(date == 1000, repl_age == 0, thresh == 1000) %>%
#     filter(date == max(date), repl_age == 0, thresh < 1000) %>%
#     # filter(date == 200, repl_age == 0) %>%
#     group_by(rep, line) %>%
#     summarize(s = ifelse(N == 0, 0L, 1L)) %>%
#     ungroup() %>%
#     spread(line, s) %>%
#     dplyr::select(-rep) %>%
#     as.matrix()
#
# colMeans(survs)
# cor(survs)
#
# pool_sims %>%
#     filter(thresh < 1000) %>%
#     # filter(thresh == 1000, date <= 500) %>%
#     ggplot(aes(date, N, group = factor(rep))) +
#     geom_line(alpha = 0.1) +
#     facet_wrap(~ factor(line))


outcomes <- pool_sims %>%
    filter(date == max(date)) %>%
    group_by(thresh, rep) %>%
    summarize(lines = list(sort(line[N > 0])),
              outcome = paste(sort(line[N > 0]), collapse = "_")) %>%
    group_by(thresh, outcome) %>%
    summarize(n = n(),
              lines = lines[1],
              reps = list(rep)) %>%
    ungroup() %>%
    dplyr::select(-outcome) %>%
    split(.$thresh) %>%
    map(~ dplyr::select(.x, -thresh))


# Colorblind palette:
pal <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1',
         '#4575b4','#313695')
# Removing the very light colors:
pal <- pal[-c(5, 6)]


set.seed(64060)
alt_states_plots <- map(1:length(outcomes),
                        function(i) {
                            pool_sims %>%
                                filter(thresh == levels(thresh)[i],
                                       rep %in% map_int(outcomes[[i]]$reps,
                                                        ~ .x[sample.int(length(.x), 1)]),
                                       date <= 1000
                                ) %>%
                                ggplot(aes(date, N, color = factor(line))) +
                                geom_line(size = 0.5) +
                                # ggplot(aes(date, fill = factor(line))) +
                                # geom_ribbon(aes(ymin = prop_start, ymax = prop_end),
                                #             color = NA) +
                                facet_wrap(~ factor(rep), ncol = 2) +
                                scale_color_manual(values = pal, guide = FALSE) +
                                # scale_color_brewer("Aphid\nline:", palette = "RdYlBu",
                                #                    guide = FALSE) +
                                # scale_fill_brewer("Aphid\nline:", palette = "Dark2",
                                #                    guide = FALSE) +
                                scale_x_continuous(breaks = seq(0, 1000, 500)) +
                                scale_y_continuous(breaks = seq(0, 3000, 1500)) +
                                # scale_y_continuous(breaks = seq(0, 1.0, 0.5)) +
                                theme(strip.text = element_blank(),
                                      axis.text = element_blank(),
                                      axis.title = element_blank(),
                                      axis.ticks = element_blank()) +
                                NULL
                        })





ggsave(filename = sprintf("figs/alt_states_%02i.pdf", 1),
       plot = alt_states_plots[[1]],
       width = 14, height = 3 * 4, units = "cm", bg = "white",
       useDingbats = FALSE)
ggsave(filename = sprintf("figs/alt_states_%02i.pdf", 2),
       plot = alt_states_plots[[2]],
       width = 14, height = 5 * 4, units = "cm", bg = "white",
       useDingbats = FALSE)



leg <- alt_states_plots %>%
    .[[1]] %>%
    {. + scale_color_brewer("Aphid\nline:", palette = "Dark2") +
            theme(legend.text = element_blank(),
                  legend.title = element_blank(),
                  legend.position = c(0.5, 0.5),
                  legend.justification = c(0.5,0.5)) +
            guides(color = guide_legend(override.aes = list(size = 1.5),
                                        direction = "horizontal", nrow = 2))
    } %>%
    ggpubr::get_legend() %>%
    ggpubr::as_ggplot() %>%
    {
        . +
            theme(plot.margin = margin(0, 0, 0, 0))
    }

ggsave(filename = "figs/alt_states_legend.pdf", plot = leg,
       width = 7.5, height = 5, units = "cm", bg = "white",
       useDingbats = FALSE)


#





# with(sim_env, ls())
#
# cor(sim_env$R, sim_env$A)
# cor(exp(sim_env$D_vec), sim_env$A)
# cor(sim_env$plant_mort_0, sim_env$A)
# cor(sim_env$plant_mort_1, sim_env$plant_mort_0)
#
# plot(sim_env$R ~ sim_env$A)
# plot(exp(sim_env$D_vec) ~ sim_env$A)
# plot(sim_env$plant_mort_0 ~ sim_env$A)
# plot(sim_env$plant_mort_1 ~ sim_env$plant_mort_0)


#






# ======================================================================================
# ======================================================================================
# ======================================================================================
# ======================================================================================
# ======================================================================================
# ======================================================================================
# ======================================================================================
# ======================================================================================
# ======================================================================================

# # Mean total number of aphids per plant:
# pool_sims_N_pp <- readr::read_rds("data-raw/pool_sims_N.rds") %>%
#     group_by(pool, rep) %>%
#     summarize(N = sum(N)) %>%
#     ungroup() %>%
#     mutate(pool_size = as.integer(floor(log10(pool)) + 1L))

pool_sims_X <- readr::read_rds("data-raw/pool_sims_X.rds")
pool_sims_Z <- readr::read_rds("data-raw/pool_sims_Z.rds")
pool_sims <- readr::read_rds("data-raw/pool_sims_N_by_linerep.rds")

# pool_sims_N <- pool_sims_N %>%
#     mutate(pool_size = as.integer(floor(log10(pool)) + 1L))
pool_sims_X <- pool_sims_X %>%
    mutate(pool_size = as.integer(floor(log10(pool)) + 1L))
pool_sims_Z <- pool_sims_Z %>%
    mutate(pool_size = as.integer(floor(log10(pool)) + 1L))
pool_sims <- pool_sims %>%
    mutate(pool_size = as.integer(floor(log10(pool)) + 1L))


sim_env <- readr::read_rds("data-raw/sim_env.rds")



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

growth <-
    load_data(impute_fxn = impute, filter_pars = NULL) %>%
    mutate(line = paste(line)) %>%
    # bind_rows(clonewars:::load_pz_data(impute_fxn = impute, filter_pars = NULL)) %>%
    mutate_at(vars(line, rep), list(factor)) %>%
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

n_lines <- growth$line %>% levels() %>% length()
n_plants <- sim_env$n_plants





# ================================================================================
# ================================================================================

# Pooled simulations

# ================================================================================
# ================================================================================



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
    mutate(line = factor(line, levels = 1:n_lines, labels = levels(growth$line))) %>%
    identity()


pool_ranks %>%
    ggplot(aes(line, rank, color = factor(line))) +
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
    mutate(line = factor(line, levels = 1:n_lines, labels = levels(growth$line))) %>%
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


df_ <- pool_sims %>%
    arrange(pool_size) %>%
    split(.$pool_size) %>%
    .[[8]]



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

        if (pool_size_ == sim_env$n_lines) M <- cbind(M[!is.na(M)])

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
                              i_ <- levels(growth$line)[as.integer(i)]
                              ind <- disp_estimates$line == i_
                              stopifnot(sum(ind) == 1)
                              pd <- 2
                              gtools::inv.logit(
                                  plant_death$after_max_mort_coefs$b0[ind] +
                                      plant_death$after_max_mort_coefs$b1[ind] * pd)
                          }),
           disp = map_dbl(line,
                          function(i) {
                              i_ <- levels(growth$line)[as.integer(i)]
                              ind <- disp_estimates$line == i_
                              stopifnot(sum(ind) == 1)
                              b0 <- exp(disp_estimates$b0[ind])
                              return(b0 * 500)
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




# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================

# Plot of certainty in dominance rankings (not necessary for poster)

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
