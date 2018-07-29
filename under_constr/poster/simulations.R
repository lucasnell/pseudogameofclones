
source("under_constr/poster/_preamble.R")



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
N_0 <- matrix(rep(48/16, n_lines * n_plants), n_plants, n_lines)
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
                                repl_times, repl_age, extinct_N, n_cores) %>%
        mutate(pair = sprintf("%i_%i", line1_, line2_), X = log(N),
               line = map_int(line, ~ lines_[.x]))
    return(sim_df)
}


# set.seed(654651985L)
# sim_df <- map_dfr(pw_combs, ~ pw_sims(.x[1], .x[2]))
# readr::write_rds(sim_df, path = "data-raw/pair_sims.rds")
# sim_df <- readr::read_rds("data-raw/pair_sims.rds")
#
# sim_by_cage <- sim_df %>%
#     mutate(line = factor(line,
#                          levels = seq_along(levels(growth$line)),
#                          labels = levels(growth$line)),
#            rep = factor(rep),
#            pair = factor(pair)) %>%
#     group_by(pair, rep, line, date) %>%
#     summarize(N = sum(N)) %>%
#     arrange(pair, rep, date, line) %>%
#     group_by(pair, rep, date) %>%
#     mutate(prop = N / sum(N),
#            prop_end = cumsum(prop),
#            prop_start = lag(prop_end, default = 0)) %>%
#     ungroup() %>%
#     mutate(X = log(N))
# readr::write_rds(sim_by_cage, path = "data-raw/pair_sims_by_cage.rds")

sim_by_cage <- readr::read_rds("data-raw/pair_sims_by_cage.rds")


{
    # This should have zero rows bc having total extinction makes no sense:
    sim_by_cage %>%
        filter(date == max_t) %>%
        group_by(pair, rep) %>%
        summarize(N = sum(N)) %>%
        ungroup() %>%
        filter(N == 0) %>%
        nrow() %>%
        `==`(0) %>%
        print()
    # This should also have zero rows bc it means aphids spontaneously appeared
    sim_by_cage %>%
        group_by(pair, rep, line) %>%
        filter(N == 0, dplyr::lead(N) > 0) %>%
        ungroup() %>%
        nrow() %>%
        `==`(0) %>%
        print()
}



survs <- sim_by_cage %>%
    filter(date == max_t) %>%
    group_by(pair, line) %>%
    summarize(n = mean(N > 0)) %>%
    ungroup()

# Rows are focal, columns are opponents
surv_mat <- matrix(0, 8, 8)
smf <- survs %>%
    group_by(pair) %>%
    mutate(i = as.integer(line),
           j = as.integer(rev(line))) %>%
    ungroup()
for (k in 1:length(pw_combs)) {
    i_ <- pw_combs[[k]][1]
    j_ <- pw_combs[[k]][2]
    surv_mat[i_, j_] <- smf %>% filter(i == i_, j == j_) %>% .[["n"]]
    surv_mat[j_, i_] <- smf %>% filter(i == j_, j == i_) %>% .[["n"]]
}

surv_mat
# Rows are focal, columns are opponents
line_mat <- matrix(NA_character_, 8, 8)
for (i in 1:n_lines) line_mat[i,] <- levels(growth$line)[i]

# Differences in survivals
diff_survs <- surv_mat - t(surv_mat)

# Now sort line_mat by surv_mat:
for (i in 1:n_lines) {
    inds <- sort(diff_survs[,i], index.return = TRUE, decreasing = TRUE)$ix
    line_mat[,i] <- line_mat[inds,i]
}


pair_df <- survs %>%
    group_by(pair) %>%
    mutate(line2 = rev(line), diff = n - rev(n)) %>%
    ungroup() %>%
    select(line1 = line, line2, diff)


pair_df %>%
    arrange(line1) %>%
    mutate(diff = case_when(
        diff < 0  ~ 0,
        diff > 0 ~ 1,
        TRUE ~ 2
    ) %>% factor(levels = 0:2, labels = c("lose", "win", "draw"))) %>%
    mutate(line1 = factor(paste(line1), levels = rev(levels(line1)))) %>%
    ggplot(aes(line2, line1)) +
    geom_point(aes(color = diff), size = 8) +
    scale_color_manual(values = c("red", "blue", "black")) +
    xlab("Opponent line") +
    ylab("Focal line")




pair_ranks <- expand.grid(opponent = 1:n_lines, line = 1:n_lines) %>%
    tbl_df() %>%
    arrange(opponent, line) %>%
    group_by(opponent) %>%
    # Descending rank:
    mutate(rank = rank(-diff_survs[,opponent[1]])) %>%
    ungroup() %>%
    mutate(opponent = factor(opponent, levels = 1:n_lines, labels = levels(growth$line)),
           line = factor(line, levels = 1:n_lines, labels = levels(growth$line)))


pair_ranksN <- sim_by_cage %>%
    group_by(pair, line) %>%
    summarize(N = mean(N)) %>%
    group_by(pair) %>%
    mutate(opponent = rev(line), diff = (N - rev(N)) / N) %>%
    ungroup() %>%
    select(opponent, line, N, diff) %>%
    add_row(opponent = sim_by_cage$line %>% unique() %>% sort(),
            line = sim_by_cage$line %>% unique() %>% sort(),
            N = NA, diff = 0) %>%
    arrange(opponent, line) %>%
    group_by(opponent) %>%
    # Descending rank:
    mutate(rank = rank(-diff)) %>%
    ungroup()



# # Dealing with tied ranks:
# pair_ranks$rank[pair_ranks$rank %% 1 != 0] <- pair_ranks %>%
#     filter(rank %% 1 != 0) %>%
#     .[["rank"]] %>%
#     `+`(rep(c(1, -1) * 0.25, 3))
#
# pair_ranks %>%
#     ggplot(aes(opponent, rank)) +
#     geom_text(aes(label = line))






# library(igraph)
# graph <- graph_from_adjacency_matrix(adjmatrix = surv_mat, diag = FALSE)
# transitivity(graph)

lines_ <- 1:n_lines
set.seed(549489)
pool_sims <- cwsims::sim_cages(n_cages, N_0[,lines_], max_t, R[lines_], A[lines_],
                               D_binom[lines_,], D_nb[lines_,], process_error,
                               plant_mort_0[lines_], plant_mort_1[lines_],
                               plant_death_age_mean, plant_death_age_sd,
                               repl_times, repl_age, extinct_N, n_cores) %>%
    mutate(X = log(N))

pool_by_cage <- pool_sims %>%
    mutate(line = factor(line,
                         levels = seq_along(levels(growth$line)),
                         labels = levels(growth$line)),
           rep = factor(rep)) %>%
    group_by(rep, line, date) %>%
    summarize(N = sum(N)) %>%
    arrange(rep, date, line) %>%
    group_by(rep, date) %>%
    mutate(prop = N / sum(N),
           prop_end = cumsum(prop),
           prop_start = lag(prop_end, default = 0)) %>%
    ungroup() %>%
    mutate(X = log(N))


# pool_ranksN <-

pool_ranks <- pool_by_cage %>%
    group_by(line) %>%
    summarize(N = mean(N)) %>%
    ungroup() %>%
    # Descending rank by mean N:
    mutate(N = rank(-N)) %>%
    # Same but by survival:
    mutate(surv = pool_by_cage %>%
               filter(date == max_t) %>%
               group_by(line) %>%
               summarize(n = mean(N > 0)) %>%
               ungroup() %>%
               .[["n"]] %>%
               {rank(-1 * .)}) %>%
    gather("method", "rank", N:surv)


pool_ranks
pair_ranks %>%
    group_by(line) %>%
    summarize(rank = mean(rank)) %>%
    .[["rank"]] %>%
    `*`(-1) %>%
    rank()
survs %>%
    group_by(line) %>%
    summarize(n = mean(n)) %>%
    .[["n"]] %>%
    `*`(-1) %>%
    rank()




rank_plot <- pair_ranks %>%
    ggplot(aes(rank, as.integer(opponent))) +
    geom_vline(xintercept = c(1, 8), linetype = 1, size = 0.5) +
    geom_hline(yintercept = 1:8, linetype = 3, size = 0.25) +
    geom_vline(data = pool_ranks %>% filter(method == "N"), aes(xintercept = rank),
               size = 2, linetype = 1, color = palette$default_primary) +
    # geom_path(aes(color = line), size = 1) +
    # geom_point(aes(color = line), size = 3) +
    geom_path(data = pair_ranksN, aes(color = line), size = 1, linetype = 1,
               color = palette$accent) +
    geom_point(data = pair_ranksN, aes(color = line), size = 3, shape = 16,
               color = palette$accent) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    scale_linetype_manual(values = c(1, 3), guide = FALSE) +
    facet_wrap(~ line, nrow = 2, scales = "free_x") +
    # reversed to make line #1 up top:
    scale_y_reverse("Opponent",
                    breaks = 1:n_lines,
                    labels = levels(pair_ranks$opponent)) +
    scale_x_continuous("Rank", breaks = 1:n_lines) +
    theme(axis.line.y.left = element_blank(), axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          strip.text = element_blank(),
          panel.spacing.y = unit(2, "cm")) +
    NULL



ggsave(filename = "figs/line_ranks.pdf", plot = rank_plot,
       width = 20, height = 10, units = "cm", bg = "white", useDingbats = FALSE)



surv_plot <- pool_by_cage %>%
    filter(date == max_t) %>%
    group_by(line) %>%
    summarize(n = mean(N > 0)) %>%
    ggplot(aes(line, color = line)) +
    geom_linerange(aes(ymin = 0, ymax = n), size = 1, color = palette$dark_primary) +
    geom_point(aes(y = n), size = 6, color = palette$accent) +
    # scale_color_brewer(palette = "Dark2", guide = FALSE) +
    scale_y_continuous("Proportion of times survived",
                       limits = c(0, 0.5), breaks = seq(0, 1, 0.25)) +
    xlab(NULL) +
    coord_flip() +
    theme(axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          strip.text = element_blank(),
          panel.spacing.y = unit(2, "cm")) +
    NULL

# pool_by_cage %>%
#     group_by(line, rep) %>%
#     summarize(N = mean(N)) %>%
#     # group_by(line) %>%
#     # summarize(n = mean(N), lower = quantile(N, 0.025), upper = quantile(N, 0.975)) %>%
#     ggplot(aes(line, color = line)) +
#     # geom_linerange(aes(ymin = lower, ymax = upper), size = 0.5) +
#     # geom_point(aes(y = n), size = 3) +
#     geom_point(aes(y = N), alpha = 0.1,
#                position = position_jitter(height = 0, width = 0.25)) +
#     stat_summary(aes(y = N), geom = "point", fun.y = mean, size = 4, color = "black") +
#     scale_color_brewer(palette = "Dark2", guide = FALSE) +
#     scale_y_continuous("Mean N over all reps", trans = "log") +
#     xlab(NULL) +
#     coord_flip() +
#     NULL
#
#
# pool_by_cage %>%
#     ggplot(aes(date, prop, color = line)) +
#     # ggplot(aes(date, gtools::logit(prop), group = rep, color = line)) +
#     geom_line(aes(group = rep), alpha = 0.1) +
#     facet_wrap(~ line, ncol = 4) +
#     # geom_smooth(method = "loess", se = FALSE, span = 0.4, color = "black", linetype = 2) +
#     # stat_summary(geom = "line", fun.y = mean, color = "black", size = 0.5) +
#     scale_color_brewer(palette = "Dark2", guide = FALSE) +
#     ylab("Relative abundance") +
#     xlab("Day") +
#     NULL


ggsave(filename = "figs/line_survs.pdf", plot = surv_plot,
       width = 20, height = 8, units = "cm", bg = "white", useDingbats = FALSE)

