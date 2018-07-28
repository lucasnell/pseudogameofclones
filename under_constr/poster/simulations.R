
source("under_constr/poster/_preamble.R")

# add_labels <- function() {}


ggplot2::theme_set(
    ggplot2::theme_get() +
        theme(axis.ticks = element_line(),
              axis.title = element_text(),
              axis.text = element_text())
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
# Correcting for weirdness from WI-L4Ø due to only having one complete rep:
D_nb[D_nb$line == "WI-L4Ø",-1] <- colMeans(D_nb[D_nb$line != "WI-L4Ø", -1])
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
# saveRDS(sim_df, file = "data-raw/pair_sims.rds")
sim_df <- readr::read_rds("data-raw/pair_sims.rds")

sim_by_cage <- sim_df %>%
    mutate(line = factor(line,
                         levels = seq_along(levels(growth$line)),
                         labels = levels(growth$line)),
           rep = factor(rep),
           pair = factor(pair)) %>%
    group_by(pair, rep, line, date) %>%
    summarize(N = sum(N)) %>%
    arrange(pair, rep, date, line) %>%
    group_by(pair, rep, date) %>%
    mutate(prop = N / sum(N),
           prop_end = cumsum(prop),
           prop_start = lag(prop_end, default = 0)) %>%
    ungroup() %>%
    mutate(X = log(N))



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
for (i in 1:n_lines) line_mat[,i] <- levels(growth$line)[i]

# Now sort by surv_mat:
for (i in 1:n_lines) {
    # Survival of line j minus survival of line i:
    surv_diffs <- numeric(n_lines)
    for (j in 1:n_lines) surv_diffs[j] <- surv_mat[j,i] - surv_mat[i,j]
    inds <- sort(surv_mat[,i] - surv_mat[i,], index.return = TRUE, decreasing = TRUE)$ix
    line_mat[i,] <- line_mat[i,inds]
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
        diff == 0 ~ 1,
        TRUE ~ 2
    ) %>% factor(levels = 0:2, labels = c("lose", "draw", "win"))) %>%
    mutate(line1 = factor(paste(line1), levels = rev(levels(line1)))) %>%
    ggplot(aes(line2, line1)) +
    geom_point(aes(color = diff), size = 8) +
    scale_color_manual(values = c("red", "black", "green")) +
    xlab("Opponent line") +
    ylab("Focal line")


# Differences in survivals
diff_survs <- surv_mat - t(surv_mat)

# library(igraph)
# graph <- graph_from_adjacency_matrix(adjmatrix = surv_mat, diag = FALSE)
# transitivity(graph)




survs %>%
    ggplot(aes(line, color = line)) +
    geom_linerange(aes(ymin = 0, ymax = n), size = 1) +
    geom_point(aes(y = n), size = 3) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    scale_y_continuous("Proportion of times survived", limits = c(0,1)) +
    xlab(NULL) +
    facet_wrap(~ pair, nrow = 7, scales = "free_x") +
    # coord_flip() +
    theme(axis.ticks = element_line(),
          axis.title = element_text(),
          axis.text = element_text()) +
    NULL


