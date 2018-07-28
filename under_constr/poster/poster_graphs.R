
source("under_constr/poster/_preamble.R")


# =======================================================================================
# =======================================================================================

#           Aphid growth plots

# =======================================================================================
# =======================================================================================


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


reps <- growth %>%
    filter(line == "WI-2016-593") %>%
    .[["rep"]] %>%
    unique() %>%
    paste() %>%
    identity()
reps <- c(reps, reps[1:3])
# reps <- reps[1]

example_ts <- lapply(reps,
       function(r) {
           growth %>%
               filter(line == "WI-2016-593", rep == r) %>%
               dplyr::select(line, rep, date, N, dD) %>%
               ggplot(aes(date)) +
               geom_line(aes(y = N), size = 1.5, color = palette$default_primary) +
               geom_line(aes(y = dD * 10), size = 1.5, color = palette$secondary_text) +
               xlab(NULL) +
               scale_y_continuous(NULL, limits = c(0, 740),
                                  sec.axis = sec_axis(~ . / 10,
                                                      name = NULL))
})


example_ts[1:7]


# for (i in 1:length(example_ts)) {
for (i in 1:1) {
    fn <- sprintf("figs/growth_ts_%02i.pdf", i)
    ggsave(filename = fn, plot = example_ts[[i]],
           width = 8, height = 8, units = "cm", bg = "white")
}








# =======================================================================================
# =======================================================================================

# Simulations

# =======================================================================================
# =======================================================================================



source(textConnection(readLines("under_constr/poster_graphs.R")[1:68]))


stan_fit <- read_rds("data-raw/stan_fit.rds")

n_plants <- 16
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
extinct_N <- 5
n_cages <- 100
n_cores <- parallel::detectCores() - 1

lines <- 1:8
set.seed(654651984L)
sim_df <- cwsims::sim_cages(n_cages, N_0[,lines], max_t, R[lines], A[lines],
                            D_binom[lines,], D_nb[lines,], process_error,
                            plant_mort_0[lines], plant_mort_1[lines],
                            plant_death_age_mean, plant_death_age_sd,
                            repl_times, repl_age, extinct_N, n_cores) %>%
    mutate(X = log(N)) %>%
    identity()

sim_by_cage <- sim_df %>%
    mutate(line = factor(line,
                         levels = seq_along(levels(growth$line)[lines]),
                         labels = levels(growth$line)[lines]),
           rep = factor(rep)) %>%
    group_by(rep, line, date) %>%
    summarize(N = sum(N)) %>%
    group_by(rep, date) %>%
    arrange(line) %>%
    mutate(prop = N / sum(N),
           prop_end = cumsum(prop),
           prop_start = lag(prop_end, default = 0)) %>%
    ungroup() %>%
    mutate(X = log(N))


example_ts <- sim_by_cage %>%
    filter(rep == 1) %>%
    ggplot(aes(date, N, color = line)) +
    geom_line(size = 1) +
    scale_color_manual(values = line_palette, guide = FALSE)

ggsave(filename = "figs/example_ts.pdf", plot = example_ts,
       width = 16, height = 8, units = "cm", bg = "white", useDingbats = FALSE)



# Example cage time series:

example_cage_ts <- sim_by_cage %>%
    filter(rep == 1, date %% 7 == 0) %>%
    ggplot(aes(date, fill = line)) +
    geom_ribbon(aes(ymin = prop_start, ymax = prop_end), color = NA) +
    geom_line(aes(y = prop_end), color = palette$divider, size = 0.25) +
    scale_fill_manual(values = line_palette, guide = FALSE)

ggsave(filename = "figs/example_cage_ts.pdf", plot = example_cage_ts,
       width = 16, height = 8, units = "cm", bg = "white", useDingbats = FALSE)

