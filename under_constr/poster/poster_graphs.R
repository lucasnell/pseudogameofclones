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

#           Alternative states plot

# =======================================================================================
# =======================================================================================


source(textConnection(readLines("under_constr/poster_graphs.R")[1:37]))

k1 = k2 = 100

a12 = a21 = 2
# Slopes for isocline lines:
m1 = -1 / a12
m2 = -a21

# Non-stable equilibrium coordinates
ns_equil <- c(x = {k1 - a12*k2} / {a12*(m2-m1)},
              y = m2 * ({k1 - a12*k2} / {a12*(m2-m1)}) + k2)
# Stable equilibrium coordinates
equil <- list(sp1 = c(x = k1, y = 0), sp2 = c(x = 0, y = k2))

# ----------
# Creating data frame for red and blue lines (the "zero isoclines")
# ----------
line_df <- data_frame(x = seq(0, 100, length.out = 11)) %>%
    mutate(sp1 = m1 * x + (k1 / a12), sp2 = m2 * x + k2) %>%
    gather(spp, y, sp1:sp2) %>%
    mutate(spp = factor(spp)) %>%
    filter(y >= 0)

# ----------
# Simulating starting points and where they'd go:
# ----------

sim_lv <- function(N01, N02, nt,
                   r = cbind(0.1, 0.1),
                   k = cbind(k1, k2),
                   a = cbind(a12, a21)) {

    N <- matrix(0, nrow = nt + 1, ncol = 2)
    N[1,] <- c(N01, N02)
    for (t in 1:nt) {
        N[t+1,] <- N[t,] + r * N[t,] * (k - {N[t,] + a * rev(N[t,])}) / k
    }

    N <- N %>%
        as_data_frame() %>%
        set_names(c("x", "y")) %>%
        identity()

    return(N)
}

# Creating geoms to add to plot

geom_lotka <- function(N01, N02, nt, arrow_n, ...) {
    lv <- sim_lv(N01 = N01, N02 = N02, nt = nt, ...)
    # geom_out <- lapply(arrow_n,
    #                    function(n_) {
    #                        geom_path(data = eval(lv[1:n_,]),
    #                                  size = 1, linejoin = "mitre",
    #                                  arrow = grid::arrow(type = "closed", angle = 20))
    #                    })
    geom_out <- list()
    geom_out <- c(geom_out,
                  list(geom_path(data = eval(lv), size = 1, linetype = 2),
                       geom_point(data = data_frame(x = N01, y = N02), size = 4,
                                  color = "black")))
    return(geom_out)
}


# ----------
# Creating ggplot object of alternative states
# ----------
as_plot <- ggplot(line_df, aes(x, y)) +
    # Red and blue lines
    geom_line(aes(color = spp), size = 2) +
    # 1:1 line:
    geom_path(data = data_frame(x = c(0, k2 + 5), y = c(0, k2 + 5)),
              size = 1, color = "gray80", linetype = 3) +
    # Hacky axes
    geom_path(data = data_frame(x = c(0, k1 + 10), y = c(0,0)), size = 1) +
    geom_path(data = data_frame(x = c(0,0), y = c(0, k2 + 10)), size = 1) +
    # Trajectories
    geom_lotka(90, 80, 200, c(5, 70)) +
    geom_lotka(90, 45, 200, 10) +
    # geom_lotka(90, 30, 200, 10) +
    geom_lotka(80, 90, 200, c(5, 70)) +
    geom_lotka(45, 90, 100, 5) +
    # geom_lotka(30, 90, 200, 10) +
    geom_lotka(10, 5,  100, 20) +
    geom_lotka(5,  10, 100, 20) +
    # Equilibrium points, first unstable, then stable
    geom_point(data = data_frame(x = ns_equil['x'], y = ns_equil['y']),
               shape = 21, color= 'black', fill = 'white', size = 8) +
    geom_point(data = data_frame(x = c(k1, 0), y = c(0, k2),
                                 species = factor(c('sp1', 'sp2'))),
               aes(color = species), size = 8, shape = 16) +
    # Aesthetics
    scale_color_manual(values = with(palette, c(default_primary, accent)),
                       guide = FALSE)


# Now simplifying it to have no axes:
as_plot <- as_plot +
    theme(axis.line.x = element_blank(),
          axis.line.y.left = element_blank())

as_plot

ggsave(filename = "figs/alt_states.pdf", plot = as_plot,
       width = 16, height = 16, units = "cm", bg = "white", useDingbats = FALSE)


n_ <- 150
as_traj_plot <- sim_lv(90, 45, n_) %>%
    bind_rows() %>%
    mutate(t = 0:n_) %>%
    gather("spp", "N", x:y, factor_key = TRUE) %>%
    ggplot(aes(t, N, color = spp)) +
    geom_line(size = 2) +
    scale_y_continuous(limits = c(0, k1 + 30)) +
    scale_color_manual(values = with(palette, c(default_primary, accent)),
                       guide = FALSE)


ggsave(filename = "figs/alt_states_traj.pdf", plot = as_traj_plot,
       width = 8, height = 8, units = "cm", bg = "white", useDingbats = FALSE)







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

