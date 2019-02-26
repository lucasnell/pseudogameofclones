
suppressPackageStartupMessages({
    library(clonewars)
})



source(".Rprofile")


# columns ending in `r` are rankings, from best to worst (high R, low A, high D are good):

data.frame(line = stan_estimates$names,
           R = stan_estimates$R,
           Rr = rank(1 / stan_estimates$R),
           A = stan_estimates$A,
           Ar = rank(stan_estimates$A),
           D = disp_estimates$b0,
           Dr = rank(1 / disp_estimates$b0))


#' Notes on parameter values:
#'   * higher `repl_threshold` and `zeta_t_thresh` gives advantage to Clover-2017-2
#'   * lower `repl_threshold` and `zeta_t_thresh` gives advantage to WI-L4Ø
#'

# 100 reps takes ~ 17 sec
# set.seed(78902356)
# sim <-
sim_reps(n_reps = 16, max_t = 1000, save_every = 1, N0 = 6,
                n_patches = 8,
                repl_times = seq(3, 1000, 3),
                repl_threshold = 10e3,
                zeta_t_thresh = 20,
                n_cores = 4) %>%  # ,
                # process_error = 0, disp_error = FALSE,
                # log_zeta_sd = 0, extinct_N = 0) %>%
    mutate_at(vars(rep, patch, line), factor) %>%
    # group_by(time, patch) %>%
    # mutate(prop = N / sum(N),
    #        lower_prop = c(0, cumsum(prop)[-n()]),
    #        upper_prop = cumsum(prop)) %>%
    # ungroup() %>%
    identity() %>%
    # filter(rep %in% 0:15) %>%
    group_by(rep, time, line) %>%
    summarize(N = sum(N)) %>%
    ungroup() %>%
    ggplot(aes(time, N)) +
    # geom_ribbon(aes(y = NULL, ymin = lower_prop, ymax = upper_prop, fill = line)) +
    geom_line(aes(color = line), size = 0.5) +
    # geom_line(aes(color = patch)) +
    facet_wrap( ~ rep, nrow = 4) +
    # facet_wrap( ~ patch, nrow = 2) +
    # theme(legend.position = "none") +
    theme(strip.text = element_blank()) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
#






# {
#     n_reps = 10
#     max_t = 100
#     n_patches = 20
#     N0 = NULL
#     R = NULL
#     A = NULL
#     D_vec = NULL
#     repl_times = seq(5, 100, 5)
#     repl_threshold = 500
#     process_error = 0
#     disp_error = FALSE
#     log_zeta_sd = 0
#     log_zeta_mean = NULL
#     mu_time = NULL
#     extinct_N = 1e-4
#     save_every = max_t %/% 100
#     by_patch = FALSE
#     n_cores = 1
#     show_progress = FALSE
#     line_names = NULL
# }
