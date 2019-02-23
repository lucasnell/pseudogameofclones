
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
           dD = ifelse(dD < 0, 0, dD),
           past_death = date - which(N == max(N))[1] # days past plant death
           ) %>%
    ungroup()





stan_fit <- read_rds("data-raw/stan_fit.rds")


library(grid)


pop_pars <-
    bind_rows(
        apply(rstan::extract(stan_fit, "R", permuted = FALSE), 3, quantile,
              probs = c(0.05, 0.95)) %>%
            t() %>%
            tbl_df() %>%
            set_names(c("lower", "upper")) %>%
            mutate(mean = apply(rstan::extract(stan_fit, "R", permuted = FALSE), 3, mean),
                   par = "R", line = sort(unique(growth$line))),
        apply(rstan::extract(stan_fit, "A", permuted = FALSE), 3, quantile,
              probs = c(0.05, 0.95)) %>%
            t() %>%
            tbl_df() %>%
            set_names(c("lower", "upper")) %>%
            mutate(mean = apply(rstan::extract(stan_fit, "A", permuted = FALSE), 3, mean),
                   par = "A", line = sort(unique(growth$line)))
    ) %>%
    mutate(par = factor(par, levels = c("R", "A")))


pop_par_plot <- pop_pars %>%
    mutate(line = factor(paste(line), levels = rev(levels(line)))) %>%
    ggplot(aes(line)) +
    geom_linerange(aes(ymin = lower, ymax = upper)) +
    geom_point(aes(y = mean), size = 3, shape = 16) +
    facet_wrap(~ par, scales = "free_x") +
    scale_y_continuous(breaks = c(0.0015, 0.002, 0.0025, 0.3, 0.33, 0.36)) +
    coord_flip() +
    theme(axis.ticks.x = element_line(),
          axis.text.x = element_text(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.text.y = element_text())



ggsave("figs/pop_ests.pdf", pop_par_plot, width = 16, height = 16, units = "cm",
       bg = "white", useDingbats = FALSE)


# =======================================================================================
# =======================================================================================

# Dispersal

# =======================================================================================
# =======================================================================================


dispersal_df <- load_data(filter_pars = NULL,
                    remove_unfinished = FALSE) %>%
    filter(date != 0) %>%
    mutate(rep = factor(rep)) %>%
    group_by(line, rep) %>%
    arrange(date) %>%
    mutate(disp = ifelse(is.na(disp), 0, disp),
           N = round(impute(N)),
           X = log(N),
           pN = N - disp,  # numbers of aphids on plant
           dD = disp - lag(disp, default = 0), # number of new dispersed aphids
           dD = ifelse(dD < 0, 0, dD),
           past_death = date - which(N == max(N))[1], # days past plant death
           disp_b = ifelse(dD == 0, 0, 1)) %>%
    ungroup() %>%
    arrange(line, rep, date) %>%
    identity()


logit <- function(p) {
    suppressWarnings({x <- log(p/(1 - p))})
    # x <- ifelse(is.nan(x), NA, x)
    return(x)
}
inv_logit <- function(x) {
    p <- 1 / (1 + exp(-x))
    p <- ifelse(is.na(p) & !is.na(x), 1, p)
    return(p)
}
predict_disp <- function(N, past_death, aphid_line) {
    inds <- map_int(aphid_line, ~ which(disp_estimates$binom$line == .x))
    # Binomial estimate (Pr(dispersal > 0)):
    b0 <- disp_estimates$binom$b0[inds]
    b1 <- disp_estimates$binom$b1[inds]
    b2 <- disp_estimates$binom$b2[inds]
    pr_disp <- inv_logit(b0 + b1 * past_death + b2 * past_death^2)
    # Negative binomial estimate (# dispersed | dispersal occurs):
    b0 <- disp_estimates$nb$b0[inds]
    n_disp <- exp(b0) * N
    return(pr_disp * n_disp)
}
# disp_plot <- dispersal_df %>%
#     mutate(pred_dD = predict_disp(N, past_death, line)) %>%
#     ggplot(aes(past_death, dD)) +
#     geom_point(alpha = 0.5, shape = 16, color = palette$default_primary) +
#     geom_point(aes(y = pred_dD), size = 0.75, color = palette$accent, shape = 1) +
#     scale_color_brewer(palette = "Set3", guide = FALSE) +
#     facet_wrap(~ line, nrow = 2) +
#     scale_y_continuous("Dispersed aphids", trans = "log1p",
#                        breaks = c(0, 4^(0:3))) +
#     scale_x_continuous("Days past plant death") +
#     NULL

disp_plot <- expand.grid(past_death = -20:5, N = 100, line = levels(dispersal_df$line)) %>%
    tbl_df() %>%
    mutate(disp_prop = predict_disp(N, past_death, line) / N) %>%
    ggplot(aes(past_death, disp_prop)) +
    geom_line(aes(color = line), size = 1) +
    scale_color_brewer(NULL, palette = "Dark2") +
    scale_y_continuous("Proportion dispersed", trans = "log",
                       breaks = 10^(-4:-2)) +
    scale_x_continuous("Days past plant death") +
    theme(legend.position = c(0.1, 1), legend.justification = c(0, 1)) +
    NULL

ggsave("figs/disp_est.pdf", disp_plot, width = 16, height = 16, units = "cm",
       bg = "white", useDingbats = FALSE)

