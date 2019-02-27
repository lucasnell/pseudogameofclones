
# This file uses peak abundance for density dependences
# and max per-capita growth to esimate intrinsic growth rates

library(clonewars)
source(".Rprofile")




# Calculate growth rate
get_r <- function(.df) {
    pcg <- .df$X[-1] - lag(.df$X)[-1] # per-capita growth
    n_t <- 4 # number of time steps to take average by
    if (length(pcg) < 4) return(NULL)
    # Max average pcg over `n_t` time steps:
    r_ <- max(sapply(1:(length(pcg) - n_t + 1), function(i) mean(pcg[i:(i+n_t-1)])))
    tibble(line = .df$line[1], rep = .df$rep[1], r = r_)
}
# Calculate density dependence
get_a <- function(.df) {
    # Inverse of max density over time series:
    a_ <- 1 / max(.df$N, na.rm = TRUE)
    tibble(line = .df$line[1], rep = .df$rep[1], a = a_)
}

# I can use all data for the growth rates:
R <- load_data(filter_pars = NULL, remove_unfinished = FALSE) %>%
    split(interaction(.$line, .$rep, drop = TRUE)) %>%
    map_dfr(~ get_r(.x))

A <- load_data(filter_pars = NULL, remove_unfinished = FALSE) %>%
    split(interaction(.$line, .$rep, drop = TRUE)) %>%
    map_dfr(~ get_a(.x)) %>%
    mutate(k = 1 / a)


library(lme4)

lmer(r ~ (1 | line), R, REML = FALSE) %>% AIC()
lm(r ~ 1, R) %>% AIC()

lmer(a ~ (1 | line), A, REML = FALSE) %>% AIC()
lm(a ~ 1, A) %>% AIC()

lmer(k ~ (1 | line), A, REML = FALSE, control = lmerControl(optimizer = "bobyqa")) %>% AIC()
lm(k ~ 1, A) %>% AIC()



R_plot <- R %>%
    ggplot(aes(line, r)) +
    geom_jitter(aes(color = line), height = 0, width = 0.25) +
    stat_summary(fun.data = "mean_cl_boot", shape = 5) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 11, color = "black")) +
    coord_flip() +
    ylab("Growth rate")
R_plot
A_plot <- A %>%
    ggplot(aes(line, k)) +
    geom_jitter(aes(color = line), height = 0, width = 0.25) +
    stat_summary(fun.data = "mean_cl_boot", shape = 5) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 11, color = "black")) +
    coord_flip() +
    ylab("Carrying capacity")
A_plot

ggsave("~/Desktop/simple_R.pdf", R_plot, width = 4, height = 4)
ggsave("~/Desktop/simple_A.pdf", A_plot, width = 4, height = 4)
