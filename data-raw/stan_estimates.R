

library(clonewars)

stan_fit <- system.file("extdata", "stan_fit.rds", package = "clonewars",
                        mustWork = TRUE)
stan_fit <- readRDS(stan_fit)

R <- apply(rstan::extract(stan_fit, "R", permuted = FALSE), 3, median)
A <- apply(rstan::extract(stan_fit, "A", permuted = FALSE), 3, median)
process_error <- median(rstan::extract(stan_fit, "sigma_epsilon", permuted = FALSE))

logZ <- log(apply(rstan::extract(stan_fit, "Z", permuted = FALSE), 3, median))
log_zeta_mean <- mean(logZ)
log_zeta_sd <- sd(logZ)

# all estimates are the same, so we can just take the first one:
mu_time <- rstan::extract(stan_fit, "mu_time")[[1]][[1]]

line_names <- load_data() %>%
    .[["line"]] %>%
    levels()

stan_estimates <- list(
    names = line_names,
    R = as.numeric(R),
    A = as.numeric(A),
    process_error = process_error,
    log_zeta_mean = log_zeta_mean,
    log_zeta_sd = log_zeta_sd,
    mu_time = mu_time
)

usethis::use_data(stan_estimates)

