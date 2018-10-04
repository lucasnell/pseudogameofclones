
# library(clonewars)
# library(testthat)
# source(".Rprofile")

#' Simulate one time series for one line.
#'
#' @noRd
#'
sim_ts <- function(r, alpha, X0, sd_epsilon, max_t) {
    X <- numeric(max_t + 1)
    X[1] <- X0
    for (t in 1:max_t) {
        X[t+1] <- X[t] + r * (1 - alpha * exp(X[t])) + rnorm(1, sd = sd_epsilon)
    }
    return(X)
}

#' Simulate multiple time series for one line.
#'
#' @noRd
#'
sim_one_line <- function(r, alpha, X0, sd_epsilon, max_t, n_ts) {
    X <- matrix(NA_real_, max_t + 1, n_ts)
    for (i in 1:n_ts) {
        X[,i] <- sim_ts(r, alpha, X0, sd_epsilon, max_t)
    }
    return(X)
}

#' Simulate multiple time series for multiple lines.
#'
#'
#' @noRd
#'
sim_lines <- function(n_ts, max_t, rs, alphas, sd_epsilon, X0s = NULL) {

    stopifnot(length(n_ts) == 1)
    stopifnot(length(max_t) == 1)

    n_lines <- length(rs)
    stopifnot(length(alphas) == n_lines)
    if (is.null(X0s)) X0s <- rep(log(2), n_lines)
    stopifnot(length(X0s) == n_lines)

    X_list <- mapply(rs, alphas, X0s, FUN = sim_one_line,
                     MoreArgs = list(sd_epsilon = sd_epsilon, max_t = max_t, n_ts = n_ts),
                     SIMPLIFY = FALSE, USE.NAMES = FALSE)

    X_mat <- lapply(1:n_lines, function(i) cbind(i, rep(1:n_ts, each = max_t + 1),
                                                 rep(1:(max_t + 1), n_ts),
                                                 as.numeric(X_list[[i]])))
    X_mat <- do.call(rbind, X_mat)
    colnames(X_mat) <- c("line", "rep", "date", "X")
    X_df <- as_data_frame(X_mat)
    mutate_at(X_df, vars(line, rep, date), as.integer)
}




X_df <- sim_lines(n_ts = 5, max_t = 20,
                  rs = seq(0.20, 0.45, 0.05),
                  alphas = seq(0.0015, 0.0040, 0.0005),
                  sd_epsilon = 0.25)
stopifnot(min(X_df$X) > 0)

X_df %>%
    ggplot(aes(date, X, color = factor(rep))) +
    geom_line() +
    facet_wrap(~ factor(line)) +
    scale_color_brewer(palette = "Dark2", guide = FALSE)


stan_fit <- fit_lines(X_df, control = list(adapt_delta = 0.85))

for (i in 1:6) {
    x <- c("sigma_epsilon", "rho", "sigma_rho", "phi", "sigma_phi_a", "sigma_phi_w")[i]
    y <- theta[seq(2, length(theta),2)][i]
    z <- apply(rstan::extract(stan_fit, x, permuted = FALSE), 3, sd)
    cat(sprintf("parameter: %s\n prior SD:     %.5g\n posterior SD: %.5g\n\n", x, y, z))
}; rm(x,y,i)


apply(rstan::extract(stan_fit, "sigma_epsilon", permuted = FALSE), 3, mean)
apply(rstan::extract(stan_fit, "sigma_epsilon", permuted = FALSE), 3, quantile,
      probs = c(0.025, 0.975))

exp(
    apply(rstan::extract(stan_fit, "rho", permuted = FALSE), 3, mean) +
        apply(rstan::extract(stan_fit, "sigma_rho", permuted = FALSE), 3, mean) *
        apply(rstan::extract(stan_fit, "Z_r", permuted = FALSE), 3, mean))
seq(0.20, 0.45, 0.05)

gtools::inv.logit(
    apply(rstan::extract(stan_fit, "phi", permuted = FALSE), 3, mean) +
        apply(rstan::extract(stan_fit, "sigma_phi_a", permuted = FALSE), 3, mean) *
        apply(rstan::extract(stan_fit, "Z_a_a", permuted = FALSE), 3, mean))
seq(0.0015, 0.0040, 0.0005)


apply(rstan::extract(stan_fit, "sigma_phi_a", permuted = FALSE), 3, mean)
apply(rstan::extract(stan_fit, "sigma_phi_w", permuted = FALSE), 3, mean)

