

#' Simulate log-transformed counts of aphid populations for one or more lines.
#'
#' @param R0 Vector of population growth rates for each line.
#' @param alpha Vector of density dependencies for each line.
#' @param n_reps Number of reps per aphid line.
#' @param nobs_ts Number of observations per time series.
#' @param sigma_process Process error. Defaults to `0.1`.
#'
#' @return A matrix with simulated log(N) through time.
#' @export
#'
sim_lines <- function(R0, alpha, n_reps, nobs_ts, sigma_process = 0.1) {

    # Number of aphid lines:
    n_lines <- length(R0)
    if (length(alpha) != n_lines) {
        stop("\nIn sim_lines, R0 and alpha must be the same length.",
             call. = FALSE)
    }
    if (length(n_reps) != n_lines) {
        stop("\nIn sim_lines, R0 and n_reps must be the same length.",
             call. = FALSE)
    }
    # Aphid line per time series
    line_ts <- c(lapply(1:length(n_reps), function(i) rep(i, n_reps[i])),
                 recursive = TRUE)
    # Number of time series
    N_ts <- sum(n_reps)

    # # Reps per aphid line
    # n_reps <- as.integer(round(runif(n_lines, 3, 7)))
    # # # observations per time series
    # nobs_ts <- as.integer(round(runif(N_ts, 10, 15)))
    # R0 <- round(abs(rnorm(n_lines, 0.29, 0.003)), 5)
    # alpha <- round(abs(rnorm(n_lines, 0.002, 0.0003)), 5)


    sigma_process <- sigma_process[1]


    X <- matrix(NA, max(nobs_ts), N_ts)
    X[1,] <- log(runif(ncol(X), 15, 25))
    for (i in 1:N_ts) {
        for (t in 1:(nobs_ts[i]-1)) {
            X[(t+1),i] <- X[t,i] + R0[line_ts[i]] *
                (1 - alpha[line_ts[i]] * exp(X[t,i])) + rnorm(1, 0, sigma_process)
        }
    }

    return(X)
}
