


#' Fit one aphid-growth time series.
#'
#' @param X Name of the object in `data` that contains the log-transformed counts
#'     through time. `NA`s are silently removed.
#' @param data An optional data frame, list, or environment that contains `X`.
#'     By default, variables are taken from the environment from which
#'     the function was called.
#' @param ... Arguments passed to `rstan::sampling` (e.g., iter, chains).
#'
#' @return A `stanfit` object containing the model fit.
#' @export
#'
#' @examples
#'
#' R0 <- 0.25
#' alpha <- 0.002
#' sigma_process <- 0.1
#'
#' X <- as.numeric(sim_lines(R0, alpha, n_reps = 1, nobs_ts = 20, process = 0.1))
#'
#' # (You would normally not change iter to 100)
#' fit <- fit_one_ts(X, iter = 100)
#'
#'
fit_one_ts <- function(X, data = sys.frame(sys.parent()), ...) {

    model_data <- eval(quote(list(X = X, N = length(X))), envir = data)
    model_data$X <- model_data$X[!is.na(model_data$X)]

    growth_fit <- rstan::sampling(stanmodels$one_ts, data = model_data, ...)

    return(growth_fit)
}



#' Fit multiple time series for a single aphid line.
#'
#' @param X Name of the matrix in `data` that contains the log-transformed counts
#'     through time (1 time series per column).
#'     Each column should contain `NA`s at the end (ONLY the end) if
#'     it wasn't observed as many times as was the time series in the matrix with
#'     the most observations.
#'     \emph{Other than at the end, missing values are not yet supported}.
#' @param data An optional data frame, list, or environment that contains `X`.
#'     By default, variables are taken from the environment from which
#'     the function was called.
#' @param ... Arguments passed to `rstan::sampling` (e.g., iter, chains).
#'
#' @return A `stanfit` object containing the model fit.
#' @export
#'
#' @examples
#'
#' R0 <- 0.29
#' alpha <- 0.002
#' process <- 0.1
#' n_reps <- 5
#' nobs_ts <- 11:15
#'
#' X <- sim_lines(R0, alpha, n_reps, nobs_ts, process)
#'
#' # (You would normally not change iter to 100)
#' fit <- fit_one_line(X, iter = 100)
#'
fit_one_line <- function(X, data = sys.frame(sys.parent()), ...) {

    X_ <- eval(X, envir = data)
    reps_ts_ <- apply(X_, 2, function(x) sum(!is.na(x)))
    X_[is.na(X_)] <- 0

    model_data <- list(N_ts = ncol(X_), X = X_,
                       reps_ts = reps_ts_, max_reps = max(reps_ts_))

    growth_fit <- rstan::sampling(stanmodels$one_line, data = model_data, ...)

    return(growth_fit)
}

