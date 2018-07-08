


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






#' Fit multiple time series for multiple aphid lines.
#'
#' @param X Name of the matrix in `data` that contains the log-transformed counts
#'     through time (1 time series per column).
#'     This option is not required if the proper object inside `data` is literally
#'     named `X` (as would be the case if you used `line_data()`).
#'     Each column should contain `NA`s at the end (ONLY the end) if
#'     it wasn't observed as many times as was the time series in the matrix with
#'     the most observations.
#'     \emph{Other than at the end, missing values are not yet supported}.
#' @param line_ts Vector (of the same length as number of columns in `X`) containing
#'     the line number for each time-series column in `X`.
#'     This option is not required if the proper object inside `data` is literally
#'     named `line_ts` (as would be the case if you used `line_data()`).
#' @param data An optional data frame, list, or environment that contains `X`.
#'     By default, variables are taken from the environment from which
#'     the function was called.
#' @param ... Arguments passed to `rstan::sampling` (e.g., iter, chains).
#'
#' @return A `stanfit` object containing the model fit.
#' @export
#'
#'
#'
fit_lines <- function(data, X, line_ts, plants = FALSE, ...) {

    if (missing(data)) data <- sys.frame(sys.parent())

    if (missing(X)) X <- quote(X)
    if (missing(line_ts)) line_ts <- quote(line_ts)
    X <- substitute(X)
    line_ts <- substitute(line_ts)

    X_ <- eval(X, envir = data)
    line_ts_ <- eval(line_ts, envir = data)
    n_lines_ <- length(unique(line_ts_))
    if (length(line_ts_) != ncol(X_)) {
        stop("\nIn `fit_lines`, `line_ts` must have the same length as number of ",
             "columns in `X`", call. = FALSE)
    }


    nobs_ts_ <- apply(X_, 2, function(x) sum(!is.na(x)))
    X_[is.na(X_)] <- 0

    model_data <- list(N_ts = ncol(X_),
                       max_reps = max(nobs_ts_),
                       n_lines = n_lines_,
                       nobs_ts = nobs_ts_,
                       line_ts = line_ts_,
                       X = X_,

                       # -----------
                       # Priors:
                       # -----------
                       w_0 = 0.15000,
                       eta = 0.50000,
                       mu_theta = -1.30700,
                       sigma_theta = 0.79220,
                       x_0 = 0.07922,
                       gamma = 1.00000,
                       mu_phi = -6.13000,
                       sigma_phi = 3.16200,
                       y_0 = 0.31620,
                       delta = 5.00000,
                       z_0 = 0.67350,
                       zeta = 7.19700)

    if (plants) {
        growth_fit <- rstan::sampling(stanmodels$all_lines_plants, data = model_data, ...)
    } else {
        growth_fit <- rstan::sampling(stanmodels$all_lines, data = model_data, ...)
    }

    return(growth_fit)
}

