


#' Fit multiple time series for multiple aphid lines.
#'
#' @param data An optional data frame, list, or environment that contains `X`.
#'     By default, variables are taken from the environment from which
#'     the function was called.
#' @param X Name of the matrix in `data` that contains the log-transformed counts
#'     through time (1 time series per column).
#'     This option is not required if the proper object inside `data` is literally
#'     named `X` (as would be the case if you used `line_data()`).
#'     Each column should contain `NA`s at the end (ONLY the end) if
#'     it wasn't observed as many times as was the time series in the matrix with
#'     the most observations.
#'     \emph{Other than at the end, missing values are not yet supported}.
#' @param L Vector (of the same length as number of columns in `X`) containing
#'     the line number for each time-series column in `X`.
#'     This option is not required if the proper object inside `data` is literally
#'     named `L` (as would be the case if you used `line_data()`).
#' @param ... Arguments passed to `rstan::sampling` (e.g., iter, chains).
#'
#' @return A `stanfit` object containing the model fit.
#' @export
#'
#'
#'
fit_lines <- function(data, X, L, ...) {

    if (missing(data)) data <- sys.frame(sys.parent())

    if (missing(X)) X <- quote(X)
    if (missing(L)) L <- quote(L)
    X <- substitute(X)
    L <- substitute(L)

    X_ <- eval(X, envir = data)
    L_ <- eval(L, envir = data)
    n_lines_ <- length(unique(L_))
    if (length(L_) != ncol(X_)) {
        stop("\nIn `fit_lines`, `L` must have the same length as number of ",
             "columns in `X`", call. = FALSE)
    }


    nobs_ts_ <- apply(X_, 2, function(x) sum(!is.na(x)))
    X_[is.na(X_)] <- 0

    model_data_ <- list(N_ts = ncol(X_),
                        max_reps = max(nobs_ts_),
                        n_lines = n_lines_,
                        nobs_ts = nobs_ts_,
                        L = L_,
                        X = X_,
                        # Priors:
                        theta = theta)

    growth_fit <- rstan::sampling(stanmodels$all_lines_plants, data = model_data_, ...)

    return(growth_fit)
}


