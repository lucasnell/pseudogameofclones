


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
fit_one_ts <- function(X, data = sys.frame(sys.parent()), ...) {

    model_data <- eval(quote(list(X = X, N = length(X))), envir = data)
    model_data$X <- model_data$X[!is.na(model_data$X)]

    growth_fit <- rstan::sampling(stanmodels$one_ts, data = model_data, ...)

    return(growth_fit)
}

