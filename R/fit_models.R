


#' Fit multiple time series for multiple aphid lines.
#'
#' @param data_df A data frame that contains columns for `line`, `rep`, and `X`.
#' @param line An optional string specifying the name for the column in `data_df` that
#'     contains which line each observation belongs to.
#'     Defaults to `"line"`.
#' @param rep An optional string specifying the name for the column in `data_df` that
#'     contains which rep each observation belongs to.
#'     Defaults to `"rep"`.
#' @param date An optional string specifying the name for the column in `data_df` that
#'     contains the date for each observation.
#'     Defaults to `"date"`.
#' @param X An optional string specifying the name for the column in `data_df` that
#'     contains the `X` value for each observation (`X = log(N)`).
#'     Defaults to `"X"`.
#' @param theta_ Optional vector of length 12 that specifies prior hyperparameters.
#' @param model_name String specifying the name of the model to fit.
#'     The options are
#'     `"full_model"`, `"no_within_alpha"`, `"no_among_alpha"`,
#'     `"one_alpha"`, `"one_r"`, `"one_r_alpha"`, `"pass_sigma_epsilon"`,
#'     or `"full_model_plant_death"`.
#'     Defaults to `"full_model"`.
#' @param sigma_epsilon Value for the SD of the process error if using a model
#'     that doesn't esimate this value itself.
#'     This argument is not used if using a model that estimates the process error SD.
#'     Defaults to `NULL`.
#' @param ... Arguments passed to `rstan::sampling` (e.g., iter, chains).
#'
#' @return A `stanfit` object containing the model fit.
#' @export
#'
#'
#'
fit_lines <- function(data_df, line, rep, date, X, theta_,
                      model_name = c("full_model", "no_within_alpha",
                                     "no_among_alpha", "one_alpha",
                                     "one_r", "one_r_alpha",
                                     "pass_sigma_epsilon",
                                     "full_model_plant_death"),
                      sigma_epsilon = NULL,
                      ...) {

    model_name <- match.arg(model_name)

    cat(sprintf("\n~~~~~~~~~~~~~~~~\nUsing \"%s model\"\n~~~~~~~~~~~~~~~~\n\n",
                model_name))

    stopifnot(inherits(data_df, "data.frame"))
    if (missing(line)) line <- quote(line)
    if (missing(rep)) rep <- quote(rep)
    if (missing(date)) date <- quote(date)
    if (missing(X)) X <- quote(X)
    # theta is already defined:
    if (missing(theta_)) theta_ <- theta

    line <- substitute(line)
    rep <- substitute(rep)
    date <- substitute(date)
    X <- substitute(X)

    data_df <- data_df %>%
        dplyr::select(!!line, !!rep, !!date, !!X) %>%
        dplyr::arrange(!!line, !!rep, !!date) %>%
        mutate_if(is.factor, as.integer) %>%
        identity()

    # number of observations for each time series:
    n_ts_ <- data_df %>%
        dplyr::group_by(!!line, !!rep) %>%
        summarize() %>%
        nrow()
    n_obs_ <- nrow(data_df)
    n_per_ <- data_df %>%
        dplyr::group_by(!!line, !!rep) %>%
        summarize(n_ = n()) %>%
        .[["n_"]] %>%
        set_names(NULL)
    stopifnot(sum(n_per_) == n_obs_)

    X_ <- data_df[[X]]

    n_lines_ <- length(unique(data_df[[line]]))
    # Line number for each time series:
    L_ <- data_df %>%
        dplyr::group_by(!!line, !!rep) %>%
        summarize() %>%
        .[["line"]] %>%
        set_names(NULL)

    model_data_ <- list(
        n_ts = n_ts_,
        n_obs = n_obs_,
        n_per = n_per_,
        X = X_,
        n_lines = n_lines_,
        L = L_,
        theta = theta_
    )

    if (model_name == "pass_sigma_epsilon") {
        err <- FALSE
        if (is.null(sigma_epsilon)) {
            err <- TRUE
        } else if (!is.numeric(sigma_epsilon) | length(sigma_epsilon) != 1) {
            err <- TRUE
        } else if (sigma_epsilon < 0) {
            err <- TRUE
        }
        if (err) {
            stop("\nIf using a model that doesn't estimate the process error SD ",
                 "then you must pass a single number >= 0 to the `sigma_epsilon` ",
                 "argument.",
                 call. = FALSE)
        }
        model_data_$sigma_epsilon <- sigma_epsilon
    }
    if (model_name == "full_model_plant_death") {
        model_data_$theta <- NULL
    }

    growth_fit <- rstan::sampling(stanmodels[[model_name]], data = model_data_, ...)

    return(growth_fit)
}


