


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
#' @param model_name String specifying the name of the model to fit.
#'     The options are
#'     `"full_model_plant_death"` or `"full_model"`.
#'     Defaults to `"full_model_plant_death"`.
#' @param ... Arguments passed to `rstan::sampling` (e.g., iter, chains).
#'
#' @return A `stanfit` object containing the model fit.
#' @export
#'
#'
#'
#'
fit_lines <- function(data_df, line, rep, date, X,
                      model_name = c("full_model_plant_death",
                                     "full_model_plant_death_R",
                                     "full_model_plant_death_RN",
                                     "full_model"),
                      ...) {

    model_name <- match.arg(model_name)

    cat(sprintf("\n~~~~~~~~~~~~~~~~\nUsing \"%s model\"\n~~~~~~~~~~~~~~~~\n\n",
                model_name))

    stopifnot(inherits(data_df, "data.frame"))
    if (missing(line)) line <- quote(line)
    if (missing(rep)) rep <- quote(rep)
    if (missing(date)) date <- quote(date)
    if (missing(X)) X <- quote(X)

    line <- substitute(line)
    rep <- substitute(rep)
    date <- substitute(date)
    X <- substitute(X)

    data_df <- data_df %>%
        dplyr::select(!!line, !!rep, !!date, !!X) %>%
        arrange(!!line, !!rep, !!date) %>%
        mutate_if(is.factor, as.integer) %>%
        identity()

    # number of observations for each time series:
    n_ts_ <- data_df %>%
        group_by(!!line, !!rep) %>%
        summarize() %>%
        nrow()
    n_obs_ <- nrow(data_df)
    n_per_ <- data_df %>%
        group_by(!!line, !!rep) %>%
        summarize(n_ = n()) %>%
        .[["n_"]] %>%
        set_names(NULL)
    stopifnot(sum(n_per_) == n_obs_)

    X_ <- data_df[[X]]

    n_lines_ <- length(unique(data_df[[line]]))
    # Line number for each time series:
    L_ <- data_df %>%
        group_by(!!line, !!rep) %>%
        summarize() %>%
        .[["line"]] %>%
        set_names(NULL)

    model_data_ <- list(
        n_ts = n_ts_,
        n_obs = n_obs_,
        n_per = n_per_,
        X = X_,
        n_lines = n_lines_,
        L = L_
    )

    growth_fit <- rstan::sampling(stanmodels[[model_name]], data = model_data_, ...)

    return(growth_fit)
}


