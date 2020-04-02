

# Note from 7 Jan 2020: Not sure any of this is useful...



# For all below, X = log(N), where N is the # aphids.


#' Function for figuring out how many aphids were mentioned in the comments.
#'
#' @param comments A vector of comments
#'
#' @noRd
#'
parse_comments <- function(comments) {
    comments_ = strsplit(comments, "")
    out <- comments_ %>%
        map_int(~ keep(.x, function(y) y %in% c(paste(0:9), " ", ",")) %>%
                    paste(collapse = "") %>%
                    trimws() %>%
                    strsplit(", | |,") %>%
                    unlist() %>%
                    as.integer() %>%
                    sum())
    out <- ifelse(is.na(out) | is.null(out), 0, out)
    return(out)
}




#' Initial read of CSV file.
#'
#'
#'
#' @noRd
#'
initial_read <- function() {

    file <- system.file("extdata", "aphid_counts_raw.csv", package = "clonewars",
                        mustWork = TRUE)

    # Aphid lines we're still using
    kept_lines <- c("R10", "WIA-5D", "WI-L4", "WI-L4Ø", "UT3", "WI-2016-593",
                    "Clover-2017-2", "Clover-2017-6")

    # Column types
    cols_ <- readr::cols(line = "c",
                         comments = "c",
                         observer = "c",
                         .default = "i")

    read_csv(file, col_types = cols_) %>%
        mutate(line = ifelse(line == 'WI-L4 (H+3)', 'WI-L4', line),
               line = ifelse(line == 'WI-L4ØA', 'WI-L4Ø', line),
               line = ifelse(line == 'WI-L4 (Ham-)', 'WI-L4Ø', line),
               date = as.Date(paste(year, month, day, sep = "-"))) %>%
        # Change any NAs to zeros:
        mutate_at(vars(matches("_juv$|_adults$")),
                  function(x) ifelse(is.na(x), 0, x)) %>%
        # Filter for lines that we still have and should keep for analyses
        filter(line %in% kept_lines) %>%
        mutate(disp = parse_comments(comments),  # <-- "dispersed" aphids
               N = stem1_juv + stem1_adults + leaf1_juv + leaf1_adults +
                   stem2_juv + stem2_adults + leaf2_juv + leaf2_adults +
                   stem3_juv + stem3_adults + leaf3_juv + leaf3_adults +
                   disp,
               # makes no sense for it to be 0, then >0 the next day:
               N = ifelse(N == 0, 1, N)) %>%
        dplyr::select(-matches("_juv$|_adults$"), -year, -month, -day) %>%
        mutate_at(vars(rep, N, disp), list(as.integer)) %>%
        mutate(X = log(N),
               # Retained lines that have *Hamiltonella defensa*:
               ham = ifelse(line %in% c("R10", "WI-L4", "UT3", "Clover-2017-2"),
                            1, 0)) %>%
        group_by(line, rep) %>%
        mutate(date = as.integer(date - min(date))) %>%
        ungroup() %>%
        arrange(line, rep, date)
}



#' Filter one data frame based on another.
#'
#' @param in_df A data frame to be filtered, with the columns `line` and `rep`.
#' @param comp_df A data frame to use for filtering, with the columns `line` and `rep`.
#' @param exclude A logical for whether to exclude rows from `in_df` based on `comp_df`.
#'
#' @return A filtered version of `in_df`.
#'
#' @noRd
#'
filter_line_rep <- function(in_df, comp_df, exclude) {
    in_comp <- map2_lgl(in_df$line, in_df$rep,
                        ~ any(map_lgl(1:nrow(comp_df),
                                      function(i) {
                                          .x == comp_df$line[i] & .y == comp_df$rep[i]
                                      })))
    if (exclude) return(in_df %>% filter(!in_comp))
    return(in_df %>% filter(in_comp))
}



#' Removing lines that aren't yet done.
#'
#' @inheritParams handle_NAs
#' @inheritParams load_data
#'
#' @noRd
#'
handle_unfinished <- function(growth, remove_unfinished) {
    if (remove_unfinished) {
        not_done <- growth %>%
            group_by(line, rep) %>%
            arrange(date) %>%
            summarize(p = tail(N, 1) / max(N),
                      three_down = all((tail(N, 4) %>% diff() %>% sign(.)) == -1)) %>%
            ungroup() %>%
            filter(p > 0.8 & !three_down) %>%
            arrange(line, rep) %>%
            dplyr::select(line, rep)
        if (nrow(not_done) > 0) {
            growth <- filter_line_rep(growth, not_done, exclude = TRUE)
        }
    }
    return(growth)
}


#' Dealing with missing values in input data frame.
#'
#' @param growth Input data frame that's been read from excel sheet.
#' @inheritParams load_data
#'
#' @noRd
#'
handle_NAs <- function(growth, allow_NA, impute_fxn) {
    missing <- growth %>%
        group_by(line, rep) %>%
        arrange(date) %>%
        summarize(dates = list(which(! min(date):max(date) %in% date) + min(date) - 1),
                  diff = length(date) - length(min(date):max(date))) %>%
        ungroup() %>%
        filter(diff != 0) %>%
        identity()
    if (any(missing$diff > 0)) {
        missing %>%
            filter(diff > 0) %>%
            print()
        stop("\nNo duplicates allowed.", call. = FALSE)
    }

    if (allow_NA & nrow(missing) > 0) {

        missing <- missing %>%
            unnest(cols = dates) %>%
            rename(date = dates)

        df_ <- growth[1:nrow(missing),] %>%
            mutate_if(~ inherits(.x, "integer"), function(x) NA_integer_) %>%
            mutate_if(~ inherits(.x, "numeric"), function(x) NA_real_) %>%
            mutate_if(~ inherits(.x, "character"), function(x) NA_character_) %>%
            mutate(line = missing$line, rep = missing$rep, date = missing$date)

        growth <- bind_rows(growth, df_) %>%
            arrange(line, rep, date)

        # Impute data if impute function provided
        if (!is.null(impute_fxn)) {
            stopifnot(inherits(impute_fxn, "function"))
            growth <- growth %>%
                group_by(line, rep) %>%
                arrange(date) %>%
                mutate(X = impute_fxn(X),
                       N = ifelse(is.na(N), exp(X), N)) %>%
                ungroup()
        }
    } else if (nrow(missing) > 0) {
        growth <- filter_line_rep(growth, missing, exclude = TRUE)
    }
    growth <- growth %>%
        mutate(date = as.integer(date))
    return(growth)
}



#' Standardize day 0 aphid counts.
#'
#' Some counts have the day we first added aphids included, while some don't.
#' This function fixes this.
#'
#' @inheritParams handle_NAs
#'
#' @noRd
#'
standardize_day0 <- function(growth) {
    growth <- growth %>%
        filter(!(date == 0 & disp == N)) %>%
        group_by(line, rep) %>%
        mutate(date = as.integer(date - min(date) + 1)) %>%
        ungroup() %>%
        split(.$line) %>%
        map_dfr(~ split(.x, .x$rep) %>%
                    map_dfr(function(y_) {
                        add_row(y_, line = y_$line[1], rep = y_$rep[1], date = 0L,
                                N = 2, X = log(2), ham = y_$ham[1]) %>%
                            arrange(date)
                    }))
    return(growth)
}





#' Function for filtering the end of a time series (for last X above a threshold).
#'
#' Note that this function is designed to be run within a single time series, meaning
#' within a line and rep combination.
#' So for a data frame, you should `dplyr::group` by line and rep, then use
#' `dplyr::summarize` when using this function.
#'
#' @param X_vec A vector of log(N) through time. This is assumed to be sorted by date!
#' @param p The proportion of the max value of `X_vec` that is kept after the max is
#'     reached.
#'     For example, if `p = 0.9` and `mX` is the max value of `X_vec`, then we would
#'     retain all points until `X_vec == mX`, PLUS we would retain any values
#'     after `X_vec == mX` where `X_vec >= p * mX` (in this case, when `X_vec`
#'     values are at least 90% of the maximum).
#'
#' @noRd
#'
end_filter <- function(X_vec, p) {
    max_X <- max(X_vec, na.rm = TRUE)
    max_ind <- which(X_vec == max_X)[1]
    end_ <- tail(which(X_vec >= p * max_X), 1)
    return(1:length(X_vec) <= end_)
}



#' Function for filtering the start of a time series (for first X above a threshold).
#'
#' Note that this function is designed to be run within a single time series, meaning
#' within a line and rep combination.
#' So for a data frame, you should `dplyr::group` by line and rep, then use
#' `dplyr::summarize` when using this function.
#'
#' @param X_vec A vector of log(N) through time. This is assumed to be sorted by date!
#' @param p The proportion of the max value of `X_vec` where points start being retained.
#'     For example, if `p = 0.5` and `mX` is the max value of `X_vec`, then we would
#'     retain all points starting with `X_vec >= p * mX`
#'     (in this case, when `X_vec` values are at least 50% of the maximum).
#'
#' @noRd
#'
start_filter <- function(X_vec, p) {
    max_X <- max(X_vec, na.rm = TRUE)
    threshold <- max_X * p
    # Index to first above the threshold
    first_ind <- which(X_vec >= threshold)[1]
    return(1:length(X_vec) >= first_ind)
}




#' Filter a time series if desired.
#'
#' @inheritParams handle_NAs
#' @inheritParams load_data
#'
#' @noRd
#'
filter_data <- function(growth, filter_pars) {
    if (!is.null(filter_pars)) {
        growth <- growth %>%
            # Filter the beginning and end of the time series:
            group_by(line, rep) %>%
            arrange(date)
        # Start filter:
        if (!is.null(filter_pars$start)) {
            growth <- growth %>%
                filter(start_filter(X, filter_pars$start))
        }
        # End filter:
        if (!is.null(filter_pars$end)) {
            growth <- growth %>%
                filter(end_filter(X, filter_pars$end))
        }
        growth <- growth %>%
            ungroup()
    }
    return(growth)
}




#' Load aphid growth data.
#'
#' @param file A filename to read from. If left empty, this will read from the
#'     default path on your DropBox folder.
#' @param filter_pars A list with the names `"start"` and `"end"`, containing single
#'     numbers with threshold for filtering the beginning and ending of time series,
#'     respectively. Set to `NULL` to avoid filtering entirely.
#'     Defaults to `list(start = 0.0, end = 0.9)`.
#' @param allow_NA Boolean for whether time series with NAs should be included.
#'     Defaults to `TRUE`.
#'
#' @return A data frame containing the aphid population-growth data.
#'
#' @export
#'
#'
#' @examples
#'
#' growth <- load_data()
#'
load_data <- function(filter_pars = list(start = 0.0, end = 0.9),
                      remove_unfinished = FALSE, allow_NA = TRUE,
                      impute_fxn = impute) {

    growth <-
        # First read through, to clean up the excel sheet:
        initial_read() %>%
        # We no longer need these columns:
        dplyr::select(-observer, -comments) %>%
        # Removing lines that aren't yet done:
        handle_unfinished(remove_unfinished) %>%
        # Dealing with missing values:
        handle_NAs(allow_NA, impute_fxn) %>%
        # Remove rows where the first two aphids were recorded and re-add it so that
        # these dates are present for all reps (they're currently not all present)
        standardize_day0() %>%
        # Filter part(s) of time series if desired:
        filter_data(filter_pars) %>%
        # Change to factor now to avoid having to drop levels, and make sure it's
        # ordered properly
        mutate(line = factor(line)) %>%
        arrange(line, rep, date)

    return(growth)
}





#' Make prediction data frame from model output and the original data frame.
#'
#' @param stan_fit The `stanfit` object containing the model fit.
#' @param orig_data The original data frame.
#' @param line An optional parameter specifying the name of the column in `orig_data`
#'     that contains info on the aphid line.
#' @param rep An optional parameter specifying the name of the column in `orig_data`
#'     that contains info on the rep within each aphid line.
#'
#' @export
#'
make_pred_df <- function(stan_fit, orig_data, line, rep, date, alpha = 0.05) {

    if (missing(line)) line <- quote(line)
    if (missing(rep)) rep <- quote(rep)
    if (missing(date)) date <- quote(date)
    line <- substitute(line)
    rep <- substitute(rep)
    date <- substitute(date)

    pred_df <- orig_data %>%
        dplyr::arrange(!!line, !!rep, !!date) %>%
        mutate(X_pred = rstan::extract(stan_fit, "X_pred", permuted = FALSE) %>%
                   apply(3, median),
               X_lower = rstan::extract(stan_fit, "X_pred", permuted = FALSE) %>%
                   apply(3, quantile, probs = alpha / 2),
               X_upper = rstan::extract(stan_fit, "X_pred", permuted = FALSE) %>%
                   apply(3, quantile, probs = 1 - alpha / 2)) %>%
        dplyr::select(line, rep, date, X, X_pred, X_lower, X_upper, everything())

    return(pred_df)
}





#' Impute missing data.
#'
#' From [Moritz et al. (2015)](https://arxiv.org/abs/1510.03924), the
#' seasonal kalman filter from the `zoo` package works well, so I'm creating this function
#' here to impute missing values below.
#' We're only missing 1 or 2 time points across a time series of >15, so this
#' shouldn't have too big an impact on our results.
#'
#' @param X Vector of log(N) through time.
#'
#' @export
#'
impute <- function(X) {
    if (any(is.na(X))) {
        X_ <- X %>%
            # Change to weekly frequency for the fit
            zoo::zooreg(frequency = 7) %>%
            zoo::na.StructTS() %>%
            as.numeric()
        return(X_)
    }
    return(X)
}
