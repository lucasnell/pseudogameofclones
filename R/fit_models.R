


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

                        # -----------
                        # Priors:
                        # -----------
                        tau = -1.8970,
                        sigma_tau = 1.0000,
                        mu_theta = -1.3070,
                        sigma_theta = 0.7922,
                        gamma = -2.5360,
                        sigma_gamma = 1.0000,
                        mu_phi = -6.1300,
                        sigma_phi = 3.1620,
                        delta = -1.1510,
                        sigma_delta = 2.0000,
                        zeta = -0.9733,
                        sigma_zeta = 2.8340)

    growth_fit <- rstan::sampling(stanmodels$all_lines_plants, data = model_data_, ...)

    return(growth_fit)
}




#' Simulate cages with the same starting conditions.
#'
#'
#'
#' @param N_0 Matrix of `# aphids` at time t=0. Should be matrix of size
#'     `n_plants` x `n_lines`.
#' @param max_t Time steps per cage. Single integer.
#' @param R Max growth rates per aphid line. Vector of length `n_lines`.
#' @param A Density dependence per aphid line. Vector of length `n_lines`.
#' @param D_slope Dispersal slope per aphid line. Vector of length `n_lines`.
#' @param D_inter Dispersal intercept per aphid line. Vector of length `n_lines`.
#' @param process_error SD of process error. Single double.
#' @param plant_mort_coefs Coefficients (b0, then b1) for aphid mortality after
#'     plant starts dying.
#' @param plant_death_age Number of days after which a plant starts dying.
#' @param plant_repl Data frame containing two columns, the first
#'     indicating the date of replacement, the second indicating the plant that
#'     gets replaced on that date. There should be no repeat plants per date,
#'     but multiple plants can be replaced on the same date.
#' @param n_reps Number of reps (i.e., cages) per chain.
#' @param n_chains Number of chains to use. Note that values >1 will result in
#'     extra info being printed to the console. You can ignore this. Defaults to `1`.
#' @param ... Other options passed to `rstan::sampling()`.
#'
#' @export
#'
sim_cages <- function(N_0, max_t, R, A, D_slope, D_inter,
                      process_error,
                      plant_mort_coefs, plant_death_age, plant_repl,
                      n_reps, n_chains = 1, ...) {

    n_plants <- nrow(N_0)
    n_lines <- ncol(N_0)

    # plant_repl <- data_frame(date = sort(sample.int(10, 20, TRUE)) * 10L) %>%
    #     group_by(date) %>%
    #     mutate(plant = sample.int(5, n(), FALSE)) %>%
    #     ungroup()


    plant_repl <- plant_repl %>%
        set_names(c("date", "plant")) %>%
        arrange(date, plant)
    # This prevents a weird error if the number of rows in plant_repl is 1
    if ((plant_repl$date %>% unique() %>% length()) == 1) {
        plant_repl <- plant_repl %>%
            add_row(date = max_t + 100, plant = 1)
    }

    # Check that there aren't duplicates within a date:
    plant_repl %>%
        group_by(date) %>%
        summarize(unq = length(unique(plant)) == n()) %>%
        .[["unq"]] %>%
        all() %>%
        {if (!.) stop("No duplicates allowed within dates.")}

    # Time points when plant replacements occur
    repl_times <- plant_repl %>%
        .[["date"]] %>%
        unique()
    # Number of replacements
    n_repl <- repl_times %>%
        length()
    # Max # plant replaced per replacement
    max_repl <- plant_repl %>%
        group_by(date) %>%
        summarize(n = n()) %>%
        .[["n"]] %>%
        max()

    # Plants to replace when replacements occur (zeros are ignored):
    #     matrix<lower=0, upper=n_plants>[n_repl, max_repl] repl_plants;
    repl_plants <- array(0, dim = c(n_repl, max_repl))
    for (i in 1:n_repl) {
        d <- repl_times[i]
        plants_ <- plant_repl %>%
            filter(date == d) %>%
            .[["plant"]]
        repl_plants[i, 1:length(plants_)] <- plants_
    }


    model_data_ <- list(
        n_plants = n_plants,
        n_lines = n_lines,
        N_0 = N_0,
        max_t = max_t,
        R = R,
        A = A,
        D_slope = D_slope,
        D_inter = D_inter,
        process_error = process_error,
        plant_mort_coefs = plant_mort_coefs,
        plant_death_age = plant_death_age,
        n_repl = n_repl,
        max_repl = max_repl,
        repl_times = repl_times,
        repl_plants = repl_plants
    )

    z <- capture.output({
        stan_sims <- rstan::sampling(stanmodels$sim_cage, data = model_data_,
                                     algorithm = "Fixed_param", chains = n_chains,
                                     iter = n_reps, warmup = 0, ...)
    })

    n_reps <- n_reps * n_chains

    N_array <- rstan::extract(stan_sims, "N_out")[[1]]

    N_mat <- lapply(1:n_reps, function(j) {
        M <- lapply(1:n_plants, function(i) N_array[j,i,,])
        M <- do.call(rbind, M)
        cbind(rep(1:n_plants, each = max_t + 1), M)
    })
    N_mat <- do.call(what = rbind, N_mat)
    N_mat <- cbind(rep(1:n_reps, each = n_plants * (max_t + 1)), N_mat)

    N_df <- N_mat %>%
        as_data_frame() %>%
        set_names(c("rep", "plant", paste0("line", 1:n_lines))) %>%
        gather("line", "N", -rep, -plant) %>%
        mutate(line = gsub("line", "", line) %>% as.integer(),
               date = rep(0:max_t, n_reps * n_plants * n_lines),
               N = ifelse(N < 0, -Inf, N)) %>%
        identity()

    return(N_df)
}
