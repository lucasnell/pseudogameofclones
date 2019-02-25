


#' Simulate multiple reps and simplify output.
#'
#'
#' @param n_reps Number of reps to simulate.
#' @param n_patches Number of patches to simulate.
#' @param max_t Max time points to simulate for each rep.
#' @param N_0 Starting abundances for each aphid line on each patch.
#'     Can be a single number if you want the same value for each line on each patch,
#'     or a matrix if you want to specify everything.
#' @param R Growth rates for each line.
#' @param A Density dependence for each line.
#' @param D_vec
#' @param process_error SD of process error. Set to 0 for no process error.
#' @param disp_error Boolean for whether to include dispersal stochasticity.
#' @param log_zeta_mean Mean of the distribution of log(zeta) values.
#' @param log_zeta_sd SD of the distribution of log(zeta) values.
#' @param mu_time Mean of time values.
#' @param repl_times Vector of times at which to replace plants.
#' @param repl_threshold Threshold above which patches are replaced.
#' @param extinct_N Threshold below which a line is considered extinct.
#' @param save_every Abundances will be stored every `save_every` time points.
#' @param by_patch Logical for whether to summarize abundances by patch, rather
#'     than separately by line and patch.
#' @param n_cores Number of cores to use. Defaults to \code{1}.
#' @param show_progress Boolean for whether to show progress bar. Defaults to
#'     \code{FALSE}.
#' @param line_names Vector of names to assign to lines.
#'
#'
#' @importFrom purrr map_dfr
#' @importFrom dplyr as_tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom dplyr arrange
#'
#' @export
#'
#'
sim_reps <- function(n_reps,
                     n_patches = 8,
                     max_t = 1000,
                     N0 = NULL,
                     R = NULL,
                     A = NULL,
                     D_vec = NULL,
                     process_error = NULL,
                     disp_error = TRUE,
                     log_zeta_mean = NULL,
                     log_zeta_sd = NULL,
                     zeta_t_thresh = 20,
                     mu_time = NULL,
                     repl_times = NULL,
                     repl_threshold = 500,
                     extinct_N = 1e-4,
                     save_every = max_t %/% 100,
                     by_patch = FALSE,
                     n_cores = 1,
                     show_progress = FALSE,
                     line_names = NULL) {

    stan_fit <- system.file("extdata", "stan_fit.rds", package = "clonewars",
                            mustWork = TRUE)
    stan_fit <- readRDS(stan_fit)

    if (is.null(R)) {
        R <- apply(rstan::extract(stan_fit, "R", permuted = FALSE), 3, median)
    }
    if (is.null(A)) {
        A <- apply(rstan::extract(stan_fit, "A", permuted = FALSE), 3, median)
    }
    n_lines <- length(A)
    if (is.null(N0)) {
        N0 <- matrix(1, n_patches, n_lines)
    } else if (is.numeric(N0)) {
        if (length(N0) != 1) stop("N0 must be a matrix or length-1 vector")
        N0 <- matrix(N0, n_patches, n_lines)
    }
    if (is.null(D_vec)) {
        D_vec <- disp_estimates$b0
    }
    if (is.null(process_error)) {
        process_error <- median(rstan::extract(stan_fit, "sigma_epsilon", permuted = FALSE))
    }
    if (is.null(disp_error)) {
        disp_error <- TRUE
    }
    if (is.null(log_zeta_mean) || is.null(log_zeta_sd)) {
        logZ <- log(apply(rstan::extract(stan_fit, "Z", permuted = FALSE), 3, median))
        if (is.null(log_zeta_mean)) log_zeta_mean <- mean(logZ)
        if (is.null(log_zeta_sd)) log_zeta_sd <- sd(logZ)
    }
    if (is.null(mu_time)) {
        # all estimates are the same, so we can just take the first one:
        mu_time <- rstan::extract(stan_fit, "mu_time")[[1]][[1]]
    }
    if (is.null(repl_times)) repl_times <- as.integer(max_t + 100)

    D_vec <- cbind(D_vec)
    D_vec <- exp(D_vec)

    if (is.null(line_names)) {
        line_names <- clonewars::load_data() %>%
            .[["line"]] %>%
            levels()
    }

    sims <- clonewars:::sim_reps_(n_reps = n_reps,
                      max_t = max_t,
                      N0 = N0,
                      R = R,
                      A = A,
                      D_vec = D_vec,
                      process_error = process_error,
                      disp_error = disp_error,
                      log_zeta_mean = log_zeta_mean,
                      log_zeta_sd = log_zeta_sd,
                      zeta_t_thresh = zeta_t_thresh,
                      mu_time = mu_time,
                      repl_times = repl_times,
                      repl_threshold = repl_threshold,
                      extinct_N = extinct_N,
                      save_every = save_every,
                      by_patch = by_patch,
                      n_cores = n_cores,
                      show_progress = show_progress)

    sims <- as.data.frame(sims)
    sims <- as_tibble(sims)

    if (!by_patch) {
        colnames(sims) <- c("rep", "time", "patch", "line", "N")
    } else colnames(sims) <- c("rep", "time", "patch", "N")

    return(sims)
}




