


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
#' @param D_vec Vector of `b0` values, where the predicted number of dispersed aphids is
#'     given by `exp(b0) * N` when `N` is the number of total aphids.
#' @param process_error SD of process error. Set to 0 for no process error.
#' @param disp_error Boolean for whether to include dispersal stochasticity.
#' @param log_zeta_mean Mean of the distribution of log(zeta) values.
#' @param log_zeta_sd SD of the distribution of log(zeta) values.
#' @param zeta_t_thresh Threshold for `exp(zeta * (t - mu_time))` that makes that
#'     patch get replaced. This is equivalent to the threshold for patch "health"
#'     that would make an experimenter replace it.
#' @param mu_time Mean of time values.
#' @param repl_times Vector of times at which to replace patches.
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

    if (is.null(R)) R <- stan_estimates$R
    if (is.null(A)) A <- stan_estimates$A
    n_lines <- length(A)
    if (is.null(N0)) N0 <- matrix(1, n_patches, n_lines)
    if (inherits(N0, "numeric") || inherits(N0, "integer")) {
        if (length(N0) != 1) stop("N0 must be a matrix or length-1 numeric vector")
        N0 <- matrix(N0, n_patches, n_lines)
    }
    if (is.null(D_vec)) D_vec <- disp_estimates$b0
    if (is.null(process_error)) process_error <- stan_estimates$process_error
    if (is.null(disp_error)) disp_error <- TRUE
    if (is.null(log_zeta_mean)) log_zeta_mean <- stan_estimates$log_zeta_mean
    if (is.null(log_zeta_sd)) log_zeta_sd <- stan_estimates$log_zeta_sd
    if (is.null(mu_time)) mu_time <- stan_estimates$mu_time
    if (is.null(repl_times)) repl_times <- as.integer(max_t + 100)
    if (is.null(line_names)) line_names <- stan_estimates$names
    if (length(R) != n_lines) stop("\nlength(R) != length(A).")
    if (ncol(N0) != n_lines) stop("\nncol(N0) != length(A).")
    if (length(D_vec) != n_lines) stop("\nlength(D_vec) != length(A).")
    if (!by_patch && length(line_names) != n_lines) {
        stop("\nlength(line_names) != length(A).")
    }

    D_vec <- cbind(D_vec)
    D_vec <- exp(D_vec)

    sims <- sim_reps_(n_reps = n_reps,
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
        sims$line <- factor(sims$line, levels = 0:(n_lines-1), labels = line_names)
    } else colnames(sims) <- c("rep", "time", "patch", "N")

    return(sims)
}




