
#' Simplify a single cage's simulation output
#'
#'
#'
#' @noRd
#'
simplify_cage <- function(cage_array, rep, N_0_, max_t_, by_cage_) {
    # Dimensions: n_plants, n_lines, n_times
    dims_ <- c(nrow(N_0_), ncol(N_0_), max_t_ + 1)
    # If this is true, then it's been aggregated by cage, so n_plants is effectively 1
    if (by_cage_) dims_[1] <- 1
    date_ = rep(1:dims_[3] - 1L, each = prod(dims_[1:2]))
    plant_ = rep(1:dims_[1], prod(dims_[2:3]))
    line_ = rep(1:dims_[2], each = dims_[1]) %>% rep(dims_[3])
    out_df <- data_frame(plant = plant_,
               line = line_,
               date = date_,
               N = as.numeric(cage_array),
               rep = rep) %>%
        arrange(rep, line, plant, date)
    if (by_cage_) out_df <- out_df %>% dplyr::select(-plant)
    return(out_df)
}




#' Simulate multiple cages and simplify output.
#'
#'
#' @param n_cages Number of cages to simulate.
#' @param N_0 Starting abundances for each aphid line on each plant.
#' @param max_t Max time points to simulate for each cage.
#' @param R Growth rates for each line.
#' @param A Density dependence for each line.
#' @param D_binom Data frame containing intercept and coefficient estimates for
#'     binomial GLM of `<whether dispersal > 0> ~ N_t + N_t^2`.
#'     There should be one set of estimates for each line regardless of whether the model
#'     estimated separate values for each line.
#' @param D_nb Data frame containing intercept and overdispersion estimates for
#'     negative binomial GLM of `<# dispersed | dispersal occurs> ~ 1`.
#'     There should be one set of estimates for each line regardless of whether the model
#'     estimated separate values for each line.
#' @param process_error SD of process error.
#' @param plant_mort_0 Intercept for plant-death-induced aphid mortality, one per line.
#' @param plant_mort_1 Coefficient for plant-death-induced aphid mortality, one per line.
#' @param plant_death_age_mean Mean of distribution of plant-death ages.
#' @param plant_death_age_sd SD of distribution of plant-death ages.
#' @param repl_times Vector of dates on which to replace plants. Use 0-based indexing!
#' @param repl_age Number of days past death when a plant will be replaced.
#' @param n_cores Number of cores to use. Defaults to \code{1}.
#' @param show_progress Boolean for whether to show progress bar. Defaults to
#'     \code{FALSE}.
#'
#' @export
#'
sim_cages <- function(n_cages, N_0, max_t, R, A, D_binom, D_nb, process_error,
                      plant_mort_0, plant_mort_1,
                      plant_death_age_mean, plant_death_age_sd,
                      repl_times, repl_age, extinct_N,
                      n_cores = 1,
                      by_cage = FALSE,
                      condense = TRUE,
                      show_progress = FALSE) {

    if (!identical(D_binom$line, D_nb$line)) {
        stop("\nline columns should be identical in both D_binom and D_nb.")
    }
    # So I don't have to do this every iteration:
    D_nb$b0 <- exp(D_nb$b0)
    # Combining D_binom and D_nb into one dispersal matrix:
    D_mat <- as.matrix(cbind(D_binom[,c("b0", "b1", "b2")],
                             D_nb[,c("b0", "theta")]))
    colnames(D_mat) <- NULL

    sims <- sim_cages_(n_cages = n_cages,
                       N_0 = N_0,
                       max_t = max_t,
                       R = R,
                       A = A,
                       D_mat,
                       process_error = process_error,
                       plant_mort_0 = plant_mort_0,
                       plant_mort_1 = plant_mort_1,
                       plant_death_age_mean = plant_death_age_mean,
                       plant_death_age_sd = plant_death_age_sd,
                       repl_times = repl_times,
                       repl_age = repl_age,
                       extinct_N = extinct_N,
                       n_cores = n_cores,
                       by_cage = by_cage,
                       show_progress = show_progress)

    # If doing many simulations, you can create data frames with too many rows,
    # that exceed R's limitations
    if (condense) {
        sims <- map2_dfr(sims, 1:length(sims),
                         ~ simplify_cage(.x, .y, N_0, max_t, by_cage))
    } else {
        sims <- map2(sims, 1:length(sims),
                     ~ simplify_cage(.x, .y, N_0, max_t, by_cage))
    }

    return(sims)
}




