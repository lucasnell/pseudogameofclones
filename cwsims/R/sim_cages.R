
#' Simplify a single cage's simulation output
#'
#'
#'
#' @noRd
#'
simplify_cage <- function(cage_array, rep, N_0_, max_t_) {
    # Dimensions: n_plants, n_lines, n_times
    dims_ <- c(nrow(N_0_), ncol(N_0_), max_t_ + 1)
    date_ = rep(1:dims_[3] - 1L, each = prod(dims_[1:2]))
    plant_ = rep(1:dims_[1], prod(dims_[2:3]))
    line_ = rep(1:dims_[2], each = dims_[1]) %>% rep(dims_[3])
    data_frame(plant = plant_,
               line = line_,
               date = date_,
               N = as.numeric(cage_array),
               rep = rep) %>%
        arrange(rep, plant, line, date)
}




#' Simulate multiple cages and simplify output.
#'
#'
#' @param n_cages Number of cages to simulate.
#' @param N_0 Starting abundances for each aphid line on each plant.
#' @param max_t Max time points to simulate for each cage.
#' @param R Growth rates for each line.
#' @param A Density dependence for each line.
#' @param D_mat Matrix containing intercept and overdispersion estimates for
#'     quasi-Poisson GLM of `<# dispersed> ~ 1 + offset(log(N))`.
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
sim_cages <- function(n_cages, N_0, max_t, R, A, D_mat, process_error,
                      plant_mort_0, plant_mort_1,
                      plant_death_age_mean, plant_death_age_sd,
                      repl_times, repl_age,
                      n_cores = 1, show_progress = FALSE) {

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
                       n_cores = n_cores,
                       show_progress = show_progress)

    sims <- map2_dfr(sims, 1:length(sims), ~ simplify_cage(.x, .y, N_0, max_t))

    return(sims)
}




