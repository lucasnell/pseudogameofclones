


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
sim_cages <- function(n_cages, N_0, max_t, R, A, D_vec, process_error,
                      plant_mort_0, plant_mort_1,
                      plant_death_age_mean, plant_death_age_sd,
                      repl_times, repl_age, extinct_N, repl_threshold,
                      disp_error = TRUE,
                      plant_death_ages = NULL,
                      n_cores = 1,
                      show_progress = FALSE,
                      line_names = NULL) {

    D_vec <- cbind(D_vec)
    D_vec <- exp(D_vec)

    if (is.null(line_names)) line_names <- 1:length(R)
    if (is.null(plant_death_ages)) plant_death_ages <- integer(0)

    sims <- sim_cages_(n_cages = n_cages,
                       N_0 = N_0,
                       max_t = max_t,
                       R = R,
                       A = A,
                       D_vec,
                       process_error = process_error,
                       disp_error = disp_error,
                       plant_mort_0 = plant_mort_0,
                       plant_mort_1 = plant_mort_1,
                       plant_death_age_mean = plant_death_age_mean,
                       plant_death_age_sd = plant_death_age_sd,
                       plant_death_ages = plant_death_ages,
                       repl_times = repl_times,
                       repl_age = repl_age,
                       extinct_N = extinct_N,
                       repl_threshold = repl_threshold,
                       n_cores = n_cores,
                       show_progress = show_progress)

    sims$N <- map_dfr(
        as.integer(1:dim(sims$N)[3]),
        function(j) {
            sims$N[,,j,drop=FALSE] %>%
                as_data_frame() %>%
                mutate(date = 1:nrow(sims$N)) %>%
                gather("line", "N", -date) %>%
                mutate(line = as.integer(gsub("V", "", line)),
                       pool = as.integer(paste(line_names, collapse = "")),
                       pool_size = as.integer(floor(log10(pool)) + 1L),
                       rep = j) %>%
                mutate(line = line_names[line]) %>%
                arrange(rep, line, date) %>%
                dplyr::select(pool, rep, line, date, N)
        })

    sims$X <- sims$X %>%
        t() %>%
        as_data_frame() %>%
        mutate(rep = 1:ncol(sims$X)) %>%
        gather("plant", "X", -rep) %>%
        mutate(plant = as.integer(gsub("V", "", plant)),
               pool = as.integer(paste(line_names, collapse = "")),
               pool_size = as.integer(floor(log10(pool)) + 1L)) %>%
        arrange(rep, plant) %>%
        dplyr::select(pool, rep, plant, X)

    sims$Z <- sims$Z %>%
        t() %>%
        as_data_frame() %>%
        mutate(rep = 1:ncol(sims$Z)) %>%
        gather("plant", "Z", -rep) %>%
        mutate(plant = as.integer(gsub("V", "", plant)),
               pool = as.integer(paste(line_names, collapse = "")),
               pool_size = as.integer(floor(log10(pool)) + 1L)) %>%
        arrange(rep, plant) %>%
        dplyr::select(pool, rep, plant, Z)


    return(sims)
}




