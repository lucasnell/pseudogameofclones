# Create line info ----

#' Organize clonal line information
#'
#' To smooth stage structure out over time,
#' not all 4th instar aphids move to adulthood immediately.
#' If adulthood starts on day `t`, then half of aphids at age `t-2` move to
#' adulthood, and half at age `t-1` do, too.
#' I adjusted age `t-2`, too, to avoid this affecting the growth rate too much.
#'
#'
#' @param name String for clonal line name.
#' @param density_0 A single number or 2-column matrix indicating starting
#'     aphid densities.
#'     If a number, there will be no starting alates and a total apterous
#'     density equal to the number provided.
#'     Stages will be calculated based on the stable age distribution.
#'     If a matrix, it must have 5 rows (rows indicate instar)
#'     or a row per aphid age (rows indicate age in days).
#'     Matrix column indicates apterous vs alate.
#' @param resistant Logical or length-2 vector of survivals of
#'     singly attacked and multiply attacked aphids.
#'     If a logical, `FALSE` is equivalent to `c(0,0)` and results in no
#'     resistance.
#'     `TRUE` results in the resistance values for a resistance line
#'     from unpublished work by Anthony Ives.
#'     Defaults to `FALSE`.
#' @param surv_juv_apterous A single number for the juvenile survival rate for
#'     apterous aphids. Defaults to `NULL`, which results in estimates
#'     from a medium-reproduction line.
#' @param surv_adult_apterous A vector of adult survival probabilities for
#'     apterous aphids. Defaults to `NULL`, which results in estimates
#'     from a medium-reproduction line.
#' @param repro_apterous A vector of fecundities for for apterous aphids.
#'     Defaults to `NULL`, which results in estimates from a
#'     medium-reproduction line.
#' @param surv_juv_alates A single number for the juvenile survival rate for
#'     alates aphids. Defaults to `"low"`, which results in estimates
#'     from a low-reproduction line.
#' @param surv_adult_alates A vector of adult survival probabilities for
#'     alates aphids. Defaults to `"low"`, which results in estimates
#'     from a low-reproduction line.
#' @param repro_alates A vector of fecundities for for alates aphids.
#'     Defaults to `"low"`, which results in estimates from a low-reproduction
#'     line.
#' @param surv_paras A single number for the juvenile survival rate for
#'     parasitized aphids.
#'     Because parasitized aphids don't make it to adulthood, this is the only
#'     survival rate necessary.
#'     Defaults to `"low"`, which results in estimates from a
#'     low-reproduction line.
#' @param temp Single string specifying `"low"` (20º C) or `"high"` (27º C)
#'     temperature. Defaults to `"low"`.
#' @param p_instar_smooth A value greater than zero here makes not all 4th
#'     instar aphids go to the adult stage at the same time.
#'     Of the 4th instars that would've moved to adulthood, `p_instar_smooth`
#'     remain as 4th instars.
#'     To keep the growth rate approximately equal, for the aphids that
#'     would've transitioned to 4th instar, `p_instar_smooth` move to
#'     adulthood instead.
#'     Defaults to `0.5`.
#'
#' @return A list with the necessary info to pass onto sim_gameofclones.
#'
#' @export
#'
clonal_line <- function(name,
                        density_0,
                        resistant = FALSE,
                        surv_juv_apterous = NULL,
                        surv_adult_apterous = NULL,
                        repro_apterous = NULL,
                        surv_juv_alates = "low",
                        surv_adult_alates = "low",
                        repro_alates = "low",
                        surv_paras = "low",
                        temp = "low",
                        p_instar_smooth = 0.5) {


    temp <- match.arg(temp, c("low", "high"))
    temp <- paste0(temp, "T")


    # --------------*
    # Construct Leslie matrices
    # --------------*

    leslie <- list(apterous = NA,
                   alates = NA,
                   paras = NA)
    # In `leslie_mat` below, items in vector are aphid lines, slices are
    # apterous/alate/parasitized.

    # Set with default values for Leslie matrix calculations
    def_L_args <- list(instar_days = dev_times$instar_days[[temp]],
                       surv_juv = mean(do.call(c, populations$surv_juv)),
                       surv_adult = colMeans(do.call(rbind,
                                                     populations$surv_adult)),
                       repro = colMeans(do.call(rbind, populations$repro)))

    # Note that we don't have values for parasitized aphid adult survival
    # or fecundity because they don't reach adulthood anyway.
    # Below, the defaults will be filled in for these, then the downstream
    # code will ignore that info.
    inputs <- list(surv_juv_apterous,
                   surv_juv_alates,
                   surv_paras,
                   surv_adult_apterous,
                   surv_adult_alates,
                   repro_apterous,
                   repro_alates)
    names(inputs) <- c("surv_juv_apterous", "surv_juv_alates", "surv_juv_paras",
                       "surv_adult_apterous", "surv_adult_alates",
                       "repro_apterous", "repro_alates")

    for (x in c("apterous", "alates", "paras")) {
        leslie_args <- def_L_args
        for (y in names(inputs)[grepl(paste0(x, "$"), names(inputs))]) {
            arg_name <- gsub(paste0("_", x), "", y)
            if (!is.null(inputs[[y]])) {
                if (is.numeric(inputs[[y]])) {
                    leslie_args[[arg_name]] <- inputs[[y]]
                } else if (inputs[[y]] %in% c("low", "high")) {
                    leslie_args[[arg_name]] <- populations[[arg_name]][[
                        inputs[[y]]]]
                } else {
                    msg <- paste0("\nERROR: input argument `", y,
                                  "` to the `clonal_line` function should be ",
                                  "NULL, a numeric vector, \"low\", or ",
                                  "\"high\"")
                    stop(msg)
                }
            }
        }
        leslie[[x]] <- do.call(leslie_matrix, leslie_args)
    }

    if (!identical(dim(leslie[[1]]), dim(leslie[[2]])) ||
        !identical(dim(leslie[[1]]), dim(leslie[[3]]))) {
        stop("\nERROR: Leslie matrices for apterous, alates, and parasitized",
             " aphids must all be of the same dimensions\n")
    }

    leslie_array <- array(do.call(c, leslie), dim = c(dim(leslie[[1]]), 3))
    ns <- nrow(leslie_array)  # number of stages; used later

    # To make 4th instars not always move to adulthood at the same time:
    if (p_instar_smooth > 0) {
        .adult <- sum(head(dev_times$instar_days[[temp]], -1)) + 1
        # I included `- 1` bc we don't to adjust the Leslie matrix for the
        # parasitized aphids
        for (j in 1:(dim(leslie_array)[3] - 1)) {
            # Of the aphids that would've moved to adulthood, make
            # `p_instar_smooth` remain as 4th instars of age `.adult-1` instead:
            .t <- .adult - 1
            surv_t <- leslie_array[.adult, .t, j]
            leslie_array[.t, .t, j] <- surv_t * p_instar_smooth
            leslie_array[.adult, .t, j] <- surv_t * (1 - p_instar_smooth)
            # Of the aphids that would've moved to age `.adult-1`, make
            # `p_instar_smooth` move to adulthood instead:
            .t <- .adult - 2
            surv_t <- leslie_array[.t+1, .t, j]
            leslie_array[.adult, .t, j] <- surv_t * p_instar_smooth
            leslie_array[.t+1, .t, j] <- surv_t * (1 - p_instar_smooth)
        }
    }



    # --------------*
    # Fill other info
    # --------------*

    attack_surv <- c(0, 0)
    if (is.logical(resistant) && resistant) attack_surv <- wasp_attack$attack_surv
    if (is.numeric(resistant) && length(resistant) == 2) attack_surv <- resistant


    d0_dbl <- is.numeric(density_0) && length(density_0) == 1
    d0_m52 <- is.numeric(density_0) && is.matrix(density_0) &&
        identical(dim(density_0), c(5L, 2L))
    d0_mn2 <- is.numeric(density_0) && is.matrix(density_0) &&
        identical(dim(density_0), c(ns, 2L))
    if (!(d0_dbl || d0_m52 || d0_mn2)) {
        stop("\nERROR: The `density_0` arg to the ",
             "`clonal_line` function must be a single number or a 5x2 or ",
             ns, "x2 numeric matrix")
    }
    if (d0_dbl) {
        # If provided with just one number, assume stable age distribution:
        sad <- as.numeric(sad_leslie(leslie_array[,,1]))
        density_0 <- cbind(sad * density_0, rep(0, length(sad)))
    } else if (d0_m52) {
        d0 <- density_0
        density_0 <- matrix(0, ns, 2)
        # Going from instar to days old, using stable age distribution to
        # calculate the proportion for each age within each instar.
        sads <- cbind(as.numeric(sad_leslie(leslie_array[,,1])),
                      as.numeric(sad_leslie(leslie_array[,,2])))
        idays <- dev_times$instar_days[[temp]]
        d0_inds <- cbind(c(1, head(cumsum(idays), -1) + 1),
                         c(head(cumsum(idays), -1), nrow(sads)))
        for (i in 1:nrow(d0)) {
            irange <- d0_inds[i,1]:d0_inds[i,2]
            for (j in 1:2) {
                sad_ij = sads[irange, j] / sum(sads[irange,j])
                density_0[irange,j] <- sad_ij * d0[i,j]
            }
        }
    }
    if (any(density_0 < 0)) {
        stop("\nERROR: The `density_0` arg to the `clonal_line` function ",
             "cannot contain values < 0.")
    }


    output <- list(name = name,
                   density_0 = density_0,
                   attack_surv = attack_surv,
                   leslie = leslie_array)
    class(output) <- "aphid"

    return(output)

}


#'
#' @export
#' @noRd
#'
print.aphid <- function(x, ...) {

    cat("< Aphid clonal line >\n")
    cat("Name: ", x$name, "\n", sep = "")
    cat("Fields:\n")
    cat("  * name <string>\n")
    cat("  * density_0 <matrix>\n")
    cat("  * attack_surv <vector>\n")
    cat("  * leslie <3D array>\n")

    invisible(x)

}

#'
#' @export
#' @noRd
#'
c.aphid <- function(...) {
    z <- list(...)
    names(z) <- sapply(z, function(x) x$name)
    class(z) <- "multiAphid"
    return(z)
}

#'
#' @export
#' @noRd
#'
print.multiAphid <- function(x, ...) {
    cat("< ", length(x), "aphid clonal lines >\n")
    cat("Lines:\n")
    for (i in 1:length(x)) {
        cat("  ", i, ". ", x[[i]]$name, sep = "")
        if (sum(x[[i]][["attack_surv"]]) > 0) cat(" (resistant)")
        cat("\n")
    }
    invisible(x)
}

# Basic type checks ----
uint_check <- function(x, n, .max = NULL, .min = NULL) {
    if (!(is.numeric(x) && length(x) == 1 && x %% 1 == 0 && x >= 0)) {
        stop(paste("\nERROR:", n, "cannot be properly cast as an",
                   "unsigned integer.\n"))
    }
    if (!is.null(.min) && x < .min) {
        stop(paste0("\nERROR: ", n, " is below the minimum allowed value (",
                    .min, ").\n"))
    }
    if (!is.null(.max) && x > .max) {
        stop(paste0("\nERROR: ", n, " is above the maximum allowed value (",
                    .max, ").\n"))
    }
}
uint_vec_check <- function(x, n, .max = NULL, .min = NULL) {
    if (!(is.numeric(x) && all(x %% 1 == 0) && all(x >= 0))) {
        stop(paste("\nERROR:", n, "cannot be properly cast as an",
                   "unsigned integer vector.\n"))
    }
    if (!is.null(.min) && any(x < .min)) {
        stop(paste0("\nERROR: ", n, " contains values below the minimum ",
                    "allowed (", .min, ").\n"))
    }
    if (!is.null(.max) && any(x > .max)) {
        stop(paste0("\nERROR: ", n, " contains values above the maximum ",
                    "allowed (", .max, ").\n"))
    }
}
dbl_check <- function(x, n, .max = NULL, .min = NULL) {
    if (!(is.numeric(x) && length(x) == 1)) {
        stop(paste("\nERROR:", n, "cannot be properly cast as a",
                   "double.\n"))
    }
    if (!is.null(.min) && x < .min) {
        stop(paste0("\nERROR: ", n, " is below the minimum allowed value (",
                   .min, ").\n"))
    }
    if (!is.null(.max) && x > .max) {
        stop(paste0("\nERROR: ", n, " is above the maximum allowed value (",
                   .max, ").\n"))
    }
}
dbl_vec_check <- function(x, n, .max = NULL, .min = NULL) {
    if (!(is.numeric(x) && is.null(dim(x)))) {
        stop(paste("\nERROR:", n, "cannot be properly cast as a",
                   "numeric vector.\n"))
    }
    if (!is.null(.min) && any(x < .min)) {
        stop(paste0("\nERROR: ", n, " contains values below the minimum ",
                    "allowed (", .min, ").\n"))
    }
    if (!is.null(.max) && any(x > .max)) {
        stop(paste0("\nERROR: ", n, " contains values above the maximum ",
                    "allowed (", .max, ").\n"))
    }
}
dbl_mat_check <- function(x, n, .max = NULL, .min = NULL) {
    if (!(is.numeric(x) && inherits(x, "matrix"))) {
        stop(paste("\nERROR:", n, "cannot be properly cast as a",
                   "numeric matrix.\n"))
    }
    if (!is.null(.min) && any(x < .min)) {
        stop(paste0("\nERROR: ", n, " contains values below the minimum ",
                    "allowed (", .min, ").\n"))
    }
    if (!is.null(.max) && any(x > .max)) {
        stop(paste0("\nERROR: ", n, " contains values above the maximum ",
                    "allowed (", .max, ").\n"))
    }
}
cube_list_check <- function(x, n) {
    if (!(inherits(x, "list") &&
          all(sapply(x, inherits, what = "array")) &&
          all(sapply(x, function(y) length(dim(y)) == 3)))) {
        stop(paste("\nERROR:", n, "cannot be properly cast as a",
                   "list of cubes.\n"))
    }
}

# full fun docs ----
#' Simulate multiple reps and simplify output - all options.
#'
#' This mainly differs from `sim_gameofclones` in that it allows for
#' stochasticity and for simulating the process of plants dying and
#' being replaced.
#' The latter process doesn't appear important except for very small scales
#' and for defenseless plants (like fava bean).
#' Stochasticity isn't necessarily important, either, and the code has
#' change in other ways such that I have largely abandoned the
#' stochastic simulations.
#' So use this at your own risk.
#'
#'
#' @param n_reps Number of reps to simulate.
#' @param clonal_lines A `multiAphid` object containing the aphid-line-specific
#'     info for all the lines in the simulations.
#'     Each object in the `multiAphid` is an `aphid` object that results from
#'     the `clonal_line` function.
#'     Combine them using `c(aphid_obj1, aphid_obj2)`.
#' @param n_fields The number of fields to simulate.
#'     Wasps operate at the field level, whereas aphids operate at the
#'     plant level.
#'     Both wasps and aphids operate separately across fields but can be
#'     connected via dispersal (see arguments `wasp_disp_p` and
#'     `alate_field_disp_p`).
#'     Defaults to `2`.
#' @param n_plants The number of plants per field.
#'     This argument is only useful if you want to simulate the process
#'     of plants or groups of plants dying.
#'     Wasps operate at the field level, whereas aphids operate at the
#'     plant level.
#'     "Plant" can indicate a group of plants where the aphids can freely
#'     disperse across them.
#'     Aphids can disperse across plants via alate production and the
#'     `alate_plant_disp_p` argument.
#'     Defaults to `1`.
#' @param max_t How many days to simulate. Defaults to `250`.
#' @param plant_check_gaps Gap(s) between when plants are check on.
#'     This is used if you want to see what will happen if you can't check
#'     on things every day in an experiment.
#'     You could use `c(3, 4)` for checking twice per week, for example.
#'     Note that this argument determines how often
#'     (1) wasps and alates disperse across fields and
#'     (2) plants are check for exceeding `max_plant_age` or `max_N`.
#'     So using this argument to simulate harvesting fields will also cause
#'     wasp and alate dispersal to not occur daily.
#'     Defaults to `1`.
#' @param max_plant_age Age at which plants are cleared.
#'     This can be useful to simulate harvesting fields.
#'     A value of `0` turns this off completely, but it also requires that
#'     `max_N > 0` because the internal code wants some threshold to look for.
#'     Defaults to `0` which effectively turns this off.
#' @param clear_surv Survival of aphids and mummies when a harvest occurs.
#'     Defaults to `0`.
#' @param max_N Maximum number of aphids at which plants are cleared.
#'     This is useful for simulating an experiment where you have a threshold
#'     at which you harvest plants.
#'     Defaults to `0`, which turns this off.
#' @param temp A string indicating which temperature to simulate.
#'     Options are `"low"` (20ºC) or `"high"` (27ºC).
#'     Defaults to `"low"`.
#' @param no_error Logical for whether to have no error at all.
#'     This being `TRUE` overrides all other `_error` arguments.
#'     Defaults to `TRUE`.
#' @param disp_error Logical for whether to have stochasticity in the
#'     dispersal process.
#'     Defaults to `FALSE`.
#' @param environ_error Logical for whether to have environmental stochasticity
#'     Defaults to `FALSE`.
#' @param plant_K_error Logical for whether to have stochasticity in aphid
#'     density dependence across plants.
#'     Defaults to `FALSE`.
#' @param wilted_effects_error Logical for whether to have stochasticity in
#'     the effects of plants dying on aphid populations.
#'     Defaults to `FALSE`.
#' @param sigma_x Standard deviation of environmental stochasticity for aphids.
#'     Defaults to the internal object `environ$sigma_x`,
#'     which is from Meisner et al. (2014).
#' @param sigma_y Standard deviation of environmental stochasticity for wasps.
#'     Defaults to the internal object `environ$sigma_y`,
#'     which is from Meisner et al. (2014).
#' @param mean_K Mean of the distribution of `K` values (affecting aphid
#'     density dependence) among plants.
#'     Defaults to `1 / 4.67e-4`, which is from Meisner et al. (2014).
#' @param sd_K Mean of the distribution of `K` values (affecting aphid
#'     density dependence) among plants. Defaults to `0`.
#' @param K_y_mult The number multiplied by `K` to get density dependence for
#'     parasitized aphids.
#'     Defaults to `1 / 1.57`, which is from Meisner et al. (2014).
#' @param wilted_prop Proportion of carrying capacity that causes the plant
#'     to become wilted. Values > 1 cause this to be ignored.
#'     Defaults to `1.1`.
#' @param rho Environmental correlation among instars.
#'     Defaults to `environ$rho`.
#' @param a Parasitoid attack rate. Defaults to the internal object
#'     `wasp_attack$a`, which is from Meisner et al. (2014).
#' @param k Aggregation parameter of the negative binomial distribution.
#'     Defaults to the internal object `wasp_attack$k`,
#'     which is from Meisner et al. (2014).
#' @param h Parasitoid handling time. Defaults to the internal object
#'     `wasp_attack$h`, which is from Meisner et al. (2014).
#' @param wasp_density_0 Starting adult wasp density.
#'     Must be a single number or a `n_fields`-length vector.
#'     Defaults to `c(3, 0)`.
#' @param wasp_delay Delay in days between when the aphids start and the
#'     wasps are added.
#'     This can be a single integer or a `n_fields`-length vector.
#'     Defaults to `8`.
#' @param wasp_disp_p Proportion of adult wasps from each field that
#'     are added to the dispersal pool.
#'     After adding wasps to the pool, they are then evenly distributed
#'     to all fields.
#'     This happens only on days indicated by `plant_check_gaps`.
#'     Defaults to `0`.
#' @param sex_ratio Sex ratio of adult wasps. Defaults to `0.5`.
#' @param s_y Daily survival rate of adult wasps.
#'     Defaults to `populations$s_y`, which is from Meisner et al. (2014).
#' @param constant_wasps Logical for whether to keep adult wasps at the same
#'     density throughout simulations.
#'     This can be a single logical or a `n_fields`-length vector.
#'     Defaults to `FALSE`.
#' @param rel_attack Relative parasitoid attack rate among instars.
#'     Defaults to `wasp_attack$rel_attack`, which is from
#'     Meisner et al. (2014).
#' @param mum_density_0 Starting mummy density. Defaults to `0`.
#' @param mum_smooth Proportion of mummies that will NOT take exactly 3 days
#'     to develop. As this value approaches 2/3, it will provide greater
#'     smoothing of wasp numbers through time.
#'     Defaults to `0.4`.
#' @param max_mum_density Maximum mummy density (ignored if zero).
#'     Used to test the effects of a potential experimental treatment.
#'     Defaults to `0`.
#' @param pred_rate Daily predation rate on aphids and mummies. Defaults to `0`.
#' @param alate_plant_disp_p Proportion of alates that emigration to new plants
#'     every day. Note that this does not depend on `plant_check_gaps`
#'     like among-field alate dispersal does.
#'     This is because we would have to replicate among-field dispersal
#'     by hand in the experiments we were planning, and we wanted to be sure
#'     that not checking daily wouldn't affect anything significantly.
#'     It didn't.
#'     This can be a single number or a vector with the same length as the
#'     number of aphid lines.
#'     Defaults to `0.1`.
#' @param disp_mort Mortality rate of alates that disperse among plants.
#'     Defaults to `0`.
#' @param alate_b0 The proportion of offspring from apterous aphids is
#'     `inv_logit(alate_b0 + alate_b1 * N)` where `N` is the total number of
#'     aphids on that plant.
#'     Defaults to `logit(0.3)`.
#' @param alate_b1 The proportion of offspring from apterous aphids is
#'     `inv_logit(alate_b0 + alate_b1 * N)` where `N` is the total number of
#'     aphids on that plant.
#'     Defaults to `0`, which makes alate production not density dependent.
#' @param alate_field_disp_p Proportion of alates from each field that
#'     are added to the dispersal pool.
#'     After adding alates to the pool, they are then evenly distributed
#'     to all fields.
#'     This happens only on days indicated by `plant_check_gaps`.
#'     Defaults to `0.1`.
#' @param shape1_wilted_mort Shape 1 for the beta distribution that generates
#'     mortality parameters for aphid populations living on a wilted plant.
#'     If `wilted_effects_error` is `FALSE`, then the wilted-plant-induced
#'     mortality is fixed at
#'     `shape1_wilted_mort / (shape1_wilted_mort + shape2_wilted_mort)`.
#'     Defaults to `3.736386`, which is based on previous work in the lab.
#' @param shape2_wilted_mort Shape 2 for the beta distribution that generates
#'     mortality parameters for aphid populations living on a wilted plant.
#'     Defaults to `5.777129`, which is based on previous work in the lab.
#' @param extinct_N Threshold below which a line is considered extinct.
#'     Defaults to `1`.
#' @param save_every Abundances will be stored every `save_every` time points.
#'     Defaults to `1`.
#' @param n_threads Number of threads to use if OpenMP is enabled
#'     (ignored otherwise).
#'     Find out whether it's enabled using `gameofclones:::using_openmp()`.
#'     Defaults to `1`.
#' @param show_progress Boolean for whether to show progress bar.
#'     Defaults to `FALSE`.
#' @param perturb Information for perturbing populations in the simulations.
#'     It should be a dataframe with 4 columns:
#'     * `when`: Integers indicating at what timepoint(s) to do the
#'       perturbations. These can be repeated if you want to perturb
#'       multiple things at the same time.
#'     * `where`: What field to do the perturbations in.
#'     * `who`: Which population to perturb.
#'       This can be a character vector where values must be the name of an
#'       aphid line, `"wasps"`, or `"mummies"`.
#'       It can also be an integer vector where, for `n` aphid lines,
#'       values `<= n` indicate an aphid line,
#'       values `== n+1` indicate mummies,
#'       and values `== n+2` indicate adult wasps.
#'       Note that perturbing the mummy population also perturbs the
#'       still-living but parasitized aphids, too.
#'     * `how`: Numbers `>= 0` that are multiplied by the desired population
#'       to cause the perturbation.
#'
#' @importFrom purrr map_dfr
#' @importFrom dplyr as_tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom dplyr arrange
#'
#' @export
#'
# full fun code ----
sim_gameofclones_full <- function(n_reps,
                               clonal_lines,
                               n_fields = 2,
                               n_plants = 1,
                               max_t = 250,
                               plant_check_gaps = 1,
                               max_plant_age = 0,
                               clear_surv = 0,
                               max_N = 0,
                               temp = "low",
                               no_error = TRUE,
                               disp_error = FALSE,
                               environ_error = FALSE,
                               plant_K_error = FALSE,
                               wilted_effects_error = FALSE,
                               sigma_x = environ$sigma_x,
                               sigma_y = environ$sigma_y,
                               mean_K = 1 / 4.67e-4,
                               sd_K = 0,
                               K_y_mult = 1 / 1.57,
                               wilted_prop = 1.1,
                               rho = environ$rho,
                               a = wasp_attack$a,
                               k = wasp_attack$k,
                               h = wasp_attack$h,
                               wasp_density_0 = c(3, 0),
                               wasp_delay = 8,
                               wasp_disp_p = 0,
                               sex_ratio = populations$sex_ratio,
                               s_y = populations$s_y,
                               constant_wasps = FALSE,
                               rel_attack = wasp_attack$rel_attack,
                               mum_density_0 = 0,
                               mum_smooth = 0.4,
                               max_mum_density = 0,
                               pred_rate = 0,
                               alate_plant_disp_p = 0.1,
                               disp_mort = 0,
                               alate_b0 = logit(0.093),
                               alate_b1 = 0,
                               alate_field_disp_p = 0.1,
                               shape1_wilted_mort = 3.736386,
                               shape2_wilted_mort = 5.777129,
                               extinct_N = 1,
                               save_every = 1,
                               n_threads = 1,
                               show_progress = FALSE,
                               perturb = NULL) {

    # This is like match.call but includes default arguments and
    # evaluates everything
    call_ <- mget(names(formals()), sys.frame(sys.nframe()))

    if (!inherits(clonal_lines, "multiAphid")) {
        if (inherits(clonal_lines, "aphid")) {
            clonal_lines <- c(clonal_lines)
        } else stop("\nERROR: `clonal_lines` must be a multiAphid object.\n")
    }

    temp <- match.arg(temp, c("low", "high"))

    n_lines <- length(clonal_lines)

    check_for_clear <- cumsum(rep(plant_check_gaps,
                                  ceiling(max_t / sum(plant_check_gaps))))
    check_for_clear <- check_for_clear[check_for_clear <= max_t]


    if (no_error) {
        disp_error <- FALSE
        environ_error <- FALSE
        plant_K_error <- FALSE
        wilted_effects_error <- FALSE
    }
    if (!environ_error) {
        demog_error <- FALSE
        sigma_x <- 0
        sigma_y <- 0
    } else demog_error <- TRUE
    if (!plant_K_error) sd_K <- 0

    if (!wilted_effects_error) {
        if (shape2_wilted_mort > 0) {
            shape1_wilted_mort <- shape1_wilted_mort /
                (shape1_wilted_mort + shape2_wilted_mort)
            shape2_wilted_mort <- 0
        }
    }

    if (length(pred_rate) == 1) pred_rate <- rep(pred_rate, n_fields)
    if (length(alate_plant_disp_p) == 1) alate_plant_disp_p <-
        rep(alate_plant_disp_p, n_lines)
    if (length(disp_mort) == 1) disp_mort <- rep(disp_mort, n_lines)
    if (length(alate_b0) == 1) alate_b0 <- rep(alate_b0, n_lines)
    if (length(alate_b1) == 1) alate_b1 <- rep(alate_b1, n_lines)
    if (length(K_y_mult) == 1) K_y_mult <- rep(K_y_mult, n_fields)
    if (length(wasp_density_0) == 1) wasp_density_0 <- rep(wasp_density_0,
                                                           n_fields)
    if (length(s_y) == 1) s_y <- rep(s_y, n_fields)
    if (length(constant_wasps) == 1) constant_wasps <- rep(constant_wasps,
                                                           n_fields)
    if (length(wasp_delay) == 1) wasp_delay <- rep(wasp_delay, n_fields)

    stopifnot(length(rel_attack) == 5)

    if (length(mum_density_0) == 1) {
        mum_density_0 <- matrix(mum_density_0, dev_times$mum_days[2], n_plants)
    }


    living_days <- rep(dev_times$mum_days[[1]], n_lines)
    disp_start <- rep(sum(head(dev_times$instar_days[[paste0(temp, "T")]], -1)),
                      n_lines)



    # ---------------*
    # Extract from `clonal_lines`
    # ---------------*
    aphid_names <- names(clonal_lines)

    densities_0 <- lapply(clonal_lines, function(x) x$density_0)
    if (!length(unique(sapply(densities_0, nrow))) == 1) {
        stop("\nERROR: All aphid lines must have the same sized density ",
             "matrices\n")
    }
    aphid_density_0 <- array(do.call(c, densities_0),
                             dim = c(dim(densities_0[[1]]), n_lines))
    aphid_density_0 <- replicate(n_plants, aphid_density_0, simplify = FALSE)

    attack_surv <- do.call(cbind, lapply(clonal_lines, `[[`, i = "attack_surv"))

    leslie_cubes <- lapply(clonal_lines, function(x) x$leslie)

    stopifnot(length(unique(sapply(leslie_cubes, nrow))) == 1)
    stopifnot(length(unique(sapply(leslie_cubes, ncol))) == 1)
    stopifnot(length(unique(sapply(leslie_cubes, function(x) dim(x)[3]))) == 1)

    if (length(rel_attack) == 5) {
        rel_attack <- rel_attack / sum(rel_attack)
        dt <- dev_times$instar_days[[paste0(temp, "T")]]
        n_adult_days <- nrow(leslie_cubes[[1]]) - sum(head(dt, -1))
        stopifnot(n_adult_days >= 0)
        if (tail(dt, 1) != n_adult_days) dt[length(dt)] <- n_adult_days
        # Commented version isn't used bc it wasn't done this way when fitting
        # the model.
        # rel_attack__ <- mapply(function(.x, .y) rep(.x, .y) / .y,
        #                        rel_attack,  dt)
        rel_attack__ <- mapply(rep, rel_attack,  dt)
        rel_attack <- do.call(c, rel_attack__)
    } else stopifnot(length(rel_attack) == nrow(leslie_cubes[[1]]))

    if (is.null(perturb)) {
        perturb_when = integer(0)
        perturb_where = integer(0)
        perturb_who = integer(0)
        perturb_how = numeric(0)
    } else {
        stopifnot(inherits(perturb, "data.frame"))
        stopifnot(all(c("when", "where", "who", "how") %in% colnames(perturb)))
        perturb <- dplyr::arrange(perturb, when)
        perturb_when <- perturb$when
        stopifnot(all(perturb$where <= n_fields & perturb$where > 0))
        perturb_where <- perturb$where - 1
        perturb_how <- perturb$how
        if (is.character(perturb$who)) {
            stopifnot(all(perturb$who %in% c(aphid_names, "mummies", "wasps")))
            perturb_who <- integer(length(perturb$who))
            perturb_who[perturb$who %in% aphid_names] <- -1 +
                match(perturb$who[perturb$who %in% aphid_names], aphid_names)
            perturb_who[perturb$who == "mummies"] <- length(aphid_names)
            perturb_who[perturb$who == "wasps"] <- length(aphid_names) + 1
        } else {
            stopifnot(all(perturb$who <= (n_lines+2) & perturb$who > 0))
            perturb_who <- perturb$who - 1
        }
    }

    uint_check(n_reps, "n_reps", .min = 1)
    uint_check(n_fields, "n_fields", .min = 1)
    uint_check(max_plant_age, "max_plant_age")
    dbl_check(max_N, "max_N")
    uint_vec_check(check_for_clear, "check_for_clear")
    dbl_check(clear_surv, "clear_surv")
    uint_check(max_t, "max_t", .min = 1)
    uint_check(save_every, "save_every", .min = 1)
    dbl_check(mean_K, "mean_K")
    dbl_check(sd_K, "sd_K")
    dbl_vec_check(K_y_mult, "K_y_mult")
    dbl_check(wilted_prop, "wilted_prop")
    dbl_check(shape1_wilted_mort, "shape1_wilted_mort")
    dbl_check(shape2_wilted_mort, "shape2_wilted_mort")
    dbl_mat_check(attack_surv, "attack_surv")
    stopifnot(inherits(disp_error, "logical") && length(disp_error) == 1)
    stopifnot(inherits(demog_error, "logical") && length(demog_error) == 1)
    dbl_check(sigma_x, "sigma_x")
    dbl_check(sigma_y, "sigma_y")
    dbl_check(rho, "rho")
    dbl_check(extinct_N, "extinct_N")
    stopifnot(inherits(aphid_names, "character"))
    cube_list_check(leslie_cubes, "leslie_cubes")
    cube_list_check(aphid_density_0, "aphid_density_0")
    stopifnot(inherits(alate_b0, "numeric"))
    stopifnot(inherits(alate_b1, "numeric"))
    dbl_check(alate_field_disp_p, "alate_field_disp_p", .min = 0, .max = 1)
    stopifnot(inherits(alate_plant_disp_p, "numeric"))
    stopifnot(inherits(disp_mort, "numeric"))
    uint_vec_check(disp_start, "disp_start")
    uint_vec_check(living_days, "living_days")
    stopifnot(inherits(pred_rate, "numeric"))
    dbl_mat_check(mum_density_0, "mum_density_0")
    dbl_check(mum_smooth, "mum_smooth", .min = 0, .max = 1)
    dbl_check(max_mum_density, "max_mum_density", .min = 0)
    stopifnot(inherits(rel_attack, "numeric"))
    dbl_check(a, "a")
    dbl_check(k, "k")
    dbl_check(h, "h")
    dbl_vec_check(wasp_density_0, "wasp_density_0")
    uint_vec_check(wasp_delay, "wasp_delay")
    dbl_check(wasp_disp_p, "wasp_disp_p", .min = 0, .max = 1)
    dbl_check(sex_ratio, "sex_ratio")
    dbl_vec_check(s_y, "s_y")
    stopifnot(inherits(constant_wasps, "logical"))
    uint_vec_check(perturb_when, "perturb_when")
    uint_vec_check(perturb_where, "perturb_where")
    uint_vec_check(perturb_who, "perturb_who")
    dbl_vec_check(perturb_how, "perturb_how", .min = 0)
    uint_check(n_threads, "n_threads", .min = 1)
    stopifnot(inherits(show_progress, "logical") && length(show_progress) == 1)


    sims <- sim_gameofclones_cpp(n_reps, n_fields, max_plant_age, max_N,
                              check_for_clear, clear_surv,
                              max_t, save_every, mean_K, sd_K, K_y_mult,
                              wilted_prop, shape1_wilted_mort,
                              shape2_wilted_mort,
                              attack_surv, disp_error, demog_error, sigma_x,
                              sigma_y, rho, extinct_N, aphid_names,
                              leslie_cubes,
                              aphid_density_0, alate_b0, alate_b1,
                              alate_field_disp_p,
                              alate_plant_disp_p, disp_mort, disp_start,
                              living_days, pred_rate,
                              mum_density_0, mum_smooth, max_mum_density,
                              rel_attack, a, k, h,
                              wasp_density_0, wasp_delay, wasp_disp_p,
                              sex_ratio, s_y, constant_wasps,
                              perturb_when, perturb_where,
                              perturb_who, perturb_how,
                              n_threads, show_progress)

    sims[["aphids"]] <- sims[["aphids"]] %>%
        mutate(across(c("rep", "time", "field", "plant"), as.integer)) %>%
        mutate(across(c("rep", "field", "plant"), ~ .x + 1L)) %>%
        mutate(line = ifelse(type == "mummy", NA_character_, line)) %>%
        as_tibble()
    sims[["wasps"]] <- sims[["wasps"]] %>%
        mutate(across(c("rep", "time", "field"), as.integer)) %>%
        mutate(across(c("rep", "field"), ~ .x + 1L)) %>%
        as_tibble()

    sims[["all_info"]] <- make_all_info(sims)

    sims[["call"]] <- call_

    class(sims) <- "cloneSims"

    return(sims)
}






#' Relatively simple simulations to set up experiments.
#'
#' Deterministic simulations of multiple fields of aphids and wasps.
#'
#'
#' @inheritParams sim_gameofclones_full
#'
#' @param K Aphid density dependence.
#'     Defaults to `12.5e3` because this caused simulations to
#'     approximately match experiments.
#' @param alate_b0 The proportion of offspring from apterous aphids is
#'     `inv_logit(alate_b0 + alate_b1 * N)` where `N` is the total number of
#'     aphids on that plant.
#'     Defaults to `-5` because we're making this weakly density dependent.
#' @param alate_b1 The proportion of offspring from apterous aphids is
#'     `inv_logit(alate_b0 + alate_b1 * N)` where `N` is the total number of
#'     aphids on that plant.
#'     Defaults to `K * 1.76e-07`, which makes alate production only mildly
#'     density dependent.
#' @param pred_rate Daily predation rate on aphids and mummies.
#'     Defaults to `0.1` to compensate for losses from plants dying.
#'
#'
#' @export
#'
# main fun code ----
sim_experiments <- function(clonal_lines,
                            n_fields = 2,
                            max_t = 250,
                            K = 12.5e3,
                            alate_b0 = -5,
                            alate_b1 = K * 1.76e-07,
                            alate_field_disp_p = 0.1,
                            K_y_mult = 1 / 1.57,
                            s_y = populations$s_y,
                            a = wasp_attack$a,
                            k = wasp_attack$k,
                            h = wasp_attack$h,
                            rel_attack = wasp_attack$rel_attack,
                            wasp_density_0 = c(3, 0),
                            wasp_delay = 8,
                            wasp_disp_p = 0,
                            constant_wasps = FALSE,
                            mum_smooth = 0.4,
                            pred_rate = 0.1,
                            extinct_N = 1,
                            save_every = 1,
                            perturb = NULL,
                            plant_check_gaps = 1,
                            max_plant_age = 0,
                            clear_surv = 0,
                            show_progress = FALSE) {


    sims <- sim_gameofclones_full(clonal_lines = clonal_lines,
                               n_fields = n_fields,
                               max_t = max_t,
                               plant_check_gaps = plant_check_gaps,
                               max_plant_age = max_plant_age,
                               clear_surv = clear_surv,
                               mean_K = K,
                               K_y_mult = K_y_mult,
                               a = a,
                               k = k,
                               h = h,
                               wasp_density_0 = wasp_density_0,
                               wasp_delay = wasp_delay,
                               wasp_disp_p = wasp_disp_p,
                               s_y = s_y,
                               constant_wasps = constant_wasps,
                               rel_attack = rel_attack,
                               mum_smooth = mum_smooth,
                               pred_rate = pred_rate,
                               alate_b0 = alate_b0,
                               alate_b1 = alate_b1,
                               alate_field_disp_p = alate_field_disp_p,
                               extinct_N = extinct_N,
                               save_every = save_every,
                               perturb = perturb,
                               show_progress = show_progress,
                               # Things not changeable in this simpler version:
                               n_reps = 1,
                               n_plants = 1,
                               max_N = 0,
                               temp = "low",
                               no_error = TRUE,
                               disp_error = FALSE,
                               environ_error = FALSE,
                               plant_K_error = FALSE,
                               wilted_effects_error = FALSE,
                               sigma_x = environ$sigma_x,
                               sigma_y = environ$sigma_y,
                               sd_K = 0,
                               wilted_prop = 1.1,
                               rho = environ$rho,
                               sex_ratio = populations$sex_ratio,
                               mum_density_0 = 0,
                               max_mum_density = 0,
                               alate_plant_disp_p = 0.1,
                               disp_mort = 0,
                               shape1_wilted_mort = 3.736386,
                               shape2_wilted_mort = 5.777129,
                               n_threads = 1)

    # Adjusting for different name from `sim_gameofclones_full`:
    sims$call[["K"]] <- sims$call[["mean_K"]]
    sims$call[["mean_K"]] <- NULL

    return(sims)

}





# Internal function to create the list of dataframes containing all info.
make_all_info <- function(sims_obj) {
    if (!inherits(sims_obj$all_info_xptr, "externalptr")) {
        stop("\nSomething has happened to the \"all_info_xptr\" field ",
             "in this `cloneSims` object. Please re-run simulation.")
    }
    all_info <- lapply(fields_to_list(sims_obj$all_info_xptr),
                                as.data.frame)
    for (i in 1:length(all_info)) {
        is_wasp <- all_info[[i]][["type"]] == "wasp"
        is_mummy <- all_info[[i]][["type"]] == "mummy"
        all_info[[i]][["plant"]][is_wasp] <- NA
        all_info[[i]][["line"]][is_wasp | is_mummy] <- NA
    }
    return(all_info)
}






#' Restart experimental simulations.
#'
#' Note that defaults for all this function's arguments except for
#' `sims_obj`, `new_starts`, `max_t`, and `save_every` are `NULL`,
#' which results in them being the same as for the original simulations.
#'
#' @param sims_obj A `cloneSims` object output from `sim_experiments`.
#' @param new_starts A dataframe or list of dataframes indicating the
#'     new starting abundances for all populations (wasps, mummies,
#'     all aphid lines) and stages.
#'     It should be the exact same format as what's in `sims_obj$all_info`.
#'     (But don't change `sims_obj$all_info` to make this be true!)
#' @param stage_ts_out Single logical for whether to output stage-structured
#'     information for all time points.
#'     If `TRUE`, the output object will contain this information in the
#'     `stage_ts` field, which will be a list of data frames.
#'     Defaults to `FALSE`.
#' @inheritParams sim_experiments
#'
#' @export
#'
restart_experiment <- function(sims_obj,
                               new_starts = NULL,
                               stage_ts_out = FALSE,
                               max_t = 250,
                               save_every = 1,
                               alate_field_disp_p = NULL,
                               K = NULL,
                               alate_b0 = NULL,
                               alate_b1 = NULL,
                               K_y_mult = NULL,
                               s_y = NULL,
                               a = NULL,
                               k = NULL,
                               h = NULL,
                               wasp_disp_p = NULL,
                               mum_smooth = NULL,
                               pred_rate = NULL,
                               plant_check_gaps = NULL,
                               max_plant_age = NULL,
                               clear_surv = NULL,
                               show_progress = NULL) {

    stopifnot(inherits(sims_obj, "cloneSims") |
                  inherits(sims_obj, "cloneSimsRestart"))
    if (!inherits(sims_obj$all_info_xptr, "externalptr")) {
        stop("\nSomething has happened to the \"all_info_xptr\" field ",
             "in this `cloneSims*` object. Please re-run simulation.")
    }

    # If `new_starts` is provided, check it carefully then adjust the
    # underlying C++ objects accordingly.
    if (!is.null(new_starts)) {
        stopifnot(inherits(new_starts, c("list", "data.frame")))
        test_against_starts <- make_all_info(sims_obj)
        n_reps <- length(test_against_starts)
        if (length(sims_obj$all_info) != n_reps) {
            stop("\nPlease do not directly edit `sims_obj$all_info`. ",
                 "You can try to fix this problem by running ",
                 "`sims_obj$all_info <- ",
                 "gameofclones:::make_all_info(sims_obj)`")
        }
        if (is.data.frame(new_starts)) {
            new_starts <- rep(list(new_starts), n_reps)
        }
        needed_cols <- c("field", "plant", "line", "type", "stage", "N")
        for (i in 1:n_reps) {
            not_present_cols <- needed_cols[!needed_cols %in%
                                                colnames(new_starts[[i]])]
            if (length(not_present_cols) > 0) {
                stop("\n`new_starts` contains at least one dataframe without ",
                     "the following needed column(s): ",
                     paste(not_present_cols, collapse = ", "))
            }
            x <- test_against_starts[[i]][, head(needed_cols, -1)]
            y <- new_starts[[i]][, head(needed_cols, -1)]
            z <- sims_obj$all_info[[i]][, head(needed_cols, -1)]
            if (!isTRUE(all.equal(x, y))) {
                stop("\n`new_starts` contains at least one dataframe that ",
                     "differs from its counterpart in ",
                     "`sims_obj$all_info`. ",
                     "Note that changing sims_obj$all_info does not fix ",
                     "this problem.")
            }
            if (!isTRUE(all.equal(x, z))) {
                stop("\nPlease do not directly edit `sims_obj$all_info`. ",
                     "You can try to fix this problem by running ",
                     "`sims_obj$all_info <- ",
                     "gameofclones:::make_all_info(sims_obj)`")
            }
        }
        N_vecs <- lapply(new_starts, function(x) x[["N"]])
    } else {
        N_vecs <- list()
    }


    # Fill in defaults / previously provided arguments:
    if (is.null(alate_field_disp_p)) alate_field_disp_p <-
            sims_obj$call[["alate_field_disp_p"]]
    if (is.null(K)) K <- sims_obj$call[["K"]]
    if (is.null(alate_b0)) alate_b0 <- sims_obj$call[["alate_b0"]]
    if (is.null(alate_b1)) alate_b1 <- sims_obj$call[["alate_b1"]]
    if (is.null(K_y_mult)) K_y_mult <- sims_obj$call[["K_y_mult"]]
    if (is.null(s_y)) s_y <- sims_obj$call[["s_y"]]
    if (is.null(a)) a <- sims_obj$call[["a"]]
    if (is.null(k)) k <- sims_obj$call[["k"]]
    if (is.null(h)) h <- sims_obj$call[["h"]]
    if (is.null(wasp_disp_p)) wasp_disp_p <- sims_obj$call[["wasp_disp_p"]]
    if (is.null(mum_smooth)) mum_smooth <- sims_obj$call[["mum_smooth"]]
    if (is.null(pred_rate)) pred_rate <- sims_obj$call[["pred_rate"]]
    if (is.null(plant_check_gaps)) plant_check_gaps <-
            sims_obj$call[["plant_check_gaps"]]
    if (is.null(max_plant_age)) max_plant_age <-
            sims_obj$call[["max_plant_age"]]
    if (is.null(clear_surv)) clear_surv <- sims_obj$call[["clear_surv"]]
    if (is.null(show_progress)) show_progress <-
            sims_obj$call[["show_progress"]]

    # Create / edit some items:
    n_fields <- sims_obj$call$n_fields
    n_lines <- length(sims_obj$call$clonal_lines)
    check_for_clear <- cumsum(rep(plant_check_gaps,
                                  ceiling(max_t / sum(plant_check_gaps))))
    check_for_clear <- check_for_clear[check_for_clear <= max_t]
    if (length(pred_rate) == 1) pred_rate <- rep(pred_rate, n_fields)
    if (length(alate_b0) == 1) alate_b0 <- rep(alate_b0, n_lines)
    if (length(alate_b1) == 1) alate_b1 <- rep(alate_b1, n_lines)
    if (length(K_y_mult) == 1) K_y_mult <- rep(K_y_mult, n_fields)
    if (length(s_y) == 1) s_y <- rep(s_y, n_fields)


    # This is like match.call but includes default arguments and
    # evaluates everything
    call_ <- mget(names(formals()), sys.frame(sys.nframe()))
    call_$n_fields <- n_fields
    call_$n_lines <- n_lines


    # Check validity of arguments:
    uint_check(max_t, "max_t", .min = 1)
    uint_check(save_every, "save_every", .min = 1)
    dbl_check(alate_field_disp_p, "alate_field_disp_p", .min = 0, .max = 1)
    dbl_check(K, "K", .min = 0)
    dbl_vec_check(alate_b0, "alate_b0")
    dbl_vec_check(alate_b1, "alate_b1")
    dbl_vec_check(K_y_mult, "K_y_mult", .min = 0)
    dbl_vec_check(s_y, "s_y", .min = 0, .max = 1)
    dbl_check(a, "a", .min = 0)
    dbl_check(k, "k", .min = 0)
    dbl_check(h, "h", .min = 0)
    dbl_check(wasp_disp_p, "wasp_disp_p", .min = 0, .max = 1)
    dbl_check(mum_smooth, "mum_smooth", .min = 0, .max = 1)
    dbl_vec_check(pred_rate, "pred_rate", .min = 0, .max = 1)
    uint_vec_check(check_for_clear, "check_for_clear")
    stopifnot(all(!duplicated(check_for_clear)))
    uint_check(max_plant_age, "max_plant_age")
    dbl_check(clear_surv, "clear_surv", .min = 0, .max = 1)
    stopifnot(inherits(show_progress, "logical") && length(show_progress) == 1)
    stopifnot(inherits(stage_ts_out, "logical") && length(stage_ts_out) == 1)


    # Make new C++ pointer to use for simulations:
    new_all_info_xptr <- restart_fill_other_pars(sims_obj$all_info_xptr,
                                                 K, alate_b0, alate_b1,
                                                 alate_field_disp_p, K_y_mult,
                                                 s_y, a, k, h, wasp_disp_p,
                                                 mum_smooth, pred_rate,
                                                 max_plant_age, clear_surv)
    # Now fill in abundances if they were input:
    if (length(N_vecs) > 0) {
        fields_from_vectors(new_all_info_xptr, N_vecs)
    }

    sims <- restart_experiments_cpp(new_all_info_xptr, max_t, save_every,
                                    check_for_clear, stage_ts_out, show_progress)

    sims[["aphids"]] <- sims[["aphids"]] %>%
        mutate(across(c("rep", "time", "field", "plant"), as.integer)) %>%
        mutate(across(c("rep", "field", "plant"), ~ .x + 1L)) %>%
        mutate(line = ifelse(type == "mummy", NA_character_, line)) %>%
        as_tibble()
    sims[["wasps"]] <- sims[["wasps"]] %>%
        mutate(across(c("rep", "time", "field"), as.integer)) %>%
        mutate(across(c("rep", "field"), ~ .x + 1L)) %>%
        as_tibble()

    sims[["all_info"]] <- make_all_info(sims)

    sims[["call"]] <- call_

    if (stage_ts_out) {
        sims[["stage_ts"]] <- lapply(sims[["stage_ts"]], as.data.frame)
    }

    class(sims) <- "cloneSimsRestart"

    return(sims)

}



# cloneSims* print methods ----

#'
#' @export
#' @noRd
#'
print.cloneSims <- function(x, ...) {
    cat("< cloneSims object >\n\n")
    cat("$aphids\n")
    print(head(x$aphids))
    cat("\n$wasps\n")
    print(head(x$wasps))
    cat("\n$all_info\n")
    print(lapply(x$all_info, head, n = 5))
    invisible(x)
}

#'
#' @export
#' @noRd
#'
print.cloneSimsRestart <- function(x, ...) {
    cat("< cloneSimsRestart object >\n\n")
    cat("$aphids\n")
    print(head(x$aphids))
    cat("\n$wasps\n")
    print(head(x$wasps))
    cat("\n$all_info\n")
    print(lapply(x$all_info, head, n = 5))
    invisible(x)
}
