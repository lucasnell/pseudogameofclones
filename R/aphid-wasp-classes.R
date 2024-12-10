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
#' @param alate_b0 The proportion of offspring from apterous aphids is
#'     `inv_logit(alate_b0 + alate_b1 * N)` where `N` is the total number of
#'     aphids in that field.
#'     Defaults to `-5`.
#' @param alate_b1 The proportion of offspring from apterous aphids is
#'     `inv_logit(alate_b0 + alate_b1 * N)` where `N` is the total number of
#'     aphids in that field.
#'     Defaults to `0.0022`, which makes alate production only mildly
#'     density dependent.
#' @param distr_0 A 2-column matrix indicating the starting stage distribution.
#'     (The total abundance of aphids is adjusted in the simulation function.)
#'     If `NULL`, there will be no starting alates and apterous densities
#'     will follow the stable age distribution.
#'     If a matrix, it must have 5 rows (rows indicate instar)
#'     or a row per aphid age (rows indicate age in days).
#'     Matrix column indicates apterous vs alate.
#'     The matrix will be coerced to sum to 1.
#'     Defaults to `NULL`.
#' @param resistant Logical or vector of survivals of
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
#' @param sigma_x Standard deviation of environmental stochasticity for aphids.
#'     This argument has no effect if used in simulations with no environmental
#'     error (i.e., when `environ_error = FALSE`).
#'     Defaults to the internal object `environ$sigma_x`,
#'     which is from Meisner et al. (2014).
#' @param rho Environmental correlation among instars.
#'     This argument has no effect if used in simulations with no environmental
#'     error (i.e., when `environ_error = FALSE`).
#'     Defaults to `environ$rho`.
#'
#' @return A list with the necessary info to pass onto `all_fields`.
#'
#' @export
#'
clonal_line <- function(name,
                        alate_b0 = -5,
                        alate_b1 = 0.0022,
                        distr_0 = NULL,
                        resistant = FALSE,
                        surv_juv_apterous = NULL,
                        surv_adult_apterous = NULL,
                        repro_apterous = NULL,
                        surv_juv_alates = "low",
                        surv_adult_alates = "low",
                        repro_alates = "low",
                        surv_paras = "low",
                        temp = "low",
                        p_instar_smooth = 0.5,
                        sigma_x = environ$sigma_x,
                        rho = environ$rho) {


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
                    msg <- paste0("\ninput argument `", y,
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
        stop("\nLeslie matrices for apterous, alates, and parasitized",
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
    # Fill starting age distributions
    # --------------*

    if (is.null(distr_0)) {
        # If provided with just one number, assume stable age distribution:
        sad <- as.numeric(sad_leslie(leslie_array[,,1]))
        distr_0 <- cbind(sad, rep(0, length(sad)))
    } else if (is.numeric(distr_0) && is.matrix(distr_0) &&
               identical(dim(distr_0), c(5L, 2L))) {
        d0 <- distr_0
        distr_0 <- matrix(0, ns, 2)
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
                distr_0[irange,j] <- sad_ij * d0[i,j]
            }
        }
    } else if (! (is.numeric(distr_0) && is.matrix(distr_0) &&
                  identical(dim(distr_0), c(ns, 2L)))) {
        stop("\nThe `distr_0` arg to the ",
             "`clonal_line` function must be a single number or a 5x2 or ",
             ns, "x2 numeric matrix")
    }
    if (any(distr_0 < 0)) {
        stop("\nThe `distr_0` arg to the `clonal_line` function ",
             "cannot contain values < 0.")
    }
    if (sum(distr_0) != 1) distr_0 <- distr_0 / sum(distr_0)


    # --------------*
    # Fill other info
    # --------------*

    attack_surv <- c(0, 0)
    if (!(is.logical(resistant) && length(resistant) == 1) &&
        !(is.numeric(resistant) && all(resistant >= 0 & resistant <= 1))) {
        stop("`resistant` must be a single logical or a numeric vector with ",
             "all elements >= 0 and <= 1.")
    }
    if (is.logical(resistant) && resistant) attack_surv <- wasp_attack$attack_surv
    if (is.numeric(resistant)) attack_surv <- resistant

    # assume only adult alates can disperse across fields:
    field_disp_start <- sum(head(dev_times$instar_days[[temp]], -1))

    living_days <- dev_times$mum_days[[1]]


    # --------------*
    # Make final output
    # --------------*

    output <- structure(list(name = name,
                             sigma_x = sigma_x,
                             rho = rho,
                             attack_surv = attack_surv,
                             leslie_mat = leslie_array,
                             distr_0 = distr_0,
                             alate_b0 = alate_b0,
                             alate_b1 = alate_b1,
                             field_disp_start = field_disp_start,
                             living_days = living_days,
                             temp = temp),
                        class =  "AphidPop")

    return(output)

}


#'
#' @export
#' @noRd
#'
print.AphidPop <- function(x, ...) {

    cat("< Aphid clonal line >\n")
    cat("Name: ", x$name, "\n", sep = "")
    cat("Fields:\n")
    cat("  * name <string>\n")
    cat("  * distr_0 <matrix>\n")
    cat("  * attack_surv <vector>\n")
    cat("  * leslie <3D array>\n")
    cat("  * temp <string>\n")

    invisible(x)

}

#'
#' @export
#' @noRd
#'
c.AphidPop <- function(...) {
    z <- list(...)
    # combine survival probabilities so they're the same length:
    last_val <- list(attack_surv = 0)
    for (n in c("attack_surv")) {
        if (length(unique(sapply(z, function(x) length(x[[n]])))) > 1) {
            max_len <- max(sapply(z, function(x) length(x[[n]])))
            for (i in 1:length(z)) {
                if (length(z[[i]][[n]]) < max_len) {
                    if (tail(z[[i]][[n]], 1) != last_val[[n]]) {
                        err_msg <- paste("Error combining aphid lines: cannot",
                                         "combine values from field %s when",
                                         "one aphid line has fewer provided",
                                         "values than another and when the",
                                         "line having fewer values doesn't",
                                         "have a last value of %.0f") |>
                            sprintf(n, last_val[[n]])
                        stop(err_msg)
                    }
                    n_more <- max_len - length(z[[i]][[n]])
                    z[[i]][[n]] <- c(z[[i]][[n]], rep(last_val[[n]], n_more))
                }
            }
        }
    }

    densities_0 <- lapply(z, function(x) x$distr_0)
    if (!length(unique(sapply(densities_0, nrow))) == 1) {
        stop("\nAll aphid lines must have the same sized density matrices\n")
    }

    leslie_cubes <- lapply(z, function(x) x$leslie)
    if (! all(apply(sapply(leslie_cubes, dim), 1,
                    function(zz) length(unique(zz)))) == 1) {
        stop("\nAll aphid lines must have the same sized Leslie matrices\n")
    }

    temp <- z[[1]][["temp"]]
    if (length(z) > 1) {
        for (i in 2:length(z)) {
            if (z[[i]][["temp"]] != temp) {
                stop("\n`temp` fields in aphid lines must match")
            }
        }
    }
    maphids <- structure(z,
                         names = sapply(z, function(x) x$name),
                         class = "multiAphid",
                         ptr = make_aphids_ptr(z))
    return(maphids)
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





# Create wasp pop info ----
#' Organize wasp population info.
#'
#'
#' @param s_y Daily survival rate of adult wasps.
#'     Defaults to `populations$s_y`, which is from Meisner et al. (2014).
#' @param sex_ratio Sex ratio of adult wasps. Defaults to `0.5`.
#' @param a Parasitoid attack rate. Defaults to the internal object
#'     `wasp_attack$a`, which is from Meisner et al. (2014).
#' @param k Aggregation parameter of the negative binomial distribution.
#'     Defaults to the internal object `wasp_attack$k`,
#'     which is from Meisner et al. (2014).
#' @param h Parasitoid handling time. Defaults to the internal object
#'     `wasp_attack$h`, which is from Meisner et al. (2014).
#' @param rel_attack Relative parasitoid attack rate among instars.
#'     Defaults to `wasp_attack$rel_attack`, which is from
#'     Meisner et al. (2014).
#' @param sigma_y Standard deviation of environmental stochasticity for wasps.
#'     This argument has no effect if `environ_error = FALSE`.
#'     Defaults to the internal object `environ$sigma_y`,
#'     which is from Meisner et al. (2014).
#' @param mum_smooth Proportion of mummies that will NOT take exactly 3 days
#'     to develop. As this value approaches 2/3, it will provide greater
#'     smoothing of wasp numbers through time.
#'     Defaults to `0.4`.
#'
#' @export
#'
#' @return List containing an external pointer of a C++ object that can be
#'     passed onto `all_fields`.
#'
wasp_pop <- function(s_y = populations$s_y,
                     sex_ratio = populations$sex_ratio,
                     a = wasp_attack$a,
                     k = wasp_attack$k,
                     h = wasp_attack$h,
                     rel_attack = wasp_attack$rel_attack,
                     sigma_y = environ$sigma_y,
                     mum_smooth = 0.4) {

    dbl_check(s_y, "s_y", .max = 1, .min = 0)
    dbl_check(sex_ratio, "sex_ratio", .min = 0, .max = 1)
    dbl_check(a, "a", .min = 0)
    dbl_check(k, "k", .min = 0)
    dbl_check(h, "h", .min = 0)
    dbl_vec_check(rel_attack, "rel_attack", .min = 0)
    if (sum(rel_attack) != 1) stop("\n`rel_attack` must sum to 1")
    dbl_check(sigma_y, "sigma_y", .min = 0)
    dbl_check(mum_smooth, "mum_smooth", .max = 1, .min = 0)


    ptr <- make_wasps_ptr(rel_attack = rel_attack,
                          a = a,
                          k = k,
                          h = h,
                          sex_ratio = sex_ratio,
                          s_y = s_y,
                          sigma_y = sigma_y,
                          mummy_smooth = mum_smooth,
                          mummy_dev_time = dev_times$mum_days[2])

    out <- structure(list("ptr" = ptr, rel_attack = rel_attack),
                     class = "WaspPop")

    return(out)

}




#' Create AllFields object containing aphid, wasp, and field parameters
#'
#' @param clonal_lines A `multiAphid` object containing the aphid-line-specific
#'     info for all the lines in the simulations.
#'     Each object in the `multiAphid` is an `aphid` object that results from
#'     the `clonal_line` function.
#'     Combine them using `c(aphid_obj1, aphid_obj2)`.
#' @param wasps A `WaspPop` object containing the parasitoid wasp info
#'     for the simulations.
#' @param n_fields The number of fields to simulate.
#'     Both wasps and aphids operate separately across fields but can be
#'     connected via dispersal (see arguments `wasp_disp_m0`, `wasp_disp_m1`,
#'     and `alate_field_disp_p`).
#'     Defaults to `2`.
#' @param environ_error Logical for whether to have environmental stochasticity.
#'     This argument applies to both aphids and wasps.
#'     Defaults to `FALSE`.
#' @param aphid_demog_error Logical for whether to have demographic
#'     stochasticity for aphids. Defaults to `FALSE`.
#' @param wasp_demog_error Logical for whether to have demographic
#'     stochasticity for wasps. Defaults to `FALSE`.
#' @param K Aphid density dependence. Either a single number for all fields,
#'     or a vector of length `n_fields` to define separate values by field.
#'     Defaults to `12.5e3` because this caused simulations to
#'     approximately match experiments.
#' @param K_y_mult The number multiplied by `K` to get density dependence for
#'     parasitized aphids.
#'     Either a single number for all fields, or a vector of length `n_fields`
#'     to define separate values by field.
#'     Defaults to `1 / 1.57`, which is from Meisner et al. (2014).
#' @param aphid_density_0 The total abundance of aphids in each field.
#'     Either a single number for all lines in all fields, or a
#'     `n_fields`-row matrix with a column for each aphid line to define
#'     separate values by line and field.
#'     Defaults to `32`.
#' @param wasp_density_0 The abundance of adult wasps in each field.
#'     Either a single number for all fields, or a vector of length `n_fields`
#'     to define separate values by field.
#'     Defaults to `c(0, 3)`.
#' @param mummy_density_0
#'     Either a single number for all fields, or a vector of length `n_fields`
#'     to define separate values by field.
#'     Defaults to `0`.
#' @param wasp_delay
#' @param wasp_disp_m0 Proportion of adult wasps from each field that
#'     are added to the dispersal pool when there are no aphids present.
#'     After adding wasps to the pool, they are then evenly distributed
#'     to all fields unless values for `wasp_field_attract` are provided.
#'     Defaults to `0.3`.
#' @param wasp_disp_m1 Effect of aphid density on wasp emigration from a patch.
#'     Emigration is `wasp_disp_m0 * exp(-wasp_disp_m1 * log(z))`, where `z` is
#'     the total number of living aphids in the patch.
#'     Note that if `wasp_disp_m0 = 0`, then this parameter doesn't change
#'     anything.
#'     Defaults to `0.34906`.
#' @param wasp_field_attract Relative attractiveness of fields to wasps.
#'     This affects the proportion of wasps that immigrate from the dispersal
#'     pool to each field.
#'     It doesn't change the number of wasps that leave fields.
#'     This can be a single numeric or a numeric vector of length `n_fields`.
#'     If a single numeric is provided, all fields are equally attractive
#'     to wasps.
#'     If a vector is provided, then the vector is divided by its sum
#'     (to make it sum to 1), then those values are used as the proportion of
#'     wasps from the dispersal pool that immigrate to each field.
#'     Note that if `wasp_disp_m0 = 0`, then this parameter doesn't change
#'     anything.
#'     Defaults to `1`.
#' @param pred_rate Daily predation rate on aphids and mummies.
#'     Defaults to `0.1`.
#' @param alate_field_disp_p Proportion of alates from each field that
#'     are added to the dispersal pool.
#'     After adding alates to the pool, they are then evenly distributed
#'     to all fields.
#'     Defaults to `0.1`.
#' @param constant_wasps Logical for whether to keep adult wasps at the same
#'     density throughout simulations.
#'     This can be a single logical or a `n_fields`-length vector.
#'     Defaults to `FALSE`.
#' @param extinct_N Threshold below which a line is considered extinct.
#'     Defaults to `1`.
#'
#'
#' @return List containing an external pointer of a C++ object that can be
#'     passed onto `sim_fields`.
#'
#' @export
#'
#' @examples
#' wp <- wasp_pop()
#' sa <- clonal_line("susceptible")
#' ra <- clonal_line("resistant", resistant = TRUE, repro_apterous = "low")
#' af <- all_fields(c(sa, ra), wp)
#'
#'
all_fields <- function(clonal_lines,
                       wasps,
                       n_fields = 2,
                       environ_error = FALSE,
                       aphid_demog_error = FALSE,
                       wasp_demog_error = FALSE,
                       K = 1 / 4.67e-4,
                       K_y_mult = 1 / 1.57,
                       aphid_density_0 = 32,
                       wasp_density_0 = c(3, 0),
                       mummy_density_0 = 0,
                       wasp_delay = 8,
                       wasp_disp_m0 = 0.3,
                       wasp_disp_m1 = 0.34906,
                       wasp_field_attract = 1,
                       pred_rate = 0,
                       alate_field_disp_p = 0.1,
                       constant_wasps = FALSE,
                       extinct_N = 1) {

    # -----------------*
    # Dimension checks and adjustments
    n_lines <- length(clonal_lines)
    if (length(aphid_density_0) == 1) {
        aphid_density_0 <- matrix(aphid_density_0, n_fields, n_lines)
    }
    if (! all(dim(aphid_density_0) == c(n_fields, n_lines))) {
        stop(paste0("\nIf provided as a matrix, `aphid_density_0` must have ",
                    "a row per field and a column per clonal line. ",
                    "Here, it should be ", n_fields, " x ", n_lines, ". ",
                    "Yours is ", nrow(aphid_density_0), " x ",
                    ncol(aphid_density_0), "."))
    }
    if (length(wasp_density_0) == 1)
        wasp_density_0 <- rep(wasp_density_0, n_fields)
    if (length(wasp_delay) == 1) wasp_delay <- rep(wasp_delay, n_fields)
    if (length(mummy_density_0) == 1)
        mummy_density_0 <- rep(mummy_density_0, n_fields)
    if (length(K) == 1) K <- rep(K, n_fields)
    if (length(K_y_mult) == 1) K_y_mult <- rep(K_y_mult, n_fields)
    if (length(pred_rate) == 1) pred_rate <- rep(pred_rate, n_fields)
    if (length(constant_wasps) == 1)
        constant_wasps <- rep(constant_wasps, n_fields)
    if (length(wasp_field_attract) == 1)
        wasp_field_attract <- rep(wasp_field_attract, n_fields)

    # Other dimensions checks are done inside `make_field_ptr`


    # -----------------*
    # Basic type checks
    if (!inherits(clonal_lines, "multiAphid")) {
        if (inherits(clonal_lines, "AphidPop")) {
            clonal_lines <- c(clonal_lines)
        } else stop("\n`clonal_lines` must be a multiAphid object.\n")
    }
    if (!inherits(attr(clonal_lines, "ptr"), "externalptr"))
        stop("\n`attr(clonal_lines, \"ptr\")` must be an externalptr object.\n")
    if (!inherits(wasps, "WaspPop"))
        stop("\n`wasps` must be a WaspPop object.\n")
    if (!inherits(wasps$ptr, "externalptr"))
        stop("\n`wasps$ptr` must be an externalptr object.\n")
    uint_check(n_fields, "n_fields", .min = 1)
    lgl_check(environ_error, "environ_error")
    lgl_check(aphid_demog_error, "aphid_demog_error")
    lgl_check(wasp_demog_error, "wasp_demog_error")
    dbl_vec_check(K, "K")
    dbl_vec_check(K_y_mult, "K_y_mult")
    dbl_mat_check(aphid_density_0, "aphid_density_0")
    dbl_vec_check(wasp_density_0, "wasp_density_0")
    dbl_vec_check(mummy_density_0, "mummy_density_0")
    uint_vec_check(wasp_delay, "wasp_delay")
    dbl_check(wasp_disp_m0, "wasp_disp_m0", .min = 0, .max = 1)
    dbl_check(wasp_disp_m1, "wasp_disp_m1")
    dbl_vec_check(wasp_field_attract, "wasp_field_attract", .min = 0)
    stopifnot(sum(wasp_field_attract) > 0)
    dbl_vec_check(pred_rate, "pred_rate")
    dbl_check(alate_field_disp_p, "alate_field_disp_p", .min = 0, .max = 1)
    lgl_vec_check(constant_wasps, "constant_wasps")
    dbl_check(extinct_N, "extinct_N")

    # Adjusting wasp relative attack vector size to match number of aphid
    # stages in Leslie matrices
    rel_attack <- wasps$rel_attack
    if (length(rel_attack) == 5) {
        rel_attack <- rel_attack / sum(rel_attack)
        cl <- clonal_lines[[1]]
        dt <- dev_times$instar_days[[cl$temp]]
        n_adult_days <- nrow(cl$leslie_mat) - sum(head(dt, -1))
        stopifnot(n_adult_days >= 0)
        if (tail(dt, 1) != n_adult_days) dt[length(dt)] <- n_adult_days
        # Commented version isn't used bc it wasn't done this way when fitting
        # the model.
        # rel_attack__ <- mapply(function(.x, .y) rep(.x, .y) / .y,
        #                        rel_attack,  dt)
        rel_attack__ <- mapply(rep, rel_attack,  dt)
        rel_attack <- do.call(c, rel_attack__)
    } else stopifnot(length(rel_attack) == nrow(leslie_cubes[[1]]))


    ptr <- make_field_ptr(aphid_demog_error = aphid_demog_error,
                          aphid_density_0 = aphid_density_0,
                          wasp_demog_error = wasp_demog_error,
                          wasp_density_0 = wasp_density_0,
                          wasp_delay = wasp_delay,
                          mummy_density_0 = mummy_density_0,
                          environ_error = environ_error,
                          aphids_ptr = attr(clonal_lines, "ptr"),
                          wasp_ptr = wasps$ptr,
                          n_fields = n_fields,
                          K = K,
                          K_y = K * K_y_mult,
                          pred_rate = pred_rate,
                          extinct_N = extinct_N,
                          constant_wasps = constant_wasps,
                          alate_field_disp_p = alate_field_disp_p,
                          wasp_disp_m0 = wasp_disp_m0,
                          wasp_disp_m1 = wasp_disp_m1,
                          wasp_field_attract = wasp_field_attract,
                          new_rel_attack = rel_attack)

    # Note: n_fields and aphid_names used in `make_perturb_list` inside
    # `sim_fields`
    out <- structure(list(n_fields = n_fields,
                          aphid_names = sapply(clonal_lines, function(x) x$name),
                          ptr = ptr),
                     class = "AllFields")

    return(out)
}
