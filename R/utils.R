
#' Calculate r_r / r_s (relative fitness for resistance vs susceptible aphids).
#'
#' This is almost directly from code in Ives et al. (2020).
#'
#' This function will be very slow if used inside an `apply` function
#' because it already uses `sapply` within.
#' Instead, `para` should be a vector.
#'
#' @param para Numeric vector of the proportion of aphids that are parasitized.
#' @param p_res The assumed proportion of aphids that are resistant to
#'     parasitism. This is the mean proportion from Ives et al. (2020).
#'
#' @return A numeric vector of r_r / r_s (relative fitness for resistance vs
#'     susceptible aphids)
#' @export
#'
#'
rel_res_fitness <- function(para, p_res = 0.48) {

    stopifnot(is.numeric(para) && is.null(dim(para)))
    stopifnot(is.numeric(p_res) && length(p_res) == 1)
    stopifnot(all(para <= 1 & para >= 0))
    stopifnot(p_res <= 1 && p_res >= 0)

    # Now we look up values of relative fitnesses based on the attack rate
    # associated with the parasitism reported (and the assumed proportion of
    # aphids that are currently resistant)
    para_lookup <- (1 - p_res) * rr_rs_lookup$ps + p_res * rr_rs_lookup$pr
    rr_rs <- sapply(para, function(p) {
        nearest <- which(abs(para_lookup - p) < 1.5e-5)

        # dealing with known discontinuity when p_res = 0.48 (default):
        if (length(nearest) == 0 &&
            p > para_lookup[rr_rs_lookup$a == 1.0581] &&
            p < para_lookup[rr_rs_lookup$a == 1.0582]) {
            nearest <- which(rr_rs_lookup$a == 1.0581 | rr_rs_lookup$a == 1.0582)
        }
        if (length(nearest) > 0) {
            return(mean(rr_rs_lookup$rr_rs[nearest], na.rm=TRUE))
        } else stop("cannot find fit for p")
    })

    return(rr_rs)

}


#' Remove tibbles from a cloneSims or cloneSimsRestart object
#'
#' @param sims_obj A cloneSims or cloneSimsRestart object
#'
#' @return
#' The same class of object as input, but with the `aphids` and `wasps` fields
#' converted to data frames instead of tibbles.
#'
#' @export
#'
rm_tibs <- function(sims_obj) {
    stopifnot(inherits(sims_obj, "cloneSims") |
                  inherits(sims_obj, "cloneSimsRestart"))
    for (n in c("aphids", "wasps")) {
        sims_obj[[n]] <- as.data.frame(sims_obj[[n]])
    }
    return(sims_obj)
}
