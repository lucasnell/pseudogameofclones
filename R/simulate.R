

# Check `perturb` dataframe and return a list that's properly formatted for use
# in cpp code.
make_perturb_list <- function(perturb, fields) {

    n_fields <- fields$n_fields
    aphid_names <- fields$aphid_names

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
        n_lines <- length(aphid_names)
        if (is.character(perturb$who)) {
            stopifnot(all(perturb$who %in% c(aphid_names, "mummies", "wasps")))
            perturb_who <- integer(length(perturb$who))
            perturb_who[perturb$who %in% aphid_names] <- -1 +
                match(perturb$who[perturb$who %in% aphid_names], aphid_names)
            perturb_who[perturb$who == "mummies"] <- n_lines
            perturb_who[perturb$who == "wasps"] <- n_lines + 1
        } else {
            stopifnot(all(perturb$who <= (n_lines+2) & perturb$who > 0))
            perturb_who <- perturb$who - 1
        }
    }

    perturb_list <- list(when = perturb_when,
                         where = perturb_where,
                         who = perturb_who,
                         how = perturb_how)
    uint_vec_check(perturb_list$when, "perturb_list$when")
    uint_vec_check(perturb_list$where, "perturb_list$where")
    uint_vec_check(perturb_list$who, "perturb_list$who")
    dbl_vec_check(perturb_list$how, "perturb_list$how", .min = 0)

    return(perturb_list)
}







# full fun docs ----
#' Simulate fields of aphids and natural enemies.
#'
#'
#' @param fields An `AllFields` object created using the
#'     `all_fields` function.
#' @param n_reps Number of reps to simulate. Defaults to `1`.
#' @param max_t How many days to simulate. Defaults to `250`.
#' @param sep_adults Single logical for whether to separate adults vs
#'     juvenile aphids in the output.
#'     Useful for counting aphids with wings (adult alates).
#'     Defaults to `FALSE`.
#' @param save_every Abundances will be stored every `save_every` time points.
#'     Defaults to `1`.
#' @param n_threads Number of threads to use if `n_reps > 1`
#'     (ignored otherwise).
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
#' @return A list of dataframes of aphids, wasps.
#'
#' @examples
#' wp <- wasp_pop()
#' sa <- clonal_line("susceptible")
#' ra <- clonal_line("resistant", resistant = TRUE, repro_apterous = "low")
#' af <- all_fields(c(sa, ra), wp)
#' sf <- sim_fields(af)
#'
#'
sim_fields <- function(fields,
                       n_reps = 1,
                       max_t = 250,
                       sep_adults = FALSE,
                       save_every = 1,
                       n_threads = 1,
                       show_progress = FALSE,
                       perturb = NULL) {

    # This is like match.call but includes default arguments and
    # evaluates everything
    call_ <- mget(names(formals()), sys.frame(sys.nframe()))

    if (!inherits(fields, "AllFields"))
        stop("\n`fields` must be an AllFields object.\n")
    if (!inherits(fields$ptr, "externalptr"))
        stop("\n`fields$ptr` must be an externalptr object.\n")

    uint_check(n_reps, "n_reps", .min = 1)
    uint_check(max_t, "max_t", .min = 1)
    lgl_check(sep_adults, "sep_adults")
    uint_check(save_every, "save_every", .min = 1)
    uint_check(n_threads, "n_threads", .min = 1)
    lgl_check(show_progress, "show_progress")

    perturb_list <- make_perturb_list(perturb, fields)
    perturb_ptr <- make_perturb_ptr(perturb_list$when, perturb_list$where,
                                    perturb_list$who, perturb_list$how)


    sims <- sim_pseudogameofclones_cpp(fields$ptr, perturb_ptr,
                                       n_reps, max_t, save_every, sep_adults,
                                       n_threads, show_progress)

    sims[["aphids"]] <- mutate(sims[["aphids"]],
                               across(c("rep", "time", "field"), as.integer))
    sims[["aphids"]] <- mutate(sims[["aphids"]],
                               across(c("rep", "field"), ~ .x + 1L))
    sims[["aphids"]] <- mutate(sims[["aphids"]],
                               line = ifelse(type == "mummy", NA, line))
    sims[["aphids"]] <- as_tibble(sims[["aphids"]])

    sims[["wasps"]] <- mutate(sims[["wasps"]],
                              across(c("rep", "time", "field"), as.integer))
    sims[["wasps"]] <- mutate(sims[["wasps"]],
                              across(c("rep", "field"), ~ .x + 1L))
    sims[["wasps"]] <- as_tibble(sims[["wasps"]])

    sims[["all_info"]] <- make_all_info(sims)

    sims[["call"]] <- call_

    class(sims) <- "cloneSims"

    return(sims)
}








# Internal function to create the list of dataframes containing all info.
make_all_info <- function(sims_obj) {
    if (!inherits(sims_obj$all_info_xptr, "externalptr")) {
        stop("\nSomething has happened to the \"all_info_xptr\" field ",
             "in this `cloneSims` object. Please re-run simulation.")
    }
    all_info <- fields_to_data_frames(sims_obj$all_info_xptr)
    for (i in 1:length(all_info)) {
        is_wasp <- all_info[[i]][["type"]] == "wasp"
        is_mummy <- all_info[[i]][["type"]] == "mummy"
        all_info[[i]][["line"]][is_wasp | is_mummy] <- NA
    }
    return(all_info)
}






#' Restart experimental simulations.
#'
#' Note for most of this function's arguments, their default value is `NULL`,
#' which results in them being the same as for the original simulations.
#' The exceptions to this are the arguments `sims_obj`, `new_starts`,
#' `stage_ts_out`, `max_t`, `save_every`, `perturb`, and `no_warns`.
#' Arguments `sims_obj`, `new_starts`, and `stage_ts_out` are new to this
#' function (see their documentation below for details).
#' Arguments `max_t` and `save_every` are set to `250` and `1`, respectively,
#' which results in a shorter time series with greater resolution, the typical
#' use case for this function.
#' Argument `perturb` has a default value of `NULL` which results in no
#' perturbations.
#' You can input a new data frame if you want to perturb these restarted
#' experiments.
#' If the initial simulations had perturbations but the new one won't,
#' this function will provide a warning indicating this.
#' You can suppress this warning by setting `no_warns = TRUE`.
#'
#'
#' @param sims_obj A `cloneSims` object output from `sim_fields`.
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
#' @inheritParams sim_fields
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
                               wasp_disp_m0 = NULL,
                               wasp_disp_m1 = NULL,
                               wasp_field_attract = NULL,
                               mum_smooth = NULL,
                               pred_rate = NULL,
                               sep_adults = NULL,
                               show_progress = NULL,
                               perturb = NULL,
                               no_warns = FALSE,
                               n_threads = 1L) {

    stopifnot(inherits(sims_obj, "cloneSims") |
                  inherits(sims_obj, "cloneSimsRestart"))
    if (!inherits(sims_obj$all_info_xptr, "externalptr")) {
        stop("\nSomething has happened to the \"all_info_xptr\" field ",
             "in this `cloneSims*` object. Please re-run simulation.")
    }

    # Extract aphid line info (the location differs depending on whether it's
    # a restart)
    if (inherits(sims_obj, "cloneSimsRestart")) {
        clonal_lines <- sims_obj$clonal_lines
    } else {
        clonal_lines <- sims_obj$call$clonal_lines
    }
    stopifnot(!is.null(clonal_lines))
    if (!inherits(clonal_lines, "multiAphid")) {
        if (inherits(clonal_lines, "aphid")) {
            clonal_lines <- c(clonal_lines)
        } else stop("\n`clonal_lines` must be a multiAphid/aphid object.")
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
                 "pseudogameofclones:::make_all_info(sims_obj)`")
        }
        if (is.data.frame(new_starts)) {
            new_starts <- rep(list(new_starts), n_reps)
        }
        needed_cols <- c("field", "line", "type", "stage", "N")
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
                     "pseudogameofclones:::make_all_info(sims_obj)`")
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
    if (is.null(wasp_disp_m0)) wasp_disp_m0 <- sims_obj$call[["wasp_disp_m0"]]
    if (is.null(wasp_disp_m1)) wasp_disp_m1 <- sims_obj$call[["wasp_disp_m1"]]
    if (is.null(wasp_field_attract)) {
        wasp_field_attract <- sims_obj$call[["wasp_field_attract"]]
    }
    if (is.null(mum_smooth)) mum_smooth <- sims_obj$call[["mum_smooth"]]
    if (is.null(pred_rate)) pred_rate <- sims_obj$call[["pred_rate"]]
    if (is.null(sep_adults)) sep_adults <- sims_obj$call[["sep_adults"]]
    if (is.null(show_progress)) show_progress <-
            sims_obj$call[["show_progress"]]

    # Error if any of these are still NULL, since this means they weren't in
    # `sims_obj$call`
    stopifnot(!is.null(alate_field_disp_p))
    stopifnot(!is.null(K))
    stopifnot(!is.null(alate_b0))
    stopifnot(!is.null(alate_b1))
    stopifnot(!is.null(K_y_mult))
    stopifnot(!is.null(s_y))
    stopifnot(!is.null(a))
    stopifnot(!is.null(k))
    stopifnot(!is.null(h))
    stopifnot(!is.null(wasp_disp_m0))
    stopifnot(!is.null(wasp_disp_m1))
    stopifnot(!is.null(wasp_field_attract))
    stopifnot(!is.null(mum_smooth))
    stopifnot(!is.null(pred_rate))
    stopifnot(!is.null(sep_adults))
    stopifnot(!is.null(show_progress))

    # Create / edit some items:
    n_fields <- sims_obj$call$n_fields
    n_lines <- length(clonal_lines)
    uint_check(n_fields, "n_fields", .min = 1)
    uint_check(n_lines, "n_lines", .min = 1)
    if (length(pred_rate) == 1) pred_rate <- rep(pred_rate, n_fields)
    if (length(alate_b0) == 1) alate_b0 <- rep(alate_b0, n_lines)
    if (length(alate_b1) == 1) alate_b1 <- rep(alate_b1, n_lines)
    if (length(K) == 1) K <- rep(K, n_fields)
    if (length(K_y_mult) == 1) K_y_mult <- rep(K_y_mult, n_fields)
    if (length(s_y) == 1) s_y <- rep(s_y, n_fields)
    stopifnot(is.numeric(wasp_field_attract))
    if (length(wasp_field_attract) == 1) {
        wasp_field_attract <- rep(1, n_fields)
    }

    stopifnot(inherits(no_warns, "logical") && length(no_warns) == 1)
    if (!no_warns && !is.null(sims_obj$call$perturb) && is.null(perturb)) {
        warning(paste("\nThe initial simulations had perturbations, but the",
                      "new one won't. Provide a dataframe to the `perturb`",
                      "argument to `restart_experiments` to add",
                      "perturbations."))
    }

    perturb_list <- make_perturb_list(perturb, fields)
    perturb_ptr <- make_perturb_ptr(perturb_list$when, perturb_list$where,
                                    perturb_list$who, perturb_list$how)


    # This is like match.call but includes default arguments and
    # evaluates everything
    call_ <- mget(names(formals()), sys.frame(sys.nframe()))
    call_$n_fields <- n_fields
    call_$n_lines <- n_lines


    # Check validity of arguments:
    uint_check(max_t, "max_t", .min = 1)
    uint_check(save_every, "save_every", .min = 1)
    uint_check(n_threads, "n_threads", .min = 1)
    dbl_check(alate_field_disp_p, "alate_field_disp_p", .min = 0, .max = 1)
    dbl_vec_check(K, "K", .min = 0)
    dbl_vec_check(alate_b0, "alate_b0")
    dbl_vec_check(alate_b1, "alate_b1")
    dbl_vec_check(K_y_mult, "K_y_mult", .min = 0)
    dbl_vec_check(s_y, "s_y", .min = 0, .max = 1)
    dbl_check(a, "a", .min = 0)
    dbl_check(k, "k", .min = 0)
    dbl_check(h, "h", .min = 0)
    dbl_check(wasp_disp_m0, "wasp_disp_m0", .min = 0, .max = 1)
    dbl_check(wasp_disp_m1, "wasp_disp_m1")
    dbl_vec_check(wasp_field_attract, "wasp_field_attract", .min = 0)
    stopifnot(sum(wasp_field_attract) > 0)
    dbl_check(mum_smooth, "mum_smooth", .min = 0, .max = 1)
    dbl_vec_check(pred_rate, "pred_rate", .min = 0, .max = 1)
    stopifnot(inherits(sep_adults, "logical") && length(sep_adults) == 1)
    stopifnot(inherits(show_progress, "logical") && length(show_progress) == 1)
    stopifnot(inherits(stage_ts_out, "logical") && length(stage_ts_out) == 1)


    # Make new C++ pointer to use for simulations:
    new_all_info_xptr <- restart_fill_other_pars(sims_obj$all_info_xptr,
                                                 K, alate_b0, alate_b1,
                                                 alate_field_disp_p,
                                                 K * K_y_mult,
                                                 s_y, a, k, h,
                                                 wasp_disp_m0, wasp_disp_m1,
                                                 wasp_field_attract,
                                                 mum_smooth, pred_rate)
    # Now fill in abundances if they were input:
    if (length(N_vecs) > 0) {
        fields_from_vectors(new_all_info_xptr, N_vecs)
    }

    sims <- restart_experiments_cpp(new_all_info_xptr, perturb_ptr,
                                    max_t, save_every,
                                    stage_ts_out, sep_adults, show_progress,
                                    n_threads)

    sims[["aphids"]] <- mutate(sims[["aphids"]],
                               across(c("rep", "time", "field"), as.integer))
    sims[["aphids"]] <- mutate(sims[["aphids"]],
                               across(c("rep", "field"), ~ .x + 1L))
    sims[["aphids"]] <- mutate(sims[["aphids"]],
                               line = ifelse(type == "mummy", NA, line))
    sims[["aphids"]] <- as_tibble(sims[["aphids"]])

    sims[["wasps"]] <- mutate(sims[["wasps"]],
                              across(c("rep", "time", "field"), as.integer))
    sims[["wasps"]] <- mutate(sims[["wasps"]],
                              across(c("rep", "field"), ~ .x + 1L))
    sims[["wasps"]] <- as_tibble(sims[["wasps"]])

    sims[["all_info"]] <- make_all_info(sims)

    sims[["call"]] <- call_

    # This is important if you want to restart output from this function.
    sims[["clonal_lines"]] <- clonal_lines

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
#' @importFrom dplyr as_tibble
#'
print.cloneSims <- function(x, ...) {
    cat("< cloneSims object >\n\n")
    cat("$aphids\n")
    print(as_tibble(x$aphids), n = 5)
    cat("\n$wasps\n")
    print(as_tibble(x$wasps), n = 5)
    cat("\n$all_info\n")
    for (i in 1:length(x$all_info)) {
        cat(sprintf("[[%i]]\n", i))
        print(as_tibble(x$all_info[[i]]), n = 5)
    }
    invisible(x)
}

#'
#' @export
#' @noRd
#'
#' @importFrom dplyr as_tibble
#'
print.cloneSimsRestart <- function(x, ...) {
    cat("< cloneSimsRestart object >\n\n")
    cat("$aphids\n")
    print(as_tibble(x$aphids), n = 5)
    cat("\n$wasps\n")
    print(as_tibble(x$wasps), n = 5)
    cat("\n$all_info\n")
    for (i in 1:length(x$all_info)) {
        cat(sprintf("[[%i]]\n", i))
        print(as_tibble(x$all_info[[i]]), n = 5)
    }
    invisible(x)
}




