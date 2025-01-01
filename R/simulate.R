

# Check `perturb` dataframe, create vectors that are properly formatted,
# then create the XPtr to an object for use in C++
make_perturb_ptr_R <- function(perturb, fields) {

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

    uint_vec_check(perturb_when, "perturb_when")
    uint_vec_check(perturb_where, "perturb_where")
    uint_vec_check(perturb_who, "perturb_who")
    dbl_vec_check(perturb_how, "perturb_how", .min = 0)

    perturb_ptr <- make_perturb_ptr(perturb_when, perturb_where,
                                    perturb_who, perturb_how)

    return(perturb_ptr)
}



#
# Convert output from `sim_fields_cpp` to a cloneSims object.
#
make_cloneSims_obj <- function(sims, call_, stage_ts_out, restarted = FALSE) {

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

    sims[["restarted"]] <- restarted

    if (stage_ts_out) {
        sims[["stage_ts"]] <- lapply(sims[["stage_ts"]], as.data.frame)
    }

    class(sims) <- "cloneSims"

    return(sims)
}



# full fun docs ----
#' Simulate fields of aphids and natural enemies.
#'
#'
#' @param fields An `AllFields` object created using the
#'     `all_fields` function.
#' @param n_reps Number of reps to simulate. Defaults to `1`.
#' @param max_t How many days to simulate. Defaults to `250`.
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
#' @param save_every Abundances will be stored every `save_every` time points.
#'     Defaults to `1`.
#' @param sep_adults Single logical for whether to separate adults vs
#'     juvenile aphids in the output.
#'     Useful for counting aphids with wings (adult alates).
#'     Defaults to `FALSE`.
#' @param stage_ts_out Single logical for whether to output stage-structured
#'     information for all time points.
#'     If `TRUE`, the output object will contain this information in the
#'     `stage_ts` field, which will be a list of data frames.
#'     Defaults to `FALSE`.
#' @param show_progress Boolean for whether to show progress bar.
#'     Defaults to `FALSE`.
#' @param n_threads Number of threads to use if `n_reps > 1`
#'     (ignored otherwise).
#'     Defaults to `1`.
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
                       perturb = NULL,
                       save_every = 1,
                       sep_adults = FALSE,
                       stage_ts_out = FALSE,
                       show_progress = FALSE,
                       n_threads = 1L) {

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
    lgl_check(stage_ts_out, "stage_ts_out")

    perturb_ptr <- make_perturb_ptr_R(perturb, fields)

    # Make XPtr<std::vector<AllFields>> that will be changed inside
    # `sim_fields_cpp`:
    all_fields_vec_ptr <- make_all_fields_vec_ptr(fields$ptr, n_reps)

    sims <- sim_fields_cpp(all_fields_vec_ptr, perturb_ptr,
                           max_t, save_every, sep_adults,
                           n_threads, show_progress, stage_ts_out)

    sims <- make_cloneSims_obj(sims, call_, stage_ts_out)

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






#
# For `restart.cloneSims, if `new_starts` is provided, check it carefully then
# return a list of N vectors (1 per rep) in the exact order necessary.
#
check_new_starts_adjust_N <- function(sims_obj, new_starts, new_all_fields_vec_ptr) {

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

    fields_from_vectors(new_all_fields_vec_ptr, N_vecs)

    invisible(NULL)

}




#' Generic restart method
#'
#' @noRd
#' @export
#'
restart <- function(sims_obj, ...) {
    UseMethod("restart")
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
#' @param sims_obj A `cloneSims` object output from `sim_fields` or
#'     from`restart.cloneSims`.
#' @param new_starts A dataframe or list of dataframes indicating the
#'     new starting abundances for all populations (wasps, mummies,
#'     all aphid lines) and stages.
#'     It should be the exact same format as what's in `sims_obj$all_info`.
#'     (But don't change `sims_obj$all_info` to make this be true!)
#' @param new_field_pars A named list of field-related arguments to change in
#'     the restarted simulations.
#'     The arguments available to change are `aphid_demog_error`,
#'     `wasp_demog_error`, `K`, `K_y`, `pred_rate`, `constant_wasps`,
#'     `alate_field_disp_p`, `wasp_disp_m0`, `wasp_disp_m1`, and
#'     `wasp_field_attract`.
#'     These arguments should follow the same specifications as for the
#'     function [all_fields()].
#'     Defaults to `NULL`, resulting in no changes.
#' @param no_warns Single logical for whether not to output a warning about
#'     having a perturbations in the initial simulation without them in this
#'     restarted simulation. Defaults to `FALSE`.
#' @inheritParams sim_fields
#'
#' @export
#'
#' @examples
#' wp <- wasp_pop()
#' sa <- clonal_line("susceptible")
#' ra <- clonal_line("resistant", resistant = TRUE, repro_apterous = "low")
#' af <- all_fields(c(sa, ra), wp)
#' sf <- sim_fields(af)
#' rf <- restart(sf)
#'
restart.cloneSims <- function(sims_obj,
                              new_starts = NULL,
                              new_field_pars = NULL,
                              max_t = 250,
                              perturb = NULL,
                              no_warns = FALSE,
                              save_every = 1,
                              sep_adults = FALSE,
                              stage_ts_out = FALSE,
                              show_progress = FALSE,
                              n_threads = 1L) {

    # This is like match.call but includes default arguments and
    # evaluates everything
    restart_call <- mget(names(formals()), sys.frame(sys.nframe()))
    call_ <- sims_obj$call
    for (n in c("max_t", "perturb", "save_every", "sep_adults", "stage_ts_out",
                "show_progress", "n_threads")) {
        if (n %in% names(restart_call)) call_[[n]] <- restart_call[[n]]
    }

    if (!inherits(sims_obj$all_info_xptr, "externalptr")) {
        stop("\nSomething has happened to the \"all_info_xptr\" field ",
             "in this `cloneSims` object. Please re-run simulation.")
    }

    n_reps <- sims_obj$call$n_reps
    if (is.null(n_reps)) n_reps <- formals(sim_fields)[["n_reps"]]

    uint_check(n_reps, "n_reps", .min = 1)
    uint_check(max_t, "max_t", .min = 1)
    lgl_check(no_warns, "no_warns")
    uint_check(save_every, "save_every", .min = 1)
    lgl_check(sep_adults, "sep_adults")
    lgl_check(stage_ts_out, "stage_ts_out")
    lgl_check(show_progress, "show_progress")
    uint_check(n_threads, "n_threads", .min = 1)

    if (!no_warns && !is.null(sims_obj$call$perturb) && is.null(perturb)) {
        warning(paste("\nThe initial simulations had perturbations, but the",
                      "new one won't. Provide a dataframe to the `perturb`",
                      "argument to `restart.cloneSims` to add",
                      "perturbations."))
    }

    perturb_ptr <- make_perturb_ptr_R(perturb, call_$fields)

    # Create pointer to a copied C++ object with new arguments if requested.
    # Start with those from original simulations:
    restart_args <- get_field_pars(sims_obj$all_info_xptr)
    if (!is.null(new_field_pars)) {
        if (!inherits(new_field_pars, "list") || is.null(names(new_field_pars))) {
            stop("\n`new_field_pars` must be a named list")
        }
        restart_args[["K_y_mult"]] <- restart_args[["K_y"]] / restart_args[["K"]]
        restart_args[["K_y"]] <- NULL
        if (!all(names(new_field_pars) %in% names(restart_args))) {
            bn <- names(new_field_pars)[!names(new_field_pars) %in% names(restart_args)]
            stop(sprintf(paste("\nThis option is unavailable for the `restart_args`",
                               "argument: '%s'.\nThese are the ones to choose from: %s."),
                         bn, paste(names(restart_args), collapse = ", ")))
        }
        for (n in names(new_field_pars)) {
            z <- new_field_pars[[n]]
            if (length(z) == 1 && length(restart_args[[n]]) > 1) {
                z <- rep(z, length(restart_args[[n]]))
            }
            if (length(z) != length(restart_args[[n]])) {
                stop(sprintf("new_field_pars[[%s]] should be length 1 or %i",
                             n, length(restart_args[[n]])))
            }
            if (class(z) != class(restart_args[[n]])) {
                stop(sprintf("new_field_pars[[%s]] should be of class %s",
                             n, class(restart_args[[n]])))
            }
            restart_args[[n]] <- z
        }
        restart_args[["K_y"]] <- restart_args[["K_y_mult"]] * restart_args[["K"]]
        restart_args[["K_y_mult"]] <- NULL
    }
    restart_args[["all_fields_in_ptr"]] <- sims_obj$all_info_xptr
    new_all_fields_vec_ptr <- do.call(restarted_field_ptr, restart_args)

    # If `new_starts` is provided, check it carefully then adjust the
    # underlying C++ objects accordingly.
    if (!is.null(new_starts)) {
        check_new_starts_adjust_N(sims_obj, new_starts, new_all_fields_vec_ptr)
    }

    sims <- sim_fields_cpp(new_all_fields_vec_ptr, perturb_ptr,
                           max_t, save_every, sep_adults,
                           n_threads, show_progress, stage_ts_out)

    sims <- make_cloneSims_obj(sims, call_, stage_ts_out, restarted = TRUE)

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
    if (x$restarted) cat("** these sims were restarted, so starting abundances",
                         "in the `call` field are not accurate\n")
    invisible(x)
}
