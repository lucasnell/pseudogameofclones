
library(tidyverse)
library(clonewars)
library(parallel)
library(patchwork)


if (file.exists(".Rprofile")) source(".Rprofile")

options(mc.cores = max(parallel::detectCores()-2L, 1L))



#' These simulations are to plan and create a priori hypotheses for the
#' eco-evo experiments.
#' I started with a model from a
#' [previous paper](http://doi.wiley.com/10.1890/13-1933.1)
#' from Tony and others, then I edited it to conform to our conditions
#' (most notably, it explicitly simulates alate (winged aphid) production
#' and plant death).
#' I won't be simulating plant death because simulations show that this only
#' adds complexity without changing outcomes.
#' I'll also have alate production be density-independent because previous
#' assays in our fields indicate density has little effect.
#'
#'
#' I'm simulating 1,000 days of interacting aphid and parasitoid wasp
#' populations.
#' There are two fields, one with wasps and another without.
#' Movement between fields must happen via alates, and
#' rates of alate production are density-independent.
#' Every day, 10% of all alates in each field disperse to the other field.
#' There are two aphid lines:
#' The susceptible line is susceptible to parasitism but has a higher
#' growth rate.
#' The resistant line is resistant to parasitism but has a lower growth rate.
#' Shifts in relative frequencies of susceptible and resistant lines
#' represents evolution in the aphid population.
#'



# Define clonal lines. Don't re-define these!
# Susceptible line: no resistance, high population growth rate
line_s <- clonal_line("susceptible",
                     density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
                     surv_juv_apterous = "high",
                     surv_adult_apterous = "high",
                     repro_apterous = "high")
# Resistant line: high resistance, low population growth rate
line_r <- clonal_line("resistant",
                      density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
                      resistant = TRUE,
                      surv_paras = 0.57,
                      surv_juv_apterous = "low",
                      surv_adult_apterous = "low",
                      repro_apterous = "low")

para_lvls <- paste(c("no parasitism", "parasitism"), "patch")

# Palette for the two clonal lines.
# Equivalent to `viridis::viridis(100)[c(70, 10)]`.
clone_pal <- c("#41BE71FF", "#482173FF")


# To maintain the same y-axis limits:
# wasp_mod <- 0.0164459  # <-- should be max(disp_wasps$wasps) / max(disp_aphids$N)
# ylims <- c(0, 2405)  # <-- max is ceiling(max(disp_aphids$N))
# y_breaks <- 0:2 * 1e3
# y_labs <- paste0(0:2, "k")
wasp_mod <- 5.078307  # <-- should be max(disp_wasps$wasps) / max(log1p(disp_aphids$N))
ylims <- c(0, 7.785)  # <-- max is ceiling(max(log(disp_aphids$N)) * 1e3) / 1e3
y_breaks <- log(10^(0:3))
y_labs <- 10^(0:3)




#'
#' Base function to do all the simulations.
#' Each section below builds on this one.
#'
do_base_sims <- function(...) {

    .args <- list(clonal_lines = c(line_s, line_r))
    .other_args <- list(...)
    if (length(.other_args) > 0) {
        stopifnot(!is.null(names(.other_args)))
        stopifnot(sum(duplicated(names(.other_args))) == 0)
        stopifnot(all(names(.other_args) %in% names(formals(sim_experiments))))
        for (n in names(.other_args)) .args[[n]] <- .other_args[[n]]
    }

    .sims <- do.call(sim_experiments, .args)

    .sims[["wasps"]] <- .sims %>%
        .[["wasps"]] %>%
        select(-rep) %>%
        mutate(field = factor(field, levels = 1:0,
                             labels = para_lvls))

    .sims[["aphids"]] <- .sims %>%
        .[["aphids"]] %>%
        select(-rep, -plant) %>%
        mutate(field = factor(field, levels = 1:0,
                             labels = para_lvls),
               line = factor(line, levels = c("resistant", "susceptible")))

    return(.sims)
}


#' Group aphids by type.
#' This never groups by line unless you're requesting mummies or all
#' parasitized aphids (which includes mummies).
#'
#' This is a separate function from `do_base_sims` so that output from
#' `do_base_sims` can be used to look at alates, etc.
#'
#' `.sim_aphids_df` should be the `"aphids"` field from `sim_clonewars` output
#'
#' `.type` takes the following options:
#'     * `"living"`: all aphids including parasitized;
#'         this is better for plotting aphid abundances
#'     * `"unparasitized"`: unparasitized aphids;
#'         this is better for community matrices
#'     * `"mummies"`: just mummies
#'     * `"parasitized"`: mummies + parasitized aphids;
#'         this is better for community matrices
#'     * `"alates"`: alates
#'     * `"apterous"`: apterous
#'
group_aphids <- function(.sim_aphids_df, .type) {

    stopifnot(is.data.frame(.sim_aphids_df))

    needed_cols <- c("field", "time", "line", "type", "N")
    stopifnot(all(needed_cols %in% colnames(.sim_aphids_df)))

    # Other columns to potentially keep:
    other_cols <- colnames(.sim_aphids_df)
    other_cols <- other_cols[!other_cols %in% needed_cols]
    keep_other_cols <- c()
    if (length(other_cols) > 0) {
        for (.c in other_cols) {
            # Only keep them if they're identical for the entire data frame
            if (length(unique(.sim_aphids_df[[.c]])) == 1) {
                keep_other_cols <- c(keep_other_cols, .c)
            }
        }
    }

    if (.type == "living") {
        .out <- .sim_aphids_df %>%
            filter(type != "mummy") %>%
            group_by(field, time, line) %>%
            summarize(N = sum(N), .groups = "drop")
    } else if (.type == "unparasitized") {
        .out <- .sim_aphids_df %>%
            filter(type != "mummy", type != "parasitized") %>%
            group_by(field, time, line) %>%
            summarize(N = sum(N), .groups = "drop")
    } else if (.type == "mummies") {
        .out <- .sim_aphids_df %>%
            filter(type == "mummy") %>%
            group_by(field, time) %>%
            summarize(N = sum(N), .groups = "drop")
    } else if (.type == "parasitized") {
        .out <- .sim_aphids_df %>%
            filter(type == "mummy" | type == "parasitized") %>%
            group_by(field, time) %>%
            summarize(N = sum(N), .groups = "drop")
    } else if (.type == "alates") {
        .out <- .sim_aphids_df %>%
            filter(type == "alate") %>%
            group_by(field, time, line) %>%
            summarize(N = sum(N), .groups = "drop")
    } else if (.type == "apterous") {
        .out <- .sim_aphids_df %>%
            filter(type == "apterous") %>%
            group_by(field, time, line) %>%
            summarize(N = sum(N), .groups = "drop")
    } else {
        stop("ERROR: strange .type requested. Options are ",
             paste(c('"living"', '"mummies"', '"unparasitized"',
                     '"parasitized"', '"alates"', '"apterous"'),
                   collapse = ", "))
    }

    for (koc in keep_other_cols) .out[[koc]] <- .sim_aphids_df[[koc]][1]

    return(.out)

}

# ============================================================================*
# ============================================================================*

# Main simulations ----

# ============================================================================*
# ============================================================================*


main_sims <- rep(list(NA), 2)
names(main_sims) <- c("no aphid dispersal", "aphid dispersal")
main_sims[["no aphid dispersal"]] <- do_base_sims(alate_field_disp_p = 0)
main_sims[["aphid dispersal"]] <- do_base_sims()
for (d in names(main_sims)) {
    for (n in c("wasps", "aphids")) {
        main_sims[[d]][[n]] <- main_sims[[d]][[n]] %>%
            mutate(disp = d)
    }
}


main_aphids <- map_dfr(main_sims, ~ group_aphids(.x[["aphids"]], "living")) %>%
    mutate(disp = factor(disp, levels = names(main_sims)))
main_mummies <- map_dfr(main_sims, ~ group_aphids(.x[["aphids"]], "mummies")) %>%
    mutate(disp = factor(disp, levels = names(main_sims)))
main_wasps <- map_dfr(main_sims, ~ .x[["wasps"]]) %>%
    mutate(disp = factor(disp, levels = names(main_sims)))


main_p_list <- map(
    levels(main_aphids$disp),
    function(d) {
        mad <- main_aphids %>%
            filter(disp == d)
        wad <- main_wasps %>%
            filter(disp == d)
        mad %>%
            mutate(N = ifelse(N == 0, NA, N),
                   N = log(N)) %>%
            ggplot(aes(time, N)) +
            ggtitle(d) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_area(data = wad %>%
                          mutate(N = wasps / wasp_mod),
                      fill = "gray60", color = NA) +
            geom_line(aes(color = line)) +
            scale_color_manual(values = clone_pal, guide = "none") +
            scale_y_continuous("Aphid abundance",
                               breaks = y_breaks,
                               limits = ylims,
                               labels = y_labs,
                               sec.axis = sec_axis(~ . * wasp_mod,
                                                   "Wasp abundance",
                                                   breaks = 0:2 * 20)) +
            scale_x_continuous("Days") +
            facet_wrap( ~ field, nrow = 1) +
            theme(strip.text = element_text(size = 10),
                  plot.title = element_text(size = 12, hjust = 0.5,
                                            face = "bold",
                                            margin = margin(0,0,0,b=6)))
    })


# main_p_list[[1]] <- main_p_list[[1]] +
#     geom_text(data = tibble(line = factor(c("susceptible", "resistant")),
#                             field = factor(para_lvls[2], levels = para_lvls),
#                             time = 250,
#                             N = log(c(250, 1700))),
#               aes(label = line, color = line),
#               size = 9 / 2.8, hjust = 1, vjust = c(0, 1)) +
#     geom_text(data = tibble(field = factor(para_lvls[2], levels = para_lvls),
#                             time = 80, N = log(1.1)),
#               aes(label = "wasps"), size = 9 / 2.8, hjust = 0, vjust = 0,
#               color = "gray50")


main_p <- wrap_plots(main_p_list, ncol = 1) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))

# main_p

# save_plot("_results/plots/sims-main.pdf", main_p, 6, 5)


# If interested in discussing proportion of aphid dispersing
# (to compare to experiments):
main_sims[["aphid dispersal"]][["aphids"]] %>%
    filter(type != "mummy") %>%
    pivot_wider(names_from = type, values_from = N) %>%
    mutate(total = alate + apterous + parasitized) %>%
    group_by(time, field, line) %>%
    summarize(p_disp = 0.1 * alate / total, .groups = "drop") %>%
    # .[["p_disp"]] %>% mean()
    ggplot(aes(time, p_disp)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_line(aes(color = line)) +
    scale_color_manual(values = clone_pal, guide = "none") +
    scale_y_continuous("Proportion dispersed") +
    scale_x_continuous("Days") +
    facet_wrap( ~ field)




# ============================================================================*
# ============================================================================*

# Starting wasp abundance ----

# ============================================================================*
# ============================================================================*


#'
#' Do a call to `sim_clonewars` that can vary by the proportion of alates
#' that move between fields every day (`alate_field_disp_p`).
#' How does dispersal affect maintenance of diversity in the resistance trait?
#'
do_wasp_sims <- function(.wasp_density_0, .max_t = 500) {

    stopifnot(is.numeric(.wasp_density_0) && length(.wasp_density_0) == 1 &&
                  .wasp_density_0 > 0)

    .sims <- do_base_sims(wasp_density_0 = c(.wasp_density_0, 0),
                          max_t = .max_t)

    for (n in c("wasps", "aphids")) .sims[[n]] <- .sims[[n]] %>%
            mutate(wasps0 = .wasp_density_0)

    return(.sims)
}



wasp_sims <- mclapply(c(1, 2, 3, 10), do_wasp_sims)


wasp_aphids <- map_dfr(wasp_sims, ~ group_aphids(.x[["aphids"]], "living")) %>%
    mutate(wasps0 = factor(wasps0, labels = sprintf("starting wasp(s) = %i",
                                                      sort(unique(wasps0)))))
wasp_mummies <- map_dfr(wasp_sims, ~ group_aphids(.x[["aphids"]], "mummies")) %>%
    mutate(wasps0 = factor(wasps0, labels = sprintf("starting wasp(s) = %i",
                                                      sort(unique(wasps0)))))
wasp_wasps <- map_dfr(wasp_sims, ~ .x[["wasps"]]) %>%
    mutate(wasps0 = factor(wasps0, labels = sprintf("starting wasp(s) = %i",
                                                      sort(unique(wasps0)))))


wasp_p_list <- map(
    levels(wasp_aphids$wasps0),
    function(w0) {
        aw0 <- wasp_aphids %>%
            filter(wasps0 == w0, time <= 250)
        ww0 <- wasp_wasps %>%
            filter(wasps0 == w0, time <= 250)
        aw0 %>%
            mutate(N = ifelse(N == 0, NA, N),
                   N = log(N)) %>%
            ggplot(aes(time, N)) +
            ggtitle(w0) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_area(data = ww0 %>%
                          mutate(N = wasps / wasp_mod),
                      fill = "gray60", color = NA) +
            geom_line(aes(color = line)) +
            scale_color_manual(values = clone_pal, guide = "none") +
            scale_y_continuous("Aphid abundance",
                               breaks = y_breaks,
                               limits = ylims,
                               labels = y_labs,
                               sec.axis = sec_axis(~ . * wasp_mod,
                                                   "Wasp abundance",
                                                   breaks = 0:2 * 20)) +
            scale_x_continuous("Days") +
            facet_wrap( ~ field, nrow = 1) +
            theme(strip.text = element_text(size = 10),
                  plot.title = element_text(size = 12, hjust = 0.5))
    })


# wasp_p_list[[1]] <- wasp_p_list[[1]] +
#     geom_text(data = tibble(line = factor(c("resistant", "susceptible")),
#                             field = factor(para_lvls[1], levels = para_lvls),
#                             time = 250, N = c(300, max(wasp_aphids$N) - 300)),
#               aes(label = line, color = line),
#               size = 9 / 2.8, hjust = 1, vjust = c(0, 1)) +
#     geom_text(data = tibble(field = factor(para_lvls[2], levels = para_lvls),
#                             time = 230, N = 700),
#               aes(label = "wasps"), size = 9 / 2.8, hjust = 1, vjust = 0.5,
#               color = "gray50")



wasp_p <- wrap_plots(wasp_p_list, ncol = 2) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))


# save_plot("_results/plots/sims-wasp.pdf", wasp_p, 10, 5)




# ============================================================================*
# ============================================================================*

# Vary dispersal ----

# ============================================================================*
# ============================================================================*


#'
#' Do a call to `sim_clonewars` that can vary by the proportion of alates
#' that move between fields every day (`alate_field_disp_p`).
#' How does dispersal affect maintenance of diversity in the resistance trait?
#'
do_disp_sims <- function(.alate_field_disp_p, .max_t = 5000) {

    .sims <- do_base_sims(alate_field_disp_p = .alate_field_disp_p,
                          max_t = .max_t)

    for (n in c("wasps", "aphids")) .sims[[n]] <- .sims[[n]] %>%
        mutate(alate_p = .alate_field_disp_p)

    return(.sims)
}



disp_sims <- mclapply(c(0.05, 0.1, 0.2, 0.4), do_disp_sims)


disp_aphids <- map_dfr(disp_sims, ~ group_aphids(.x[["aphids"]], "living")) %>%
    mutate(alate_p = factor(alate_p, labels = sprintf("%.0f%% dispersal",
                                                    sort(unique(alate_p)) * 100)))
disp_mummies <- map_dfr(disp_sims, ~ group_aphids(.x[["aphids"]], "mummies")) %>%
    mutate(alate_p = factor(alate_p, labels = sprintf("%.0f%% dispersal",
                                                    sort(unique(alate_p)) * 100)))
disp_wasps <- map_dfr(disp_sims, ~ .x[["wasps"]]) %>%
    mutate(alate_p = factor(alate_p, labels = sprintf("%.0f%% dispersal",
                                                    sort(unique(alate_p)) * 100)))


disp_p_list <- map(
    levels(disp_aphids$alate_p),
    function(ap) {
        aap <- disp_aphids %>%
            filter(alate_p == ap, time <= 500) %>%
            mutate(N = ifelse(N == 0, NA, N),
                   N = log(N))
        wap <- disp_wasps %>%
            filter(alate_p == ap, time <= 500)
        p <- aap %>%
            filter(!is.na(N)) %>%
            ggplot(aes(time, N)) +
            ggtitle(ap) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_area(data = wap %>%
                          mutate(N = wasps / wasp_mod),
                      fill = "gray60", color = NA) +
            geom_line(aes(color = line)) +
            scale_color_manual(values = clone_pal, guide = "none") +
            scale_y_continuous("Aphid abundance",
                               breaks = y_breaks,
                               limits = ylims,
                               labels = y_labs,
                               sec.axis = sec_axis(~ . * wasp_mod,
                                                   "Wasp abundance",
                                                   breaks = 0:2 * 20)) +
            scale_x_continuous("Days") +
            facet_wrap( ~ field, nrow = 1) +
            theme(strip.text = element_text(size = 10),
                  plot.title = element_text(size = 12, hjust = 0.5))
        if (any(is.na(aap$N))) {
            ext <- filter(aap, is.na(N)) %>%
                group_by(field) %>%
                filter(time == min(time)) %>%
                ungroup() %>%
                mutate(N = 0)
            p <- p +
                geom_point(data = ext, aes(color = line), shape = 4, size = 2)
        }
        return(p)
    })


# disp_p_list[[1]] <- disp_p_list[[1]] +
#     geom_text(data = tibble(line = factor(c("resistant", "susceptible")),
#                             field = factor(para_lvls[2], levels = para_lvls),
#                             time = c(500, 500), N = c(300, max(disp_aphids$N) - 300)),
#               aes(label = line, color = line),
#               size = 9 / 2.8, hjust = 1, vjust = c(0, 1))


# disp_p_list[[4]] <- disp_p_list[[4]] +
#     geom_text(data = tibble(field = factor(para_lvls[2], levels = para_lvls),
#                             time = 400, N = 1500),
#               aes(label = "wasps"), size = 9 / 2.8, hjust = 1, vjust = 0,
#               color = "gray50")



disp_p <- wrap_plots(disp_p_list, ncol = 2) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))


# save_plot("_results/plots/sims-disp.pdf", disp_p, 10, 5)






# ============================================================================*
# ============================================================================*

# Stability - starts ----

# ============================================================================*
# ============================================================================*


#'
#' Vary the proportion of starting aphids that are resistant to parasitism.
#' How stable are results to starting conditions?
#'
do_stable_start_sims <- function(.p_res, .max_t = 5000) {

    # .p_res = 0.25
    # rm(.p_res, .n, .line_s, .line_r)

    if (round(sum(line_s$density_0[line_s$density_0 > 0]), 10) !=
        round(sum(line_r$density_0[line_r$density_0 > 0]), 10)) {
        stop("ERROR: starting densities in `line_s` and `line_r` should be ",
             "the same. (Let the later functions change this.)")
    }

    .n <- 2 * sum(line_s$density_0[line_s$density_0 > 0])

    .line_s <- line_s
    .line_r <- line_r

    # Make both densities sum to 1
    .line_s$density_0 <- .line_s$density_0 / sum(.line_s$density_0)
    .line_r$density_0 <- .line_r$density_0 /  sum(.line_r$density_0)

    # Now make them sum to correct abundances:
    .line_s$density_0 <- .line_s$density_0 * .n * (1 - .p_res)
    .line_r$density_0 <- .line_r$density_0 * .n * .p_res

    .sims <- do_base_sims(clonal_lines = c(.line_r, .line_s),
                          max_t = .max_t)

    for (n in c("wasps", "aphids")) .sims[[n]] <- .sims[[n]] %>%
        mutate(p_res = .p_res)

    return(.sims)
}




p_res_lvls <- c(0.35, 0.4, 0.8, 0.85)

stable_start_sims <- mclapply(p_res_lvls, do_stable_start_sims)


stable_start_aphids <- map_dfr(stable_start_sims,
                               ~ group_aphids(.x[["aphids"]], "living")) %>%
    mutate(p_res = factor(p_res, labels = sprintf("%.0f%% starting resistance",
                                                  sort(unique(p_res)) * 100)))
stable_start_mummies <- map_dfr(stable_start_sims,
                                ~ group_aphids(.x[["aphids"]], "mummies")) %>%
    mutate(p_res = factor(p_res, labels = sprintf("%.0f%% starting resistance",
                                                  sort(unique(p_res)) * 100)))
stable_start_wasps <- map_dfr(stable_start_sims, ~ .x[["wasps"]]) %>%
    mutate(p_res = factor(p_res, labels = sprintf("%.0f%% starting resistance",
                                                  sort(unique(p_res)) * 100)))



stable_start_p_list <- map(
    levels(stable_start_aphids$p_res),
    function(pr) {
        ssa <- stable_start_aphids %>%
            filter(p_res == pr, time <= 500) %>%
            mutate(N = ifelse(N == 0, NA, N),
                   N = log(N))
        ssw <- stable_start_wasps %>%
            filter(p_res == pr, time <= 500)
        p <- ssa %>%
            filter(!is.na(N)) %>%
            ggplot(aes(time, N)) +
            ggtitle(pr) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_area(data = ssw %>%
                          mutate(N = wasps / wasp_mod),
                      fill = "gray60", color = NA) +
            geom_line(aes(color = line)) +
            scale_color_manual(values = clone_pal, guide = "none") +
            scale_y_continuous("Aphid abundance",
                               breaks = y_breaks,
                               limits = ylims,
                               labels = y_labs,
                               sec.axis = sec_axis(~ . * wasp_mod,
                                                   "Wasp abundance",
                                                   breaks = 0:2 * 20)) +
            scale_x_continuous("Days") +
            facet_wrap(~ field, nrow = 1) +
            theme(strip.text = element_text(size = 10),
                  plot.title = element_text(size = 12, hjust = 0.5))
        if (any(is.na(ssa$N))) {
            ext <- filter(ssa, is.na(N)) %>%
                group_by(field) %>%
                filter(time == min(time)) %>%
                ungroup() %>%
                mutate(N = 0)
            p <- p +
                geom_point(data = ext, aes(color = line), shape = 4, size = 2)
        }
        return(p)
    })

# stable_start_p_list[[1]] <- stable_start_p_list[[1]] +
#     geom_text(data = tibble(line = factor(c("resistant", "susceptible")),
#                             field = factor(para_lvls[1], levels = para_lvls),
#                             time = c(500, 500),
#                             N = c(300, max(stable_start_aphids$N) - 300)),
#               aes(label = line, color = line),
#               size = 9 / 2.8, hjust = 1, vjust = c(0, 1)) +
#     geom_text(data = tibble(field = factor(para_lvls[2], levels = para_lvls),
#                             time = 400, N = 800),
#               aes(label = "wasps"), size = 9 / 2.8, hjust = 0.5, vjust = 0.5,
#               color = "gray50")

stable_start_p <- wrap_plots(stable_start_p_list, ncol = 2) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))

# save_plot("_results/plots/sims-stable-start.pdf", stable_start_p, 10, 5)











# ============================================================================*
# ============================================================================*

# Stability - perturbs ----

# ============================================================================*
# ============================================================================*




#'
#' Vary the proportion of starting aphids that are resistant to parasitism.
#' How stable are results to starting conditions?
#'
do_stable_perturb_sims <- function(.perturb_df) {

    .max_t = max(.perturb_df$when) + 5e3

    .sims <- do_base_sims(max_t = .max_t,
                          perturb = .perturb_df)

    for (n in c("wasps", "aphids")) .sims[[n]] <- .sims[[n]] %>%
        mutate(who = paste(unique(.perturb_df$who), collapse = " & "))

    return(.sims)
}

stable_perturb_sims <- list(tibble(when = 10e3, where = 1:2,
                                   who = "resistant", how = 0.5),
                            tibble(when = 10e3, where = 1:2,
                                   who = "susceptible", how = 0.5),
                            tibble(when = 10e3, where = 1:2,
                                   who = "wasps", how = 0.5),
                            tibble(when = 10e3, where = 1:2,
                                   who = "mummies", how = 0.5)) %>%
    mclapply(do_stable_perturb_sims)

pert_t_range <- c(9.8, 10.5) * 1e3

stable_perturb_aphids <- map_dfr(stable_perturb_sims,
                                 ~ group_aphids(.x[["aphids"]], "living")) %>%
    mutate(who = factor(who, levels = c("resistant", "susceptible",
                                        "wasps", "mummies"),
                        labels = c("resistant aphids", "susceptible aphids",
                                   "adult wasps",
                                   "mummies & parasitized aphids"))) %>%
    filter(time >= pert_t_range[1], time <= pert_t_range[2])
stable_perturb_mummies <- map_dfr(stable_perturb_sims,
                                  ~ group_aphids(.x[["aphids"]], "mummies")) %>%
    mutate(who = factor(who, levels = c("resistant", "susceptible",
                                        "wasps", "mummies"),
                        labels = c("resistant aphids", "susceptible aphids",
                                   "adult wasps",
                                   "mummies & parasitized aphids"))) %>%
    filter(time >= pert_t_range[1], time <= pert_t_range[2])
stable_perturb_wasps <- map_dfr(stable_perturb_sims, ~ .x[["wasps"]]) %>%
    mutate(who = factor(who, levels = c("resistant", "susceptible",
                                        "wasps", "mummies"),
                        labels = c("resistant aphids", "susceptible aphids",
                                   "adult wasps",
                                   "mummies & parasitized aphids"))) %>%
    filter(time >= pert_t_range[1], time <= pert_t_range[2])





stable_perturb_p_list <- map(
    levels(stable_perturb_aphids$who),
    function(w) {
        ssa <- stable_perturb_aphids %>%
            filter(who == w)
        ssw <- stable_perturb_wasps %>%
            filter(who == w)
        ssa %>%
            mutate(N = ifelse(N == 0, NA, N),
                   N = log(N)) %>%
            ggplot(aes(time, N)) +
            ggtitle(w) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_area(data = ssw %>%
                          mutate(N = wasps / wasp_mod),
                      fill = "gray60", color = NA) +
            geom_line(aes(color = line)) +
            scale_color_manual(values = clone_pal, guide = "none") +
            scale_y_continuous("Aphid abundance",
                               breaks = y_breaks,
                               limits = ylims,
                               labels = y_labs,
                               sec.axis = sec_axis(~ . * wasp_mod,
                                                   "Wasp abundance",
                                                   breaks = 0:2 * 20)) +
            scale_x_continuous("Days") +
            facet_wrap(~ field, nrow = 1) +
            theme(strip.text = element_text(size = 10),
                  plot.title = element_text(size = 12, hjust = 0.5))
    })

# stable_perturb_p_list[[1]] <- stable_perturb_p_list[[1]] +
#     geom_text(data = tibble(line = factor(c("resistant", "susceptible")),
#                             field = factor(para_lvls[1], levels = para_lvls),
#                             time = 10.5e3,
#                             N = c(420, max(stable_perturb_aphids$N) - 300)),
#               aes(label = line, color = line),
#               size = 9 / 2.8, hjust = 1, vjust = c(0, 1)) +
#     geom_text(data = tibble(field = factor(para_lvls[2], levels = para_lvls),
#                             time = 10200, N = 110),
#               aes(label = "wasps"), size = 9 / 2.8, hjust = 0, vjust = 0.5,
#               color = "gray50")


stable_perturb_p <- wrap_plots(stable_perturb_p_list, ncol = 2) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))

# save_plot("_results/plots/sims-stable-perturb.pdf", stable_perturb_p, 10, 5)







# ============================================================================*
# ============================================================================*

# Community matrix ----

# ============================================================================*
# ============================================================================*


#'
#' Calculate community matrix numerically.
#'
.step_size = 1e-5
.perturb_time = 100e3
do_comm_mat_perturb_sims <- function(.who, .where) {

    # .who = "resistant"; .where = 1

    # Always perturb all wasps and future wasps
    ..who <- .who
    if (.who == "parasitized") ..who <- c("wasps", "mummies")

    .perturb_df <- tibble(when = .perturb_time, where = .where,
                          who = ..who,
                          how = 1 + .step_size)

    .max_t = max(.perturb_df$when) + 100

    .sims <- do_base_sims(max_t = .max_t,
                          perturb = .perturb_df)


    for (n in c("wasps", "aphids")) .sims[[n]] <- .sims[[n]] %>%
        filter(time > (min(.perturb_df$when) - 100)) %>%
        mutate(perturb_who = paste(.who, collapse = " & "),
               perturb_where = .where)

    return(.sims)
}

# Takes ~3 min
comm_mat_perturb_sims <- crossing(.who = c("resistant", "susceptible",
                                           "parasitized"), .where = 1:2) %>%
    # no need to perturb the non-existent parasitized aphids in the
    # no parasitism field
    filter(!(.who == "parasitized" &
                 .where == which(grepl("^no", para_lvls)))) %>%
    as.list() %>%
    c(list(FUN = do_comm_mat_perturb_sims, SIMPLIFY = FALSE)) %>%
    do.call(what = mcmapply)


comm_mat_perturb_df <- list(
    # non-parasitized aphids:
    map_dfr(comm_mat_perturb_sims,
            ~ group_aphids(.x[["aphids"]], "unparasitized")) %>%
        mutate(line = paste(line)) %>%
        rename(type = line),
    # wasps + parasitized + mummies:
    map_dfr(
        seq_along(comm_mat_perturb_sims),
        function(i) {
            a <- group_aphids(comm_mat_perturb_sims[[i]][["aphids"]],
                              "parasitized")
            b <- comm_mat_perturb_sims[[i]][["wasps"]]
            full_join(a, b, by = c("field", "time", "perturb_who",
                                   "perturb_where")) %>%
                mutate(N = N + wasps, type = "parasitized") %>%
                select(-wasps)
        })) %>%
    bind_rows() %>%
    mutate(perturb_where = factor(perturb_where, levels = 2:1,
                                  labels = para_lvls)) %>%
    rename(who = type, where = field)

# Order of rows and columns:
cm_order <- tibble(
    who = rep(c("parasitized", "resistant", "susceptible"), each = 2),
    where = rep(para_lvls, 3)) %>%
    # removing this because it's always zero:
    filter(!(who == "parasitized" & where == "no parasitism"))

#'
#' One column for the community matrix: ∂ y / ∂ x for one x and all y
#'
one_cm_col <- function(x_who, x_where) {
    .df <- comm_mat_perturb_df %>%
        filter(perturb_who == x_who, perturb_where == x_where) %>%
        select(-perturb_who, -perturb_where) %>%
        filter(time %in% (.perturb_time - 1:0)) %>%
        arrange(where, who, time)
    #' The step used was multiplicative, where new values were calculated at
    #' `f( x * (s+1) ) == f( x + x * s )`, so we calculate `x * s` here.
    dx <- filter(.df, who == x_who, where == x_where)[["N"]][1] * .step_size
    dy <- map2_dbl(cm_order$who, cm_order$where,
                   function(.x, .y) {
                       filter(.df, who == .x, where == .y)[["N"]] %>%
                           diff()
                   })
    dy_dx <- cbind(dy / dx)
    return(dy_dx)
}


comm_mat <- map2(cm_order$who, cm_order$where, one_cm_col) %>%
    do.call(what = cbind)

comm_mat

eigen(comm_mat)$values

# # To share with Tony:
# saveRDS(comm_mat, "~/Desktop/eco-evo_comm_mat.rds")


