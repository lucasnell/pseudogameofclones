
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
                     density_0 = matrix(c(rep(0, 3), 32, rep(0, 6)), 5, 2),
                     surv_juv_apterous = "high",
                     surv_adult_apterous = "high",
                     repro_apterous = "high")
# Resistant line: high resistance, low population growth rate
line_r <- clonal_line("resistant",
                      density_0 = matrix(c(rep(0, 3), 32, rep(0, 6)), 5, 2),
                      resistant = TRUE,
                      surv_paras = 0.57,
                      surv_juv_apterous = "low",
                      surv_adult_apterous = "low",
                      repro_apterous = "low")

para_lvls <- c("parasitism", "no parasitism")



#'
#' Base function to do all the simulations.
#' Each section below builds on this one.
#'
do_base_sims <- function(clonal_lines, ...) {

    # clonal_lines = c(line_s, line_r)
    # rm(clonal_lines, .sims, .args)

    .args <- list(clonal_lines = clonal_lines)
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
        mutate(field = factor(field, levels = 0:max(field),
                             labels = para_lvls))

    .sims[["aphids"]] <- .sims %>%
        .[["aphids"]] %>%
        select(-rep, -plant) %>%
        mutate(field = factor(field, levels = 0:max(field),
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

# Match experiments ----

# ============================================================================*
# ============================================================================*




plot_test <- function(.sims) {

    # aphids <- .sims[["aphids"]] %>%
    #     filter(type != "mummy") %>%
    #     mutate(type = factor(type, levels = c("apterous", "alate", "parasitized")))
    aphids <- group_aphids(.sims[["aphids"]], "living")
    mummies <- .sims[["aphids"]] %>% filter(type == "mummy")
    wasps <- .sims[["wasps"]]

    mod <- max(wasps$wasps) / max(aphids$N)

    aphids %>%
        ggplot(aes(time, N / 1e3)) +
        geom_hline(yintercept = 0, color = "gray70") +
        geom_area(data = wasps %>%
                      mutate(N = wasps / mod),
                  fill = "gray80", color = NA) +
        # geom_line(aes(color = line, linetype = type)) +
        geom_line(aes(color = line), size = 1) +
        scale_color_manual(values = c("chartreuse3", "firebrick"), guide = "none") +
        scale_linetype_manual(values = c("solid", "33", "11")) +
        scale_y_continuous(expression("Aphid abundance (" %*% 1000 * ")"),
                           sec.axis = sec_axis(~ . * mod * 1000,
                                               "Wasp abundance")) +
        scale_x_continuous("Days") +
        facet_wrap(~ field, ncol = 1) +
        theme(strip.text = element_text(size = 10),
              legend.position = "top") +  # ,
              # panel.grid.major.x = element_line(size = 0.25, color = "gray90"),
              # panel.grid.major.y = element_line(size = 0.25, color = "gray90")) +
        # coord_cartesian(ylim = c(0, 3)) +
        NULL

}


#'
#' At day 8 of prelim assays, mean(N_res / N_sus) was 0.219.
#' At day 12 of prelim assays, mean(N_res / N_sus) was 0.527.
#' The overall average was 0.373.
#'
#' Average aphids for 3 fava bean plants inside a single 8" pot at 8 and 12 days for
#' each line:
#'
#'    date line          num
#'   <int> <fct>       <dbl>
#' 1     8 resistant    11.5
#' 2     8 susceptible  47.8
#' 3    12 resistant    23
#' 4    12 susceptible  76.8
#'


line_s2 <- clonal_line(name = "susceptible",
                       # >>>>>>>>>>>>> things to potentially change:
                       # surv_paras = 0.8,
                       # <<<<<<<<<<<<<
                       surv_juv_apterous = "high",
                       surv_adult_apterous = "high",
                       repro_apterous = "high",
                       # density_0 = matrix(c(rep(0, 3), 32, rep(0, 6)), 5, 2))
                       density_0 = matrix(c(rep(0, 3), 250, rep(0, 6)), 5, 2))
                       # density_0 = matrix(c(rep(0, 3), 410, rep(0, 6)), 5, 2))
# Resistant line: high resistance, low population growth rate
line_r2 <- clonal_line(name = "resistant",
                       resistant = TRUE,
                       surv_paras = 0.57,
                       surv_juv_apterous = "low",
                       surv_adult_apterous = "low",
                       repro_apterous = "low",
                       # <<<<<<<<<<<<<
                       # density_0 = matrix(c(rep(0, 3), 32, rep(0, 6)), 5, 2))
                       density_0 = matrix(c(rep(0, 3), 60, rep(0, 6)), 5, 2))
                       # density_0 = matrix(c(rep(0, 3), 123, rep(0, 6)), 5, 2))




.sims <- do_base_sims(
    # alate_field_disp_p = 0.2,
    # pred_rate = 0.1,
    # K = 11e3,
    # s_y = populations$s_y * 0.8,
    # clear_surv = 0.01,
    wasp_density_0 = c(10, 0),
    wasp_delay = 0,
    # extinct_N = 0,
    # max_plant_age = 0,
    max_t = 10e3,
    clonal_lines = c(line_r2, line_s2))
.sims %>% plot_test()

mt <- max(.sims$aphids$time)

zz <- ((.sims$aphids %>% filter(time == mt))[["N"]] -
           (.sims$aphids %>% filter(time == mt-1))[["N"]])
((.sims$wasps %>% filter(time == mt))[["wasps"]] - (.sims$wasps %>% filter(time == mt-1))[["wasps"]])

.sims$aphids %>% filter(time == mt) %>% mutate(d = zz)


# # A tibble: 14 × 6
# time field         line        type              N       d
# <int> <fct>         <fct>       <chr>         <dbl>   <dbl>
#     1 10000 parasitism    resistant   alate        805.    1.10
# 2 10000 parasitism    resistant   apterous     782.    3.46
# 3 10000 parasitism    resistant   parasitized   20.8   0.0196
# 4 10000 parasitism    susceptible alate        136.   -1.40
# 5 10000 parasitism    susceptible apterous     157.   -4.11
# 6 10000 parasitism    susceptible parasitized   52.3  -0.155
# 7 10000 parasitism    NA          mummy          5.16  0
# 8 10000 no parasitism resistant   alate        213.   -1.06
# 9 10000 no parasitism resistant   apterous     239.   -3.39
# 10 10000 no parasitism resistant   parasitized    0     0
# 11 10000 no parasitism susceptible alate       1003.    1.34
# 12 10000 no parasitism susceptible apterous     826.    4.17
# 13 10000 no parasitism susceptible parasitized    0     0
# 14 10000 no parasitism NA          mummy          0     0



mod <- max(.sims$wasps$wasps) / max(.sims$aphids$N)

.sims$aphids %>%
    filter(type != "mummy") %>%
    filter(time > 99000) %>%
    filter(field == "parasitism") %>%
    filter(line == "resistant") %>%
    filter(type == "alate") %>%
    ggplot(aes(time, N / 1e3)) +
    # geom_hline(yintercept = 0, color = "gray70") +
    geom_line(aes(color = line, linetype = type)) +
    # geom_line(aes(color = line), size = 1) +
    scale_color_manual(values = c("chartreuse3", "firebrick"), guide = "none") +
    scale_linetype_manual(values = c("solid", "33", "11")) +
    scale_y_continuous(expression("Aphid abundance (" %*% 1000 * ")")) +
    scale_x_continuous("Days") +
    facet_wrap(~ field, ncol = 1) +
    theme(strip.text = element_text(size = 10),
          legend.position = "top") +  # ,
    # panel.grid.major.x = element_line(size = 0.25, color = "gray90"),
    # panel.grid.major.y = element_line(size = 0.25, color = "gray90")) +
    # coord_cartesian(ylim = c(0, 3)) +
    NULL




ad <- .sims$all_info

inds <- ad[[1]]$field == 1 & ad[[1]]$line == "resistant" & ad[[1]]$type != "parasitized"
ad[[1]][inds, "N"] <- ad[[1]][inds, "N"] * 1.1

z <- restart_experiment(.sims, new_starts = ad, max_t = 10, save_every = 1)


.sims2 <- do_base_sims(
    # alate_field_disp_p = 0.2,
    # pred_rate = 0.1,
    # K = 11e3,
    # s_y = populations$s_y * 0.8,
    wasp_density_0 = c(10, 0),
    wasp_delay = 0,
    # extinct_N = 0,
    max_t = 10e3+1,
    clonal_lines = c(line_r2, line_s2),
    perturb = tibble(when = 10e3+1, who = "resistant", how = 1.1, where = 1))

.sims$aphids %>% filter(time == max(time), field == para_lvls[1])
.sims2$aphids %>% filter(time == max(time), field == para_lvls[1])
z$aphids %>% filter(time == 1, field == min(field))


(.sims2$aphids %>% filter(time == max(time), field == para_lvls[1]))[["N"]] -
(z$aphids %>% filter(time == 1, field == min(field)))[["N"]]



.sims$wasps %>% filter(time == max(time), field == para_lvls[1])
.sims2$wasps %>% filter(time == max(time), field == para_lvls[1])
z$wasps %>% filter(time == 1, field == min(field))



# sad_dist <- sapply(1:length(dev_times$instar_days$lowT),
#        function(i) {
#            v <- dev_times$instar_days$lowT
#            m <- sad_leslie(line_r2$leslie[,,1])
#            j0 <- ifelse(i == 1, 1, sum(v[1:(i-1)]) + 1)
#            j1 <- ifelse(i == length(v), nrow(m), sum(v[1:i]))
#            sum(m[j0:j1])
#        })
#
#
# # # This will vary in function:
# # Y_m = 40  # adult wasps
# # x = 1000  # non-parasitized aphids
# # attack_surv = c(1, 0.6) # resistance to single and multiple oviposits
#
# wasp_surv <- function(Y_m, x, attack_surv, km = 1, hm = 1, am = 1) {
#     stopifnot(length(attack_surv) > 1)
#     a = wasp_attack$a * am
#     k = wasp_attack$k * km
#     h = wasp_attack$h * hm
#     rel_attack = wasp_attack$rel_attack
#     A_ <- (a * rel_attack * Y_m) / (h * x + 1)
#     AA <- (1 + A_ / k);
#     as_mults <- lapply(1:(length(attack_surv)-1), function(z) A_ * AA^(-k-z))
#     # A <-AA^(-k) + attack_surv[1] * A_ * AA^(-k-1) +
#     #     attack_surv[2] * { 1 - ( AA^(-k) + A_ * AA^(-k-1) ) }
#     A <- AA^(-k)
#     for (i in 1:(length(attack_surv)-1)) A <- A + attack_surv[i] * as_mults[[i]]
#     A <- A + tail(attack_surv, 1) * { 1 - ( AA^(-k) + Reduce(`+`, as_mults) ) }
#     return(A)
# }
#
# f <- function(z, .as, .km) {
#     map_dbl(z, ~ sum(sad_dist * wasp_surv(Y_m = .x, x = 1000, attack_surv = .as, km = .km)))
# }
#
# curve(f(x, wasp_attack$attack_surv, .km = 1), 0, 80, ylim = c(0, 1))
# curve(f(x, c(1, 0.2), .km = 1), 0, 80, add = TRUE, col = "red")
# curve(f(x, wasp_attack$attack_surv, .km = 10), 0, 80, add = TRUE, lty = 2)
# curve(f(x, c(1, 0.2), .km = 10), 0, 80, add = TRUE, col = "red", lty = 2)
# curve(f(x, wasp_attack$attack_surv, .km = 100), 0, 80, add = TRUE, lty = 3)
# curve(f(x, c(1, 0.2), .km = 100), 0, 80, add = TRUE, col = "red", lty = 3)









# ============================================================================*
# ============================================================================*

# Vary dispersal ----

# ============================================================================*
# ============================================================================*


#'
#' Do a call to `sim_clonewars` that can vary by the proportion of alates
#' that move between fields every day (`alate_p`).
#' How does dispersal affect maintenance of diversity in the resistance trait?
#'
do_disp_sims <- function(.alate_disp_prop, .max_t = 5000) {

    .sims <- do_base_sims(.clonal_lines = c(line_r, line_s),
                          .alate_disp_prop = .alate_disp_prop,
                          .max_t = .max_t, .perturb = NULL)

    for (n in names(.sims)) .sims[[n]] <- .sims[[n]] %>%
        mutate(alate_p = .alate_disp_prop)

    return(.sims)
}



# Takes ~8 sec
disp_sims <- mclapply(c(0.01, 0.05, 0.1), do_disp_sims)


disp_aphids <- map_dfr(disp_sims, ~ group_aphids(.x[["aphids"]], "living")) %>%
    mutate(alate_p = factor(alate_p, labels = sprintf("%.0f%% dispersal",
                                                    sort(unique(alate_p)) * 100)))
disp_mummies <- map_dfr(disp_sims, ~ group_aphids(.x[["aphids"]], "mummies")) %>%
    mutate(alate_p = factor(alate_p, labels = sprintf("%.0f%% dispersal",
                                                    sort(unique(alate_p)) * 100)))
disp_wasps <- map_dfr(disp_sims, ~ .x[["wasps"]]) %>%
    mutate(alate_p = factor(alate_p, labels = sprintf("%.0f%% dispersal",
                                                    sort(unique(alate_p)) * 100)))


disp_mod <- max(disp_wasps$wasps) / max(disp_aphids$N)

disp_p <- disp_aphids %>%
    ggplot(aes(time, N / 1e3)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_area(data = disp_wasps %>%
                  filter(field == para_lvls[1]) %>%
                  mutate(N = wasps / disp_mod),
              fill = "gray80", color = NA) +
    geom_line(aes(color = line, linetype = field, size = field)) +
    geom_text(data = tibble(line = factor(c("resistant", "susceptible")),
                            field = factor(para_lvls[1], levels = para_lvls),
                            alate_p = disp_aphids$alate_p %>% unique() %>% sort() %>% .[1],
                            time = c(500, 500), N = c(300, max(disp_aphids$N) - 300)),
              aes(label = line, color = line),
              size = 9 / 2.8, hjust = 0, vjust = c(0, 1)) +
    geom_text(data = tibble(field = factor(para_lvls[1], levels = para_lvls),
                            alate_p = disp_aphids$alate_p %>% unique() %>% sort() %>% .[3],
                            time = 500, N = 1200),
              aes(label = "wasps"), size = 9 / 2.8, hjust = 0, vjust = 0,
              color = "gray50") +
    scale_color_manual(values = c("chartreuse3", "firebrick"), guide = "none") +
    scale_linetype_manual("Patch:", values = c("solid", "22")) +
    scale_size_manual("Patch:", values = c(0.5, 0.75)) +
    scale_y_continuous(expression("Aphid abundance (" %*% 1000 * ")"),
                       sec.axis = sec_axis(~ . * disp_mod * 1000,
                                           "Wasp abundance")) +
    scale_x_continuous("Days") +
    facet_wrap(~ alate_p, ncol = 1) +
    theme(strip.text = element_text(size = 10),
          legend.position = "top") +
    coord_cartesian(xlim = c(0, 2000))


# disp_p
disp_p + coord_cartesian(xlim = c(0, 250))



# save_plot("_results/plots/disp_sims.pdf", disp_p, 4, 5)






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

    .n <- 2 * line_s$density_0[line_s$density_0 > 0]

    .line_s <- line_s
    .line_r <- line_r

    if (!identical(.line_s$density_0, .line_r$density_0)) {
        stop("ERROR: starting densities in `line_s` and `line_r` should be ",
             "the same. (Let the later functions change this.)")
    }
    if (sum(.line_s$density_0 > 0) != 1) {
        stop("ERROR: starting densities in `line_s` and `line_r` should only ",
             "consist of a single stage.")
    }
    .line_s$density_0[.line_s$density_0 > 0] <- .n * (1 - .p_res)
    .line_r$density_0[.line_r$density_0 > 0] <- .n * .p_res

    .sims <- do_base_sims(.clonal_lines = c(.line_r, .line_s),
                          .alate_disp_prop = 0.05, # middle val from `disp_sims`
                          .max_t = .max_t, .perturb = NULL)

    for (n in names(.sims)) .sims[[n]] <- .sims[[n]] %>%
        mutate(p_res = .p_res)

    return(.sims)
}




p_res_lvls <- c(0.01, 0.1, 0.75, 0.8)

# Takes ~9 sec
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

stable_start_mod <- max(stable_start_wasps$wasps) / max(stable_start_aphids$N)


stable_start_p <- stable_start_aphids %>%
    ggplot(aes(time, N / 1e3)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_area(data = stable_start_wasps %>%
                  filter(field == para_lvls[1]) %>%
                  mutate(N = wasps / stable_start_mod),
              fill = "gray80", color = NA) +
    geom_line(aes(color = line, linetype = field, size = field)) +
    scale_color_manual(values = c("chartreuse3", "firebrick"), guide = "none") +
    scale_linetype_manual("Patch:", values = c("solid", "22")) +
    scale_size_manual("Patch:", values = c(0.5, 0.75)) +
    scale_y_continuous(expression("Aphid abundance (" %*% 1000 * ")"),
                       sec.axis = sec_axis(~ . * stable_start_mod * 1000,
                                           "Wasp abundance")) +
    scale_x_continuous("Days") +
    facet_wrap(~ p_res, ncol = 1) +
    theme(strip.text = element_text(size = 10),
          legend.position = "top") +
    coord_cartesian(xlim = c(0, 2000))


# stable_start_p

# save_plot("_results/plots/stable_start_sims.pdf", stable_start_p, 4, 5)











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

    .sims <- do_base_sims(.clonal_lines = c(line_r, line_s),
                          .alate_disp_prop = 0.05, # middle val from `disp_sims`
                          .max_t = .max_t,
                          .perturb = .perturb_df)

    for (n in names(.sims)) .sims[[n]] <- .sims[[n]] %>%
        mutate(who = paste(unique(.perturb_df$who), collapse = " & "))

    return(.sims)
}

# Takes ~30 sec
stable_perturb_sims <- list(tibble(when = 10e3, where = 1:2,
                                   who = "resistant", how = 0.5),
                            tibble(when = 10e3, where = 1:2,
                                   who = "susceptible", how = 0.5),
                            tibble(when = 10e3, where = rep(1:2, 2),
                                   who = rep(c("wasps", "mummies"), each = 2),
                                   how = 0.5)) %>%
    mclapply(do_stable_perturb_sims)

pert_t_range <- c(9.8, 12) * 1e3

stable_perturb_aphids <- map_dfr(stable_perturb_sims,
                                 ~ group_aphids(.x[["aphids"]], "living")) %>%
    mutate(who = factor(who, levels = c("resistant", "susceptible",
                                        "wasps & mummies"),
                        labels = c("resistant aphids", "susceptible aphids",
                                   "wasps & parasitized aphids"))) %>%
    filter(time > pert_t_range[1], time < pert_t_range[2])
stable_perturb_mummies <- map_dfr(stable_perturb_sims,
                                  ~ group_aphids(.x[["aphids"]], "mummies")) %>%
    mutate(who = factor(who, levels = c("resistant", "susceptible",
                                        "wasps & mummies"),
                        labels = c("resistant aphids", "susceptible aphids",
                                   "wasps & parasitized aphids"))) %>%
    filter(time > pert_t_range[1], time < pert_t_range[2])
stable_perturb_wasps <- map_dfr(stable_perturb_sims, ~ .x[["wasps"]]) %>%
    mutate(who = factor(who, levels = c("resistant", "susceptible",
                                        "wasps & mummies"),
                        labels = c("resistant aphids", "susceptible aphids",
                                   "wasps & parasitized aphids"))) %>%
    filter(time > pert_t_range[1], time < pert_t_range[2])



stable_perturb_mod <- max(stable_perturb_wasps$wasps) / max(stable_perturb_aphids$N)



stable_perturb_p <- stable_perturb_aphids %>%
    ggplot(aes(time, N / 1e3)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_area(data = stable_perturb_wasps %>%
                  filter(field == para_lvls[1]) %>%
                  mutate(N = wasps / stable_perturb_mod),
              fill = "gray80", color = NA) +
    geom_line(aes(color = line, linetype = field, size = field)) +
    scale_color_manual(values = c("chartreuse3", "firebrick"), guide = "none") +
    scale_linetype_manual("Patch:", values = c("solid", "22")) +
    scale_size_manual("Patch:", values = c(0.5, 0.75)) +
    scale_y_continuous(expression("Aphid abundance (" %*% 1000 * ")"),
                       sec.axis = sec_axis(~ . * stable_perturb_mod * 1000,
                                           "Wasp abundance")) +
    scale_x_continuous("Days") +
    facet_wrap(~ who, ncol = 1) +
    theme(strip.text = element_text(size = 10),
          legend.position = "top")


# stable_perturb_p


# save_plot("_results/plots/stable_perturb_sims.pdf", stable_perturb_p, 4, 5)






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

    .sims <- do_base_sims(.clonal_lines = c(line_r, line_s),
                          .alate_disp_prop = 0.05, # middle val from `disp_sims`
                          .max_t = .max_t,
                          .perturb = .perturb_df)


    for (n in names(.sims)) .sims[[n]] <- .sims[[n]] %>%
        filter(time > (min(.perturb_df$when) - 100)) %>%
        mutate(perturb_who = paste(.who, collapse = " & "),
               perturb_where = .where)

    return(.sims)
}

# # Takes ~3 min
# comm_mat_perturb_sims <- crossing(.who = c("resistant", "susceptible",
#                                            "parasitized"), .where = 1:2) %>%
#     # no need to perturb the non-existent parasitized aphids in the
#     # no parasitism field
#     filter(!(.who == "parasitized" &
#                  .where == which(grepl("^no", para_lvls)))) %>%
#     as.list() %>%
#     c(list(FUN = do_comm_mat_perturb_sims, SIMPLIFY = FALSE)) %>%
#     do.call(what = mcmapply)
# saveRDS(comm_mat_perturb_sims, "_results/comm_mat_perturb_sims.rds")

comm_mat_perturb_sims <- readRDS("_results/comm_mat_perturb_sims.rds")



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
    mutate(perturb_where = factor(perturb_where, levels = 1:2,
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




# ============================================================================*
# ============================================================================*

# Low parasitism? ----

#'
#' Trying to simulate with low parasitism instead of no parasitism.
#' Can't really get this working!
#' The wasps in the low-parasitism field always go extinct.
#'

# ============================================================================*
# ============================================================================*


#'
#' Trying to get low parasitism field not to crash.
#'
#'
#' * `kymm` affects the low parasitism field's value for the `K_y_mult` arg
#'   to `sim_clonewars` (it is multiplied by the default value) and
#'   affects the density dependence of adult wasps;
#'   higher values create *weaker* density dependence.
#' * `wdp` changes the `wasp_dispersal_p` arg to `sim_clonewars` and
#'   is the proportion of adult wasps that are contributed to the wasp
#'   dispersal pool every day;
#'   greater values mean greater dispersal
#'
do_low_sims <- function(kymm, wdp) {

    .sims <- do_base_sims(.clonal_lines = c(line_r, line_s),
                          .alate_disp_prop = 0.05, .max_t = 1000,
                          .perturb = NULL, .extinct_N = 0,
                          .K_y_mult = (1 / 1.57) * c(1, kymm),
                          .wasp_dispersal_p = wdp)

    for (n in names(.sims)) .sims[[n]] <- .sims[[n]] %>%
        mutate(field = fct_recode(field, low = "no parasitism",
                                 high = "parasitism"),
               k = kymm, w = wdp)

    return(.sims)
}



# Takes ~4 sec
low_sims <- crossing(kymm = c(0.2, 0.4, 0.6),
                     wdp = c(0.01, 0.02, 0.04)) %>%
    as.list() %>%
    c(list(FUN = do_low_sims, SIMPLIFY = FALSE)) %>%
    do.call(what = mcmapply)



low_aphids <- map_dfr(low_sims, ~ group_aphids(.x[["aphids"]], "living")) %>%
    mutate(k = factor(k), w = factor(w))
low_mummies <- map_dfr(low_sims, ~ group_aphids(.x[["aphids"]], "mummies")) %>%
    mutate(k = factor(k), w = factor(w))
low_wasps <- map_dfr(low_sims, ~ .x[["wasps"]]) %>%
    mutate(k = factor(k), w = factor(w))

low_mod <- max(low_wasps$wasps) / max(low_aphids$N)


low_p_list <- low_aphids %>%
    split(interaction(.$k, .$w, drop = TRUE)) %>%
    map(function(.d) {
        .k <- .d$k[[1]] %>% paste() %>% as.numeric()
        .w <- .d$w[[1]] %>% paste() %>% as.numeric()
        .d %>%
            ggplot(aes(time, N / 1e3)) +
            ggtitle(sprintf("k = %.1f | w = %.2f", .k, .w)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_area(data = low_wasps %>%
                          filter(k == .k, w == .w) %>%
                          mutate(N = wasps / low_mod),
                      fill = "gray80", color = NA) +
            geom_line(aes(color = line)) +
            scale_color_manual(values = c("chartreuse3", "firebrick"),
                               guide = "none") +
            scale_linetype_manual("Patch:", values = c("solid", "22")) +
            scale_size_manual("Patch:", values = c(0.5, 0.75)) +
            scale_y_continuous(expression("Aphid abundance (" %*% 1000 * ")"),
                               sec.axis = sec_axis(~ . * low_mod * 1000,
                                                   "Wasp abundance")) +
            scale_x_continuous("Days") +
            facet_wrap(~ field, ncol = 1) +
            theme(strip.text = element_text(size = 8),
                  axis.title = element_text(size = 8),
                  axis.text = element_text(size = 6),
                  plot.title = element_text(size = 9, hjust = 0.5))
    })



wrap_plots(low_p_list, nrow = 3, ncol = 3)






