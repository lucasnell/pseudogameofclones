
library(tidyverse)
library(clonewars)
library(parallel)
library(patchwork)

options(mc.cores = max(parallel::detectCores()-2L, 1L))

save_plot <- function(fn, p, w, h, ...) {
    fn_dir <- dirname(fn)
    if (!dir.exists(fn_dir)) stop("ERROR: `", fn_dir, "` doesn't exist")
    cairo_pdf(filename = fn, width = w, height = h, ...)
    plot(p)
    dev.off()
    invisible(NULL)
}


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
#' assays in our cages indicate density has little effect.
#'
#'
#' I'm simulating 1,000 days of interacting aphid and parasitoid wasp
#' populations.
#' There are two patches, one with wasps and another without.
#' Movement between patches must happen via alates, and
#' rates of alate production are density-independent.
#' Every day, 10% of all alates in each cage disperse to the other cage.
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
                      resistant = c(0.95, 0.8),
                      surv_juv_apterous = "low",
                      surv_adult_apterous = "low",
                      repro_apterous = "low")

para_lvls <- c("parasitism", "no parasitism")


#'
#' Base function to do all the simulations.
#' Each section below builds on this one.
#'
do_base_sims <- function(.clonal_lines, .alate_disp_prop, .max_t, .perturb) {

    .sims <- sim_clonewars(n_reps = 1, n_patches = 1, n_cages = 2,
                           max_N = 0,
                           max_plant_age = 1e9,
                           mean_K = formals(sim_clonewars)$mean_K * 8,
                           no_error = TRUE,
                           wasp_delay = 0,
                           wasp_density_0 = c(3, 0),
                           plant_check_gaps = 1,  # makes alates disperse daily
                           alate_b0 = logit(0.3),
                           clonal_lines = .clonal_lines,
                           alate_disp_prop = .alate_disp_prop,
                           max_t = .max_t,
                           perturb = .perturb)

    .sims[["mummies"]] <- .sims %>%
        .[["aphids"]] %>%
        filter(type == "mummy") %>%
        group_by(cage, time) %>%
        summarize(N = sum(N), .groups = "drop") %>%
        mutate(cage = factor(cage, levels = 0:max(cage),
                             labels = para_lvls))
    .sims[["aphids"]] <- .sims %>%
        .[["aphids"]] %>%
        filter(type != "mummy") %>%
        group_by(cage, time, line) %>%
        summarize(N = sum(N), .groups = "drop") %>%
        mutate(cage = factor(cage, levels = 0:max(cage),
                             labels = para_lvls),
               line = factor(line, levels = c("resistant", "susceptible")))
    .sims[["wasps"]] <- .sims %>%
        .[["wasps"]] %>%
        select(-rep) %>%
        mutate(cage = factor(cage, levels = 0:max(cage),
                             labels = para_lvls))

    return(.sims)
}


# ============================================================================*
# ============================================================================*

# Vary dispersal ----

# ============================================================================*
# ============================================================================*


#'
#' Do a call to `sim_clonewars` that can vary by the proportion of alates
#' that move between cages every day (`alate_p`).
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


disp_aphids <- map_dfr(disp_sims, ~ .x[["aphids"]]) %>%
    mutate(alate_p = factor(alate_p, labels = sprintf("%.0f%% dispersal",
                                                    sort(unique(alate_p)) * 100)))
disp_mummies <- map_dfr(disp_sims, ~ .x[["mummies"]]) %>%
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
                  filter(cage == para_lvls[1]) %>%
                  mutate(N = wasps / disp_mod),
              fill = "gray80", color = NA) +
    geom_line(aes(color = line, linetype = cage, size = cage)) +
    geom_text(data = tibble(line = factor(c("resistant", "susceptible")),
                            cage = factor(para_lvls[1], levels = para_lvls),
                            alate_p = disp_aphids$alate_p %>% unique() %>% sort() %>% .[1],
                            time = c(500, 500), N = c(300, max(disp_aphids$N) - 300)),
              aes(label = line, color = line),
              size = 9 / 2.8, hjust = 0, vjust = c(0, 1)) +
    geom_text(data = tibble(cage = factor(para_lvls[1], levels = para_lvls),
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


stable_start_aphids <- map_dfr(stable_start_sims, ~ .x[["aphids"]]) %>%
    mutate(p_res = factor(p_res, labels = sprintf("%.0f%% starting resistance",
                                                  sort(unique(p_res)) * 100)))
stable_start_mummies <- map_dfr(stable_start_sims, ~ .x[["mummies"]]) %>%
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
                  filter(cage == para_lvls[1]) %>%
                  mutate(N = wasps / stable_start_mod),
              fill = "gray80", color = NA) +
    geom_line(aes(color = line, linetype = cage, size = cage)) +
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

stable_perturb_aphids <- map_dfr(stable_perturb_sims, ~ .x[["aphids"]]) %>%
    mutate(who = factor(who, levels = c("resistant", "susceptible",
                                        "wasps & mummies"),
                        labels = c("resistant aphids", "susceptible aphids",
                                   "wasps & parasitized aphids"))) %>%
    filter(time > pert_t_range[1], time < pert_t_range[2])
stable_perturb_mummies <- map_dfr(stable_perturb_sims, ~ .x[["mummies"]]) %>%
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
                  filter(cage == para_lvls[1]) %>%
                  mutate(N = wasps / stable_perturb_mod),
              fill = "gray80", color = NA) +
    geom_line(aes(color = line, linetype = cage, size = cage)) +
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



# save_plot("_results/plots/stable_perturb_sims.pdf", stable_perturb_p, 4, 5)











#'
#' ### Stabilizing mechanism
#'
#' The vector of abundances for the aphid stages at time $t+1$ ($\mathbf{X}(t+1)$) is
#'
#' $$
#' \mathbf{X}(t+1) = S(z(t)) \cdot \mathbf{A}(x(t), Y_m(t)) \cdot (\mathbf{L} \, \mathbf{X}(t))
#' $$
#'
#' where $S$ gives survival of aphids as it relates to the the total number of aphids $z$,
#' $\mathbf{A}$ is a vector of the probabilities of the aphid stages surviving parasitism,
#' $x$ is the total number of unparasitized aphids,
#' $Y_m$ is the total number of adult wasps,
#' and
#' $\mathbf{L}$ is the aphid line's Leslie matrix.
#'
#'
#'
#' The probability of not being attacked by wasps
#' for a given aphid stage ($A_i$) is
#'
#' $$
#' A_i = \left( 1 + \frac{a p_i Y_m}{k (h x + 1)} \right)^{-k}
#' \textrm{.}
#' $$
#'
#' Here, $a$ changes the overall attack rate,
#' $p_i$ is the relative rate on stage $i$ aphids,
#' $h$ is the handling time, and
#' $k$ is the aggregation parameter of the negative binomial distribution.
#'
#'
#' If we focus on one of the sub-panels from the stability simulations and
#' plot the attack probabilities through time, we can see that it
#'


a <- wasp_attack$a
k <- wasp_attack$k
h <- wasp_attack$h
p_i <- wasp_attack$rel_attack[3] / 2
K <- formals(sim_clonewars)$mean_K * 8
A <- function(x, Y_m) (1 + (a * p_i * Y_m) / (k * (h * x + 1)))^(-k)
S <- function(z) 1 / (1 + z / K)

stab_mech_df <- stable_perturb_aphids %>%
    filter(who == "susceptible aphids", cage == "parasitism") %>%
    select(time, line, N)
# stab_mech_wasps_df <-
stab_mech_df %>%
    pivot_wider(names_from = line, values_from = N) %>%
    mutate(total = resistant + susceptible) %>%
    left_join(stable_perturb_wasps %>%
                  filter(who == "susceptible aphids", cage == "parasitism") %>%
                  select(time, wasps),
              by = c("time")) %>%
    # Benefits conferred to susceptible line by resistant via reduced parasitism:
    mutate(res_benefits = A(total, wasps) / A(susceptible, wasps)) %>%
    # Costs conferred to susceptible line by resistant via competition:
    mutate(res_costs = S(total) / S(susceptible)) %>%
    mutate(w_res = S(total) * A(total, wasps),
           wo_res = S(susceptible) * A(susceptible, wasps),
           res_effect = w_res / wo_res) %>%
    .[["res_effect"]] %>% range()




# surv_p <-

stab_mech_df %>%
    ggplot(aes(time, N)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_line(data = stab_mech_wasps_df %>%
                  mutate(res_benefits = res_benefits * max(stab_mech_df$N)),
              aes(time, res_benefits),
              color = "dodgerblue", size = 0.5) +
    geom_line(aes(color = line), size = 0.5) +
    scale_color_manual(NULL, values = c("chartreuse3", "firebrick")) +
    scale_linetype_manual(values = c(1, 2)) +
    scale_y_continuous("Aphid abundance", #breaks = 0:2 * 4000,
                       sec.axis = sec_axis(~ . / max(stab_mech_df$N),
                                           "Wasp survival")) +
    theme(legend.position = "none")

# surv_p + coord_cartesian(xlim = c(0, 100))
surv_p


ggsave("~/Desktop/wasp_survival.pdf", surv_p, width = 6, height = 4)


curve(S(x) * A(x, 10), 1, 10e3, ylim = c(0, 0.6), ylab = "survival", xlab = "# aphids")
curve(S(x) * A(x, 100), 1, 10e3, add = TRUE, col = "red")

.n <- 10000
curve((S(.n) * A(.n, 10)) - (S(.n*x) * A(.n*x, 10)), 0, 1,
      ylim = c(-0.2, 0.2), ylab = "survival", xlab = "proportion resistant")
curve((S(.n) * A(.n, 100)) - (S(.n*x) * A(.n*x, 100)), 0, 1,
      add = TRUE, col = "red")
abline(h = 0, lty = 2, col = "gray70")





stable_aphids %>%
    pivot_wider(names_from = line, values_from = N) %>%
    left_join(stable_wasps, by = c("s_y", "p_res", "time")) %>%
    # Benefits conferred to susceptible line by resistant:
    mutate(res_benefits = A(resistant + susceptible, wasps) -
               A(susceptible, wasps)) %>%
    mutate(res_benefits = res_benefits / max(res_benefits) *
                      max(max(resistant), max(susceptible))) %>%
    ggplot(aes(time)) +
    geom_line(aes(y = res_benefits), color = "gray70") +
    geom_line(aes(y = resistant), color = "chartreuse3") +
    geom_line(aes(y = susceptible), color = "firebrick") +
    # coord_cartesian(xlim = c(0, 100)) +
    facet_grid(p_res ~ s_y, labeller = label_parsed)


# ppp <-
stable_aphids %>%
    pivot_wider(names_from = line, values_from = N) %>%
    left_join(stable_wasps, by = c("s_y", "p_res", "time")) %>%
    filter(wasps > 0) %>%
    # Benefits conferred to susceptible line by resistant:
    mutate(total = resistant + susceptible,
           res_benefits = A(total, median(wasps)) - A(susceptible, median(wasps)),
           p_res = resistant / total) %>%
    filter(total > 3000) %>%
    ggplot(aes(p_res, res_benefits)) +
    geom_point(shape = 1, alpha = 0.4) +
    ylab("Reduction in wasp mortality\nfor susceptible line") +
    xlab("Proportion of aphids that are resistant") +
    NULL

ppp

# ggsave("~/Desktop/benefits_from_resistant.pdf", ppp, width = 5, height = 3)








# # Playing around with individual combinations:
# z <- do_stable_sims(.p_res = 0, .s_y = populations$s_y / 1000,  max_t = 1000)
# z$mod <- max(z$wasps$wasps) / max(z$aphids$N)
#
#
# z$aphids %>%
#     mutate(line = factor(line, levels = c("resistant", "susceptible"))) %>%
#     ggplot(aes(time, N)) +
#     geom_line(data = z$wasps %>%
#                   mutate(N = wasps / z$mod),
#               color = "gray70", size = 0.25) +
#     geom_line(aes(color = line), size = 0.5) +
#     scale_color_manual(NULL, values = c("chartreuse3", "firebrick")) +
#     scale_linetype_manual(values = c(1, 2)) +
#     scale_y_continuous("Aphid abundance", breaks = 0:2 * 4000,
#                        sec.axis = sec_axis(~ . * z$mod,
#                                            "Adult wasp abundance")) +
#     theme(legend.position = "none")

