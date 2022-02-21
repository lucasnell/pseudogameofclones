
suppressPackageStartupMessages({
    library(tidyverse)
    library(clonewars)
    library(parallel)
    library(patchwork)
})
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


save_plot("_results/plots/disp_sims.pdf", disp_p, 4, 5)






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
#' ### Stability to perturbations
#'
#' Below shows that the coexistence results can be relatively stable
#' to various perturbations, especially as adult wasp survival decreases.
#' We perturbed the susceptible line, resistant line, and adult wasps.
#' See captions for details.
#'




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
        mutate(who = paste(.perturb_df$who, collapse = " & "))

    return(.sims)
}

# Takes ~30 sec
stable_perturb_sims <- list(tibble(when = 10e3, who = "resistant", how = 0.5),
                            tibble(when = 10e3, who = "susceptible", how = 0.5),
                            tibble(when = 10e3, who = c("wasps", "mummies"),
                                   how = 0.5)) %>%
    mclapply(do_stable_perturb_sims)

pert_t_range <- c(9.8, 12) * 1e3

stable_perturb_aphids <- map_dfr(stable_perturb_sims, ~ .x[["aphids"]]) %>%
    mutate(who = factor(who, levels = c("resistant", "susceptible", "wasps"))) %>%
    filter(time > pert_t_range[1], time < pert_t_range[2])
stable_perturb_mummies <- map_dfr(stable_perturb_sims, ~ .x[["mummies"]]) %>%
    mutate(who = factor(who, levels = c("resistant", "susceptible", "wasps"))) %>%
    filter(time > pert_t_range[1], time < pert_t_range[2])
stable_perturb_wasps <- map_dfr(stable_perturb_sims, ~ .x[["wasps"]]) %>%
    mutate(who = factor(who, levels = c("resistant", "susceptible", "wasps"))) %>%
    filter(time > pert_t_range[1], time < pert_t_range[2])



stable_perturb_mod <- max(stable_perturb_wasps$wasps) / max(stable_perturb_aphids$N)



# stable_perturb_p <-
stable_perturb_aphids %>%
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










# fig4_cap <- paste("Aphids and wasps over time, with varying adult",
#                   "wasp survival and perturbation types.",
#                   "Sub-panel columns separate the adult wasp survival rates used,",
#                   "and labels above indicate how they relate",
#                   "to the default value ($\\hat{s}_y$).",
#                   "Sub-panel rows separate the following perturbation types:",
#                   "'none' indicates no perturbation,",
#                   "'susceptible' indicates the susceptible line was reduced by 50\\%,",
#                   "'resistant' indicates the resistant line was reduced by 50\\%,",
#                   "and",
#                   "'wasps' indicates the wasps were reduced by 90\\%.",
#                   "The black arrows indicate when the perturbations occurred.",
#                   "Gray lines are adult wasps, green are resistant aphids,",
#                   "and red are susceptible aphids.")
# fig5_cap <- paste("Aphids and wasps over time, with varying perturbation types.",
#                   "Sub-panels separate the following perturbation types:",
#                   "'susceptible' indicates the susceptible line was reduced by 50\\%,",
#                   "and",
#                   "'resistant' indicates the resistant line was reduced by 50\\%.",
#                   "The black arrows indicate when the perturbations occurred.",
#                   "Gray lines are adult wasps, green are resistant aphids,",
#                   "and red are susceptible aphids.",
#                   "All used an adult survival rate of $\\hat{s}_y / 8$.")


# Each takes ~10 sec
# p_res_ <- 0.5
# pert_sims <- list()
# pert_sims[["none"]] <- crossing(.p_res = p_res_) %>%
#     pmap(do_stable_sims, max_t = 1000)
# pert_sims[["resistant"]] <- crossing(.p_res = p_res_) %>%
#     pmap(do_stable_sims, max_t = 1000,
#          perturb = tibble(when = 500, who = "resistant", how = 0.7))
# pert_sims[["susceptible"]] <- crossing(.p_res = p_res_) %>%
#     pmap(do_stable_sims, max_t = 1000,
#          perturb = tibble(when = 500, who = "susceptible", how = 0.7))
# pert_sims[["wasps"]] <- crossing(.p_res = p_res_) %>%
#     pmap(do_stable_sims, max_t = 1000,
#          perturb = tibble(when = 500, who = "wasps", how = 0.1))
#
# saveRDS(pert_sims, "under_constr/pert_sims.rds")

pert_sims <- readRDS("under_constr/pert_sims.rds")


pert_aphids <- map_dfr(names(pert_sims),
                       function(.z) {
                           map_dfr(pert_sims[[.z]], ~ .x[["aphids"]]) %>%
                               select(-p_res) %>%
                               mutate(pert = .z) %>%
                               stable_fct()
                       })
# pert_mummies <- map_dfr(names(pert_sims),
#                         function(.z) {
#                             map_dfr(pert_sims[[.z]], ~ .x[["mummies"]]) %>%
#                                 select(-p_res) %>%
#                                 mutate(pert = .z) %>%
#                                 stable_fct()
#                         })
pert_wasps <- map_dfr(names(pert_sims),
                      function(.z) {
                          map_dfr(pert_sims[[.z]], ~ .x[["wasps"]]) %>%
                              select(-p_res) %>%
                              mutate(pert = .z) %>%
                              stable_fct()
                      })


pert_mod <- max(pert_wasps$wasps) / max(pert_aphids$N)

pert_p <- pert_aphids %>%
    ggplot(aes(time, N)) +
    geom_hline(yintercept = 0) +
    geom_segment(data = pert_wasps %>%
                   filter(pert != "none") %>%
                   distinct(cage, pert),
               aes(x = 500, xend = 500,
                   y = max(pert_aphids$N) + 1000, yend = max(pert_aphids$N)),
               arrow = arrow(length = unit(3, "pt"))) +
    geom_area(data = pert_wasps %>%
                  mutate(N = wasps / pert_mod),
              fill = "gray80", color = NA) +
    geom_line(aes(color = line), size = 0.25) +
    facet_grid(pert ~ cage, labeller = label_parsed) +
    scale_color_manual(NULL, values = c("chartreuse3", "firebrick")) +
    scale_linetype_manual(values = c(1, 2)) +
    scale_y_continuous("Aphid abundance", #breaks = 0:2 * 4000,
                       sec.axis = sec_axis(~ . * pert_mod,
                                           "Adult wasp abundance")) +
    theme(legend.position = "none")

# pert_p + coord_cartesian(xlim = c(0, 100))

pert_p

# ggsave("~/Desktop/pert_plots.pdf", pert_p, width = 6.5, height = 5)






#'
#' In the one case where a perturbation caused exclusion of the wasps and
#' resistant line, this appears to be a case of timing of the perturbation.
#' Below, we used only $\hat{s}_y / 5$ for adult wasp survival, then simulated
#' the perturbations of the susceptible and resistant lines again, but 100 days
#' later than in the previous figure.
#' We chose 100 days because that puts the perturbation in the opposite
#' phase of the fluctuations, when the resistant line is increasing and the
#' susceptible is decreasing.
#'




.p_time <- 400

# # Takes ~10 sec
# pert_sims2 <- map(c("resistant", "susceptible"),
#                   function(.x) {
#                       do_stable_sims(.p_res = 0.4, .s_y = populations$s_y / 8,
#                                      max_t = 1000,
#                                      perturb = tibble(when = .p_time,
#                                                       who = .x, how = 0.5))
#                   })
# names(pert_sims2) <- c("resistant", "susceptible")
# saveRDS(pert_sims2, "under_constr/pert_sims2.rds")

pert_sims2 <- readRDS("under_constr/pert_sims2.rds")



pert_aphids2 <- map_dfr(names(pert_sims2),
                       function(.z) {
                           pert_sims2[[.z]][["aphids"]] %>%
                               select(-rep, -p_res, -s_y) %>%
                               mutate(pert = .z) %>%
                               stable_fct()
                       })
pert_wasps2 <- map_dfr(names(pert_sims2),
                      function(.z) {
                          pert_sims2[[.z]][["wasps"]] %>%
                              select(-rep, -p_res, -s_y) %>%
                              mutate(pert = .z) %>%
                              stable_fct()
                      })


pert_mod2 <- max(pert_wasps2$wasps) / max(pert_aphids2$N)

pert_p2 <- pert_aphids2 %>%
    ggplot(aes(time, N)) +
    geom_segment(x = .p_time, xend = .p_time,
                 y = max(pert_aphids2$N), yend = max(pert_aphids2$N) - 2000,
               arrow = arrow(length = unit(3, "pt"))) +
    geom_line(aes(color = line), size = 0.25) +
    geom_line(data = pert_wasps2 %>%
                  mutate(N = wasps / pert_mod2),
              color = "gray70", size = 0.25) +
    facet_wrap(~ pert, ncol = 1) +
    scale_color_manual(NULL, values = c("chartreuse3", "firebrick")) +
    scale_y_continuous("Aphid abundance",
                       sec.axis = sec_axis(~ . * pert_mod2,
                                           "Adult wasp abundance")) +
    theme(legend.position = "none")

pert_p2


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
K <- formals(sim_clonewars)$mean_K * 4
A <- function(x, Y_m) (1 + (a * p_i * Y_m) / (k * (h * x + 1)))^(-k)
S <- function(z) 1 / (1 + z / K)

surv_p <- stable_aphids %>%
    filter(s_y == "hat(s)[y] / 6", p_res == "p[res] == 0.4") %>%
    ggplot(aes(time, N)) +
    geom_hline(yintercept = c(0, max(stable_aphids$N)), color = "gray70") +
    geom_line(data = stable_aphids %>%
                  filter(s_y == "hat(s)[y] / 6", p_res == "p[res] == 0.4") %>%
                  group_by(s_y, p_res, time) %>%
                  summarize(N = sum(N), .groups = "drop") %>%
                  left_join(stable_wasps, by = c("s_y", "p_res", "time")) %>%
                  mutate(wasp_surv = A(N, wasps)) %>%
                  mutate(wasp_surv = wasp_surv * max(stable_aphids$N)),
              aes(time, wasp_surv),
              color = "dodgerblue", size = 0.25) +
    geom_line(aes(color = line), size = 0.25) +
    scale_color_manual(NULL, values = c("chartreuse3", "firebrick")) +
    scale_linetype_manual(values = c(1, 2)) +
    scale_y_continuous("Aphid abundance", #breaks = 0:2 * 4000,
                       sec.axis = sec_axis(~ . / max(stable_aphids$N),
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

