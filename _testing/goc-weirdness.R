

# library(gameofclones)
library(pseudogameofclones)
library(tidyverse)


col_pal <- c(resistant = "#218F8DFF", susceptible = "#DDE318FF")

# # Define clonal lines. Don't re-define these!
# # Susceptible line: no resistance, high population growth rate
# line_s <- clonal_line("susceptible",
#                       density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
#                       surv_juv_apterous = "high",
#                       surv_adult_apterous = "high",
#                       repro_apterous = "high")
# # Resistant line: high resistance, low population growth rate
# line_r <- clonal_line("resistant",
#                       density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
#                       resistant = TRUE,
#                       surv_paras = 0.57,
#                       surv_juv_apterous = "low",
#                       surv_adult_apterous = "low",
#                       repro_apterous = "low")

n_fields <- 28L
wasp_density_0 <- rep(0.2, n_fields)

delta.list <- exp(0.5*(-10:10))

# wasp dispersal
wasp.var.list <- 0.1*(0:20)
pick.wasp.var <- 11L

set.seed(0)
#sd.wasp.field.attract <- 0.5844
sd.wasp.field.attract <- 0.5844 * wasp.var.list[pick.wasp.var]
wasp_field_attract <- sort(exp(-rnorm(0, sd = sd.wasp.field.attract, n = n_fields)))
wasp_field_attract <- wasp_field_attract/sum(wasp_field_attract)
wasp_field_attract

wasp_disp_m0 <- 0.3
#wasp_disp_m1 <- 0.34906
wasp_disp_m1 <- 0.34906 * wasp.var.list[pick.wasp.var]

harvesting.length <- 28
day.interval <- harvesting.length/n_fields
# max_t needs to be a multiple of the harvesting length
max_t <- harvesting.length*500

extinct_N <- 1e-10

min.surv <- .01
max.surv <- .04
n.events <- round(max_t/day.interval)

perturb <- data.frame(when=rep(day.interval*(1:n.events), each=3),
                      where=rep(1:n_fields, each=3),
                      who=c("resistant","susceptible","mummies"),
                      how=runif(3*n.events*n_fields, min=min.surv, max=max.surv)) |>
    as_tibble() |>
    # kill all mummies
    mutate(how = ifelse(who == "mummies", 0, how)) |>
    slice(1:(3*n.events))

# sim <- sim_experiments(clonal_lines = c(line_s, line_r),
#                        n_fields = n_fields,
#                        wasp_density_0 = wasp_density_0,
#                        wasp_disp_m0 = wasp_disp_m0,
#                        wasp_disp_m1 = wasp_disp_m1,
#                        perturb = perturb,
#                        extinct_N = 1,
#                        wasp_field_attract = wasp_field_attract,
#                        max_t = max_t)
#
# sim$aphids |>
#     filter(time == 1, is.nan(N)) |> distinct(field)
#
#
#
# sim0 <- sim_experiments(clonal_lines = c(line_s, line_r),
#                        n_fields = n_fields,
#                        wasp_density_0 = wasp_density_0,
#                        wasp_disp_m0 = 0,
#                        wasp_disp_m1 = 0,
#                        perturb = perturb,
#                        extinct_N = 1e-2,
#                        wasp_field_attract = 1,
#                        max_t = max_t)
#
# sim0$aphids |>
#     filter(time == max(time), !is.na(line)) |>
#     group_by(line) |>
#     summarize(N = sum(N))
# sim0$wasps |>
#     filter(time == max(time)) |>
#     summarize(wasps = sum(wasps))
#
#
# sim$aphids |>
#     filter(time == max(time), !is.na(line)) |>
#     group_by(line) |>
#     summarize(N = sum(N))
# sim$wasps |>
#     filter(time == max(time)) |>
#     summarize(wasps = sum(wasps))
#
#
#
#
# safe_pals <- list(main = c("#000000", "#2271B2", "#3DB7E9", "#F748A5",
#                            "#359B73", "#d55e00", "#e69f00", "#f0e442"),
#                   alt = c("#000000", "#AA0DB4", "#FF54ED", "#00B19F",
#                           "#EB057A", "#F8071D", "#FF8D1A", "#9EFF37"))
# # pt_dyn_df <- sim$aphids |>
# pt_dyn_df <- sim0$aphids |>
#     filter(!is.na(line),
#            # field %in% round(seq(1, 28, length.out = 3)),
#            time > harvesting.length*400,
#            time < harvesting.length*500) |>
#     group_by(time, field, line) |>
#     summarize(N = sum(N), .groups = "drop") |>
#     pivot_wider(names_from = "line", values_from = "N") |>
#     mutate(time = time - min(time),
#            prop = resistant / (resistant + susceptible),
#            field = factor(field))
# pt_pts_df <- pt_dyn_df |>
#     group_by(field) |>
#     filter(prop %in% range(prop)) |>
#     mutate(type = ifelse(prop == min(prop), "trough", "peak")) |>
#     group_by(field, prop) |>
#     # Make sure there's only one per field and proportion:
#     filter(time == min(time)) |>
#     ungroup() |>
#     select(time, field, type, prop)
# pt_summ_df <- pt_pts_df |>
#     group_by(type) |>
#     summarize(hi = max(prop), prop = min(prop), .groups = "drop") |>
#     mutate(time = max(pt_dyn_df$time) * 1.05)
#
# pt_dyn_df |>
#     ggplot(aes(time, prop)) +
#     geom_hline(yintercept = 0, color = "gray70") +
#     geom_hline(yintercept = 1, color = "gray70") +
#     geom_line(aes(color = field), linewidth = 0.75, alpha = 0.25) +
#     geom_point(data = pt_pts_df, aes(color = field, shape = type),
#                size = 3) +
#     geom_errorbar(data = pt_summ_df, aes(ymin = prop, ymax = hi),
#                   width = 40, color = safe_pals$main[2], linewidth = 1) +
#     geom_text(data = pt_summ_df, aes(y = (prop + hi) / 2, label = type),
#               color = safe_pals$main[2], hjust = 0, nudge_x = 50,
#               fontface = "bold", size = 12 / 2.8) +
#     scale_color_viridis_d(guide = "none", option = "A") +
#     scale_shape_manual(values = c(17, 15), guide = "none") +
#     scale_y_continuous("Proportion resistant", limits = c(0, 1)) +
#     scale_x_continuous("Days", limits = c(0, max(pt_dyn_df$time) * 1.125),
#                        breaks = seq(0, max(pt_dyn_df$time) %/% 500 * 500, 500))


set.seed(1957298725)
perturb10 <- tibble(when = rep(1:max_t, each = 3),
                    where = when %% 10,
                    who = rep(c("resistant", "susceptible", "mummies"), max_t),
                    how = runif(max_t*3, min=min.surv, max=max.surv)) |>
    mutate(where = ifelse(where == 0, 10, where),
           how = ifelse(who == "mummies", 0, how))


# LEFT OFF ----
#' For `sf2` below, why are aphids in field 1 `NaN` starting at time 9?
#' This only happens when `wasp_disp_m1 != 0` and when a perturbation
#' makes aphids in field 1 turn to zero before the wasps are added.
#' Similarly, a wasp delay <= 7 removes this issue, which is still
#' after field 1 aphids go to zero. Very weird.

wp <- wasp_pop()
line_s <- clonal_line("susceptible",
                      distr_0 = cbind(c(0,0,0,0,1), rep(0, 5)),
                      surv_juv_apterous = "high",
                      surv_adult_apterous = "high",
                      repro_apterous = "high")
line_r <- clonal_line("resistant",
                      distr_0 = cbind(c(0,0,0,0,1), rep(0, 5)),
                      resistant = TRUE,
                      surv_paras = 0.57,
                      surv_juv_apterous = "low",
                      surv_adult_apterous = "low",
                      repro_apterous = "low")
af <- all_fields(c(line_s, line_r), wp, n_fields = 10, K = 12500,
                 aphid_density_0 = 32,
                 wasp_density_0 = 3,
                 wasp_delay = 8,
                 wasp_disp_m0 = 0,
                 wasp_disp_m1 = 0,
                 extinct_N = 1,
                 wasp_field_attract = 1,
                 pred_rate = 0.1)
sf <- sim_fields(af,
                 perturb = perturb10,
                 max_t = 100) ##max_t)

af2 <- all_fields(c(line_s, line_r), wp, n_fields = 10, K = 12500,
                 aphid_density_0 = 32,
                 wasp_density_0 = 3,
                 wasp_delay = 8,
                 wasp_disp_m0 = wasp_disp_m0,
                 wasp_disp_m1 = wasp_disp_m1,
                 extinct_N = 1,
                 wasp_field_attract = wasp_field_attract[1:10],
                 pred_rate = 0.1)
sf2 <- sim_fields(af2,
                  perturb = perturb10,
                  max_t = 100) ##max_t)

sf2$aphids |>
    filter(!is.na(line)) |>
    group_by(time, field, line) |>
    summarize(N = sum(N), .groups = "drop") |>
    filter(is.nan(N))







sf2$aphids |>
    filter(is.nan(N))
sf2$aphids |>
    filter(time == 9 & field == 1 & !is.nan(N))
sf2$aphids |>
    filter(time == 8 & field == 1)




sf$aphids |>
    filter(!is.na(line)) |>
    group_by(time, line) |>
    summarize(N = sum(N), .groups = "drop") |>
    ggplot(aes(time, N, color = line)) +
    geom_line() +
    # geom_line(data = sf$wasps, aes(y = wasps * 5e3), color = "#440154FF") +
    scale_color_manual(values = col_pal)
sf2$aphids |>
    filter(!is.na(line)) |>
    group_by(time, line) |>
    summarize(N = sum(N), .groups = "drop") |>
    ggplot(aes(time, N, color = line)) +
    geom_line() +
    # geom_line(data = sf$wasps, aes(y = wasps * 5e3), color = "#440154FF") +
    scale_color_manual(values = col_pal)



rprop_df <- sf$aphids |>
    filter(!is.na(line)) |>
    group_by(time, field, line) |>
    summarize(N = sum(N), .groups = "drop") |>
    pivot_wider(names_from = "line", values_from = "N") |>
    mutate(time = time - min(time),
           prop = resistant / (resistant + susceptible),
           field = factor(field))


rprop_df |>
    ggplot(aes(time, prop)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_hline(yintercept = 1, color = "gray70") +
    geom_line(aes(color = field), linewidth = 0.75, alpha = 0.25) +
    scale_color_viridis_d(guide = "none", option = "A") +
    scale_y_continuous("Proportion resistant", limits = c(0, 1)) +
    scale_x_continuous("Days")
