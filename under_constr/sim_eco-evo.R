
library(tidyverse)
library(clonewars)


ut3 <- clonal_line("UT3",
                   density_0 = matrix(c(rep(0, 3), 16, rep(0, 6)), 5, 2),
                   resistant = TRUE,
                   # resistant = c(1, 1),
                   surv_juv_paras = 0,
                   surv_juv_apterous = "low",
                   surv_adult_apterous = "low",
                   repro_apterous = "low")
wia5d <- clonal_line("WIA-5D",
                     density_0 = matrix(c(rep(0, 3), 16, rep(0, 6)), 5, 2),
                     surv_juv_apterous = "high",
                     surv_adult_apterous = "high",
                     repro_apterous = "high")

sims <- sim_clonewars(1, c(ut3, wia5d),
                      max_t = 1000,
                      alate_b0 = -4, alate_b1 = 8/2000,
                      # alate_b0 = -1000, alate_b1 = 0,
                      max_N = 0, max_plant_age = 25,
                      # death_prop = 0.4,
                      # shape1_death_mort = 0.5,
                      # shape2_death_mort = 0,
                      clear_surv = 0.5,
                      mean_K = formals(sim_clonewars)$mean_K * 4,
                      sd_K = formals(sim_clonewars)$sd_K * 4,
                      s_y = populations$s_y / 10,
                      no_error = TRUE,
                      # environ_error = TRUE,
                      # sigma_x = 0,
                      # sigma_y = 0,
                      sex_ratio = 0.66,
                      wasp_delay = 5,
                      wasp_density_0 = 4)


mod <- max(max(sims$wasps$wasps),
           max(sims$aphids$N[sims$aphids$type == "mummy"])) /
    max(sims$aphids$N[sims$aphids$type != "mummy"])

sims %>%
    .[["aphids"]] %>%
    filter(type != "mummy") %>%
    # filter(rep == 0) %>%
    # filter(time < 10) %>%
    mutate(line = factor(line, levels = c("UT3", "WIA-5D")),
           type = factor(type, levels = c("apterous", "alate"))) %>%
    group_by(rep, patch, time, line) %>%
    summarize(N = sum(N), .groups = "drop") %>%
    ggplot(aes(time, N)) +
    geom_area(data = sims %>%
                  .[["aphids"]] %>%
                  filter(type == "mummy") %>%
                  # filter(rep == 0) %>%
                  # filter(time < 10) %>%
                  select(-line, -type) %>%
                  group_by(rep, patch, time) %>%
                  summarize(N = sum(N), .groups = "drop") %>%
                  mutate(N = N / mod),
              fill = "gray80", color = "gray50", size = 0.1) +
    geom_line(data = sims %>%
                  .[["wasps"]] %>%
                  # filter(rep == 0) %>%
                  # filter(time < 10) %>%
                  mutate(rep = factor(rep),
                         N = wasps / mod),
              size = 1, color = "gray50") +
    geom_line(aes(color = line)) +
    # facet_grid(rep ~ patch) +
    facet_wrap( ~ patch) +
    scale_color_manual(values = c("chartreuse3", "firebrick")) +
    scale_linetype_manual(values = c(1, 2)) +
    scale_y_continuous("Aphid abundance",
                       sec.axis = sec_axis(~ . * mod,
                                           "Wasp & mummy abundance"))


sims %>%
    .[["aphids"]] %>%
    filter(type != "mummy") %>%
    mutate(line = factor(line, levels = c("UT3", "WIA-5D")),
           type = factor(type, levels = c("apterous", "alate"))) %>%
    group_by(rep, time, line) %>%
    summarize(N = sum(N), .groups = "drop") %>%
    ggplot(aes(time, N)) +
    geom_area(data = sims %>%
                  .[["aphids"]] %>%
                  filter(type == "mummy") %>%
                  group_by(rep, time) %>%
                  summarize(N = sum(N), .groups = "drop") %>%
                  mutate(N = N / mod),
              fill = "gray80", color = "gray50", size = 0.1) +
    geom_line(data = sims %>%
                  .[["wasps"]] %>%
                  mutate(rep = factor(rep),
                         N = wasps / mod),
              size = 1, color = "gray50") +
    geom_line(aes(color = line)) +
    scale_color_manual(values = c("chartreuse3", "firebrick")) +
    scale_linetype_manual(values = c(1, 2)) +
    scale_y_continuous("Aphid abundance",
                       sec.axis = sec_axis(~ . * mod,
                                           "Wasp & mummy abundance"))




sims %>%
    .[["aphids"]] %>%
    filter(type == "mummy") %>%
    group_by(rep, time) %>%
    summarize(N = sum(N), .groups = "drop") %>%
    filter(rep == 0) %>%
    ggplot(aes(time, N)) +
    geom_line() +
    geom_line(data = sims %>%
                  .[["wasps"]] %>%
                  filter(rep == 0) %>%
                  # filter(time < 10) %>%
                  mutate(rep = factor(rep),
                         N = wasps), color = "blue")

sims %>%
    .[["aphids"]] %>%
    filter(type == "mummy") %>%
    filter(rep == 0) %>%
    # filter(time < 10) %>%
    select(-line, -type) %>%
    .[["N"]] %>% max()

sims %>%
    .[["aphids"]] %>%
    filter(type != "mummy") %>%
    group_by(time, patch) %>%
    summarize(N = sum(N), .groups = "drop") %>%
    ggplot(aes(time, N)) +
    geom_line() +
    facet_wrap(~ patch)

