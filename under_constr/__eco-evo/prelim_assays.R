


#'
#' This file looks for differences between lines WIA-5D and UT3 so they
#' can (hopefully) be used for eco-evo experiments.
#' We're assessing whether WIA-5D is a better competitor (in the absence
#' of parasitoids) and whether UT3 is more resistant to parasitoids.
#' Data were collected in fall 2020.
#'


suppressPackageStartupMessages({
    library(tidyverse)
    library(readxl)
    library(clonewars)
})



# ========================================================================
# ========================================================================

# Pop. growth ----

# ========================================================================
# ========================================================================



pop_df <- read_excel("~/Box Sync/2020/aphids/eco-evo-prelims.xlsx",
                     sheet = "population-growth") %>%
    mutate(date = as.Date(paste(year, month, day, sep = "-")),
           rep = factor(rep)) %>%
    group_by(rep, line) %>%
    mutate(date = difftime(date, min(date), units = "days") %>% as.integer()) %>%
    select(rep, line, date, num) %>%
    group_by(rep, date) %>%
    mutate(extinct = factor(sum(num) == 0)) %>%
    ungroup()

sc_pop_p <- pop_df %>%
    group_by(rep) %>%
    mutate(num = (num - mean(num)) / sd(num)) %>%
    ungroup() %>%
    pivot_wider(names_from = "line", values_from = "num") %>%
    mutate(diff_z = `WIA-5D` - `UT3`) %>%
    ggplot(aes(date, diff_z)) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
    geom_line(aes(group = rep), color = "gray50") +
    geom_point(aes(fill = diff_z, shape = extinct), color = "gray40", size = 2) +
    scale_fill_gradient2(low = "firebrick", high = "chartreuse", mid = "gray60",
                         midpoint = 0, guide = FALSE) +
    scale_shape_manual(values = c(21, 4), guide = FALSE) +
    scale_y_continuous(expression("Scaled green" - "red"), breaks = c(-2, 0, 2)) +
    scale_x_continuous("Days after start") +
    NULL

# sc_pop_p

ggsave("~/Desktop/scaled_pop_comp.pdf", sc_pop_p, width = 5, height = 3)


pop_p <- pop_df %>%
    mutate(line = factor(line, levels = c("UT3", "WIA-5D")),
           zero = factor(num == 0)) %>%
    ggplot(aes(date, num, color = line)) +
    geom_line() +
    geom_point(aes(shape = zero), size = 2) +
    facet_wrap(~ rep, nrow = 2) +
    scale_color_manual(values = c("chartreuse", "firebrick")) +
    scale_shape_manual(values = c(19, 4), guide = FALSE) +
    scale_y_continuous("Number of aphids") +
    scale_x_continuous("Days after start") +
    NULL

# pop_p

ggsave("~/Desktop/pop_comp.pdf", pop_p, width = 5, height = 3)



sum_diffs <- pop_df %>%
    group_by(rep, line) %>%
    summarize(num = sum(num), .groups = "drop")


get_avg_diff <- function(.data, .trans = identity, .col = "num") {
    mean(.trans(.data[[.col]][.data$line == "WIA-5D"])) -
        mean(.trans(.data[[.col]][.data$line == "UT3"]))
}

obs_diff <- get_avg_diff(sum_diffs)

set.seed(654067829)
diff_perms <- replicate(2000, {
    ..sum_diffs <- sum_diffs
    ..sum_diffs[["line"]] <- sample(..sum_diffs[["line"]])
    get_avg_diff(..sum_diffs)
})

mean(diff_perms > obs_diff)

hist(diff_perms, xlim = c(min(diff_perms),
                          max(max(diff_perms), obs_diff)) * 1.05)
abline(v = obs_diff, lty = 2, col = "red")





# ========================================================================
# ========================================================================

# Wasp resistance ----

# ========================================================================
# ========================================================================


wasp_df <- read_excel("~/Box Sync/2020/aphids/eco-evo-prelims.xlsx",
                      sheet = "wasp-resistance_1") %>%
    # mutate(start = as.Date(paste(`start-year`, `start-month`, `start-day`, sep = "-")),
    #        end = as.Date(paste(`end-year`, `end-month`, `end-day`, sep = "-")),
    #        rep = factor(rep)) %>%
    # select(-starts_with("start-"), -starts_with("end-")) %>%
    select(rep, line, starts_with(c("juv","adult-", "mumm"))) %>%
    mutate(surv = pmin(`juv-assayed` - mummies, `adult-end`) / `juv-assayed`)


# wasp_df %>%
#     select(-juvenile) %>%
#     pivot_longer(starts_with(c("juv-","adult-", "mumm")),
#                  names_to = "id", values_to = "num") %>%
#     mutate(stage = case_when(id == "juv-exposed" ~ 0L,
#                              id == "juv-assayed" ~ 1L,
#                              id == "adult-end" ~ 2L,
#                              id == "mummies" ~ 2L,
#                              TRUE ~ NA_integer_) %>%
#                factor(levels = 0:2,
#                       labels = c("exposure", "assay\nstart", "assay\nend")),
#            type = as.integer(id == "mummies") %>%
#                factor(levels = 0:1, labels = c("alive", "mummy"))) %>%
#     ggplot(aes(stage, num, color = type)) +
#     geom_point() +
#     # geom_line() +
#     facet_wrap(~ rep, nrow = 4)


mummy_p <- wasp_df %>%
    mutate(p_mummies = mummies / `juv-assayed`,
           line = factor(line, levels = c("UT3", "WIA-5D")),
           zero = factor(p_mummies == 0)) %>%
    ggplot(aes(line, p_mummies, color = line, shape = zero)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_jitter(height = 0, width = 0.25) +
    scale_color_manual(values = c("chartreuse", "firebrick"), guide = FALSE) +
    scale_shape_manual(values = c(19, 4), guide = FALSE) +
    ylab("Mummy proportion") +
    NULL

ggsave("~/Desktop/mummy.pdf", mummy_p, width = 3, height = 3)


juv_p <- wasp_df %>%
    mutate(line = factor(line, levels = c("UT3", "WIA-5D")),
           zero = factor(juvenile == 0)) %>%
    ggplot(aes(line, juvenile, color = line, shape = zero)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_jitter(height = 0, width = 0.25) +
    scale_color_manual(values = c("chartreuse", "firebrick"), guide = FALSE) +
    scale_shape_manual(values = c(19, 4), guide = FALSE) +
    ylab("Number of juveniles") +
    NULL

ggsave("~/Desktop/juveniles.pdf", juv_p, width = 3, height = 3)



surv_p <- wasp_df %>%
    mutate(line = factor(line, levels = c("UT3", "WIA-5D")),
           zero = factor(surv == 0)) %>%
    ggplot(aes(line, surv, color = line, shape = zero)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_jitter(height = 0, width = 0.25) +
    scale_color_manual(values = c("chartreuse", "firebrick"), guide = FALSE) +
    scale_shape_manual(values = c(19, 4), guide = FALSE) +
    ylab("Survival proportion") +
    NULL


ggsave("~/Desktop/survival.pdf", surv_p, width = 3, height = 3)



obs_juv <- get_avg_diff(wasp_df, .col = "juvenile")
obs_surv <- get_avg_diff(wasp_df, .col = "surv")

set.seed(678530)
juv_perms <- replicate(2000, {
    ..wasp_df <- wasp_df
    ..wasp_df[["line"]] <- sample(..wasp_df[["line"]])
    get_avg_diff(..wasp_df, .col = "juvenile")
})
surv_perms <- replicate(2000, {
    ..wasp_df <- wasp_df
    ..wasp_df[["line"]] <- sample(..wasp_df[["line"]])
    get_avg_diff(..wasp_df, .col = "surv")
})

mean(juv_perms < obs_juv)
mean(surv_perms < obs_surv)


