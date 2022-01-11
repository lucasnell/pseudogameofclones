
library(tidyverse)
library(clonewars)


source(".Rprofile")


# Date of most recent downloaded datasheet:
.date <- "2022-01-06"

#'
#' Reps to exclude from plot.
#' - Rep 7 failed bc wasps got into the cage
#'
#'

.exclude_reps <- c(7)

exp_df <- read_csv(paste0("~/Box Sync/eco-evo_experiments/results_csv/",
                          .date, "_eco-evo_datasheet_round-2.csv"),
         col_types = cols()) %>%
    filter(! rep %in% .exclude_reps) %>%
    mutate(date = as.Date(date, format = "%m/%d/%Y"),
           start_date = as.Date(start_date, format = "%m/%d/%Y"),
           days = difftime(date, start_date, units = "days") %>%
               as.integer(),
           # Adjusting for Calvin's different way of counting alate OUT
           # before "2021-08-05"
           alates_total_red = ifelse(date < as.Date("2021-08-05") &
                                       observer == "CES",
                                   alates_total_red * 2, alates_total_red),
           alates_total_green = ifelse(date < as.Date("2021-08-05") &
                                         observer == "CES",
                                   alates_total_green * 2,
                                   alates_total_green)) %>%
    # Adjusting for these not being input for isolated treatment:
    mutate(across(starts_with("alates"), ~ ifelse(treatment == "ISOLATED" &
                                                      is.na(.), 0, .))) %>%
    mutate(treatment = treatment %>%
               tolower() %>%
               factor(levels = c("dispersal", "isolated")),
           cage = cage %>%
               tolower() %>%
               factor(levels = c("wasp", "no wasp")),
           rep = factor(rep, levels = sort(unique(rep))),
           cage_id = interaction(treatment, rep, cage, drop = TRUE, sep = "__",
                                 lex.order = TRUE),
           trt_id = interaction(treatment, rep, drop = TRUE, sep = ", rep ",
                                lex.order = TRUE)) %>%
    select(cage_id, trt_id, everything(), -start_date, -red_line,
           -green_line) %>%
    # To show when an experimental cage was terminated:
    group_by(cage_id) %>%
    mutate(terminated = date == max(date[plant1_red >= 0])) %>%
    ungroup() %>%
    # Using 2nd latest date to prevent all cages that weren't sampled on the
    # latest date from returning TRUE here
    mutate(terminated = (terminated & date < sort(unique(date),
                                                  decreasing = TRUE)[2]))



# Should be zero!
exp_df %>%
    filter(if_any(starts_with("plant"), is.na)) %>%
    nrow()

# Should be TRUE!
nrow(distinct(exp_df, treatment, rep)) == nrow(distinct(exp_df, rep))





# Cage-level wasp/mummy data:
wasp_cage_df <- exp_df %>%
    select(cage_id, trt_id, treatment, rep, cage, days, wasps, mummies,
           wasps_removed) %>%
    filter(wasps >= 0)



# Cage-level aphid abundance data:
aphid_cage_df <- exp_df %>%
    mutate(red = rowSums(select(., starts_with("plant") & ends_with("_red"))),
           green = rowSums(select(., starts_with("plant") &
                                      ends_with("_green")))) %>%
    select(-contains("plant"), -starts_with("wasps"), -mummies,
           -observer, -notes, -starts_with("alates")) %>%
    pivot_longer(all_of(c("red", "green")), names_to = "line",
                 values_to = "aphids") %>%
    mutate(line = case_when(line == "red" ~ "susceptible",
                            line == "green" ~ "resistant") %>%
               factor(levels = c("resistant", "susceptible"))) %>%
    filter(aphids >= 0) %>%
    mutate(log_aphids = log1p(aphids)) %>%
    select(cage_id, trt_id, treatment, rep, cage, days, line, aphids,
           log_aphids, terminated)



# Cage-level alate and replaced plant data:
alate_cage_df <- exp_df %>%
    mutate(n_replaced = replaced_plants %>%
               str_split(",") %>%
               map_int(length)) %>%
    select(cage_id, trt_id, treatment, rep, cage, days,
           starts_with("alates_"), n_replaced) %>%
    pivot_longer(starts_with("alates_total_"), names_to = "line",
                 values_to = "alates_total") %>%
    mutate(line = str_remove(line, "alates_total_"),
           alates_in = ifelse(line == "red", alates_in_red, alates_in_green),
           line = case_when(line == "red" ~ "susceptible",
                            line == "green" ~ "resistant") %>%
               factor(levels = c("resistant", "susceptible"))) %>%
    select(cage_id, trt_id, treatment, rep, cage, days, line,
           alates_total, alates_in, n_replaced)


aphid_y <- "aphids"

wasp_mod <- max(aphid_cage_df[[aphid_y]], na.rm = TRUE) /
    max(max(wasp_cage_df$wasps), max(alate_cage_df$alates_in, na.rm = TRUE))

# Use col `Y` when facet scales are free, `Z` when they're not
lab_df <- aphid_cage_df %>%
    split(.$trt_id, drop = TRUE) %>%
    map_dfr(function(.d) {
        max_aphids <- max(.d[[aphid_y]], na.rm = TRUE)
        .dw <- wasp_cage_df %>%
            filter(trt_id == .d$trt_id[1])
        .da <- alate_cage_df %>%
            filter(trt_id == .d$trt_id[1])
        max_y2 <- max(max(.dw$wasps, na.rm = TRUE),
                      max(.da$alates_in, na.rm = TRUE)) * wasp_mod
        .d %>%
            filter(cage == "wasp", days == 0) %>%
            distinct(trt_id, cage, days) %>%
            mutate(Y = max(max_aphids, max_y2),
                   Z = max(aphid_cage_df[[aphid_y]], na.rm = TRUE))
    })



aphid_cage_p <- aphid_cage_df %>%
    ggplot(aes(days)) +
    geom_area(data = wasp_cage_df %>%
                  mutate(wasps = wasps * wasp_mod),
              aes(y = wasps), fill = "gray80", color = NA) +
    geom_hline(yintercept = 0, color = "gray60") +
    geom_vline(data = aphid_cage_df %>%
                   filter(rep %in% 6:10) %>%
                   distinct(trt_id) %>%
                   mutate(x = difftime(as.Date("2021-08-23"),
                                       as.Date("2021-06-14"),
                                       units = "days") %>% as.integer()),
               aes(xintercept = x), linetype = 3, color = "gray70") +
    geom_vline(data = aphid_cage_df %>%
                   filter(rep %in% 6:10) %>%
                   distinct(trt_id) %>%
                   mutate(x = difftime(as.Date("2021-07-26"),
                                       as.Date("2021-06-14"),
                                       units = "days") %>% as.integer()),
               aes(xintercept = x), linetype = 2, color = "gray70") +
    geom_point(data = aphid_cage_df %>%
                     filter(rep == 8, cage == "wasp") %>%
                     distinct(trt_id, cage) %>%
                     mutate(days = difftime(as.Date("2021-09-16"),
                                            as.Date("2021-06-14"),
                                            units = "days") %>% as.integer()),
               aes(y = max(aphid_cage_df[[aphid_y]]) * 0.9),
               shape = 8, size = 3) +
    geom_point(data = aphid_cage_df %>%
                     filter(rep == 8, cage == "wasp") %>%
                     distinct(trt_id, cage) %>%
                     mutate(days = difftime(as.Date("2021-12-13"),
                                            as.Date("2021-06-14"),
                                            units = "days") %>% as.integer()),
               aes(y = max(aphid_cage_df[[aphid_y]]) * 0.9),
               shape = 3, size = 3) +
    # # Points for number of plants replaced:
    # geom_point(data = alate_cage_df %>%
    #                filter(line == "resistant", n_replaced > 0),
    #            aes(y = n_replaced * wasp_mod),
    #           size = 0.75, color = "gray40") +
    # Points for number of alates input:
    geom_point(data = alate_cage_df %>% filter(alates_in > 0),
               aes(y = alates_in * wasp_mod, color = line), size = 0.75) +
    geom_line(aes_(y = as.name(aphid_y), color = ~line), size = 0.75) +
    # # Points for when no aphids of a line present:
    # geom_point(data = aphid_cage_df %>%
    #                filter(aphids == 0) %>%
    #                mutate(Y = max(aphid_cage_df[[aphid_y]]) *
    #                           ifelse(line == "resistant", 0.1, 0)),
    #           aes_(y = ~Y, color = ~line), shape = 4, size = 3) +
    #
    # Points for when experiment terminated:
    geom_point(data = aphid_cage_df %>% filter(terminated),
              aes_(y = as.name(aphid_y), color = ~line), shape = 4, size = 3) +
    geom_text(data = lab_df,
              aes(y = Z, label = trt_id),
              # # use this when you have new reps:
              # aes(y = Y, label = trt_id),
              size = 10 / 2.8, hjust = 0, vjust = 1) +
    facet_grid(trt_id ~ cage,
               scales = "fixed") +
               # # use this when you have new reps:
               # scales = "free_y") +
    scale_color_manual(NULL, values = c("chartreuse3", "firebrick")) +
    scale_y_continuous(ifelse(aphid_y == "log_aphids",
                              "log(aphid abundance + 1)",
                              "Aphid abundance"),
                       sec.axis = sec_axis(~ . / wasp_mod,
                                           paste0("Wasps (gray shading)\n",
                                                  "Alates input (points)"))) +
    scale_x_continuous("Day of experiment") +
    theme(strip.text.y = element_blank(), legend.position = "top") +
    NULL

aphid_cage_p


if (!file.exists(sprintf("~/Desktop/eco-evo_cages__%s.pdf", .date))) {
    ggsave(sprintf("~/Desktop/eco-evo_cages__%s.pdf", .date),
           aphid_cage_p, width = 6.5, height = 8)
}



# Checking specific situations:
aphid_cage_df %>%
    filter(days %in% sort(unique(days), decreasing = TRUE)[1:2]) %>%
    filter(rep %in% c(9,10), line == "resistant", cage == "no wasp") %>%
    select(treatment, rep, days, cage, line, aphids) %>%
    arrange(treatment, rep, days, cage, line)



#  alates ----

# Alate production decreasing?
# alate_p <-
aphid_cage_df %>%
    filter(treatment != "isolated", cage == "no wasp") %>%
    # filter(days > 7) %>%
    group_by(trt_id, treatment, rep, days, line) %>%
    summarize(across(c(aphids, alates_total), sum), .groups = "drop") %>%
    mutate(P = ifelse(aphids > 0, alates_total / aphids, NA)) %>%
    {stopifnot(all(.$P <= 1)); .} %>%
    {amod <<- max(.$P, na.rm = TRUE) / (max(.$aphids) * 1.2); .} %>%
    ggplot(aes(days, P, color = line)) +
    stat_summary(aes(y = aphids * amod, group = trt_id), geom = "line",
                 fun = sum, color = "gray70") +
    geom_line(aes(y = aphids * amod)) +
    geom_text(data = aphid_cage_df %>%
                  filter(trt_id == "dispersal, rep 10") %>%
                  filter(aphids == max(aphids)),
              aes(days + 5, aphids * amod * 1.1, label = "both\nlines"),
              hjust = 0, color = "gray50", size = 9 / 2.8, lineheight = 0.75) +
    geom_point() +
    facet_wrap(~ trt_id) +
    scale_color_manual("aphid line", values = c("chartreuse3", "firebrick")) +
    scale_y_continuous("Counted alates / aphids\n(points)",
                       sec.axis = sec_axis(~ . / amod,
                                           "Total aphids\n(lines)")) +
    scale_x_continuous("Day of experiment", breaks = seq(0, 50, 25))

alate_p

ggsave(sprintf("~/Desktop/eco-evo_alates__%s.pdf", .date),
       alate_p, width = 6.5, height = 4)






# observer differences ----



# Expected vs observed alates IN:
exp_df %>%
    filter(treatment != "isolated") %>%
    select(observer, rep, cage, days, starts_with("alates")) %>%
    arrange(rep, days, cage) %>%
    group_by(rep, days) %>%
    summarize(observer = paste(unique(observer), collapse = "__"),
              in_pred = list(c(rev(alates_total_red) / 2,
                               rev(alates_total_green) / 2)),
              in_obs = list(c(alates_in_red, alates_in_green)),
              line = list(rep(c("red", "green"), each = 2)),
              cage = list(rep(unique(cage), 2)),
              .groups  = "drop") %>%
    filter(!grepl("__", observer)) %>%
    unnest(in_pred:cage) %>%
    # problematic day:
    # filter(in_pred > 7, in_obs < 2, observer == "CES")
    ggplot(aes(in_pred, in_obs, color = cage)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray70") +
    geom_point() +
    facet_wrap(~ observer, nrow = 1) +
    coord_equal() +
    NULL


# Differences in plant replacement by observer:
exp_df %>%
    mutate(n_replaced = replaced_plants %>%
               str_split(",") %>%
               map_int(length)) %>%
    select(observer, n_replaced) %>%
    ggplot(aes(observer, n_replaced)) +
    geom_jitter() +
    stat_summary(fun.data = "mean_cl_boot", color = "red")








# =============================================================================*
# =============================================================================*

# For Schoville lab meeting -  9 Dec 2021

# =============================================================================*
# =============================================================================*



bounds_df <- aphid_cage_df %>%
    summarize(days = list(range(days)),
              !!as.name(aphid_y) := list(range(!!as.name(aphid_y)))) %>%
    unnest(cols = c(days, !!as.name(aphid_y)))



trt <- "dispersal"
# trt <- "isolated"


p <- aphid_cage_df %>%
    filter(treatment == trt) %>%
    ggplot(aes(days)) +
    # To set bounds independently of facet_* call:
    geom_blank(data = bounds_df) +
    geom_area(data = wasp_cage_df %>%
                  filter(treatment == trt) %>%
                  mutate(wasps = wasps * wasp_mod),
              aes(y = wasps), fill = "gray80", color = NA) +
    geom_hline(yintercept = 0, color = "gray60") +
    geom_vline(data = aphid_cage_df %>%
                   filter(treatment == trt) %>%
                   filter(rep %in% 6:10) %>%
                   distinct(trt_id) %>%
                   mutate(x = difftime(as.Date("2021-08-23"),
                                       as.Date("2021-06-14"),
                                       units = "days") %>% as.integer()),
               aes(xintercept = x), linetype = 3, color = "gray70") +
    geom_vline(data = aphid_cage_df %>%
                   filter(treatment == trt) %>%
                   filter(rep %in% 6:10) %>%
                   distinct(trt_id) %>%
                   mutate(x = difftime(as.Date("2021-07-26"),
                                       as.Date("2021-06-14"),
                                       units = "days") %>% as.integer()),
               aes(xintercept = x), linetype = 2, color = "gray70") +
    geom_line(aes_(y = as.name(aphid_y), color = ~line), size = 0.75) +
    # Points for when experiment terminated:
    geom_point(data = aphid_cage_df %>%
                   filter(treatment == trt) %>%
                   filter(terminated),
               aes_(y = as.name(aphid_y), color = ~line), shape = 4, size = 3) +
    facet_grid(trt_id ~ cage,
               scales = "fixed") +
    scale_color_manual(NULL, values = c("chartreuse3", "firebrick")) +
    scale_y_continuous(ifelse(aphid_y == "log_aphids",
                              "log(aphid abundance + 1)",
                              "Aphid abundance"),
                       sec.axis = sec_axis(~ . / wasp_mod,
                                           paste0("Wasps (gray shading)\n",
                                                  "Alates input (points)"))) +
    scale_x_continuous("Day of experiment") +
    theme(strip.text.y = element_blank(), legend.position = "right") +
    NULL

if (trt == "dispersal") {
    p <- p +
        geom_point(data = aphid_cage_df %>%
                       filter(treatment == trt) %>%
                       filter(rep == 8, cage == "wasp") %>%
                       distinct(trt_id, cage) %>%
                       mutate(days = difftime(as.Date("2021-09-16"),
                                              as.Date("2021-06-14"),
                                              units = "days") %>% as.integer()),
                   aes(y = max(aphid_cage_df[[aphid_y]]) * 0.9),
                   shape = 8, size = 3)
}

ggsave(paste0("~/Desktop/eco-evo-", trt, ".pdf"),
       p, width = 6.5, height = ifelse(trt == "dispersal", 4, 3))

