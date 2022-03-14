
#'
#' This is the "working" edition of the experiment plot, to be used while
#' it's still running.
#'

library(tidyverse)
library(clonewars)

source(".Rprofile")


# # Palette for the two clonal lines.
# # Equivalent to `viridis::viridis(100)[c(70, 10)]`.
# clone_pal <- c("#41BE71FF", "#482173FF")

# Above is used for paper, below is used here
clone_pal <- c("chartreuse3", "firebrick")


# Date of most recent downloaded datasheet:
.date <- "~/Box Sync/eco-evo_experiments/results_csv/" %>%
    list.files("_eco-evo_datasheet_round-2.csv") %>%
    str_split("_") %>%
    map_chr(~ .x[[1]]) %>%
    as.Date() %>%
    max()



exp_df <- read_csv(paste0("~/Box Sync/eco-evo_experiments/results_csv/",
                          .date, "_eco-evo_datasheet_round-2.csv"),
         col_types = cols()) %>%
    #'
    #' Rep 7 failed bc wasps got into the no-wasp cage, and we didn't
    #' adequately respond to this.
    #'
    filter(rep != 7) %>%
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
           trt_id = interaction(treatment, rep, drop = TRUE, sep = ", rep ",
                                lex.order = TRUE)) %>%
    select(trt_id, everything(), -start_date, -red_line,
           -green_line) %>%
    # To show when an experimental cage was terminated:
    group_by(treatment, rep, cage) %>%
    mutate(terminated = date == max(date[plant1_red >= 0])) %>%
    ungroup() %>%
    # Using `- 4` to prevent all cages that weren't sampled on the
    # latest date from returning TRUE here
    mutate(terminated = (terminated &
                             date < (max(date) - 4) &
                             days < (max(days) - 4)))



# Should be zero!
exp_df %>%
    filter(if_any(starts_with("plant"), is.na)) %>%
    nrow()

# Should be TRUE!
nrow(distinct(exp_df, treatment, rep)) == nrow(distinct(exp_df, rep))



# "Pesky" wasps - wasps that made it into no-wasp cages:

# Get date for latest pesky wasp datasheet:
.pesky_date <- "~/Box Sync/eco-evo_experiments/results_csv/" %>%
    list.files("pesky-wasps.csv") %>%
    str_split("_") %>%
    map_chr(~ .x[[1]]) %>%
    as.Date() %>%
    max()


pesky_df <- read_csv(paste0("~/Box Sync/eco-evo_experiments/results_csv/",
                            .pesky_date, "_pesky-wasps.csv"),
                     col_types = cols()) %>%
    select(rep, date, mummies, starts_with("adult")) %>%
    rename(females = `adult females`, males = `adult males`, unsexed = `adults unk.`) %>%
    mutate(wasps = females + males + unsexed,
           date = as.Date(date, "%d-%b-%y")) %>%
    select(rep, date, mummies, wasps, everything()) %>%
    group_by(date, rep) %>%
    summarize_all(sum) %>%
    ungroup() %>%
    # These are changed to make it join easier with the other wasp data:
    mutate(rep = factor(rep, levels = levels(exp_df$rep)),
           cage = factor("no wasp", levels = levels(exp_df$cage)),
           treatment = ifelse(rep %in% c(6,8,10,11),"dispersal","isolated") %>%
               factor(levels = levels(exp_df$treatment)),
           trt_id = paste(treatment, rep, sep = ", rep ") %>%
               factor(levels = levels(exp_df$trt_id)),
           start_date = map_chr(rep, ~ filter(exp_df, rep == .x) %>%
                                    .[["date"]] %>% min() %>% paste()) %>%
               as.Date(),
           days = difftime(date, start_date, units = "days") %>%
               as.integer()) %>%
    select(-start_date)




# ============================================================================*
# ============================================================================*

# Summarize by cage ----

# ============================================================================*
# ============================================================================*

# Cage-level wasp/mummy data:
wasp_cage_df <- exp_df %>%
    # We don't want to include dates from pesky wasps bc it's redundant:
    filter(!interaction(rep, date, cage, drop = TRUE) %in%
               interaction(pesky_df$rep, pesky_df$date, pesky_df$cage,
                           drop = TRUE)) %>%
    select(trt_id, treatment, rep, cage, days, date, wasps, mummies, wasps_rm, mumm_rm) %>%
    filter(wasps >= 0) %>%
    mutate(pesky = 0) %>%
    bind_rows(pesky_df %>%
                  select(trt_id, treatment, rep, cage, days, date, wasps, mummies) %>%
                  mutate(wasps_rm = wasps, mumm_rm = mummies,
                         pesky = 1)) %>%
    mutate(pesky = factor(pesky, levels = 0:1))







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
    select(trt_id, treatment, rep, cage, days, line, aphids,
           log_aphids, terminated)



# Cage-level alate and replaced plant data:
alate_cage_df <- exp_df %>%
    mutate(n_replaced = replaced_plants %>%
               str_split(",") %>%
               map_int(length)) %>%
    select(trt_id, treatment, rep, cage, days,
           starts_with("alates_"), n_replaced) %>%
    pivot_longer(starts_with("alates_total_"), names_to = "line",
                 values_to = "alates_total") %>%
    mutate(line = str_remove(line, "alates_total_"),
           alates_in = ifelse(line == "red", alates_in_red, alates_in_green),
           line = case_when(line == "red" ~ "susceptible",
                            line == "green" ~ "resistant") %>%
               factor(levels = c("resistant", "susceptible"))) %>%
    select(trt_id, treatment, rep, cage, days, line,
           alates_total, alates_in, n_replaced)


# ============================================================================*
# ============================================================================*

# Main plot ----

# ============================================================================*
# ============================================================================*


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
               shape = 8, size = 2.5) +
    geom_point(data = aphid_cage_df %>%
                     filter(rep == 8, cage == "wasp") %>%
                     distinct(trt_id, cage) %>%
                     mutate(days = difftime(as.Date("2021-12-13"),
                                            as.Date("2021-06-14"),
                                            units = "days") %>% as.integer()),
               aes(y = max(aphid_cage_df[[aphid_y]]) * 0.9),
               shape = 3, size = 2.5) +
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
    scale_color_manual(NULL, values = clone_pal) +
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
           aphid_cage_p, width = 7, height = 8)
}




# ============================================================================*
# ============================================================================*

#  Dispersal proportion ----

# ============================================================================*
# ============================================================================*



left_join(aphid_cage_df, alate_cage_df,
          by = c("treatment", "rep", "cage", "days", "line", "trt_id")) %>%
    # Filter out isolated treatments:
    filter(treatment != "isolated") %>%
    # Filter out observations from before we were dispersing alates
    # on the tops of plants:
    filter(!(rep %in% 6:10 &
                 days < as.integer(difftime(as.Date("2021-08-23"),
                                            as.Date("2021-06-14"),
                                            units = "days")))) %>%
    # Filter out the beginnings of each period because it takes some time
    # for alates to show up
    filter(days > 20) %>%
    mutate(dispersed = alates_total / 2) %>%
    # Control for rare scenario where there were lots of alates on the
    # sides/tops but few on plants, so alates_total > aphids
    mutate(aphids = ifelse(alates_total > aphids, alates_total + aphids, aphids)) %>%
    filter(aphids > 0) %>%
    mutate(p_dispersed = dispersed / aphids) %>%
    group_by(cage, line) %>%
    summarize(p_dispersed = mean(p_dispersed), nobs = n(), .groups = "drop")

