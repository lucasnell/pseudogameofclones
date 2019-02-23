

library(tidyverse)
library(readxl)


source(".Rprofile")

comp_df <- read_excel(paste0("~/Dropbox/Aphid Project 2017/Lauren_competition/",
                             "competition_data_entry.xlsx")) %>%
    rename_all(~ tolower(.)) %>%
    filter(!is.na(green_line)) %>%
    mutate(green_line = ifelse(green_line == "Clover-2017-6 (reg+",
                               "Clover-2017-6 (reg+)", green_line),
           red_line = ifelse(red_line == "WI-L4(Ham-)", "WI-L4 (Ham-)", red_line)) %>%
    mutate_at(vars(green_line, red_line), ~ gsub(" \\(Ham\\-\\)", "Ø",
                                                 gsub(" \\(reg\\+\\)| \\(Ham\\+\\)", "",
                                                      .x))) %>%
    mutate(total_red = plant_red + dispersed_red,
           total_green = plant_green + dispersed_green) %>%
    mutate(date = as.Date(paste(year, month, day, sep = "-"))) %>%
    select(-year, -month, -day, -observer, -comments,
           -starts_with("dispersed_"), -starts_with("plant_")) %>%
    select(date, everything()) %>%
    mutate(combo = paste0(green_line, red_line) %>%
               factor() %>% as.integer() %>% factor()) %>%
    # gather("color", "line", green_line, red_line) %>%
    # gather("color2", "total", total_green, total_red) %>%
    # mutate(color = gsub("_line", "", color),
    #        color2 = gsub("total_", "", color2)) %>%
    # filter(color == color2) %>%
    # select(-color2) %>%
    group_by(treatment, combo, rep) %>%
    mutate(date = as.integer(date - min(date))) %>%
    ungroup() %>%
    rename(green = green_line, red = red_line) %>%
    mutate(diff = total_green - total_red)


comp_df %>%
    filter(grepl("^27", treatment)) %>%
    # filter(grepl("^20", treatment)) %>%
    split(interaction(.$combo, .$rep)) %>%
    map_dfr(function(.x) {
        .mean <- mean(c(.x$total_red, .x$total_green))
        .sd <- sd(c(.x$total_red, .x$total_green))
        mutate(.x,
               total_red = (total_red - .mean) / .sd,
               total_green = (total_green - .mean) / .sd,
               diff_z = total_green - total_red)
    }) %>%
    ggplot(aes(date)) +
    # geom_line(aes(y = total_green, linetype = factor(rep))) +
    # geom_line(aes(y = total_red, linetype = factor(rep)), color = "red") +
    geom_hline(yintercept = 0, linetype = 3, color = "gray80") +
    geom_line(aes(y = diff_z, group = factor(rep)), color = "gray50") +
    geom_point(aes(y = diff_z, color = diff_z), shape = 16) +
    facet_grid(green ~ red, labeller = label_both) +
    theme_classic() +
    # scale_color_brewer(palette = "Dark2") +
    scale_color_gradient2(low = "#d95f02", high = "#1b9e77", mid = "gray80", midpoint = 0)

