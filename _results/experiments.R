


library(tidyverse)
library(gameofclones)
library(patchwork)
library(grid)

source(".Rprofile")


# Palette for the two clonal lines.
# Equivalent to `viridis::viridis(100)[c(85, 10)]`.
clone_pal <- c("#99D83DFF", "#482173FF")

# What to name wasp and no wasp cages for plot:
cage_lvls <- c("parasitism cage" = "wasp", "no parasitism cage" = "no wasp")

# Date of most recent downloaded datasheet:
.date <- "~/Box Sync/gameofclones/results_csv/" %>%
    list.files("_eco-evo_datasheet_round-2.csv") %>%
    str_split("_") %>%
    map_chr(~ .x[[1]]) %>%
    as.Date() %>%
    max()


exp_df <- read_csv(paste0("~/Box Sync/gameofclones/results_csv/",
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
           # Adjusting for CES's different way of counting alate OUT
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
               factor(levels = c("no wasp", "wasp")),
           rep = factor(rep, levels = sort(unique(rep)))) %>%
    select(everything(), -start_date, -red_line,
           -green_line) %>%
    # To show when an experimental cage was terminated bc of an extinction:
    group_by(treatment, rep, cage) %>%
    mutate(terminated = date == max(date[plant1_red >= 0])) %>%
    ungroup() %>%
    # Using `- 4` to prevent all cages that weren't sampled on the
    # latest date from returning TRUE here
    mutate(terminated = (terminated &
                             date < (max(date) - 4) &
                             days < (max(days) - 4)))



# "Pesky" wasps - wasps that made it into no-wasp cages:

# Get date for latest pesky wasp datasheet:
.pesky_date <- "~/Box Sync/gameofclones/results_csv/" %>%
    list.files("pesky-wasps.csv") %>%
    str_split("_") %>%
    map_chr(~ .x[[1]]) %>%
    as.Date() %>%
    max()


pesky_df <- read_csv(paste0("~/Box Sync/gameofclones/results_csv/",
                            .pesky_date, "_pesky-wasps.csv"),
                     col_types = cols()) %>%
    select(rep, date, mummies, starts_with("adult")) %>%
    #' On 3/16 and 3/22 I removed wasps/mummies but didn't record the number.
    filter(!is.na(mummies)) %>%
    rename(females = `adult females`,
           males = `adult males`,
           unsexed = `adults unk.`) %>%
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
           start_date = map_chr(rep, ~ filter(exp_df, rep == .x) %>%
                                    .[["date"]] %>% min() %>% paste()) %>%
               as.Date(),
           days = difftime(date, start_date, units = "days") %>%
               as.integer()) %>%
    select(-start_date) %>%
    filter(date <= max(exp_df$date))




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
    select(treatment, rep, cage, days, date, wasps, mummies, wasps_rm, mumm_rm) %>%
    filter(wasps >= 0) %>%
    mutate(pesky = 0) %>%
    bind_rows(pesky_df %>%
                  select(treatment, rep, cage, days, date, wasps, mummies) %>%
                  mutate(wasps_rm = wasps, mumm_rm = mummies,
                         pesky = 1)) %>%
    mutate(pesky = factor(pesky, levels = 0:1),
           cage = fct_recode(cage, !!!cage_lvls))




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
    mutate(log_aphids = log1p(aphids),
           cage = fct_recode(cage, !!!cage_lvls)) %>%
    select(treatment, rep, cage, days, line, aphids,
           log_aphids, terminated)



# Cage-level alate and replaced plant data:
alate_cage_df <- exp_df %>%
    mutate(n_replaced = replaced_plants %>%
               str_split(",") %>%
               map_int(length)) %>%
    select(treatment, rep, cage, days,
           starts_with("alates_"), n_replaced) %>%
    pivot_longer(starts_with("alates_total_"), names_to = "line",
                 values_to = "alates_total") %>%
    mutate(line = str_remove(line, "alates_total_"),
           alates_in = ifelse(line == "red", alates_in_red, alates_in_green),
           line = case_when(line == "red" ~ "susceptible",
                            line == "green" ~ "resistant") %>%
               factor(levels = c("resistant", "susceptible")),
           cage = fct_recode(cage, !!!cage_lvls)) %>%
    select(treatment, rep, cage, days, line,
           alates_total, alates_in, n_replaced)


# ============================================================================*
# ============================================================================*

# Main plot ----

# ============================================================================*
# ============================================================================*


# Adjust this if you'd prefer to use log-transformed or raw scale:
aphid_y <- c("aphids", "log_aphids")[2]

wasp_mod <- max(max(wasp_cage_df$wasps),
                max(alate_cage_df$alates_in, na.rm = TRUE)) /
    max(aphid_cage_df[[aphid_y]], na.rm = TRUE)

# maximum N used to define y axis limits and items that are near the top
# of plots:
max_N <- max(aphid_cage_df[[aphid_y]]) / 0.9


exp_p_list <- aphid_cage_df %>%
    distinct(treatment, rep) %>%
    arrange(treatment, rep) %>%
    .[["rep"]] %>%
    paste() %>%
    set_names() %>%
    map(
    function(r) {
        # r = levels(aphid_cage_df$rep)[1]
        # rm(r, acd, lcd, wcd, trt_title)
        acd <- aphid_cage_df %>%
            filter(rep == r) %>%
            rename(N = !!sym(aphid_y))
        lcd <- alate_cage_df %>%
            filter(rep == r, alates_in > 0) %>%
            mutate(N = alates_in / wasp_mod)
        wcd <- wasp_cage_df %>%
            filter(rep == r) %>%
            mutate(N = wasps / wasp_mod)
        #' Y breaks and labels will differ depending on whether we're
        #' using log1p(aphids) or just aphids.
        if (aphid_y == "aphids") {
            y_breaks <- 0:2 * 2e3
            y_labs <- paste0(0:2 * 2, "k")
        } else {
            y_breaks <- c(0, log1p(4 * 10^c(1,3)))
            y_labs <- c(0, 4 * 10^c(1,3))
        }
        p <- acd %>%
            ggplot(aes(days, N, color = line))
        p <- p +
            # Vertical line(s) for early termination:
            geom_vline(data = acd %>% filter(terminated),
                       aes(xintercept = days),
                       size = 0.5, linetype = "22", color = "gray60") +
            geom_area(data = wcd, fill = "gray60", color = NA) +
            geom_hline(yintercept = 0, color = "gray70") +
            # Main abundance lines:
            geom_line() +
            scale_color_manual(values = viridis::viridis(100)[c(85, 10)],
                               guide = "none") +
            scale_y_continuous(sec.axis = sec_axis(~ . * wasp_mod,
                                                   breaks = 0:2 * 40),
                               breaks = y_breaks, labels = y_labs) +
            scale_x_continuous("Days", limits = c(0, 250),
                               breaks = 0:5 * 50) +
            facet_grid( ~ cage, scales = "fixed") +
            theme(strip.text = element_blank(),
                  axis.title = element_blank(),
                  axis.text.x = element_blank(),
                  panel.background = element_rect(fill = "transparent"),
                  plot.background = element_rect(fill = "transparent",
                                                 color = NA),
                  legend.background = element_rect(fill = "transparent"),
                  legend.box.background = element_rect(fill = "transparent")) +
            coord_cartesian(clip = FALSE, ylim = c(0, max_N))
        #'
        #' EDIT: not doing either one in main text.
        #' Providing both in supplement.
        #'
        #' #' On day 94 of rep 8 (2021-09-16), we added 3 female wasps to the
        #' #' wasp cage.
        #' #' On day 182 of rep 8 (2021-12-13), we added 3 UT3 alates to the
        #' #' wasp cage.
        #' #' I'm only adding an arrow for the former, and I'll mention the
        #' #' latter in the text.
        #' #'
        #' if (r == "8") {
        #'     p <- p +
        #'         geom_segment(data = tibble(cage = factor(names(cage_lvls)[1]),
        #'                                    days = 94,
        #'                                    N = max_N,
        #'                                    N2 = N - (max_N * 0.075)),
        #'                      aes(xend = days, yend = N2),
        #'                      size = 0.3, linejoin = "mitre", color = "black",
        #'                      arrow = arrow(length = unit(0.1, "lines"),
        #'                                    type = "closed"))
        #' }
        return(p)
    })



# wrap_plots(exp_p_list, ncol = 1)


# for (i in 1:length(exp_p_list)) {
#     fn <- sprintf("_results/plots/experiments/%i-rep%s.pdf",
#                   i, names(exp_p_list)[i])
#     save_plot(fn, exp_p_list[[i]], 4, 1.25)
# }





# ============================================================================*
# ============================================================================*

# Wasp culling effects ----

# ============================================================================*
# ============================================================================*


#' We started culling adult wasps 2x / week on 26 July 2021.
#' We reduced this to 1x / week on 23 Aug 2021.
#' This only affected reps 6, 8, and 10 because the only other rep going
#' at that time (rep 9) had already lost all its aphids by the time we
#' started culling 2x / week.


cull_df <- exp_df %>%
    filter(date <= as.Date("2021-08-23"), rep %in% c(6,8,10),
           cage == "wasp", days >= 7) %>%
    select(rep, date, wasps, mummies)

#' Differences in peak wasp and mummy abundances before and after we
#' started culling:
cull_df %>%
    mutate(cull = date >= as.Date("2021-07-26")) %>%
    group_by(cull, rep) %>%
    summarize(wasps = max(wasps),
              mummies = max(mummies),
              .groups = "drop") %>%
    group_by(cull) %>%
    summarize(wasps = mean(wasps),
              mummies = mean(mummies))



cull_df %>%
    ggplot(aes(date, wasps)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ rep, ncol = 1)








# ============================================================================*
# ============================================================================*

# Plants replaced ----

# ============================================================================*
# ============================================================================*




rplants_mod <- max(alate_cage_df$n_replaced) /
    max(aphid_cage_df$aphids, na.rm = TRUE)


exp_rplants_p_list <- aphid_cage_df %>%
    distinct(treatment, rep) %>%
    arrange(treatment, rep) %>%
    .[["rep"]] %>%
    paste() %>%
    set_names() %>%
    map(
        function(r) {
            # r = levels(aphid_cage_df$rep)[1]
            # rm(r, acd, prd)
            acd <- aphid_cage_df %>%
                filter(rep == r) %>%
                rename(N = aphids)
            prd <- alate_cage_df %>%
                filter(rep == r, line == "resistant") %>%
                mutate(N = n_replaced / rplants_mod)
            p <- acd %>%
                ggplot(aes(days, N / 1e3, color = line)) +
                geom_hline(yintercept = 0, color = "gray70") +
                # Lines for number of plants replaced:
                geom_area(data = prd, fill = "gray60", color = NA) +
                # geom_line(data = prd, size = 0.5, color = "gray60") +
                # Main abundance lines:
                geom_line() +
                # Points for early termination:
                geom_point(data = acd %>% filter(terminated), shape = 4, size = 3) +
                scale_color_manual(values = clone_pal, guide = "none") +
                scale_y_continuous(sec.axis = sec_axis(~ . * rplants_mod * 1000,
                                                       breaks = 0:2 * 5),
                                   limits = c(0, max_N / 1000),
                                   breaks = 0:2 * 2) +
                scale_x_continuous("Days", limits = c(0, 250),
                                   breaks = 0:5 * 50) +
                facet_grid( ~ cage, scales = "fixed") +
                theme(strip.text.x = element_blank(),
                      strip.text.y = element_blank(),
                      axis.title.y = element_blank()) +
                coord_cartesian(clip = FALSE)
            if (r != "13") {
                p <- p +
                    theme(axis.title.x = element_blank(),
                          axis.text.x = element_blank())
            }
            return(p)
        })


exp_rplants_lm_p <- aphid_cage_df %>%
    mutate(id = interaction(treatment, rep, cage, line, drop = TRUE)) %>%
    split(.$id, drop = TRUE) %>%
    map_dfr(function(dd) {
        left_join(select(dd, days, log_aphids),
                  alate_cage_df %>%
                      filter(treatment == dd$treatment[1], rep == dd$rep[1],
                             cage == dd$cage[1], line == dd$line[1]) %>%
                      select(days, n_replaced),
                  by = "days") %>%
            mutate(pcg = log_aphids - lag(log_aphids)) %>%
            filter(!is.na(pcg)) %>%
            mutate(treatment = dd$treatment[1], rep = dd$rep[1],
                   cage = dd$cage[1], line = dd$line[1])
    }) %>%
    ggplot(aes(n_replaced, pcg)) +
    geom_hline(yintercept = 0, color = "gray60") +
    geom_point(alpha = 0.25) +
    stat_smooth(formula = y ~ x, method = "lm", se = TRUE) +
    xlab(expression("Number of plants replaced on day " * italic(t))) +
    ylab(expression("Per-capita growth rate (log[" * italic(N[t]) /
                        italic(N[t-1]) * "])")) +
    theme(strip.text = element_text(size = 10),
          axis.title = element_text(size = 9))


exp_rplants_p <- wrap_elements(grid::textGrob(expression("Aphid abundance ("
                                                         %*% 1000 * ")"),
                             x = 0, vjust = 1, rot = 90)) +
    wrap_plots(exp_rplants_p_list, ncol = 1) +
    wrap_elements(grid::textGrob("Plants replaced (gray shading)",
                                 x = 1, vjust = 1, rot = -90)) +
    exp_rplants_lm_p +
    plot_layout(widths = c(0.06, 1, 0.06), heights = c(1, 0.3),
                nrow = 2, ncol = 3, design = "123
                                              #4#") +
    plot_annotation(tag_levels = list(c("", LETTERS[1:7], "", LETTERS[8]))) &
    theme(plot.tag = element_text(size = 14, face = "bold"))






# save_plot("_results/plots/repl-plants.pdf", exp_rplants_p, 6, 12)





# ============================================================================*
# ============================================================================*

#  Dispersal proportion ----

# ============================================================================*
# ============================================================================*



left_join(aphid_cage_df, alate_cage_df,
          by = c("treatment", "rep", "cage", "days", "line")) %>%
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






# ============================================================================*
# ============================================================================*

# Stats on things the techs counted ----

# ============================================================================*
# ============================================================================*


base_df <- read_csv(paste0("~/Box Sync/gameofclones/results_csv/",
                           .date, "_eco-evo_datasheet_round-2.csv"),
                    col_types = cols()) %>%
    filter(observer != "LAN") %>%
    rowwise() %>%
    mutate(n_aphids = sum(c_across(starts_with("plant"))),
           n_alates_in = sum(c_across(starts_with("alates_in_")))) %>%
    mutate(n_repl_plants = replaced_plants %>%
               str_split(",") %>%
               map_int(length)) %>%
    select(-starts_with("plant"), -starts_with("alates_in_"), -replaced_plants)

# Total aphids counted, alates dispersed, plants replaced:
for (n in c("n_aphids", "n_alates_in", "n_repl_plants")) {
    x <- base_df %>%
        filter(!!sym(n) >= 0) %>%
        .[[n]] %>%
        sum()
    cat(sprintf("Total %s = %i\n", n, x))
}; rm(n, x)

# Total wasps and mummies counted:
for (n in c("wasps", "mummies")) {
    x <- base_df %>%
        filter(rep == 7) %>%
        filter(!!sym(n) >= 0) %>%
        filter(!is.na(!!sym(n))) %>%
        .[[n]] %>%
        sum()
    y <- wasp_cage_df %>%
        filter(!!sym(n) >= 0) %>%
        filter(!is.na(!!sym(n))) %>%
        .[[n]] %>%
        sum()
    cat(sprintf("Total %s = %i\n", n, x + y))
}; rm(n, x, y)

