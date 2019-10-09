

#'
#' This version is for when alates are also counted. They were done summer 2019.
#'

# library(clonewars)
library(tidyverse)
library(readxl)


source(".Rprofile")

fn <- "~/Dropbox/Aphid Project 2017/Lauren_competition/competition_data_entry24July19.xlsx"


# comp_df <-
read_excel(fn) %>%
    rename_all(~ tolower(.)) %>%
    filter(!is.na(green_line)) %>%
    mutate(green_line = ifelse(green_line == "Clover-2017-6 (reg+",
                               "Clover-2017-6 (reg+)", green_line),
           red_line = ifelse(red_line == "WI-L4(Ham-)", "WI-L4 (Ham-)", red_line)) %>%
    .[["green_line"]] %>% unique()
    mutate_at(vars(green_line, red_line), ~ gsub(" \\(Ham\\-\\)", "Ø",
                                                 gsub(" \\(reg\\+\\)| \\(Ham\\+\\)", "",
                                                      .x))) %>%
    mutate(total_red = apterous_red + alate_red,
           total_green = apterous_green + alate_green) %>%
    mutate(date = as.Date(paste(year, month, day, sep = "-"), format = "%Y-%B-%d")) %>%
    select(-year, -month, -day, -observer, -comments) %>%
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


# Filter out when all aphids die out quickly
comp_df <- comp_df %>%
    group_by(green, red, treatment, rep) %>%
    mutate(N = sum({(total_green + total_red) > 0})) %>%
    ungroup() %>%
    filter(N >= 5) %>%
    select(-N)

# Filter out when there is never an increase in either red or green aphids
comp_df <- comp_df %>%
    group_by(green, red, treatment, rep) %>%
    filter(! { all(diff(total_red) <= 0) | all(diff(total_green) <= 0) }) %>%
    ungroup()




comp_plot <- function(.temp) {
    comp_df %>%
        filter(grepl(sprintf("^%s", paste(.temp)), treatment)) %>%
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
        geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
        geom_line(aes(y = diff_z, group = factor(rep)), color = "gray30", alpha = 0.5) +
        geom_point(aes(y = diff_z, fill = diff_z), shape = 21, color = "gray20", size = 2,
                   alpha = 0.5) +
        facet_grid(red ~ green) + #, labeller = label_both) +
        # scale_color_gradient2(low = "#d95f02", high = "#1b9e77", mid = "gray80",
        scale_fill_gradient2(low = "firebrick", high = "chartreuse", mid = "gray60",
                             midpoint = 0, guide = FALSE) +
        scale_y_continuous(expression("Scaled green" - "red"), breaks = c(-2, 0, 2)) +
        scale_x_continuous("Date", breaks = seq(0, 30, 10)) +
        theme(strip.text = element_text(size = 10)) +
        theme_classic()
}

# comp_plot(20) %>%
#     ggsave(filename = "~/Desktop/comp_20.pdf", height = 6, width = 6)
comp_plot(27)
# %>%
#     ggsave(filename = "~/Desktop/comp_27.pdf", height = 6, width = 6)




# Below is to look at each rep separately, to check for things to filter out
for (.rep in unique(comp_df$rep)) {
    .plot <- comp_df %>%
        filter(rep == .rep) %>%
        ggplot(aes(date)) +
        ggtitle(paste("rep #", .rep)) +
        geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
        geom_hline(yintercept = log(2), linetype = 2, color = "gray80") +
        geom_line(aes(y = log(total_green)), color = "green") +
        geom_line(aes(y = log(total_red)), color = "red") +
        geom_point(aes(y = log(total_green)), color = "green") +
        geom_point(aes(y = log(total_red)), color = "red") +
        facet_grid(red ~ green) +
        scale_y_continuous("log(N)", limits = c(-1, NA)) +
        scale_x_continuous("Date", breaks = seq(0, 30, 10), limits = c(0, 30)) +
        theme(strip.text = element_text(size = 10)) +
        theme_classic()
    ggsave(filename = sprintf("~/Desktop/comp/rep_%i.pdf", .rep),
           plot = .plot, height = 6, width = 6)
}; rm(.rep, .plot)



logit <- function(p) log(p / (1 - p))
asinsqrt <- function(p) asin(sqrt(p))
inv_asinsqrt <- function(x) sin(x)^2



alate_df <- comp_df %>%
    gather("color", "line", red:green, factor_key = TRUE) %>%
    mutate(apterous = ifelse(color == "red", apterous_red, apterous_green),
           alate = ifelse(color == "red", alate_red, alate_green)) %>%
    mutate(pr_alate = alate / (alate + apterous),
           pr_alate = ifelse(is.nan(pr_alate), 0, pr_alate),
           total = alate_red + alate_green + apterous_red + apterous_green) %>%
    select(-starts_with("alate_"), -starts_with("apterous_"), -starts_with("total_"),
           -treatment)

alate_df %>%
    group_by(rep, combo, line, color) %>%
    summarize(pr_alate = sum(alate) / sum(alate + apterous)) %>%
    ungroup() %>%
    # ggplot(aes(line, logit(pr_alate))) +
    ggplot(aes(line, asinsqrt(pr_alate))) +
    # ggplot(aes(line, pr_alate)) +
    geom_point(aes(color = color), shape = 1, size = 2,
               position = position_jitter(width = 0.25)) +
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, size = 0.75) +
    stat_summary(fun.y = "mean", geom = "point", size = 3, shape = 124) +
    scale_color_manual(values = c("red", "darkgreen"), guide = FALSE) +
    scale_x_discrete(NULL) +
    scale_y_continuous("Proportion of alates over time series",
                       breaks = asinsqrt(seq(0, 0.12, 0.04)),
                       labels = seq(0, 0.12, 0.04)) +
                       # breaks = ggplot2::waiver()) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10, color = "black")) +
    coord_flip() +
    NULL



alate_df %>%
    # ggplot(aes(apterous, logit(pr_alate))) +
    # ggplot(aes(date, asinsqrt(pr_alate))) +
    ggplot(aes(total, pr_alate)) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
    # geom_line(aes(group = factor(rep)), color = "gray30", alpha = 0.5) +
    geom_point(aes(color = color), shape = 16, size = 2, alpha = 0.5) +
    facet_wrap( ~ line, nrow = 2) +
    scale_color_manual(values = c("red", "green"), guide = FALSE) +
    scale_y_continuous("Alate proportion") +
    # scale_x_continuous("Date", breaks = seq(0, 30, 10)) +
    scale_x_continuous("Total aphids") +
    theme(strip.text = element_text(size = 10)) +
    theme_classic()



alate_df

library(nlme)


nlme(alate ~ exp(total + line), alate_df,
     fixed = total + line ~ 1,
     random = line ~ 1,
     correlation = corAR1(form = ~ date | combo))

alate_df$combooo <- factor(paste(alate_df$combo, alate_df$rep, alate_df$color))
gls(log1p(alate) ~ total + line + 0, alate_df,
     correlation = corAR1(form = ~ date | combooo))






# comp_diffs <- comp_df %>%
#     filter(grepl("^27", treatment)) %>%
#     split(interaction(.$combo, .$rep)) %>%
#     map_dfr(function(.x) {
#         .mean <- mean(c(.x$total_red, .x$total_green))
#         .sd <- sd(c(.x$total_red, .x$total_green))
#         mutate(.x,
#                total_red = (total_red - .mean) / .sd,
#                total_green = (total_green - .mean) / .sd,
#                diff_z = total_green - total_red)
#     }) %>%
#     select(date, green, red, rep, combo, diff_z)
#
#
#
#
#
# library(clonewars)
# source(".Rprofile")
#
# {
#     ggplot(data = tibble(x = 0:10), aes(x)) +
#         stat_function(fun = sin, size = 1) +
#         stat_function(fun = function(x) sin(x - (pi/1.5)), color = "red", size = 1) +
#         coord_cartesian(ylim = c(-1,1)) +
#         theme(axis.title = element_blank(),
#               axis.text = element_blank(),
#               axis.ticks = element_blank())
# } %>%
#     ggsave(filename = "~/Desktop/wa.pdf", width = 4.75, height = 2.25)
#
# {
#     ggplot(data = tibble(x = 0:20), aes(x)) +
#         stat_function(fun = sin, size = 1) +
#         stat_function(fun = function(x) sin(x - (pi/2)), color = "red", size = 1) +
#         coord_cartesian(ylim = c(-1,2)) +
#         theme(axis.title = element_blank(),
#               axis.text = element_blank(),
#               axis.ticks = element_blank())
#     } %>%
#     ggsave(filename = "~/Desktop/wa_beetles.pdf", width = 4.75, height = 2.25)
#
# {
#     ggplot(data = tibble(x = 0:20), aes(x)) +
#         stat_function(fun = function(x) {
#             z <- x
#             z[x <= 1.5 * pi] <- (1 + sin(x[x <= 1.5 * pi])) / 2
#             z[x > 1.5 * pi] <- 1.5 * (1 + sin(x[x > 1.5 * pi])) / 2
#             z[x > 2.5 * pi] <- 0.5 + (1 + sin(x[x > 2.5 * pi])) / 2
#             return(z)
#         },
#         size = 1) +
#         stat_function(fun = function(x) {
#             ifelse(x <= 2 * pi, 1, 0.5) * (1 + sin(x - (pi/2))) / 2
#         },
#         color = "red", size = 1) +
#         theme(axis.title = element_blank(),
#               axis.text = element_blank(),
#               axis.ticks = element_blank()) +
#         NULL
#     } %>%
#     ggsave(filename = "~/Desktop/wa_beetles2.pdf", width = 4.75, height = 2.25)
#
#
#
# #

















