

#'
#' This version is for when alates are also counted. They were done summer 2019.
#'

library(tidyverse)
library(readxl)


source(".Rprofile")
source("under_constr/theme_black.R")


# Lines in order they should appear in plots:
lines_ <- c("R10", "WIA-5D", "WI-L4", "WI-L4Ø", "UT3", "WI-2016-593", "Clover-2017-2",
            "Clover-2017-6")


fn <- paste0("~/Dropbox/Aphid Project 2017/Lauren_competition/pre-summer_2019/",
             "competition_data_entry.xlsx")


old_comp_df <- read_excel(fn) %>%
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
    group_by(treatment, combo, rep) %>%
    mutate(date = as.integer(date - min(date))) %>%
    ungroup() %>%
    rename(green = green_line, red = red_line) %>%
    mutate(diff = total_green - total_red,
           treatment = gsub("°", " ", treatment))


# Filter out when all aphids die out quickly
old_comp_df <- old_comp_df %>%
    group_by(green, red, treatment, rep) %>%
    mutate(N = sum({(total_green + total_red) > 0})) %>%
    ungroup() %>%
    filter(N >= 5) %>%
    select(-N)

# Filter out when there is never an increase in either red or green aphids
old_comp_df <- old_comp_df %>%
    group_by(green, red, treatment, rep) %>%
    filter(! { all(diff(total_red) <= 0) | all(diff(total_green) <= 0) }) %>%
    ungroup()



fn <- paste0("~/Dropbox/Aphid Project 2017/Lauren_competition/",
             "competition_data_entry24July19.xlsx")



comp_df <- read_excel(fn) %>%
    rename_all(~ tolower(.)) %>%
    filter(!is.na(green_line)) %>%
    mutate(red_line = ifelse(red_line == "WI-L4 Ham-", "WI-L4Ø", red_line)) %>%
    mutate_at(vars(green_line, red_line), ~ gsub(" Reg\\+| Ham\\+", "",
                                                      .x)) %>%
    mutate(total_red = apterous_red + alate_red,
           total_green = apterous_green + alate_green) %>%
    mutate(date = as.Date(paste(year, month, day, sep = "-"), format = "%Y-%B-%d")) %>%
    select(-year, -month, -day, -observer, -comments) %>%
    select(date, everything()) %>%
    mutate(combo = paste0(green_line, red_line) %>%
               factor() %>% as.integer() %>% factor()) %>%
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




comp_plot <- function(.temp, .df = comp_df) {
    .df %>%
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
        scale_x_continuous("Day", breaks = seq(0, 30, 10)) +
        theme(strip.text = element_text(size = 10)) +
        theme_classic() +
        NULL
}

# comp_plot(20)
# comp_plot(20) %>%
#     ggsave(filename = "~/Desktop/comp_20.pdf", height = 6, width = 6)
comp_plot(27)
# %>%
#     ggsave(filename = "~/Desktop/comp_27.pdf", height = 6, width = 6)


bind_rows(old_comp_df %>%
              filter(treatment == "27 C") %>%
              select(date, green, red, rep, total_red, total_green, treatment),
          comp_df %>%
              filter(treatment == "27 C") %>%
              select(date, green, red, rep, total_red, total_green, treatment) %>%
              mutate(rep = rep + 100)) %>%
    mutate(combo = paste0(green, red) %>%
               factor() %>% as.integer() %>% factor()) %>%
    comp_plot(.temp = 27) %>%
    ggsave(filename = "~/Desktop/comp_27.pdf", height = 5, width = 6)



# comp_plot(27) %>%
#     ggsave(filename = "~/Desktop/comp_27.pdf", height = 5, width = 6)







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
           total = alate_red + alate_green + apterous_red + apterous_green,
           apt_total = apterous_red + apterous_green) %>%
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



alate_ts_p <- alate_df %>%
    mutate(line = factor(line, levels = lines_)) %>%
    # ggplot(aes(apterous, logit(pr_alate))) +
    # ggplot(aes(date, asinsqrt(pr_alate))) +
    ggplot(aes(date, alate)) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
    # geom_line(aes(group = factor(rep)), color = "gray30", alpha = 0.5) +
    geom_point(aes(color = line), shape = 16, size = 2, alpha = 0.75) +
    facet_wrap( ~ line, nrow = 2) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    scale_y_continuous("Number of alates") +
    scale_x_continuous("Day", breaks = seq(0, 30, 10)) +
    # scale_x_continuous("Total apterous aphids", breaks = c(0, 250, 500)) +
    theme(strip.text = element_text(size = 10)) +
    theme_black()

ggsave("~/Desktop/alate_ts.pdf", alate_ts_p, width = 6, height = 4)





library(lme4)

z_trans <- function(x) (x - mean(x)) / sd(x)
# z_inv_trans <- function(x, m, s) (x * s) + m


adf <- alate_df %>%
    group_by(rep, combo, line) %>%
    summarize(alate = sum(alate),
              apterous = sum(apterous),
              total = sum(total),
              days = n()) %>%
    ungroup() %>%
    # mutate(log_apterous = log(apterous),
    #        log_total = log(total),
    #        log_days = log(days))
    mutate(z_apterous = z_trans(apterous),
           z_total = z_trans(total),
           z_days = z_trans(days))


# Formulas for models:
forms <- list(alate ~ z_total + (1 + z_total | line),
              alate ~ z_total + (1 | line),
              alate ~ z_total,
              alate ~ z_days + (1 + z_days | line),
              alate ~ z_days + (1 | line),
              alate ~ z_days,
              alate ~ z_apterous + (1 + z_apterous | line),
              alate ~ z_apterous + (1 | line),
              alate ~ z_apterous)


fit_mod <- function(.f) {
    if (grepl("\\|", deparse(.f))) {
        m <- glmer(.f, adf, family = poisson)
    } else {
        m <- glm(.f, adf, family = poisson)
    }
    return(as.numeric(logLik(m)))
}

# Log-likelihoods
LLs <- sapply(forms, fit_mod)

# Which one's best?
form <- forms[[which(LLs == max(LLs))]]


mod <- glmer(form, adf, family = poisson)

plot(mod)


nd <- crossing(z_total = seq((0 - mean(adf$total)) / sd(adf$total),
                             (1000 - mean(adf$total)) / sd(adf$total),
                             length.out = 100),
               line = unique(adf$line))
nd$alate <- predict(mod, newdata = nd)



nd %>%
    ggplot() +
    geom_line(aes(z_total, alate, group = line, color = line)) +
    scale_y_continuous("Total alates produced",
                       breaks = log(2 * 4^(-1:2)),
                       labels = 2 * 4^(-1:2)) +
    scale_x_continuous("Total aphids on plant", limits = c(NA, max(nd$z_total) * 1.25),
                       breaks = (seq(250, 1000, 250) - mean(adf$total)) / sd(adf$total),
                       labels = seq(250, 1000, 250)) +
    theme_classic() +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    geom_text(data = nd %>%
                  filter(z_total == max(nd$z_total)) %>%
                  arrange(desc(alate)) %>%
                  mutate(alate = alate + c(0, 0.1, rep(0, 6))),
              aes(z_total + 0.05, alate, label = line, color = line),
              hjust = 0) +
    NULL





adf2 <- alate_df %>%
    group_by(rep, combo, line) %>%
    summarize(alate = max(alate),
              total = max(total)) %>%
    ungroup() %>%
    mutate(z_total = z_trans(total),
           log_total = log(total),
           log_alate = log(alate))


mod2 <- glmer(alate ~ log_total + (1 + log_total | line), adf2, family = poisson)




nd2 <- crossing(log_total = log(seq(25, 1000, length.out = 100)),
               line = unique(adf$line))
nd2$alate <- exp(predict(mod2, newdata = nd2))


nd2 %>%
    ggplot(aes(log_total, alate, color = line)) +
    geom_point(data = adf2) +
    geom_line(size = 1) +
    scale_y_continuous("Max alates produced") +
    scale_x_continuous("Max aphids on plant",
                       breaks = log(10^(0:3)),
                       labels = 10^(0:3)) +
    # geom_text(data = nd2 %>%
    #               filter(log_total == max(nd2$log_total)) %>%
    #               arrange(desc(alate)) %>%
    #               mutate(alate = alate + c(0, 0, 0.1, 0, 0, -0.1, 0, 0)),
    #           aes(log_total + 0.1, alate, label = line, color = line),
    #           hjust = 0) +
    facet_wrap(~ line, nrow = 2) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    theme_classic() +
    NULL




ranefs <- function(.m) {
    R <- as.matrix(ranef(.m)[["line"]])
    F_ <- matrix(fixef(.m), nrow(R), 2, byrow = TRUE)
    as.numeric(R + F_)
}


# # Takes ~ 13 min
# rboots <- bootMer(mod2, ranefs, 2000, use.u = FALSE, type = "parametric",
#                   seed = 89165456, ncpus = 4)
# saveRDS(rboots, "rboots.rds")
#
# rboots <- readRDS("rboots.rds")
# coefs <- rboots$t

# # Takes ~ 8.5 min
# rboots <- bootMer(mod2, ranefs, 2000, use.u = TRUE, type = "parametric",
#                   seed = 89165456, ncpus = 4)
# saveRDS(rboots, "rboots_useuTRUE.rds")

rboots <- readRDS("rboots_useuTRUE.rds")


# rboots %>%
#     .[["t"]] %>%
#     {tibble(low = apply(., 2, quantile, probs = 0.25),
#        med = apply(., 2, median),
#        high = apply(., 2, quantile, probs = 0.75),
#        est = ranefs(mod2),
#        par = rep(c("intercept", "slope"), each = 8),
#        line = rep(rownames(ranef(mod2)[["line"]]), 2))} %>%
#     ggplot(aes(line)) +
#     ggtitle("use.u = FALSE") +
#     geom_pointrange(aes(y = est, ymin = low, ymax = high, color = line)) +
#     facet_wrap(~ par, scales = "free_x") +
#     scale_color_brewer(palette = "Dark2", guide = FALSE) +
#     theme_classic() +
#     theme(axis.text.x = element_text(color = "black"),
#           axis.title.x = element_blank()) +
#     coord_flip() +
#     NULL


boot_ests <- rboots %>%
    .[["t"]] %>%
    {tibble(low = apply(., 2, quantile, probs = 0.16),
       med = apply(., 2, median),  # median est. from bootstrapping
       high = apply(., 2, quantile, probs = 0.84),
       est = as.numeric(unlist(coef(mod2)[["line"]])), # est. from model
       par = rep(c("intercept", "slope"), each = 8),
       line = rep(rownames(ranef(mod2)[["line"]]), 2))}



est_p <- boot_ests %>%
    mutate(line = factor(line, levels = rev(lines_))) %>%
    ggplot(aes(line)) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
    geom_pointrange(aes(y = med, ymin = low, ymax = high, color = line)) +
    # geom_point(aes(y = est, color = line), shape = 18, size = 3) +
    facet_wrap(~ par, scales = "free_x") +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    theme_black() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 12)) +
    coord_flip() +
    NULL


ggsave("~/Desktop/estimates.pdf", est_p, width = 6, height = 4)







#'
#' Similar to above, but using median bootstrapped estimates instead of estimates
#' from model
#'

pred_p <- boot_ests %>%
    select(line, med, par) %>%
    spread(par, med) %>%
    crossing(log_total = log(seq(25, 1000, length.out = 100))) %>%
    mutate(alate = exp(intercept + slope * log_total),
           line = factor(line, levels = lines_)) %>%
    ggplot(aes(log_total, alate, color = line)) +
    geom_point(data = mutate(adf2, line = factor(line, levels = lines_)), shape = 1) +
    geom_line(size = 1) +
    scale_y_continuous("Max alates produced") +
    scale_x_continuous("Max aphids on plant",
                       breaks = log(5 * 10^(1:2)),
                       labels = 5 * 10^(1:2)) +
    facet_wrap(~ line, nrow = 2) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    theme_black() +
    # theme(strip.background = element_blank(),
    #       strip.text = element_text(size = 11)) +
    NULL


ggsave("~/Desktop/predictions.pdf", pred_p, width = 6, height = 4)
