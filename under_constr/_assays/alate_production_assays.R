

#'
#' This file uses the alate data from when we were looking at competition among lines.
#' They were done summer 2019.
#'

library(tidyverse)
library(readxl)


source(".Rprofile")


# Lines in order they should appear in plots:
lines_ <- c("R10", "WIA-5D", "WI-L4", "WI-L4Ø", "UT3", "WI-2016-593", "Clover-2017-2",
            "Clover-2017-6")


fn <- paste0("~/Dropbox/Aphid Project 2017/Lauren_competition/",
             "competition_data_entry24July19.xlsx")

alate_df <- read_excel(fn) %>%
    rename_all(~ tolower(.)) %>%
    filter(!is.na(green_line)) %>%
    mutate(date = as.Date(paste(year, month, day, sep = "-"), format = "%Y-%B-%d")) %>%
    select(-year, -month, -day, -observer, -comments) %>%
    mutate(red_line = ifelse(red_line == "WI-L4 Ham-", "WI-L4Ø", red_line)) %>%
    mutate_at(vars(green_line, red_line), ~ gsub(" Reg\\+| Ham\\+", "",
                                                      .x)) %>%
    mutate(id = paste(green_line, red_line, rep, sep = "__") %>%
               factor()) %>%
               # factor() %>% as.integer() %>% factor()) %>%
    select(id, date, everything(), -rep, -treatment) %>%  # they were all 27ºC
    # Impossible numbers during this date:
    filter(!(id == "WI-2016-593__WI-L4Ø__1" & date == "2019-05-29")) %>%
    arrange(id, date) %>%
    group_by(id) %>%
    mutate(N = apterous_green + apterous_red + alate_green + alate_red,
           lag_N = lag(N),
           pcg = log(N / lag_N),
           lag_pcg = log(lag(N, 1) / lag(N, 2))) %>%
    group_by(date) %>%
    mutate(tN = sum(N)) %>%
    ungroup() %>%
    gather("color", "line", green_line, red_line) %>%
    mutate(apterous = ifelse(grepl("^green", color), apterous_green, apterous_red),
           alates = ifelse(grepl("^green", color), alate_green, alate_red)) %>%
    arrange(id, date, color) %>%
    select(id, date, line, apterous, alates, tN, N, lag_N, pcg, lag_pcg)

# Remove short time series:
alate_df <- alate_df %>%
    group_by(id) %>%
    # filter(date <= min(date[N == max(N)])) %>%  # <-- to remove dates after peak N
    mutate(n_obs = n()) %>%
    ungroup() %>%
    filter(n_obs > 4, N > 0) %>%
    select(-n_obs) %>%
    arrange(id, line, date)



alate_df %>%
    ggplot(aes(lag_N, alates)) +
    geom_point(aes(color = line), na.rm = TRUE) +
    facet_wrap(~ line, nrow = 2) +
    scale_x_continuous(trans = "log1p") +
    scale_y_continuous(trans = "log1p") +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    theme_classic() +
    theme(strip.background = element_blank()) +
    # coord_cartesian(ylim = c(0, NA)) +
    NULL






library(lme4)

z_trans <- function(x) {
    .sd <- sd(x, na.rm = TRUE)
    if (.sd == 0 || is.na(.sd)) .sd <- 1
    (x - mean(x, na.rm = TRUE)) / .sd
}
# z_inv_trans <- function(x, m, s) (x * s) + m




alate_mod_df <- alate_df %>%
    group_by(id, line) %>%
    mutate(z_N = z_trans(N),
           p_N = N / (max(N) + 1),
           lag_z_N = lag(z_N),
           lag_p_N = lag(p_N),
           new_alates = alates - lag(alates),
           new_alates = ifelse(new_alates < 0, 0, new_alates),
           new_apterous = apterous - lag(apterous),
           new_apterous = ifelse(new_apterous < 0, 0, new_apterous),
           lag_apterous = lag(apterous)) %>%
    group_by(id) %>%
    mutate(z_apterous = z_trans(apterous),
           z_alate = z_trans(alates)) %>%
    ungroup() %>%
    filter(!(is.na(new_alates) | is.na(lag_apterous) | is.na(p_N))) %>%
    mutate(line = factor(line, levels = alate_df$line %>% unique() %>% sort()))

alate_mod <- glmer(cbind(new_alates, new_apterous) ~ p_N + (p_N | line),
                   alate_mod_df, family = binomial)

summary(alate_mod)
boot_fun <- function(.x) {
    .dd <- ranef(.x)[['line']]
    as.numeric(c(fixef(.x), .dd[[1]], .dd[[2]]))
}

# # Takes ~22 min
# alate_boots <- bootMer(alate_mod, boot_fun, nsim = 2000, seed = 1086406336,
#                        .progress = "txt")
#
# saveRDS(alate_boots, "under_constr/_assays/alate_mod_boots.rds")

alate_boots <- readRDS("under_constr/_assays/alate_mod_boots.rds")

cbind(
apply(alate_boots$t, 2, median),
apply(alate_boots$t, 2, quantile, probs = 0.025),
apply(alate_boots$t, 2, quantile, probs = 0.975)
)




alate_mod <- glmer(cbind(new_alates, new_apterous) ~ lag_p_N + (1 + lag_p_N | line),
                    alate_mod_df, family = binomial)


alate_mod_df %>%
    ggplot(aes(lag_z_N, new_alates / (new_alates + new_apterous))) +
    geom_point(aes(color = line), na.rm = TRUE) +
    facet_wrap(~ line, nrow = 2) +
    stat_smooth(method = "loess", se = FALSE, color = "black") +
    # scale_x_continuous(trans = "logit") +
    # scale_y_continuous(trans = "asn") +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    theme_classic() +
    theme(strip.background = element_blank()) +
    # coord_cartesian(ylim = c(0, NA)) +
    NULL

alate_mod2 <- lm(log1p(alates) ~ lag_apterous * p_N, data = alate_mod_df)

summary(alate_mod2)



nd <- crossing(lag_p_N = seq(0, 1, length.out = 100),
               line = alate_mod_df$line %>% unique() %>% sort())


nd %>%
    mutate(y = predict(alate_mod2, newdata = nd)) %>%
    ggplot(aes(x = lag_p_N)) +
    geom_line(aes(y = gtools::inv.logit(y))) +
    stat_smooth(data = alate_mod_df, na.rm = TRUE,
                aes(y = new_alates / (new_alates + new_apterous)),
                method = "loess", se = FALSE) +
    geom_point(data = alate_mod_df, na.rm = TRUE,
               aes(y = new_alates / (new_alates + new_apterous), color = line)) +
    facet_wrap(~ line) +
    scale_color_brewer(palette = "Dark2") +
    theme_classic()





