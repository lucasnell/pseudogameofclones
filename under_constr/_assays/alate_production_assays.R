

#'
#' This file uses the alate data from when we were looking at competition among lines.
#' They were done summer 2019.
#'


suppressPackageStartupMessages({
    library(tidyverse)
    library(readxl)
    library(clonewars)
    library(lme4)
})


source(".Rprofile")

z_trans <- function(x) {
    .sd <- sd(x, na.rm = TRUE)
    if (.sd == 0 || is.na(.sd)) .sd <- 1
    (x - mean(x, na.rm = TRUE)) / .sd
}
# z_inv_trans <- function(x, m, s) (x * s) + m



alate_df <- paste0("~/Dropbox/Aphid Project 2017/Lauren_competition/",
                   "competition_data_entry24July19.xlsx") %>%
    read_excel() %>%
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
    filter(date <= min(date[N == max(N)])) %>%  # <-- to remove dates after peak N
    mutate(date = difftime(date, min(date), units = "days") %>% as.integer()) %>%
    mutate(n_obs = date %>% unique() %>% length()) %>%
    ungroup() %>%
    filter(n_obs > 2, N > 0) %>%
    select(-n_obs) %>%
    arrange(id, line, date)



# Including a bunch of variables for models:

alate_mod_df <- alate_df %>%
    group_by(id, line) %>%
    mutate(z_N = z_trans(N),
           p_N = N / (max(N) + 1),
           lag_z_N = lag(z_N),
           lag_p_N = lag(p_N),
           new_alates = alates - lag(alates),
           new_alates = ifelse(new_alates < 0, 0, new_alates),
           new_alate_rate = new_alates / (date - lag(date)),
           new_apterous = apterous - lag(apterous),
           new_apterous = ifelse(new_apterous < 0, 0, new_apterous),
           new_apterous_rate = new_apterous / (date - lag(date)),
           new_all = new_apterous + new_alates,
           lag_apterous = lag(apterous),
           is_alate = ifelse(new_alates == 0, 0, 1),
           time = as.numeric(date - lag(date))) %>%
    group_by(id) %>%
    mutate(z_apterous = z_trans(apterous),
           z_alate = z_trans(alates)) %>%
    ungroup() %>%
    filter(!(is.na(new_alates) | is.na(lag_apterous) | is.na(p_N))) %>%
    # Because I'm only interested in the ratio of new_alates to new_apterous, when they're
    # both zero, this is uninformative:
    filter(!(new_alates == 0 & new_apterous == 0)) %>%
    mutate(line = factor(line, levels = alate_df$line %>% unique() %>% sort()))




# It's not obvious that any of these variables actually affect anything:
alate_mod_df %>%
    ggplot(aes(lag_N, log(new_alates / new_all))) +
    geom_jitter(aes(color = line), height = 0, width = 0.25) +
    # stat_summary(fun.data = "mean_cl_boot", color = "dodgerblue") +
    facet_wrap(~ line, nrow = 2) +
    NULL



alate_mod <- glmer(cbind(new_alates, new_apterous) ~ (1 | line), data = alate_mod_df,
                    family = binomial)
fixef(alate_mod)
ranef(alate_mod)

predict(alate_mod, newdata = tibble(line = "WI-L4"), type = "response")
inv_logit(-4.601386 + 0.3964972)


# Seems to work okay:
alate_mod_df %>%
    mutate(obs = new_alates) %>%
    select(line, new_all, obs) %>%
    {mutate(., mod = new_all * predict(alate_mod, newdata = ., type = "response"))} %>%
    group_by(line) %>%
    summarize(obs = sum(obs), mod = sum(mod))

