

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

asn_trans <- function(x) 2 * asin(sqrt(x))


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
    mutate(id = paste(green_line, red_line, rep, sep = "__")) %>%
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
    select(id, date, line, apterous, alates, tN, N, lag_N, pcg, lag_pcg) %>%
    # --------
    # Remove short time series:
    # --------
    group_by(id) %>%
    # filter(date <= min(date[N == max(N)])) %>%  # <-- to remove dates after peak N
    mutate(date = difftime(date, min(date), units = "days") %>% as.integer()) %>%
    mutate(n_obs = date %>% unique() %>% length()) %>%
    ungroup() %>%
    filter(n_obs > 2) %>%
    filter(N > 0) %>%
    select(-n_obs) %>%
    arrange(id, line, date)




# Impute N for all dates from the start of an assay rep
impute_all_days <- function(N, date) {


    N <- c(N, rep(NA, sum(!(0:max(date) %in% date))))
    date <- c(date, (0:max(date))[!(0:max(date) %in% date)])

    N <- N[order(date)]
    date <- date[order(date)]

    # Fixed version of `zoo::na.StructTS`
    # with an extra checks for length 1 and for all-identical non-NAs
    zoo__na.StructTS <- function (object) {

        na.rm = FALSE; maxgap = Inf

        if (length(object) == 1) return(object)

        na.StructTS.0 <- function(y) {
            yf <- y
            isna <- is.na(y)
            if (length(unique(y[!isna])) == 1) {
                yf[isna] <- yf[!isna][1]
                return(yf)
            }
            yf[isna] <- rowSums(tsSmooth(StructTS(y))[, -2, drop=FALSE])[isna]
            zoo:::.fill_short_gaps(y, yf, maxgap = maxgap)
        }
        object[] <- if (length(dim(object)) == 0)
            na.StructTS.0(object)
        else apply(object, 2, na.StructTS.0)
        if (na.rm)
            na.trim(object, is.na = "all")
        else object
    }

    part_a <- N[1:(which(N == max(N, na.rm=TRUE))[1])] %>%
        as.ts() %>%
        zoo__na.StructTS() %>%
        as.numeric()

    part_b <- N[(which(N == max(N, na.rm=TRUE))[1]):length(N)] %>%
        as.ts() %>%
        zoo__na.StructTS() %>%
        as.numeric()

    tibble(N = c(part_a, part_b[-1]), date = date)

}




# Find N on the day that new adults were born
# Default `delay` is that between birth and adulthood (i.e., when alates are observed)
born_N <- function(N, date,
                   .delay = sum(head(dev_times$instar_days$lowT, -1)) - 1) {

    N_all <- impute_all_days(N, date)

    delay_lookup <- function(.d) {
        if (.d < .delay) return(NA_real_)
        N_all$N[N_all$date == (.d - .delay)]
    }

    map_dbl(date, delay_lookup)

}


# Including a bunch of variables for models:

alate_mod_df <- alate_df %>%
    group_by(id, line) %>%
    mutate(b_N = born_N(N, date),
           pb_N = b_N / (max(N) + 1),
           pb_nc_N = pb_N - 0.5,
           new_alates = alates - lag(alates),
           new_alates = ifelse(new_alates < 0, 0, new_alates),
           new_apterous = apterous - lag(apterous),
           new_apterous = ifelse(new_apterous < 0, 0, new_apterous),
           new_all = new_apterous + new_alates,
           lag_apterous = lag(apterous),
           is_alate = ifelse(new_alates == 0, 0, 1)) %>%
    ungroup() %>%
    # Because I'm only interested in the ratio of new_alates to new_apterous, when they're
    # both zero, this is uninformative:
    filter(!(new_alates == 0 & new_apterous == 0)) %>%
    filter(!is.na(b_N)) %>%
    mutate(line = factor(line, levels = alate_df$line %>% unique() %>% sort()),
           plant = factor(id) %>% as.integer() %>% factor())




# b_N seems to work okay:
alate_mod_df %>%
    ggplot(aes(b_N, asn_trans(new_alates / new_all))) +
    geom_point(aes(color = line), shape = 1, alpha = 0.5) +
    stat_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    facet_wrap(~ line, nrow = 2) +
    NULL


# Full model (`~ b_N + (1 | line) + (b_N | id)`): 345.7641

forms <- list(cbind(new_alates, new_apterous) ~ b_N + (b_N | line) + (b_N | id),
              cbind(new_alates, new_apterous) ~ b_N + (1 | line) + (b_N | id),
              cbind(new_alates, new_apterous) ~ b_N + (b_N | line) + (1 | id),
              cbind(new_alates, new_apterous) ~ b_N + (1 | line) + (1 | id),
              cbind(new_alates, new_apterous) ~ b_N + (b_N | line),
              cbind(new_alates, new_apterous) ~ b_N + (1 | line),
              cbind(new_alates, new_apterous) ~ b_N + (b_N | id),
              cbind(new_alates, new_apterous) ~ b_N + (1 | id))
mods <- map(forms, ~ glmer(.x, data = alate_mod_df, family = binomial))

# These didn't fit properly, so we remove them:
forms[map_lgl(mods, ~ isSingular(.x))]
mods <- mods[map_lgl(mods, ~ !isSingular(.x))]

# AICs for all non-singular model fits:
map_dbl(mods, ~ AIC(.x))

# This one's best:
mods[[which(map_dbl(mods, ~ AIC(.x)) == min(map_dbl(mods, ~ AIC(.x))))]]


# (If the anything changes, make sure the one below is the same as above)
alate_mod <- glmer(cbind(new_alates, new_apterous) ~ b_N + (1 | line) +
                       (b_N | plant),
                   data = alate_mod_df, family = binomial)

glmer(cbind(new_alates, new_apterous) ~ (1 | line) + (b_N | plant),
      data = alate_mod_df, family = binomial)


fixef(alate_mod)
ranef(alate_mod)



predict(alate_mod, newdata = tibble(line = "WI-L4", b_N = 10, plant = 1),
        type = "response")

b0 <- fixef(alate_mod)[["(Intercept)"]] +
    ranef(alate_mod)[["line"]]["WI-L4","(Intercept)"] +
    ranef(alate_mod)[["plant"]]["1","(Intercept)"]
b1 <- fixef(alate_mod)[["b_N"]] +
    ranef(alate_mod)[["plant"]]["1","b_N"]
inv_logit(b0 + 10 * b1)



# Seems to work okay:
alate_mod_df %>%
    mutate(obs = new_alates) %>%
    select(plant, line, new_all, b_N, obs) %>%
    {mutate(., mod = new_all * predict(alate_mod, newdata = ., type = "response"))} %>%
    group_by(line) %>%
    summarize(obs = sum(obs), mod = sum(mod))


alate_prob_coefs <- alate_mod %>%
    ranef() %>%
    .[["line"]] %>%
    rownames_to_column("line") %>%
    rename(b0 = `(Intercept)`) %>%
    mutate(b0 = b0 + fixef(alate_mod)[["(Intercept)"]],
           b1 = fixef(alate_mod)[["b_N"]])

alate_prob_coefs


# How to use this (we'll ignore plant variation in simulations, since
# there's already variation in the K parameter):
.line <- "WI-L4"
.b_N <- 80
b0 <- alate_prob_coefs %>% filter(line == .line) %>% .[["b0"]]
b1 <- alate_prob_coefs %>% filter(line == .line) %>% .[["b1"]]
inv_logit(b0 + b1 * .b_N)
predict(alate_mod, newdata = tibble(line = "WI-L4", b_N = .b_N),
        re.form = ~ (1 | line), type = "response")



