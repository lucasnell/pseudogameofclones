
source("load_data.R")
library(nlme)
# suppressPackageStartupMessages({
#     library(rstan)
# })

w_ham <- c("R10", "WI-L4", "UT3", "Clover-2017-2")

growth

d <- growth %>%
    group_by(line, rep) %>%
    arrange(date) %>%
    mutate(N_t = lag(N)) %>%
    ungroup %>%
    filter(!is.na(N_t)) %>%
    mutate(X = log(N), X_t = log(N_t),
           ham = ifelse(line %in% w_ham, 1, 0)) %>%
    arrange(line, rep, date)



ricker <- deriv(~ X_t + R0 * (1 - alpha * exp(X_t)),
                c("R0", "alpha"),
                function(X_t, R0, alpha){})

# ricker <- function(N_t, R0, alpha) N_t * exp(R0 * (1 - alpha * N_t))

mod <- nlme(model = X ~ ricker(X_t, R0, alpha),
            data = d,
            fixed = c(R0 ~ ham, alpha ~ ham),
            random = R0 + alpha ~ 1 | line,
            correlation = corAR1(form = ~ date | line/rep),
            start = c(1e-2, 1e-2, 1e-2, 1e-2))
mod %>% summary

mod2 <- nlme(model = X ~ ricker(X_t, R0, alpha),
            data = d,
            fixed = c(R0 ~ ham, alpha ~ ham),
            random = R0 ~ 1 | line,
            correlation = corAR1(form = ~ date | line/rep),
            start = c(1e-2, 1e-2, 1e-2, 1e-2))


mod3 <- nlme(model = X ~ ricker(X_t, R0, alpha),
            data = d,
            fixed = c(R0 ~ ham, alpha ~ ham),
            random = alpha ~ 1 | line,
            correlation = corAR1(form = ~ date | line/rep),
            start = c(1e-2, 1e-2, 1e-2, 1e-2))


anova(mod, mod2)
anova(mod, mod3)

d %>%
    mutate(pred = as.numeric(predict(mod))) %>%
    ggplot(aes(X, pred)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(shape = 1, color = 'dodgerblue') +
    theme_classic()




d %>%
    mutate(N_pred = exp(as.numeric(predict(mod)))) %>%
    nest(date:N_pred) %>%
    {.[sample.int(nrow(.), 12), ]} %>%
    unnest() %>%
    mutate(linerep = paste(line, rep, sep = "_")) %>%
    ggplot(aes(date, N)) +
    geom_line(aes(y = N_pred)) +
    geom_point() +
    theme_classic() +
    facet_wrap(~ linerep, nrow = 3)


