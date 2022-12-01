
#'
#' This file looks for affects of clonal line and density on the
#' proportion of aphids that are alates, using observations on colonies
#' we maintain in the lab.
#' Data were collected in summer 2020.
#'


library(lme4)
library(tidyverse)
library(gameofclones)
library(here)

options(boot.ncpus = max(parallel::detectCores()-2L, 1L))



col_counts <- here("_results/_data/alate-counts.csv") |>
    read_csv(col_types = cols()) |>
    mutate(sample_date = as.Date(sample_date, "%m/%d/%Y"),
           date = as.Date(date, "%m/%d/%Y"),
           age = as.numeric(difftime(sample_date, date, units = "days"))) |>
    # These were mixed-clone samples and not useful for single-line estimates:
    filter(!grepl("^Sample", line)) |>
    # Filter out very new colonies bc we typically start with apterous aphids,
    # which would bias these estimates:
    filter(age >= 10, apterous + alate > 20) |>
    mutate(total = apterous + alate,
           pot = as.factor(paste0(line, "_", date)),
           obs = 1:n())


col_counts |>
    mutate(p = alate / total) |>
    ggplot(aes(total, p)) +
    geom_point(aes(color = line)) +
    scale_color_manual(guide = "none",
                       values = viridisLite::viridis(10)[
                           do.call(c, map(1:5, ~ .x + c(0,5)))]) +
    scale_y_continuous("Proportion alates") +
    xlab("Total aphids")



# --------------------------------------------------------*
# Is there an effect of the total number of aphids on alate proportion?
# --------------------------------------------------------*

total_mod <- glmer(cbind(alate, apterous) ~ total + (1 | obs) + (1 | pot),
                   family = binomial, data = col_counts)
summary(total_mod)

# Bootstrap `total` fixed effect estimate:
boot_total_fixef <- function(x) fixef(x)[["total"]]

# Takes ~14 sec on my machine using 6 threads
total_mod_boot <- bootMer(total_mod, boot_total_fixef, nsim = 2000,
                       seed = 96123392, parallel = "multicore")

# 95% CI and median estimate from parametric bootstrapping:
quantile(total_mod_boot$t, c(0.025, 0.5, 0.975))
#        2.5%         50%       97.5%
# -0.01029315  0.01634713  0.04274965


# --------------------------------------------------------*
# Is there an effect of aphid line on alate proportion?
# --------------------------------------------------------*


line_mod <- glmer(cbind(alate, apterous) ~ (1 | obs) + (1 | pot) + (1 | line),
                  family = binomial, data = col_counts)
summary(line_mod)

# Bootstrap random effect of line:
boot_line_ranef <- function(x) {
    b0 <- coef(x)[["line"]][,"(Intercept)"]
    names(b0) <- rownames(ranef(x)[["line"]])
    return(b0)
}

# Takes ~14 sec on my machine using 6 threads
line_mod_boot <- bootMer(line_mod, boot_line_ranef, nsim = 2000,
                         seed = 1784416833, parallel = "multicore",
                         re.form = ~ (1 | line))

# 95% CI and median estimates for each line from parametric bootstrapping:
apply(line_mod_boot$t, 2, quantile, probs = c(0.025, 0.5, 0.975))
#       Clover-2017-2 Clover-2017-6       R10       UT3 WI-2016-593 WI-2016-746
# 2.5%      -2.478109     -2.956589 -3.183403 -3.342695   -2.946504   -2.362257
# 50%       -1.928381     -2.289194 -2.367808 -2.574542   -2.249657   -1.713488
# 97.5%     -1.252620     -1.726391 -1.798317 -2.014520   -1.654194   -0.806024
#
#           WI-48      WI-L4    WI-L40    WIA-5D
# 2.5%  -3.122752 -2.3469877 -2.398196 -2.909406
# 50%   -2.363722 -1.7048844 -1.836794 -2.239044
# 97.5% -1.866954 -0.9159304 -1.114809 -1.678777


# To visualize this:

apply(line_mod_boot$t, 2, quantile, probs = c(0.025, 0.5, 0.975)) |>
    t() |>
    as.data.frame() |>
    rownames_to_column() |>
    set_names(c("line", "low", "med", "high")) |>
    as_tibble() |>
    ggplot(aes(line)) +
    geom_point(aes(y = med)) +
    geom_linerange(aes(ymin = low, ymax = high))

# If you want to see the entire distribution of the bootstrap resamples:
line_mod_boot |>
    as_tibble() |>
    pivot_longer(everything(), names_to = "line", values_to = "b0") |>
    ggplot(aes(b0)) +
    geom_vline(xintercept = median(line_mod_boot$t), linetype = 1, color = "gray70") +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_density() +
    facet_wrap(~ line, ncol = 1) +
    theme(strip.text = element_text(size = 8))


# And the 95% CI of the difference between the two lines of interest:
quantile(line_mod_boot$t[,"WIA-5D"] - line_mod_boot$t[,"UT3"], c(0.025, 0.5, 0.975))
#       2.5%        50%      97.5%
# -0.3421796  0.3062545  1.2050211




# --------------------------------------------------------*
# What number should we use then?
# --------------------------------------------------------*

simple_mod <- glmer(cbind(alate, apterous) ~ (1 | obs) + (1 | pot),
                  family = binomial, data = col_counts)
summary(simple_mod)

# Generalized linear mixed model fit by maximum likelihood (Laplace  Approximation)
#  [glmerMod]
#  Family: binomial  ( logit )
# Formula: cbind(alate, apterous) ~ (1 | obs) + (1 | pot)
#    Data: col_counts
#
#      AIC      BIC   logLik deviance df.resid
#    490.9    498.2   -242.4    484.9       81
#
# Scaled residuals:
#      Min       1Q   Median       3Q      Max
# -1.16916 -0.31780 -0.00212  0.18392  0.90547
#
# Random effects:
#  Groups Name        Variance Std.Dev.
#  obs    (Intercept) 0.4495   0.6704
#  pot    (Intercept) 1.7470   1.3217
# Number of obs: 84, groups:  obs, 84; pot, 81
#
# Fixed effects:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -2.2809     0.1903  -11.98   <2e-16 ***
#     ---
#     Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

inv_logit(fixef(simple_mod)[["(Intercept)"]])
# [1] 0.09271424



