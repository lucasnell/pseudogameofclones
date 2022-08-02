
#'
#' This file looks for affects of clonal line and density on the
#' proportion of aphids that are alates, using observations on colonies
#' we maintain in the lab.
#' Data were collected in summer 2020.
#'


library(lme4)
library(tidyverse)
library(gameofclones)


col_counts <- paste0("~/Box Sync/eco-evo_experiments/prelim_assays/alates/",
                     "Colony alate counts DATA 25Oct2020.csv") %>%
    read_csv(col_types = cols()) %>%
    mutate(sample_date = as.Date(sample_date, "%m/%d/%Y"),
           date = as.Date(date, "%m/%d/%Y"),
           age = as.numeric(difftime(sample_date, date, units = "days"))) %>%
    # These were mixed-clone samples and not useful for single-line estimates:
    filter(!grepl("^Sample", line)) %>%
    # Filter out very new colonies bc we typically start with apterous aphids,
    # which would bias these estimates:
    filter(age >= 10, apterous + alate > 20) %>%
    mutate(total = apterous + alate,
           pot = as.factor(paste0(line, "_", date)),
           obs = 1:n())


col_counts %>%
    mutate(p = alate / total) %>%
    ggplot(aes(total, p)) +
    geom_point(aes(color = line)) +
    scale_color_manual(guide = "none",
                       values = viridisLite::viridis(10) %>%
                           .[do.call(c, map(1:5, ~ .x + c(0,5)))]) +
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
                       seed = 96123392, parallel = "multicore",
                       ncpus = max(1, parallel::detectCores()-2))

# 95% CI and median estimate from parametric bootstrapping:
quantile(total_mod_boot$t, c(0.025, 0.5, 0.975))



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
                         ncpus = max(1, parallel::detectCores()-2),
                         re.form = ~ (1 | line))

# 95% CI and median estimates for each line from parametric bootstrapping:
apply(line_mod_boot$t, 2, quantile, probs = c(0.025, 0.5, 0.975))

# To visualize this:

apply(line_mod_boot$t, 2, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    set_names(c("line", "low", "med", "high")) %>%
    as_tibble() %>%
    ggplot(aes(line)) +
    geom_point(aes(y = med)) +
    geom_linerange(aes(ymin = low, ymax = high))

# If you want to see the entire distribution of the bootstrap resamples:
line_mod_boot %>%
    as_tibble() %>%
    pivot_longer(everything(), names_to = "line", values_to = "b0") %>%
    ggplot(aes(b0)) +
    geom_vline(xintercept = median(line_mod_boot$t), linetype = 1, color = "gray70") +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_density() +
    facet_wrap(~ line, ncol = 1) +
    theme(strip.text = element_text(size = 8))


# And the 95% CI of the difference between the two lines of interest:
quantile(line_mod_boot$t[,"WIA-5D"] - line_mod_boot$t[,"UT3"], c(0.025, 0.5, 0.975))





# --------------------------------------------------------*
# What number should we use then?
# --------------------------------------------------------*

simple_mod <- glmer(cbind(alate, apterous) ~ (1 | obs) + (1 | pot),
                  family = binomial, data = col_counts)
summary(simple_mod)

inv_logit(fixef(simple_mod)[["(Intercept)"]])



