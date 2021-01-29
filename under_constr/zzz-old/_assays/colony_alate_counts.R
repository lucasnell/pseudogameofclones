


#'
#' This file looks for differences between lines in numbers of alates
#' present in colonies.
#' Data were collected in summer 2020.
#'


suppressPackageStartupMessages({
    library(tidyverse)
    library(lme4)
})

inv_logit <- clonewars::inv_logit
logit <- clonewars::logit

if (file.exists(".Rprofile")) source(".Rprofile")
# ggplot theme:
theme_set(theme_classic() %+replace%
              theme(strip.background = element_blank(),
                    strip.text = element_text(size = 11),
                    legend.background = element_blank(),
                    plot.title = element_text(size = 14, hjust = 0)))


col_counts <- paste0("~/Box Sync/Aphid Project/Colony alate counts/",
                     "Colony alate counts DATA 25Oct2020.csv") %>%
    read_csv(col_types = cols()) %>%
    mutate(sample_date = as.Date(sample_date, "%m/%d/%Y"),
           date = as.Date(date, "%m/%d/%Y"),
           age = as.numeric(difftime(sample_date, date, units = "days"))) %>%
    filter(age >= 10, apterous + alate > 20) %>%
    # filter(apterous + alate > 0) %>%
    mutate(z_age = (age - mean(age)) / sd(age),
           line2 = factor(line),
           line = factor(ifelse(grepl("^Sample", line), "mixed", line)),
           total = apterous + alate,
           z_total = (total - mean(total)) / sd(total),
           z_date = as.numeric(date),
           z_date = (z_date - mean(z_date)) / sd(z_date),
           pot = as.factor(paste0(line, "_", date)),
           obs = 1:n(),
           thrips = date < as.Date("2020-08-14"))




col_counts %>%
    mutate(p = alate / total) %>%
    ggplot(aes(age, p)) +
    geom_point(aes(color = line)) +
    facet_wrap(~ line) +
    scale_color_manual(guide = FALSE,
                       values = viridisLite::viridis(10) %>%
                           .[do.call(c, map(1:5, ~ .x + c(0,5)))]) +
    scale_y_continuous("Proportion alates") +
    xlab("Pot age")

# # Looks like there's an effect of date (due to thrips early on):
# col_counts %>%
#     ggplot(aes(date, alate / total, color = line)) +
#     # geom_vline(xintercept = as.Date("2020-08-14"), linetype = 2) +
#     geom_point() +
#     stat_smooth(method = "lm", formula = y ~ x, se = FALSE) +
#     facet_wrap(~ line) +
#     theme(legend.position = "none")





# # # Convergence issues:
# # z <- glmer(cbind(alate, apterous) ~ line + age + (1|pot) + (1|obs),
# #            family = binomial, data = col_counts)
# # z0 <- glmer(cbind(alate, apterous) ~ age + (1|pot) + (1|obs),
# #             family = binomial, data = col_counts)
# # summary(z)
# # anova(z,z0)
#
# # # Convergence issues:
# # z <- glmer(cbind(alate, apterous) ~ line + total + (1|pot) + (1|obs),
# #            family = binomial, data = col_counts)
# # z0 <- glmer(cbind(alate, apterous) ~ total + (1|pot) + (1|obs),
# #             family = binomial, data = col_counts)
# # summary(z)
# # anova(z,z0)
#
# # # Convergence issues:
# # z <- glmer(cbind(alate, apterous) ~ total + (1|obs) + (1|line),
# #            family = binomial, data = col_counts)
# # z0 <- glmer(cbind(alate, apterous) ~ total + (1|obs),
# #             family = binomial, data = col_counts)
# # summary(z)
# # anova(z,z0)


z <- glm(cbind(alate, apterous) ~ line + age, data = col_counts, family = quasibinomial)
z0 <- glm(cbind(alate, apterous) ~ age, data = col_counts, family = quasibinomial)
summary(z)





# "z_date", "(1 | obs)", "z_age", "z_total", "(1 | line)"



# Line is non-significant here
# I can't include it as a random effect bc I get `boundary (singular) fit`

al_mod <- glmer(cbind(alate, apterous) ~ (1 | obs) + z_age + z_date + (1 | line),
                control = glmerControl(optimizer = "bobyqa",
                                       optCtrl = list(maxfun = 100e3)),
                family = binomial, data = col_counts)

al_mod0 <- glmer(cbind(alate, apterous) ~ (1 | obs) + z_age + z_date,
                 control = glmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 100e3)),
                 family = binomial, data = col_counts)
anova(al_mod, al_mod0)


# To bootstrap predicted Pr(alate):
pred_data <- map_dfr(levels(col_counts$line),
                     function(.line) {
                         .df <- filter(col_counts, line == .line)
                         tibble(z_age = seq(min(.df$z_age), max(.df$z_age),
                                            length.out = 11),
                                line = .line)
                     }) %>%
    mutate(z_date = 0, obs = 1)

boot_fun <- function(x) {
    cf <- coef(x)[["line"]]
    b0 <- cf[,"(Intercept)"]
    b1 <- cf[,"z_age"]
    z <- c(b0, b1)
    .lines <- rownames(ranef(x)[["line"]])
    names(z) <- c(paste0("int_", .lines), paste0("z_age_", .lines))
    return(z)
}


# # Takes ~1.3 min on my machine using 3 threads
# al_mod_boot <- bootMer(al_mod, boot_fun, 2000, seed = 96123392,
#                        parallel = "multicore", ncpus = 3,
#                        re.form = ~ (1 | line))
# saveRDS(al_mod_boot, "al_mod_boot.rds")


al_mod_boot <- readRDS("al_mod_boot.rds")


# fortify.merMod(al_mod) %>%
#     ggplot(aes(.fitted, .resid)) +
#     geom_point()




al_boot_df <- map_dfr(1:nrow(al_mod_boot$t),
                      function(i) {
                          .row <- al_mod_boot$t[i,]
                          mutate(pred_data,
                                 rep = i,
                                 p = .row[paste0("z_age_", line)] * z_age +
                                     .row[paste0("int_", line)])
                      }) %>%
    mutate(rep = factor(rep, levels = 1:2000),
           line_color = ifelse(line %in% c("R10", "WIA-5D", "WI-L4", "WI-L40",
                                           "WI-2016-746"),
                               "red", "green") %>%
               factor(levels = c("red", "green")),
           age = z_age * sd(col_counts$age) + mean(col_counts$age))



pred_data %>%
    mutate(p = predict(al_mod, type = "response", newdata = pred_data,
                       re.form = ~ (1 | line)),
           age = z_age * sd(col_counts$age) + mean(col_counts$age)) %>%
    ggplot(aes(age, p)) +
    geom_line(data = al_boot_df %>% mutate(p = inv_logit(p)),
              alpha = 0.02, aes(group = rep, color = line_color)) +
    geom_line(size = 1, color = "black") +
    geom_point(data = col_counts %>%
                   mutate(p = alate / total)) +
    facet_wrap(~ line) +
    scale_color_manual(values = c(red = "#e41a1c", green = "#4daf4a"),
                       guide = FALSE) +
    scale_y_continuous("Proportion alates") +
    xlab("Pot age")

# Logit-transformed
pred_data %>%
    mutate(p = predict(al_mod, type = "link", newdata = pred_data,
                       re.form = ~ (1 | line)),
           age = z_age * sd(col_counts$age) + mean(col_counts$age)) %>%
    ggplot(aes(age, p)) +
    geom_line(data = al_boot_df,
              alpha = 0.02, aes(group = rep, color = line_color)) +
    geom_line(size = 1, color = "black") +
    geom_point(data = col_counts %>%
                   mutate(p = logit(alate / total))) +
    facet_wrap(~ line) +
    scale_color_manual(values = c(red = "#e41a1c", green = "#4daf4a"),
                       guide = FALSE) +
    scale_y_continuous("Proportion alates") +
    xlab("Pot age")


