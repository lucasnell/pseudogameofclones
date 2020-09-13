


#'
#' This file looks for differences between lines in numbers of alates
#' present in colonies.
#' Data were collected in summer 2020.
#'


suppressPackageStartupMessages({
    library(tidyverse)
    library(lme4)
    library(clonewars)
})

source(".Rprofile")


col_counts <- read_csv("~/Desktop/Colony alate counts DATA 08Sep2020.csv") %>%
    mutate(sample_date = as.Date(sample_date, "%m/%d/%Y"),
           date = as.Date(date, "%m/%d/%Y"),
           age = difftime(sample_date, date, units = "days"),
           id = interaction(line, date, sep = "__") %>% factor(),
           line = factor(line),
           date_fct = factor(date),
           date_int = difftime(date, mean(date), units = "days"),
           total = apterous + alate) %>%
    mutate(across(c("apterous", "alate", "age", "total", "date_int"), as.integer)) %>%
    filter(total > 0)


col_counts %>%
    ggplot(aes(date, alate / total)) +
    geom_point(aes(color = line))



al_mod <- glmer(cbind(alate, apterous) ~ date_int + (total | line),
                family = binomial, data = col_counts)
al_mod %>% summary()

pred_data <- crossing(total = seq(1, 80, length.out = 11),
                      line = levels(col_counts$line),
                      date_int = 0)

boot_fun <- function(x) {
    # z <- fixef(x)[["(Intercept)"]]
    # c(z + ranef(x)[["line"]][["(Intercept)"]],
    #   ranef(x)[["line"]][["total"]])
    predict(x, newdata = pred_data, type = "response")
}





# # Takes ~5 min
# al_mod_boot <- bootMer(al_mod, boot_fun, 2000, seed = 96123392,
#                        parallel = "multicore", ncpus = 3, use.u = TRUE)
# saveRDS(al_mod_boot, "al_mod_boot.rds")

al_mod_boot <- readRDS("al_mod_boot.rds")

al_boot_df <- al_mod_boot$t %>%
    {map_dfr(1:nrow(.),
             function(i) {
                 mutate(pred_data, rep = i, p = .[i,])
             })} %>%
    select(rep, total, line, p) %>%
    mutate(rep = factor(rep, levels = 1:2000))




#
# al_boot_ci_df <- al_boot_df %>%
#     group_by(total, line) %>%
#     summarize(lo = quantile(p, 0.05),
#               mid = median(p),
#               hi = quantile(p, 0.95)) %>%
#     ungroup()


# red: #e41a1c
# green: #4daf4a


pred_data %>%
    filter(!line %in% c("WI-48", "WI-2016-746")) %>%
    mutate(p = predict(al_mod, type = "response", newdata = pred_data %>%
                           filter(!line %in% c("WI-48", "WI-2016-746")))) %>%
    ggplot(aes(total, p)) +
    geom_line(data = al_boot_df %>% filter(!line %in% c("WI-48", "WI-2016-746")),
              alpha = 0.02, aes(group = rep), color = "gray60") +
    geom_line(size = 1, color = "dodgerblue") +
    geom_point(data = col_counts %>%
                   filter(!line %in% c("WI-48", "WI-2016-746")) %>%
                   mutate(p = inv_logit(logit(alate / total) - date_int *
                                            fixef(al_mod)[["date_int"]]))) +
                   # mutate(p = alate / total)) +
    facet_wrap(~ line) +
    # scale_color_viridis_d(option = "D", end = 0.9, guide = FALSE) +
    scale_y_continuous("Proportion alates") +
    xlab("Total aphids on plant")





# col_counts %>%
#     mutate(p = alate / total) %>%
#     ggplot(aes(line, p)) +
#     geom_point() +
#     geom_point(data = tibble(line = col_counts$line %>% unique() %>% sort(),
#                              date = col_counts$date[1],
#                              ) %>%
#                    mutate(p = predict(al_mod, newdata = ., type = "response")),
#                color = "red")


