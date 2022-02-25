
#'
#' This file looks for differences between aphid clonal lines named
#' "WIA-5D" and "UT3" that are ostensibly susceptible and resistant,
#' respectively.
#' We're assessing whether WIA-5D is a better competitor (in the absence
#' of parasitoids) and whether UT3 is more resistant to parasitoids.
#' Data were collected in fall 2020.
#'


library(grid)
library(tidyverse)
library(readxl)
library(patchwork)
library(clonewars)
library(lme4)
library(nlme)

assay_excel_file <- paste0("~/Box Sync/eco-evo_experiments/prelim_assays/",
                           "eco-evo-prelims.xlsx")


# ========================================================================*
# ========================================================================*

# Plots ====

# ========================================================================*
# ========================================================================*


# -----------------------------------------------`
# __pop. growth ----
# -----------------------------------------------`


pop_df <- read_excel(assay_excel_file,
                     sheet = "population-growth") %>%
    mutate(date = as.Date(paste(year, month, day, sep = "-")),
           rep = factor(rep)) %>%
    group_by(rep, line) %>%
    mutate(date = difftime(date, min(date), units = "days") %>% as.integer(),
           line = factor(line, levels = c("UT3", "WIA-5D"),
                         labels = c("resistant", "susceptible"))) %>%
    select(rep, line, date, num) %>%
    group_by(rep, date) %>%
    mutate(extinct = factor(sum(num) == 0)) %>%
    # Remove once they start declining bc this is due to plant death:
    group_by(rep) %>%
    filter(date <= date[which(num == max(num))]) %>%
    ungroup() %>%
    arrange(rep, date, line)




pop_p <- pop_df %>%
    mutate(zero = factor(num == 0)) %>%
    ggplot(aes(date, num, color = line)) +
    geom_line() +
    geom_point(aes(shape = zero), size = 2) +
    facet_wrap(~ rep, nrow = 2) +
    scale_color_manual(values = c("chartreuse3", "firebrick")) +
    scale_shape_manual(values = c(19, 4), guide = "none") +
    scale_y_continuous("Number of aphids") +
    scale_x_continuous("Days after start") +
    theme(legend.title = element_blank(),
          strip.text = element_blank()) +
    NULL

pop_p

# save_plot("_results/plots/assays_competition.pdf", pop_p, 5, 3)





# ========================================================================*
# ========================================================================*

# -----------------------------------------------`
# __resistance ----
# -----------------------------------------------`

# ========================================================================*
# ========================================================================*

#'
#' In the first set of wasp-resistance assays, we added 10 juveniles
#' of each aphid line into a deli container with two fava bean leaves.
#' We then added 3 female wasps to the containers for about 2 hours.
#' These wasps were mated but had not been exposed to aphids because
#' successful attacks on aphids can change wasps’ future aphid-color
#' preferences (Langley et al. 2006).
#' After wasp exposures, we put the aphids into petri dishes separated by
#' clonal line, each dish containing two fava bean leaves with bases
#' inserted into agar gel.
#' We kept these at 20ºC for 10 days, after which we counted the number
#' of mummies and living aphids.
#'
#' In the second set, we first exposed 10 juvenile aphids of a single line to
#' wasps for about 3 hours, then exposed 10 juveniles of the other line to
#' the same wasps for the same duration.
#' We did this for both aphid lines being the first to be exposed.
#'

wasp_df <- bind_rows(
    read_excel(assay_excel_file, sheet = "wasp-resistance_1") %>%
        select(wasp_group, line, starts_with(c("juv","adult-", "mumm"))) %>%
        mutate(set = 1L),
    read_excel(assay_excel_file, sheet = "wasp-resistance_2") %>%
        select(wasp_group, round, line,
               starts_with(c("juv", "adult-", "mumm"))) %>%
        mutate(set = 2L)) %>%
    select(set, round, everything()) %>%
    rename_with(function(x) gsub("-", "_", x)) %>%
    mutate(survived = pmin(juv_assayed - mummy_end, adult_end),
           p_survived = survived / juv_assayed,
           non_survived = juv_assayed - survived,
           p_mummy = mummy_end / juv_assayed,
           non_mummy = juv_assayed - mummy_end,
           line = factor(line, levels = c("UT3", "WIA-5D"),
                         labels = c("resistant", "susceptible")),
           set = factor(set, levels = 1:2,
                        labels = c("set 1", "set 2")),
           across(c(round, wasp_group), factor),
           id = case_when(set == 1 ~ 0L,
                          round == 1 ~ 1L,
                          TRUE ~ 2L) %>%
               factor(),
           obs = 1:n())




mummy_p <- wasp_df %>%
    ggplot(aes(line, mummy_end, color = line)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_jitter(height = 0, width = 0.25, shape = 1) +
    stat_summary(fun.data = "mean_cl_boot", color = "black") +
    scale_color_manual(values = c("chartreuse3", "firebrick"), guide = "none") +
    ylab("Mummy proportion") +
    NULL

# mummy_p



juv_p <- wasp_df %>%
    ggplot(aes(line, juvenile_end, color = line)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_jitter(height = 0, width = 0.25, shape = 1) +
    stat_summary(fun.data = "mean_cl_boot", color = "black") +
    scale_color_manual(values = c("chartreuse3", "firebrick"), guide = "none") +
    scale_y_continuous("Number of juveniles") +
    NULL

# juv_p




surv_p <- wasp_df %>%
    ggplot(aes(line, p_survived, color = line)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_jitter(height = 0, width = 0.25, shape = 1) +
    stat_summary(fun.data = "mean_cl_boot", color = "black") +
    scale_color_manual(values = c("chartreuse3", "firebrick"), guide = "none") +
    ylab("Survival proportion") +
    NULL

# surv_p



wasp_p <- mummy_p + surv_p + juv_p +
    plot_annotation(tag_levels = "A") +
    plot_layout(nrow = 1) &
    theme(plot.tag = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 9, color = "black", angle = 45,
                                     vjust = 1, hjust = 1),
          axis.title.x = element_blank())

# save_plot("_results/plots/assays_wasps.pdf", wasp_p, 6.5, 3.5, seed = 45670)






# ========================================================================*
# ========================================================================*

# Statistics ====

# ========================================================================*
# ========================================================================*



#'
#' Get average difference of a given column by aphid line.
#' Used for permutation tests.
#'
get_avg_diff <- function(.data, .trans = identity, .col = "num") {
    lvls <- levels(.data$line)
    stopifnot(length(lvls) == 2)
    mean(.trans(.data[[.col]][.data$line == lvls[1]])) -
        mean(.trans(.data[[.col]][.data$line == lvls[2]]))
}


# -----------------------------------------------`
# __pop. growth ----
# -----------------------------------------------`


#'
#' Simple permutation test of sum of total aphids across each time series.
#'
sum_diffs <- pop_df %>%
    group_by(rep, line) %>%
    summarize(num = sum(num), .groups = "drop")

obs_diff <- get_avg_diff(sum_diffs)

set.seed(654067829)
diff_perms <- replicate(2000, {
    ..sum_diffs <- sum_diffs
    ..sum_diffs[["line"]] <- sample(..sum_diffs[["line"]])
    get_avg_diff(..sum_diffs)
})

mean(abs(diff_perms) >= abs(obs_diff))





pop_mod <- gls(log1p(num) ~ line, pop_df, corAR1(form = ~ 1 | rep))

#'
#' Blocked bootstrap for line estimates, where I bootstrap the number of aphids
#' within each replicate and aphid line.
#'
set.seed(24357986)
pop_mod_boots <- replicate(2000, {
    .pop_df <- pop_df
    for (r in levels(pop_df$rep)) {
        for (l in levels(pop_df$line)) {
            rl_inds <- which(pop_df$rep == r & pop_df$line == l)
            smpl_inds <- sample(rl_inds, replace = TRUE)
            .pop_df$num[rl_inds] <- pop_df$num[smpl_inds]
        }
    }
    .pop_mod <- gls(log1p(num) ~ line, .pop_df, corAR1(form = ~ 1 | rep))
    return(coef(.pop_mod)[["linesusceptible"]])
})

quantile(pop_mod_boots, c(0.025, 0.5, 0.975))





# -----------------------------------------------`
# __resistance ----
# -----------------------------------------------`


obs_mum <- get_avg_diff(wasp_df, .col = "mummy-end")
obs_juv <- get_avg_diff(wasp_df, .col = "juvenile-end")
obs_surv <- get_avg_diff(wasp_df, .col = "surv")

set.seed(678530)
mum_perms <- replicate(2000, {
    ..wasp_df <- wasp_df
    ..wasp_df[["line"]] <- sample(..wasp_df[["line"]])
    get_avg_diff(..wasp_df, .col = "mummy-end")
})
juv_perms <- replicate(2000, {
    ..wasp_df <- wasp_df
    ..wasp_df[["line"]] <- sample(..wasp_df[["line"]])
    get_avg_diff(..wasp_df, .col = "juvenile-end")
})
surv_perms <- replicate(2000, {
    ..wasp_df <- wasp_df
    ..wasp_df[["line"]] <- sample(..wasp_df[["line"]])
    get_avg_diff(..wasp_df, .col = "surv")
})

# Values of 0 actually indicate < 0.0005
mean(abs(mum_perms) >= abs(obs_mum))
mean(abs(juv_perms) >= abs(obs_juv))
mean(abs(surv_perms) >= abs(obs_surv))



#'
#' For survival and number of juveniles, there were some differences caused by
#' assay set / round.
#' I made a factor variable called `id` with three levels: (1) assay set 1,
#' (2) assay set 2 round 1, and (3) assay set 2 round 2.
#' I then did some regressions below to show that this doesn't affect our
#' conclusions based on the permutation tests.
#' If anything, the permutations appear to be conservative, so they are
#' what I present in text.
#'

surv_mod <- glmer(cbind(survived, non_survived) ~ line + id +
                      (1 | wasp_group), wasp_df, family = binomial)
surv_mod %>% summary()

set.seed(9876345)
surv_mod_boot <- bootMer(surv_mod, function(x) fixef(x)[["linesusceptible"]],
                         nsim = 2000)
# hist(surv_mod_boot$t)
quantile(surv_mod_boot$t, c(0.025, 0.5, 0.975))

juv_mod <- glmer(juvenile_end ~ line + id + (1 | wasp_group),
                 wasp_df, family = poisson)
juv_mod %>% summary()

set.seed(14253796)
juv_mod_boot <- bootMer(juv_mod, function(x) fixef(x)[["linesusceptible"]],
                         nsim = 2000)
# hist(juv_mod_boot$t)
quantile(juv_mod_boot$t, c(0.025, 0.5, 0.975))
