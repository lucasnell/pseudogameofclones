
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
library(patchwork)
library(gameofclones)
library(lme4)
library(here)




# ========================================================================*
# ========================================================================*

# Plots ====

# ========================================================================*
# ========================================================================*

# Palette for the two clonal lines.
# Equivalent to `viridis::viridis(100)[c(85, 10)]`.
clone_pal <- c("#99D83DFF", "#482173FF")

# -----------------------------------------------`
# __pop. growth ----
# -----------------------------------------------`



pop_df <- here("_results/_data/assays-population_growth.csv") |>
    read_csv() |>
    mutate(date = as.Date(paste(year, month, day, sep = "-")),
           rep = factor(rep)) |>
    group_by(rep, line) |>
    mutate(date = difftime(date, min(date), units = "days") |> as.integer(),
           line = factor(line, levels = c("UT3", "WIA-5D"),
                         labels = c("resistant", "susceptible"))) |>
    select(rep, line, date, num) |>
    mutate(extinct = factor(sum(num) == 0)) |>
    # Remove once they start declining bc this is due to plant death:
    group_by(rep) |>
    filter(date <= date[which(num == max(num))]) |>
    ungroup() |>
    arrange(rep, date, line)






pop_p <- pop_df |>
    mutate(zero = factor(num == 0)) |>
    ggplot(aes(date, log1p(num), color = line)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_line() +
    geom_point(aes(shape = zero), size = 2) +
    facet_wrap(~ rep, nrow = 2) +
    scale_color_manual(values = clone_pal) +
    scale_shape_manual(values = c(19, 4), guide = "none") +
    scale_y_continuous("Aphid abundance",
                       breaks = log1p(c(0, 3 * 10^(0:2))),
                       labels = c(0, 3 * 10^(0:2))) +
    scale_x_continuous("Days after start") +
    theme(legend.title = element_blank(),
          strip.text = element_blank()) +
    NULL

# pop_p

# save_plot("_results/plots/assays-competition.pdf", pop_p, 5, 3)





# ========================================================================*
# ========================================================================*

# -----------------------------------------------`
# __resistance ----
# -----------------------------------------------`

# ========================================================================*
# ========================================================================*

#'
#' We did two sets of wasp-resistance assays: choice and no-choice.
#'
#' In the choice assays, we added 10 juveniles
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
#' In the no-choice assays, we first exposed 10 juvenile aphids of a single
#' line to wasps for about 3 hours, then exposed 10 juveniles of the other
#' line to the same wasps for the same duration.
#' We did this for both aphid lines being the first to be exposed.
#'







wasp_df <- bind_rows(
    here("_results/_data/assays-wasp_resistance_choice.csv") |>
        read_csv(col_types = cols()) |>
        select(wasp_group, line, starts_with(c("juv","adult-", "mumm"))) |>
        mutate(set = 1L),
    here("_results/_data/assays-wasp_resistance_no_choice.csv") |>
        read_csv(col_types = cols()) |>
        select(wasp_group, round, line,
               starts_with(c("juv", "adult-", "mumm"))) |>
        mutate(set = 2L,
               #' Because there were 10 groups in the first set of assays,
               #' and these groups differ from those used in set 1
               wasp_group = wasp_group + 10)) |>
    select(set, round, everything()) |>
    rename_with(function(x) gsub("-", "_", x)) |>
    mutate(survived = pmin(juv_assayed - mummy_end, adult_end),
           p_survived = survived / juv_assayed,
           non_survived = juv_assayed - survived,
           p_mummy = mummy_end / juv_assayed,
           non_mummy = juv_assayed - mummy_end,
           line = factor(line, levels = c("UT3", "WIA-5D"),
                         labels = c("resistant", "susceptible")),
           set = factor(set, levels = 1:2,
                        labels = c("choice", "no-choice")),
           across(c(round, wasp_group), factor),
           id = case_when(set == "choice" ~ 0L,
                          round == 1 ~ 1L,
                          TRUE ~ 2L) |>
               factor(),
           obs = 1:n())


set.seed(1257460453)
wasp_boots <- lapply(1:2000, function(i) {
    si <- sample(which(wasp_df$line == "susceptible"), replace = TRUE)
    ri <- sample(which(wasp_df$line == "resistant"), replace = TRUE)
    od <- wasp_df[c(si[1], ri[1]), "line"]
    od[["rep"]] <- i
    for (x in c("p_mummy", "juvenile_end", "p_survived")) {
        od[[x]] <- c(mean(wasp_df[[x]][si]), mean(wasp_df[[x]][ri]))
    }
    return(od)
}) |>
    do.call(what = bind_rows)

wasp_boot_ci <- wasp_boots |>
    group_by(line) |>
    summarize(across(p_mummy:p_survived, list(~ quantile(.x, 0.025),
                                              ~ quantile(.x, 0.5),
                                              ~ quantile(.x, 0.975)))) |>
    rename_with(~ gsub("_2", "", .x), ends_with("_2"))



mummy_p <- wasp_df |>
    ggplot(aes(line, p_mummy, color = line)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_jitter(aes(shape = set), height = 0, width = 0.25, alpha = 0.5) +
    stat_summary(geom = "point", fun = mean, size = 3) +
    geom_linerange(data = wasp_boot_ci,
                    aes(ymin = p_mummy_1, ymax = p_mummy_3)) +
    scale_color_manual(values = clone_pal, guide = "none") +
    scale_shape_manual(NULL, values = c(15, 17)) +
    ylab("Mummy proportion") +
    NULL

# mummy_p




juv_p <- wasp_df |>
    ggplot(aes(line, juvenile_end, color = line)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_jitter(aes(shape = set), height = 0, width = 0.25, alpha = 0.5) +
    stat_summary(geom = "point", fun = mean, size = 3) +
    geom_linerange(data = wasp_boot_ci,
                    aes(ymin = juvenile_end_1,
                        ymax = juvenile_end_3)) +
    scale_color_manual(values = clone_pal, guide = "none") +
    scale_shape_manual(NULL, values = c(15, 17)) +
    scale_y_continuous("Number of juveniles") +
    NULL

# juv_p




surv_p <- wasp_df |>
    ggplot(aes(line, p_survived, color = line)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_jitter(aes(shape = set), height = 0, width = 0.25, alpha = 0.5) +
    stat_summary(geom = "point", fun = mean, size = 3) +
    geom_linerange(data = wasp_boot_ci,
                    aes(ymin = p_survived_1,
                        ymax = p_survived_3)) +
    scale_color_manual(values = clone_pal, guide = "none") +
    scale_shape_manual(NULL, values = c(15, 17)) +
    ylab("Survival proportion") +
    NULL

# surv_p



wasp_p <- mummy_p + surv_p + juv_p +
    plot_annotation(tag_levels = "A") +
    plot_layout(nrow = 1, guides = "collect") &
    theme(plot.tag = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 9, color = "black", angle = 45,
                                     vjust = 1, hjust = 1),
          axis.title.x = element_blank())

# save_plot("_results/plots/assays-wasps.pdf", wasp_p, 6.5, 3.5, seed = 45670)






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
#' This should be a conservative assessment.
#'
sum_diffs <- pop_df |>
    group_by(rep, line) |>
    summarize(num = sum(num), .groups = "drop")

obs_diff <- get_avg_diff(sum_diffs)

set.seed(654067829)
diff_perms <- replicate(2000, {
    ..sum_diffs <- sum_diffs
    ..sum_diffs[["line"]] <- sample(..sum_diffs[["line"]])
    get_avg_diff(..sum_diffs)
})

mean(abs(diff_perms) >= abs(obs_diff))
# [1] 0.0305






# -----------------------------------------------`
# __resistance ----
# -----------------------------------------------`


obs_mum <- get_avg_diff(wasp_df, .col = "mummy_end")
obs_juv <- get_avg_diff(wasp_df, .col = "juvenile_end")
obs_surv <- get_avg_diff(wasp_df, .col = "p_survived")

set.seed(678530)
mum_perms <- replicate(2000, {
    ..wasp_df <- wasp_df
    ..wasp_df[["line"]] <- sample(..wasp_df[["line"]])
    get_avg_diff(..wasp_df, .col = "mummy_end")
})
juv_perms <- replicate(2000, {
    ..wasp_df <- wasp_df
    ..wasp_df[["line"]] <- sample(..wasp_df[["line"]])
    get_avg_diff(..wasp_df, .col = "juvenile_end")
})
surv_perms <- replicate(2000, {
    ..wasp_df <- wasp_df
    ..wasp_df[["line"]] <- sample(..wasp_df[["line"]])
    get_avg_diff(..wasp_df, .col = "p_survived")
})

# Values of 0 actually indicate < 0.0005
mean(abs(mum_perms) >= abs(obs_mum))
# [1] 0
mean(abs(juv_perms) >= abs(obs_juv))
# [1] 0.0015
mean(abs(surv_perms) >= abs(obs_surv))
# [1] 0



#'
#' For survival and number of juveniles, there were some differences caused by
#' assay set / round.
#' Because the `wasp_group` factor differs by assay type / round, we included
#' that as a random effect in the regressions below to show that this
#' doesn't affect our conclusions based on the permutation tests.
#' If anything, the permutations appear to be conservative, so they are
#' what I present in the main text.
#'


surv_mod <- glmer(cbind(survived, non_survived) ~ line + id + (1 | wasp_group),
                  wasp_df, family = binomial)
surv_mod |> summary()

set.seed(9876345)
surv_mod_boot <- bootMer(surv_mod, function(x) fixef(x)[["linesusceptible"]],
                         nsim = 2000)
# hist(surv_mod_boot$t)
quantile(surv_mod_boot$t, c(0.025, 0.5, 0.975))
#      2.5%       50%     97.5%
# -2.469087 -1.883217 -1.362103

juv_mod <- glmer(juvenile_end ~ line + id + (1 | wasp_group),
                 wasp_df, family = poisson)
juv_mod |> summary()

set.seed(14253796)
juv_mod_boot <- bootMer(juv_mod, function(x) fixef(x)[["linesusceptible"]],
                         nsim = 2000)
# hist(juv_mod_boot$t)
quantile(juv_mod_boot$t, c(0.025, 0.5, 0.975))
#       2.5%        50%      97.5%
# -1.0873879 -0.9627950 -0.8421847

