
#'
#' Install `gameofclones` package from GitHub:
#'
if (!require("gameofclones")) {
    if (!require("remotes")) install.packages("remotes")
    remotes::install_github("lucasnell/gameofclones")
    library(gameofclones)
}

#'
#' WARNING: Some of the code below contains tibbles.
#' To remove these abominations, you can use this function on the output
#' from `sim_experiments` or `restart_experiment`:
#'
rm_tibs <- function(.sims) {
    for (n in c("aphids", "wasps")) {
        .sims[[n]] <- as.data.frame(.sims[[n]])
    }
    return(.sims)
}



#'
#' Define aphid line information.
#' Both lines start with 32 adult aphids.
#'

# Susceptible line: no resistance, high population growth rate
line_s <- clonal_line("susceptible",
                      density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
                      surv_juv_apterous = "high",
                      surv_adult_apterous = "high",
                      repro_apterous = "high")
# Resistant line: high resistance, low parasitized-aphid survival rate,
#                 low population growth rate
line_r <- clonal_line("resistant",
                      density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
                      resistant = TRUE,
                      surv_paras = 0.57,
                      surv_juv_apterous = "low",
                      surv_adult_apterous = "low",
                      repro_apterous = "low")


#' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#'
#' Testing `perturb` argument...
z <- sim_experiments(clonal_lines = c(line_s, line_r),
                     alate_field_disp_p = 0,
                     perturb = tibble(when = 50, where = 1, who = "resistant", how = 10),
                     max_t = 100, save_every = 1) %>%
    .[["aphids"]] %>%
    filter(type != "mummy") %>%
    group_by(time, field, line) %>%
    summarize(N = sum(N), .groups = "drop")
x <- sim_experiments(clonal_lines = c(line_s, line_r),
                     alate_field_disp_p = 0,
                     max_t = 100, save_every = 1) %>%
    .[["aphids"]] %>%
    filter(type != "mummy") %>%
    group_by(time, field, line) %>%
    summarize(N = sum(N), .groups = "drop")

z %>%
    ggplot(aes(time, N)) +
    geom_line(aes(color = line)) +
    geom_line(data = x, aes(group = line), color = "black", linetype = 2) +
    facet_wrap(~ field)

#' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<





#'
#' The resulting `aphid` objects have the following fields:
#'
#' - `name`       : string indicate the name of the line ("susceptible" or
#'                  "resistant")
#' - `density_0`  : 29x2 matrix with starting abundances of each stage,
#'                  for apterous aphids (col 1) and alates (col 2)
#' - `attack_surv`: length-2 vector of wasp attack survivals for singly and
#'                  multiply parasitized aphids
#' - `leslie`     : 29x29x3 array with the 3 leslie matrices for apterous
#'                  aphids, alates, and parasitized aphids, respectively.
#'                  The latter is largely ignored as in the original model,
#'                  but is included for compatibility with the C++ code.
#'


#'
#' The defaults for the `sim_experiments` function are designed to
#' approximately replicate the experimental results.
#' The two things you'll definitely need to provide it are the clonal lines
#' we've defined and the maximum time to run the simulations.
#' The default for the latter is 250 days, so that won't work for reaching
#' equilibrium.
#'
#' This takes ~ 5 sec on my computer:
#'
sims <- sim_experiments(clonal_lines = c(line_s, line_r),
                        max_t = 1e6, save_every = 1000)
sims <- rm_tibs(sims)
sims

#'
#' The resulting `cloneSims` object contains the following fields
#' (the last two don't show up when printing `sims` bc they shouldn't be
#'  messed with):
#'
#' - `aphids`       : A `tibble` of aphid abundances through time and by cage.
#'                    These data are not stage structured.
#'                    Note that in the simulations, the `field` column is
#'                    equivalent to the experimental cage.
#' - `wasps`        : A `tibble` of wasp abundances through time and by cage.
#' - `all_info`     : A list of length 1 containing a data frame with the
#'                    abundances at the end of the simulations.
#'                    These data are stage structured and should be used
#'                    as the starting point for perturbations if using
#'                    the `restart_experiment` function described below.
#' - `all_info_xptr`: A pointer to the C++ object used to run the simulations
#'                    that is used to restart the simulations for doing
#'                    perturbations. This should not be touched.
#' - `call`         : A long list containing the call information for the
#'                    simulations. This can and should be ignored.
#'

#'
#' To perturb and restart experimental simulations.
#'
#' Data frame of stage-structured ending abundances:
ss_df <- sims$all_info[[1]]

head(ss_df)

#' The data frame should have the following columns:
#'
#' - `field`: This is equivalent to experimental cage.
#'            Field 1 has wasps, field 2 doesn't.
#' - `plant`: This is a holdover from when I was simulating individual plants
#'            withering. It can be ignored.
#' - `line` : The name of the associated aphid line. This is an empty string
#'            for wasps and mummies.
#' - `type` : The organism type. This can be `"wasp"`, `"mummy"`, `"apterous"`,
#'            `"alate"`, or `"parasitized"`.
#' - `stage`: The stage (in days). There should be 29 stages for alate and
#'            apterous aphids, 7 for parasitized aphids, 4 for mummies, and
#'            1 for adult wasps.
#' - `N`    : Abundance.
#'

#'
#' How to perturb abundances and restart the simulations.
#' Three key points for perturbations:
#' 1) You shouldn't change the `all_info` field directly from the original
#'    simulation object (`sims` in this case).
#'    Always copy to a new object before editing.
#' 2) You should only change the `N` column and leave everything else the same.
#'
#'
#' Here are some examples:
#'

# Perturb wasps in the wasp cage:
inds <- ss_df$field == 1 & ss_df$type == "wasp"
ss_df[inds, "N"] <- ss_df[inds, "N"] * (1 + 1e-5)

# Perturb 10-day-old resistant alates in the no-wasp cage:
inds <- ss_df$field == 2 & ss_df$type == "alate" & ss_df$line == "resistant" & ss_df$stage == 10
ss_df[inds, "N"] <- ss_df[inds, "N"] * (1 + 1e-5)

# Do new simulation for 1 time step:
new_sims <- restart_experiment(sims, new_starts = ss_df, max_t = 1)
new_sims <- rm_tibs(new_sims)
new_sims

# Full stage-structured data from new sims:
new_ss_df <- new_sims$all_info[[1]]
head(new_ss_df)


#'
#' You can also change model parameters to see how they change results.
#'
#' See `?restart_experiment` for the arguments that you can change with
#' this function.
#' The defaults reported in this documentation are those from the
#' `sim_experiments` function, so those will be used here unless you
#' change them in the original call to `sim_experiments`.
#'
#' The main argument you may want to change here is the proportion of
#' alates that move to another "field" (i.e., cage) every day.
#' The default is 0.1.
#'
new_sims2 <- restart_experiment(sims, max_t = 1,
                                alate_field_disp_p = 0.1 * (1 + 1e-5))
new_sims2 <- rm_tibs(new_sims2)
new_ss_df2 <- new_sims2$all_info[[1]]
