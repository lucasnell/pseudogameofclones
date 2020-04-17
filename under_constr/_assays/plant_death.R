
suppressPackageStartupMessages({
    library(tidyverse)
    library(clonewars)
    library(viridisLite)
    library(MASS)
})

source(".Rprofile")


# Age at which plants start dying:
death_days <- load_data() %>%
    group_by(line, rep) %>%
    summarize(md = max(date)) %>%
    .[["md"]]

death_days %>%
    {print(mean(.)); print(sd(.)); .} %>%
    hist()

# Use this for lognormal distribution:
# (Lognormal worked better than normal or Poisson; negative binomial failed to fit)
fitdistr(death_days, "lognormal") %>%
    .[["estimate"]]


# Calculate growth rate
get_r <- function(.df) {
    .df <- .df[(which(.df$N > 20)[1]):nrow(.df),]
    pcg <- .df$X[-1] - lag(.df$X)[-1] # per-capita growth
    n_t <- 4 # number of time steps to take average by
    if (length(pcg) < 4) return(NULL)
    # Max average pcg over `n_t` time steps:
    r_ <- max(sapply(1:(length(pcg) - n_t + 1), function(i) mean(pcg[i:(i+n_t-1)])))
    tibble(line = .df$line[1], rep = .df$rep[1], r = r_)
}

R <- load_data(filter_pars = NULL, remove_unfinished = FALSE) %>%
    split(interaction(.$line, .$rep, drop = TRUE)) %>%
    map_dfr(get_r)

# Then look at what happens after death starts:
load_data(filter_pars = NULL) %>%
    group_by(line, rep) %>%
    filter(date >= min(date[N == max(N)])) %>%
    mutate(date = date - min(date),
           n_obs = n()) %>%
    ungroup() %>%
    filter(n_obs >= 3) %>%
    ggplot(aes(date, X)) +
    geom_line(aes(group = rep), alpha = 0.5) +
    facet_wrap(~ line, nrow = 2)

# Calculate growth rate after plant death
get_death_r <- function(.df) {
    if (any(.df$N == 0)) .df <- .df[1:(which(.df$N == 0)[1]),]
    pcg <- .df$X[-1] - lag(.df$X)[-1] # per-capita growth
    r_ <- mean(pcg)
    tibble(line = .df$line[1], rep = .df$rep[1], r = r_)
}


death_r <- load_data(filter_pars = NULL) %>%
    group_by(line, rep) %>%
    filter(date >= min(date[N == max(N)])) %>%
    mutate(date = date - min(date),
           n_obs = n()) %>%
    ungroup() %>%
    filter(n_obs >= 3) %>%
    split(interaction(.$line, .$rep, drop = TRUE)) %>%
    map_dfr(get_death_r) %>%
    # How much plant death modifies per-capita growth rate
    set_names(c(".line", ".rep", ".r")) %>%
    mutate(mod = pmap_dbl(list(.line, .rep, .r),
                          function(.line, .rep, .r) {
                              r <- filter(R, line == .line, rep == .rep) %>% .[["r"]]
                              return(.r - r)
                          })) %>%
    rename_all(~ gsub("^\\.", "", .x))

death_r %>%
    ggplot(aes(line, mod, color = line)) +
    geom_jitter(height = 0, width = 0.25) +
    scale_color_brewer(palette = "Dark2")


death_rates <- death_r %>%
    .[["mod"]] %>%
    exp()

death_rates %>%
    {print(mean(.)); print(sd(.)); .} %>%
    hist()

# Use these for Beta distribution:
fitdistr(death_rates, "beta", list(shape1 = 1, shape2 = 1)) %>%
    .[["estimate"]]
