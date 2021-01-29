
source("under_constr/poster/_preamble.R")


# =======================================================================================
# =======================================================================================

#           Aphid growth plots

# =======================================================================================
# =======================================================================================


growth <-
    load_data(impute_fxn = impute, filter_pars = NULL) %>%
    mutate(line = paste(line)) %>%
    bind_rows(clonewars:::load_pz_data(impute_fxn = impute, filter_pars = NULL)) %>%
    mutate_at(vars(line, rep), list(factor)) %>%
    # Now filter out early part of each time series, before N > 6
    # N <= 6 is when the stochasticity associated with only starting with 2 adults
    # appears to be strongest
    group_by(line, rep) %>%
    filter(1:n() >= which(N > 6)[1]) %>%
    mutate(date = date - min(date),
           disp = ifelse(is.na(disp), 0, disp),
           N = round(N),
           X = log(N),
           pN = N - disp,  # numbers of aphids on plant
           dD = disp - lag(disp, default = 0), # number of new dispersed aphids
           dD = ifelse(dD < 0, 0, dD)) %>%
    ungroup()


reps <- growth %>%
    filter(line == "WI-2016-593") %>%
    .[["rep"]] %>%
    unique() %>%
    paste() %>%
    identity()
reps <- c(reps, reps[1:3])
# reps <- reps[1]



example_ts <- growth %>%
    filter(line == "WI-2016-593") %>%
    dplyr::select(line, rep, date, pN, dD) %>%
    ggplot(aes(date, group = rep)) +
    geom_line(aes(y = pN), alpha = 0.5, size = 1.5, color = palette$default_primary) +
    geom_line(aes(y = dD * 10), alpha = 0.75, size = 1.5, color = palette$secondary_text) +
    xlab(NULL)


# for (i in 1:length(example_ts)) {
ggsave(filename = "figs/growth_ts.pdf", plot = example_ts,
       width = 8, height = 8, units = "cm", bg = "white", useDingbats = FALSE)




