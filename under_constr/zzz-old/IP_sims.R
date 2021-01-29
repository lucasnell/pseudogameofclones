
# This file is for use with IPs, Spring 2019

# This file uses peak abundance for density dependences
# and max per-capita growth to esimate intrinsic growth rates

library(clonewars)



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
# Calculate density dependence
get_a <- function(.df) {
    # Inverse of max density over time series:
    a_ <- 1 / max(.df$N, na.rm = TRUE)
    tibble(line = .df$line[1], rep = .df$rep[1], a = a_)
}

# I can use all data for the growth rates:
R <- load_data(filter_pars = NULL, remove_unfinished = FALSE) %>%
    split(interaction(.$line, .$rep, drop = TRUE)) %>%
    map_dfr(~ get_r(.x))

A <- load_data(filter_pars = NULL, remove_unfinished = TRUE) %>%
    split(interaction(.$line, .$rep, drop = TRUE)) %>%
    map_dfr(~ get_a(.x)) %>%
    mutate(k = 1 / a)


# ======================================================================================
# ======================================================================================

#           Simulations using these numbers

# ======================================================================================
# ======================================================================================


# A %>% group_by(line) %>% summarize(k = mean(k)) %>%
#     mutate(line = factor(line, levels = line[order(k)])) %>%
#     ggplot() +
#     geom_segment(aes(line, min(k), xend = line, yend = k)) +
#     geom_point(aes(line, k, color = line)) +
#     ylab("Carrying capacity") +
#     coord_flip() +
#     theme(axis.title.y = element_blank(),
#           axis.text.y = element_blank()) +
#     scale_color_brewer(palette = "Dark2", guide = FALSE, direction = -1)



R_byline <- R %>% group_by(line) %>% summarize(Rm = mean(r))
#
# R_byline %>%
#     mutate(line = factor(line, levels = line[order(Rm)])) %>%
#     ggplot() +
#     geom_segment(aes(line, min(Rm), xend = line, yend = Rm)) +
#     geom_point(aes(line, Rm, color = line)) +
#     ylab("Growth rate") +
#     coord_flip() +
#     theme(axis.title.y = element_blank(),
#           axis.text.y = element_blank()) +
#     scale_color_brewer(palette = "Dark2", guide = FALSE, direction = -1)


d <- c(-3.75, 0.00, -2.25, -0.25, -2, -2.5, -2, -3.5)

# tibble(line = levels(R_byline$line), disp = d) %>%
#     write_csv(path = "~/Desktop/disp_rates.csv")


set.seed(78123456)
sim_df <- sim_reps(n_reps = 100, max_t = 500, save_every = 1, N0 = 6,
         n_patches = 8,
         R = R %>% group_by(line) %>% summarize(Rm = mean(r)) %>% .[["Rm"]],
         A = A %>% group_by(line) %>% summarize(Am = mean(a)) %>% .[["Am"]] %>%
             mean() %>% rep(8),
         D_vec = d,
         repl_times = seq(3, 1000, 3),
         repl_threshold = 10e3,
         zeta_t_thresh = 20,
         # log_zeta_sd = 0, process_error = FALSE, disp_error = FALSE,
         n_cores = 4) %>%  # ,
    # process_error = 0, disp_error = FALSE,
    # log_zeta_sd = 0, extinct_N = 0) %>%
    mutate_at(vars(rep, patch), factor) %>%
    identity() %>%
    # filter(rep %in% 0:15) %>%
    group_by(rep, time, line) %>%
    summarize(N = sum(N)) %>%
    group_by(rep, time) %>%
    mutate(prop = N / sum(N),
           lower_prop = c(0, cumsum(prop)[-n()]),
           upper_prop = cumsum(prop)) %>%
    ungroup() %>%
    mutate(line = factor(line, levels = rev(paste(R_byline$line[order(R_byline$Rm)]))))


# Deterministic plot:
# sim_df %>%
#     filter(rep == 0) %>%
#     ggplot(aes(time, N)) +
#     # geom_ribbon(aes(y = NULL, ymin = lower_prop, ymax = upper_prop, fill = line)) +
#     geom_line(aes(color = line), size = 0.5) +
#     # geom_line(aes(color = patch)) +
#     # facet_wrap( ~ rep, nrow = 10) +
#     # facet_wrap( ~ patch, nrow = 2) +
#     # theme(legend.position = "none") +
#     theme(strip.text = element_blank()) +
#     scale_color_brewer(palette = "Dark2") +
#     scale_fill_brewer(palette = "Dark2") +
#     theme(axis.text = element_blank(), axis.ticks = element_blank()) +
#     guides(color = guide_legend(override.aes = list(size = 1)))




{
    sim_df %>%
        filter(as.numeric(paste(rep)) < 4) %>%  # <-- all reps is too much
        ggplot(aes(time, log(N))) +
        # geom_vline(xintercept = 365, linetype = 2) +
        geom_line(aes(color = line), size = 0.3) +
        facet_wrap(~ rep, nrow = 2) +
        scale_color_brewer(palette = "Dark2") +
        theme(strip.text = element_blank()) +
        xlab("Time (days)") +
        scale_y_continuous("log(Aphid abundance)") +
        coord_cartesian(ylim = c(0, 8.78777)) +
        # theme(axis.title = element_blank(),
        #       axis.text = element_blank(),
        #       axis.ticks = element_blank()) +
        guides(color = guide_legend(override.aes = list(size = 1)))
    } %>%
    ggsave(filename = "~/Desktop/time_series.pdf", width = 6, height = 4)

{
    sim_df %>%
        filter(time == max(time)) %>%
        group_by(rep) %>%
        summarize(N = sum(N > 0)) %>%
        ungroup() %>%
        {table(.$N)} %>%
        as.data.frame() %>%
        as_tibble() %>%
        ggplot(aes(Var1, Freq / 100)) +
        geom_bar(stat = "identity", fill = "dodgerblue", color = "gray20", size = 0.25) +
        xlab("Number of lines surviving") +
        ylab("Proportion of simulations")
    } %>%
    ggsave(filename = "~/Desktop/coexistence.pdf", width = 2.5, height = 4)

sim_df %>%
    filter(time == max(time)) %>%
    group_by(rep) %>%
    summarize(N = sum(N > 0)) %>%
    ungroup() %>%
    {table(.$N)} %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    set_names(c("n_spp", "prob_times")) %>%
    mutate(n_spp = as.integer(n_spp), prob_times = prob_times / 100) %>%
    write_csv(path = "~/Desktop/coexistence.csv")




n_won <- sim_df %>%
    filter(time == max(time)) %>%
    group_by(line) %>%
    summarize(N = sum(N > 0)) %>%
    ungroup() %>%
    mutate(line = factor(line, levels = paste(R_byline$line[order(R_byline$Rm)])))

n_won %>%
    mutate(prob_times = N / 100) %>%
    select(-N) %>%
    write_csv(path = "~/Desktop/prob_survived.csv")


{
    n_won %>%
        ggplot(aes(line, N / 100)) +
        geom_hline(yintercept = 1, linetype = 2, color = "gray60") +
        geom_segment(aes(xend = line, yend = 0)) +
        # geom_point(data = R_byline %>%
        #                mutate(Rm = ((Rm - min(Rm)) / (max(Rm) - min(Rm))) *
        #                           (diff(range(n_won$N))) + min(n_won$N)),
        #           aes(line, Rm), shape = 1) +
        geom_point(aes(color = line), size = 3) +
        scale_color_brewer(palette = "Dark2", guide = FALSE, direction = -1) +
        scale_y_continuous(breaks = c(0, 0.5, 1)) +
        coord_flip(ylim = c(0, 1.05)) +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_text(size = 11, color = "black")) +
        ylab("Proportion of\ntimes survived")
    } %>%
    ggsave(filename = "~/Desktop/prop_survived.pdf", width = 3, height = 5)



sim_df %>%
    group_by(rep, time) %>%
    summarize(nlines = sum(N>0)) %>%
    group_by(rep) %>%
    summarize(tl = time[nlines == tail(nlines, 1)][1]) %>%
    ggplot(aes(tl)) +
    geom_histogram(bins = 25, fill = "dodgerblue") +
    xlab("Time to final line number")
