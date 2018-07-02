suppressPackageStartupMessages({
    library(tidyverse)
    library(readxl)
})

theme_set(theme_classic() %+replace%
              theme(strip.background = element_blank(),
                    strip.text = element_text(size = 11),
                    legend.background = element_blank(),
                    plot.title = element_text(size = 14, hjust = 0)))

# Function for figuring out how many aphids were mentioned in the comments
parse_comments <- function(comments) {
    comments = strsplit(comments, ' ')
    out <- map(comments, ~ suppressWarnings(as.numeric(.x))) %>%
        map(~ if(all(is.na(.x))) {0} else {sum(.x[!is.na(.x)])}) %>%
        unlist
    out <- ifelse(is.na(out) | is.null(out), 0, out)
    return(out)
}
# Function for filtering out decreasing N
dec_filter <- function(X_vec, p) {
    max_N <- max(X_vec)
    max_ind <- which(X_vec == max_N)[1]
    inds <- unique(c(1:max_ind, which(X_vec > (1 - p) * max_N)))
    return(1:length(X_vec) %in% inds)
}

# Function for filtering out N above a threshold
threshold_filter <- function(X_vec, p) {
    max_N <- max(X_vec)
    threshold <- max_N * p
    # Index to first above the threshold
    first_ind <- which(X_vec >= threshold)[1]
    return(1:length(X_vec) >= first_ind)
}

# Lines that we still have and should keep for analyses
lines_to_keep <- c("R10", "WIA-5D", "WI-L4", "WI-L4Ø", "UT3", "WI-2016-593",
                   "Clover-2017-2", "Clover-2017-6")

w_ham <- c("R10", "WI-L4", "UT3", "Clover-2017-2")


growth <- read_excel(paste0('~/Dropbox/Aphid Project 2017/Lucas_traits/',
                            'traits_data_entry.xlsx')) %>%
    mutate(line = ifelse(line == 'WI-L4 (H+3)', 'WI-L4', line),
           line = ifelse(line == 'WI-L4ØA', 'WI-L4Ø', line),
           date = as.Date(paste(year, month, day, sep = "-"))) %>%
    # Change any NAs to zeros:
    mutate_at(vars(matches("_juv$|_adults$")), function(x) ifelse(is.na(x), 0, x)) %>%
    filter(line %in% lines_to_keep) %>%
    mutate(comments = parse_comments(comments),
           N = stem1_juv + stem1_adults + leaf1_juv + leaf1_adults +
               stem2_juv + stem2_adults + leaf2_juv + leaf2_adults +
               stem3_juv + stem3_adults + leaf3_juv + leaf3_adults +
               comments,
           disp = comments,  # <-- "disp" is for dispersed aphids, ones not on the plant
           # makes no sense for it to be 0, then >0 the next day:
           N = ifelse(N == 0, 1, N)) %>%
    select(line, rep, date, N, disp) %>%
    mutate(rep = as.integer(rep), N = as.integer(N), disp = as.integer(disp),
           line = factor(line)) %>%
    group_by(line, rep) %>%
    mutate(date = as.integer(date - min(date)),
           r = log(N / lag(N)) / (date - lag(date))) %>%
    arrange(date) %>%
    mutate(N_t = lag(N)) %>%
    ungroup %>%
    mutate(X = log(N), X_t = log(N_t),
           ham = ifelse(line %in% w_ham, 1, 0)) %>%
    arrange(line, rep, date) %>%
    filter(!is.na(N_t)) %>%
    # Filter the beginning and end of the time series:
    group_by(line, rep) %>%
    # Beginning filter:
    filter(threshold_filter(X, 0.50)) %>%
    # End filter:
    filter(dec_filter(N, 0.0)) %>%
    ungroup %>%
    identity()



rm(lines_to_keep, dec_filter, parse_comments, threshold_filter, w_ham)


# How I decided on the beginning and ending filtering:

z_trans <- function(x) (x - mean(x)) / sd(x)

source(".Rprofile")

no_filter <- growth %>%
    ggplot(aes(date, X, color = factor(rep))) +
    geom_line(size = 0.75) +
    facet_wrap(~ line, nrow = 2) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    ggtitle("no filter")

no_filter +
    geom_hline(yintercept = log(40), linetype = 3) +
    geom_vline(xintercept = 5, linetype = 3) +
    NULL

back_filter <- growth %>%
    group_by(line, rep) %>%
    filter(dec_filter(N, 0.0)) %>%
    ungroup %>%
    ggplot(aes(date, X, color = factor(rep))) +
    geom_line(size = 0.75) +
    facet_wrap(~ line, nrow = 2) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    ggtitle("back filtered")

growth %>%
    group_by(line, rep) %>%
    # filter(threshold_filter(X, 0.50)) %>%
    filter(date >= 5) %>%
    filter(dec_filter(X, 0.0)) %>%
    ungroup %>%
    ggplot(aes(date, X, color = factor(rep))) +
    geom_line(size = 0.75) +
    facet_wrap(~ line, nrow = 2) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    ggtitle("filtered")

no_filter


growth %>%
    group_by(line, rep) %>%
    mutate(r = z_trans(r)) %>%
    # Back filter:
    filter(dec_filter(X, 0.0)) %>%
    # # Front filter:
    # filter(date >= 5) %>%
    filter(threshold_filter(X, 0.50)) %>%
    ungroup %>%
    ggplot(aes(date, r, color = factor(rep))) +
    geom_vline(xintercept = 5, linetype = 2) +
    # geom_hline(yintercept = 0, linetype = 2) +
    geom_hline(yintercept = c(-1, 1), linetype = 2) +
    geom_line(size = 0.75) +
    facet_wrap(~ line, nrow = 2) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    ggtitle("r")

back_filter +
    geom_vline(xintercept = 5, linetype = 2)

