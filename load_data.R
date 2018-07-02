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
    mutate_at(vars(rep, N, disp), funs(as.integer)) %>%
    mutate(line = factor(line),
           X = log(N),
           ham = ifelse(line %in% w_ham, 1, 0)) %>%
    group_by(line, rep) %>%
    mutate(date = as.integer(date - min(date)),
           r = log(N / lag(N)) / (date - lag(date))) %>%
    arrange(date) %>%
    ungroup %>%
    arrange(line, rep, date) %>%
    # Filter the beginning and end of the time series:
    group_by(line, rep) %>%
    # Beginning filter:
    filter(threshold_filter(X, 0.50)) %>%
    # End filter:
    filter(dec_filter(X, 0.0)) %>%
    ungroup %>%
    identity()



rm(lines_to_keep, dec_filter, parse_comments, threshold_filter, w_ham)




