suppressPackageStartupMessages({
    library(tidyverse)
    library(readxl)
})

# Function for figuring out how many aphids were mentioned in the comments
parse_comments <- function(comments) {
    comments = strsplit(comments, ' ')
    out <- map(comments, ~ suppressWarnings(as.numeric(.x))) %>% 
        map(~ if(all(is.na(.x))) {0} else {sum(.x[!is.na(.x)])}) %>% 
        unlist
    return(out)
}
# Function for filtering out decreasing N
dec_filter <- function(N_vec, p = 0.2) {
    # We can assume blank ones are zeros
    N_vec <- ifelse(is.na(N_vec), 0, N_vec)
    max_N <- max(N_vec)
    max_ind <- which(N_vec == max_N)[1]
    inds <- unique(c(1:max_ind, which(N_vec > (1 - p) * max_N)))
    return(1:length(N_vec) %in% inds)
}

# Function for filtering out N above a threshold
threshold_filter <- function(N_vec, threshold = 50) {
    # We can assume blank ones are zeros
    N_vec <- ifelse(is.na(N_vec), 0, N_vec)
    # Index to first above the threshold
    first_ind <- which(N_vec >= threshold)[1]
    return(1:length(N_vec) >= first_ind)
}

# Lines that we still have and should keep for analyses
lines_to_keep <- c("R10", "WIA-5D", "WI-L4", "WI-L4Ø", "UT3", "WI-2016-593",
                   "Clover-2017-2", "Clover-2017-6", "Clover-2017-9")

growth <- read_excel(paste0('~/Dropbox/Aphid Project 2017/Lucas_traits/',
                            'traits_data_entry.xlsx')) %>%
    mutate(line = ifelse(line == 'WI-L4 (H+3)', 'WI-L4', line),
           line = ifelse(line == 'WI-L4ØA', 'WI-L4Ø', line)) %>% 
    filter(line %in% lines_to_keep) %>% 
    mutate(date = ifelse(date == '1/31/018', '1/31/2018', date),
           date = as.Date(date, '%m/%d/%Y'), 
           comments = parse_comments(comments), 
           N = stem1_juv + stem1_adults + leaf1_juv + leaf1_adults + 
               stem2_juv + stem2_adults + leaf2_juv + leaf2_adults + 
               stem3_juv + stem3_adults + leaf3_juv + leaf3_adults + 
               comments,
           # makes no sense for it to be 0, then >0 the next day:
           N = ifelse(N == 0, 1, N)) %>% 
    select(line, rep, date, N) %>% 
    mutate(rep = as.integer(rep), N = as.integer(N),
           line = factor(line)) %>%
    group_by(line, rep) %>%
    filter(dec_filter(N)) %>%
    filter(threshold_filter(N)) %>%
    mutate(date = as.integer(date - min(date)),
           r = log(lead(N) / N) / (lead(date) - date)) %>%
    ungroup

rm(lines_to_keep, dec_filter, parse_comments, threshold_filter)


