
source("load_data.R")

suppressPackageStartupMessages({
    library(rstan)
})


# Removing reps with missing data for now:
growth <- growth %>%
    filter(!(rep == 2 & line %in% c("R10", "UT3", "WI-L4")))

growth




