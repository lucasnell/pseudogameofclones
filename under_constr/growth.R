
source("load_data.R")

suppressPackageStartupMessages({
    library(rstan)
})


# Removing reps with missing data for now:
growth <- growth %>%
    filter(!(rep == 2 & line %in% c("R10", "UT3", "WI-L4")))

X <- {growth %>%
        filter(line == "WIA-5D", rep == 2)}$X
N <- length(X)

growth_fit <- stan(file = "aphid_growth.stan",
                   data = c("X", "N"),
                   iter = 1000, chains = 4)
print(growth_fit)
plot(growth_fit)

# -Wunknown-pragmas -Wc++11-inline-namespace -Wunused-function -Wunneeded-internal-declaration -Wmacro-redefined
#




# growth %>%
#     group_by(line, rep) %>%
#     summarize(r = mean(r, na.rm = TRUE)) %>%
#     group_by(line) %>%
#     summarize(r = mean(r, na.rm = TRUE)) %>%
#     ungroup() %>%
#     summarize(r_mean = mean(r, na.rm = TRUE),
#               r_sd = sd(r, na.rm = TRUE)) %>%
#     identity()
#
# growth %>%
#     group_by(line, rep) %>%
#     summarize(N = max(N, na.rm = TRUE)) %>%
#     group_by(line) %>%
#     summarize(N = mean(N, na.rm = TRUE)) %>%
#     ungroup() %>%
#     summarize(N_mean = mean(N, na.rm = TRUE),
#               N_sd = sd(N, na.rm = TRUE)) %>%
#     identity()



