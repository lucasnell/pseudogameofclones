

suppressPackageStartupMessages({
    library(aphidreps)
    library(tidyverse)
})
source(".Rprofile")

set.seed(9)

# Number of aphid lines:
n_lines <- 8L
# Reps per aphid line
n_reps <- as.integer(round(runif(n_lines, 3, 7)))
# Aphid line per time series
line_ts <- c(lapply(1:length(n_reps), function(i) rep(i, n_reps[i])), recursive = TRUE)
# # observations per time series
nobs_ts <- as.integer(round(runif(sum(n_reps), 10, 15)))

# R0 <- round(rnorm(n_lines, 0.29, 0.003), 5)
# alpha <- round(rnorm(n_lines, 0.002, 0.0003), 5)
R0 <- c(0.308, 0.321, 0.297, 0.301, 0.302, 0.268, 0.267, 0.258)
alpha <- c(0.00177, 0.00195, 0.00301, 0.00237, 0.00243, 0.00229, 0.00239, 0.00249)
sigma_process <- 0.1

X <- sim_lines(R0, alpha, n_reps, nobs_ts, sigma_process)


matplot(X, type = "l")
abline(h = log(1 / alpha), lty = 2, col = "red")


fit <- fit_lines(X, line_ts, control = list(adapt_delta = 0.90))


# fit <- fit_one_line(X, control = list(adapt_delta = 0.85))
# adapt_delta = 0.80: 198 divergent transitions after warmup
# adapt_delta = 0.85: 2 divergent transitions after warmup
# adapt_delta = 0.90: 8 divergent transitions after warmup
# adapt_delta = 0.95: 1 divergent transitions after warmup

print(fit, digits = 6, pars = "R0"); cat(crayon::bgRed(R0))
print(fit, digits = 6, pars = "alpha"); cat(crayon::bgRed(alpha))
print(fit, digits = 6, pars = "process")

rank(apply(rstan::extract(fit, "R0", permuted = FALSE), 3, mean)); rank(R0)
apply(rstan::extract(fit, "R0", permuted = FALSE), 3, mean); R0



posterior <- as.array(fit)

dim(posterior)

alpha_array <- posterior[,,which(dimnames(posterior)$parameters == "alpha[1]")]


running_means_alpha <- lapply(1:nrow(alpha_array),
                              function(n) {
                                  sapply(1:4, function(i) mean(log(alpha_array[1:n,i])))
                              })
running_means_alpha <- do.call(rbind, running_means_alpha)
running_means_alpha <- running_means_alpha %>%
    tbl_df() %>%
    setNames(c("chain1", "chain2", "chain3", "chain4")) %>%
    mutate(iteration = 1:n()) %>%
    gather("chain", "log_alpha", chain1:chain4, factor_key = TRUE)

running_means_alpha %>%
    ggplot(aes(iteration, log_alpha, color = chain)) +
    geom_hline(yintercept = log(alpha[1]), linetype = 2) +
    geom_line() +
    NULL




# plot(x = fit, pars = c("R0", "alpha", "sigma_process"))
mcmc_trace(posterior, pars = c("alpha")) +
    scale_y_continuous(trans = "log")


source("under_constr/load_data.R")

# Removing reps with missing data for now:
growth <- growth %>%
    filter(!(rep == 2 & line %in% c("R10", "UT3", "WI-L4")))

source(".Rprofile")


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
#     summarize(alpha = 1 / max(N, na.rm = TRUE)) %>%
#     group_by(line) %>%
#     summarize(alpha = mean(alpha, na.rm = TRUE)) %>%
#     ungroup() %>%
#     summarize(alpha_mean = mean(alpha, na.rm = TRUE),
#               alpha_sd = sd(alpha, na.rm = TRUE)) %>%
#     identity()
#


