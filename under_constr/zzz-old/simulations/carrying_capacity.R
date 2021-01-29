
#'
#' The parameter `K` in the model isn't an exact match to the carrying capacity.
#' The actual carrying capacity depends on `K` and the Leslie matrix.
#'


library(clonewars)

sim_CC <- function(K, L = NULL) {

    S <- function(z) 1 / (1 + z / K)

    if (is.null(L)) {
        L <- leslie_matrix(dev_times$instar_days$lowT, populations$surv_juv$high,
                           populations$surv_adult$high, populations$repro$high)
    }

    N <- numeric(1+1000)
    X <- sad_leslie(L)

    N[1] <- sum(X)

    for (t in 2:length(N)) {
        X <- S(N[t-1]) * L %*% X
        N[t] <- sum(X)
    }

    # plot(N, type = "l")

    return(max(N))
}

Ks <- sapply(seq(500, 10e3, 500), foo)

tibble(K = seq(500, 10e3, 500),
       cc = Ks) %>%
    lm(formula = cc ~ K)


# Predicted carrying capacity for Leslie matrix `L` and parameter `K`:
CC <- function(K) {
    (Re(eigen(L)$values[1]) - 1) * K * x
}


