## code to prepare `rr_rs_lookup` internal dataset goes here

#' This dataset is a lookup table for short-term fitness for resistant versus
#' susceptible aphids as a function of the parasitoid attack rate.
#' It is based on Ives et al. (2020) code.

# Set threads this way if you're on unix:
# options(mc.cores = THREADS_TO_USE)

growthcost = 6.5266912e-01
resist = 2.7982036e-01
Leslie_s <- matrix(c(0.5, 0, 0, 0, 2.55, 0.5, 0.5, 0, 0, 0, 0, 0.5, 0.5,
                     0, 0, 0, 0, 0.5, 0.5, 0, 0, 0, 0, 0.5, 0.8),
                   byrow = TRUE, nrow=5)
Leslie_r <- Leslie_s
Leslie_r[1,5] <- growthcost * Leslie_r[1,5]
relatt <- c(0.1200, 0.2700, 0.3900, 0.1600, 0.0600)
kk <- 0.3475
yLeslie <- matrix(c(0.5, 0, 0, 0, 0,
                    0.5, 0.5, 0, 0, 0,
                    0, 0.5, 0.5, 0, 0,
                    0, 0, 0.5, 0.667, 0,
                    0, 0, 0, 0.333, 0.95),
                  byrow = TRUE, nrow = 5)
one_calc <- function(a){
    A <- (1 + exp(a) * relatt / kk)^(-kk)
    rs <- max(abs(eigen(diag(A) %*% Leslie_s)$values))
    rr <- max(abs(eigen(diag(1 - resist * (1-A)) %*% Leslie_r)$values))
    B11 <- diag(A) %*% Leslie_s
    B21 <- matrix(0,nrow=5, ncol=5)
    B21[1,] <- 1-A
    B12 <- matrix(0,nrow=5, ncol=5)
    B22 <- yLeslie
    B <- rbind(cbind(B11, B12), cbind(B21, B22))
    SADA <- eigen(B)$vectors[,1]
    ps <- Re(sum(SADA[8:9])/sum(SADA[c(4:5,8:9)]))
    if (is.nan(ps)) ps <- 1
    B11 <- diag(1 - resist*(1-A)) %*% Leslie_r
    B21 <- matrix(0,nrow=5, ncol=5)
    B21[1,] <- resist*(1-A)
    B12 <- matrix(0,nrow=5, ncol=5)
    B22 <- yLeslie
    B <- rbind(cbind(B11, B12), cbind(B21, B22))
    SADA <- eigen(B)$vectors[,1]
    pr <- Re(sum(SADA[8:9])/sum(SADA[c(4:5,8:9)]))
    rr_rs <- rr / rs
    return(c(a, rs, rr, ps, pr, rr_rs))
}

# Takes ~15 sec on my mac using 6 threads.
if (.Platform$OS.type == "unix") {
    rr_rs_list <- parallel::mclapply(round(seq(-20, 5, 1e-4), 4), one_calc)
} else {
    rr_rs_list <- lapply(round(seq(-20, 5, 1e-4), 4), one_calc)
}
rr_rs_lookup <- as.data.frame(do.call(rbind, rr_rs_list))
colnames(rr_rs_lookup) <- c("a", "rs", "rr", "ps", "pr", "rr_rs")

usethis::use_data(rr_rs_lookup, overwrite = TRUE, internal = TRUE, compress = "xz")
