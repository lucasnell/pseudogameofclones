
#'
#' Create LaTeX table of stage-structured parameters.
#'

library(gameofclones)

line_s <- clonal_line("susceptible",
                      density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
                      surv_juv_apterous = "high",
                      surv_adult_apterous = "high",
                      repro_apterous = "high")
# Resistant line: high resistance, low population growth rate
line_r <- clonal_line("resistant",
                      density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
                      resistant = TRUE,
                      surv_paras = 0.57,
                      surv_juv_apterous = "low",
                      surv_adult_apterous = "low",
                      repro_apterous = "low")


# apterous, alates, parasitized
sl <- line_s$leslie
rl <- line_r$leslie

#' Combine fecundities and survivals from two clones into one string, and
#' add an extra hyphen to account for fact that it only goes to n-1
combine_fs <- function(sl_list){
    z <- mapply(function(x, y) sprintf("$%.4f~|~%.4f$", x, y),
                sl_list[["r"]], sl_list[["s"]])
    return(c(z, "$-$"))
}


# colSums in si and ui below is to account for smoothing parameter
si <- list(r = colSums(rl[-1,-ncol(rl),1]), s = colSums(sl[-1,-ncol(sl),1])) |>
    combine_fs()
fi <- list(r = rl[1,-1,1], s = sl[1,-1,1]) |>
    combine_fs()
ui <- list(r = colSums(rl[-1,-ncol(rl),2]), s = colSums(sl[-1,-ncol(sl),2])) |>
    combine_fs()
gi <- list(r = rl[1,-1,2], s = sl[1,-1,2]) |>
    combine_fs()

Ri <- rep(wasp_attack$rel_attack, c(2, 2, 2, 2, 21)) |>
    (function(x) x / sum(x))()
Rin <- Ri
# Because adult alates are not attacked:
Rin[(8 + 1):length(Rin)] <- 0
# Because they together sum to 1:
Rsum <- sum(Rin) + sum(Ri)
Ri  <- Ri / Rsum
Rin <- Rin / Rsum
# Convert to strings:
Ri <- sprintf("$%.4f$", Ri)
Rin <- sprintf("$%.4f$", Rin)

par_df <- data.frame(i = sprintf("$%-2s$", 1:29),
                     si, fi, ui, gi, Ri, Rin)

for (i in 1:nrow(par_df)) {
    cat(sprintf("%s \\\\\n", paste(par_df[i,], collapse = " & ")))
}


