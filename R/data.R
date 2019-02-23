

#' Data to generate dispersal predictions for a given line and log(number of aphids).
#'
#' @format An 8 x 3 data frame with columns for intercept (`inter`) and slope (`slope`)
#' of the Poisson GLM of `dispersal ~ log(N)` for each aphid line.
#' (It was actually done using a GLMM, but for purposes of prediction, you can think
#' of it as separate GLMs.)
#'
#'
#'@examples
#' The following function would make predictions:
#' ```
#' predict_disp <- function(X, aphid_line) {
#'     disp_estimates_ <- clonewars::disp_estimates
#'     m <- map_dbl(aphid_line, ~ filter(disp_estimates_, line == .x)[["slope"]])
#'     b <- map_dbl(aphid_line, ~ filter(disp_estimates_, line == .x)[["inter"]])
#'     y <- exp(m * X + b)
#'     return(y)
#' }
#' ```
#'
#'
"disp_estimates"

