

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



#' Estimates of the stan model of population parameters.
#'
#' @format A list with the following items:
#' \describe{
#' \item{names}{A length-8 character vector of names for each aphid line.}
#' \item{R}{A length-8 numeric vector of growth rates for each aphid line.}
#' \item{A}{A length-8 numeric vector of density dependences for each aphid line.}
#' \item{process_error}{SD of process error.}
#' \item{log_zeta_mean}{Mean of log(zeta) distribution.}
#' \item{log_zeta_sd}{SD of log(zeta) distribution.}
#' \item{mu_time}{Mean time among all time series. Used to calculate the effect of
#'     plant deterioration through time.}
#' }
#'
#'
"stan_estimates"

