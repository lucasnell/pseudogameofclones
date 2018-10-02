
#' Leslie matrices for fast and slow aphid populations.
#'
#' @format A list with two 27 x 27 matrices.
#'
#' @source \url{https://onlinelibrary.wiley.com/doi/abs/10.1890/13-1933.1}
#'
#'
"leslie"



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




#' Information for simulating plant death.
#'
#' @format A list with two fields, each containing a single number:
#' * `days_to_max`: An integer specifying the number of days after which the plant
#'     begins to die.
#' * `after_max_mortality`: The per-day mortality to add to the aphid populations on
#'     a plant after the plant begins to die.
#'
"plant_death"



#' Priors for use in estimating population parameters.
#'
#' @format A numeric vector of length 12. Each number is a hyperparameter value.
#'     See `vignette("01-choose_priors")` for more information.
#'
"theta"

