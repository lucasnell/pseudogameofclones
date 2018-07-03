

#' Visually compare normal, student_t, cauchy, laplace, and product_normal.
#'
#' From \url{http://mc-stan.org/rstanarm/reference/priors.html#examples}
#'
#' @param scale Scale for all
#' @param df_t Degrees of freedom for `student_t`.
#' @param xlim X limit for the plot.
#'
#' @export
#'
compare_priors <- function(scale = 1, df_t = 2, xlim = c(-10, 10)) {
    dt_loc_scale <- function(x, df, location, scale) {
        1/scale * dt((x - location)/scale, df)
    }
    dlaplace <- function(x, location, scale) {
        0.5 / scale * exp(-abs(x - location) / scale)
    }
    dproduct_normal <- function(x, scale) {
        besselK(abs(x) / scale ^ 2, nu = 0) / (scale ^ 2 * pi)
    }
    stat_dist <- function(dist, ...) {
        ggplot2::stat_function(ggplot2::aes_(color = dist), ...)
    }
    ggplot2::ggplot(data.frame(x = xlim), ggplot2::aes(x)) +
        stat_dist("normal", size = 0.75, fun = dnorm,
                  args = list(mean = 0, sd = scale)) +
        stat_dist("student_t", size = 0.75, fun = dt_loc_scale,
                  args = list(df = df_t, location = 0, scale = scale)) +
        stat_dist("cauchy", size = 0.75, linetype = 2, fun = dcauchy,
                  args = list(location = 0, scale = scale)) +
        stat_dist("laplace", size = 0.75, linetype = 2, fun = dlaplace,
                  args = list(location = 0, scale = scale)) +
        stat_dist("product_normal", size = 0.75, linetype = 2, fun = dproduct_normal,
                  args = list(scale = 1))
}


