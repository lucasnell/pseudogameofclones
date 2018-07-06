

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
compare_priors <- function(location = 0, scale = 1, df_t = 2, xlim = c(-10, 10),
                           beta = FALSE) {
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
    outp <- ggplot2::ggplot(data.frame(x = xlim), ggplot2::aes(x)) +
        stat_dist("normal", size = 0.75, fun = dnorm,
                  args = list(mean = location, sd = scale)) +
        stat_dist("student_t", size = 0.75, fun = dt_loc_scale,
                  args = list(df = df_t, location = location, scale = scale)) +
        stat_dist("cauchy", size = 0.75, linetype = 2, fun = dcauchy,
                  args = list(location = location, scale = scale)) +
        stat_dist("laplace", size = 0.75, linetype = 2, fun = dlaplace,
                  args = list(location = location, scale = scale)) +
        # stat_dist("product_normal", size = 0.75, linetype = 2, fun = dproduct_normal,
        #           args = list(scale = 1))
        ggplot2::scale_color_brewer("distribution:", palette = "Dark2") +
        NULL
    if (beta) {
        mu_ <- location
        var_ <- scale^2
        alpha_ <- ((1 - mu_) / var_ - 1 / mu_) * mu_ ^ 2
        beta_ <- alpha_ * (1 / mu_ - 1)
        outp <- outp +
            stat_dist("beta", size = 0.75, linetype = 2, fun = dbeta,
                      args = list(shape1 = alpha_, shape2 = beta_))
    }
    return(outp)
}


