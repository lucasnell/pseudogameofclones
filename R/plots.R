

#' Visually compare normal, student_t, cauchy, laplace, and product_normal.
#'
#' From \url{http://mc-stan.org/rstanarm/reference/priors.html#examples}
#'
#' @param distributions Vector of which distribution(s) to plot.
#'     Can also be a single function that takes the `location` and `scale` arguments.
#' @param location Location for all distributions
#' @param scale Scale for all
#' @param df_t Degrees of freedom for `"student_t"` distribution. Defaults to 2.
#' @param xlim X limit for the plot.
#'
#' @export
#'
compare_priors <- function(distributions, location = 0, scale = 1, df_t = 2,
                           xlim = c(-10, 10), guide = FALSE) {

    if (inherits(guide, "logical") & guide == TRUE) guide <- "legend"

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

    outp <- ggplot2::ggplot(data.frame(x = xlim), ggplot2::aes(x))


    if (inherits(distributions, "character")) {
        if ("normal" %in% distributions) {
            outp <- outp +
                stat_dist("normal", size = 0.75, fun = dnorm,
                      args = list(mean = location, sd = scale))
        }
        if ("student_t" %in% distributions) {
        outp <- outp +
            stat_dist("student_t", size = 0.75, fun = dt_loc_scale,
                      args = list(df = df_t, location = location, scale = scale))
        }
        if ("cauchy" %in% distributions) {
            outp <- outp +
                stat_dist("cauchy", size = 0.75, fun = dcauchy,
                          args = list(location = location, scale = scale))
        }
        if ("laplace" %in% distributions) {
            outp <- outp +
                stat_dist("laplace", size = 0.75, fun = dlaplace,
                          args = list(location = location, scale = scale))
        }
        if ("product_normal" %in% distributions) {
            outp <- outp +
                stat_dist("product_normal", size = 0.75, fun = dproduct_normal,
                          args = list(scale = 1))
        }
        if ("beta" %in% distributions) {
            mu_ <- location
            var_ <- scale^2
            alpha_ <- ((1 - mu_) / var_ - 1 / mu_) * mu_ ^ 2
            beta_ <- alpha_ * (1 / mu_ - 1)
            outp <- outp +
                stat_dist("beta", size = 0.75, fun = dbeta,
                          args = list(shape1 = alpha_, shape2 = beta_))
        }
    } else if (inherits(distributions, "function")) {
        outp <- outp +
            stat_dist("custom", size = 0.75, fun = distributions,
                      args = list(location = location, scale = scale))
    }

    outp <- outp +
        ggplot2::scale_color_brewer("distribution:", palette = "Dark2", guide = guide) +
        ggplot2::ylab("density") +
        ggplot2::xlab("value") +
        NULL

    return(outp)
}


