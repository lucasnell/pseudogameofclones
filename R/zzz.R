.onLoad <- function(libname, pkgname) {
    # ggplot theme:
    ggplot2::theme_set(ggplot2::`%+replace%`(ggplot2::theme_classic(),
                       ggplot2::theme(strip.background = ggplot2::element_blank(),
                        strip.text = ggplot2::element_text(size = 11),
                        legend.background = ggplot2::element_blank(),
                        plot.title = ggplot2::element_text(size = 14, hjust = 0))))
}
