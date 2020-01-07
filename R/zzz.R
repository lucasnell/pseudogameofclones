.onLoad <- function(libname, pkgname) {
    # ggplot theme:
    theme_set(theme_classic() %+replace%
                  theme(strip.background = element_blank(),
                        strip.text = element_text(size = 11),
                        legend.background = element_blank(),
                        plot.title = element_text(size = 14, hjust = 0)))
}
