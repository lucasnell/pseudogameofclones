.onLoad <- function(libname, pkgname) {
    modules <- paste0("stan_fit4", names(stanmodels), "_mod")
    if (length(stanmodels) > 0) {## for testing purposes
        for (m in modules) {
            loadModule(m, what = TRUE)
        }
    }
  # This allows it to use all 4 cores on the machine:
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  # ggplot theme:
  theme_set(theme_classic() %+replace%
                theme(strip.background = element_blank(),
                      strip.text = element_text(size = 11),
                      legend.background = element_blank(),
                      plot.title = element_text(size = 14, hjust = 0)))
}
