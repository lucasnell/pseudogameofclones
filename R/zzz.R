.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
  # This allows it to use all 4 cores on the machine:
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
}
