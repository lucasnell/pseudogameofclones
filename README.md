
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/134013089.svg)](https://zenodo.org/badge/latestdoi/134013089)


# gameofclones <img src="logo.svg" align="right" height="150" />


This R package simulates aphid–parasitoid eco-evolutionary dynamics.


# Installation

To install this package:

```r
remotes::install_github("lucasnell/gameofclones")
```


# Package functions

Each exported function has its own documentation, but here are some basic
descriptions of each:

* `clonal_line`: Create an object to store information about an aphid 
  clonal line. Outputs an `aphid` object, which you can combine using `c()`
  to create a `multiAphid` object.
* `inv_logit` and `logit`: Convenience functions for proportion data
* `rel_res_fitness`: Calculate relative fitness for resistance vs susceptible 
  aphids using calculations from 
  [Ives et al. (2020)](https://doi.org/10.1038/s41559-020-1155-0).
* `restart_experiment`: Restart simulations from `sim_experiments` or
  `restart_experiment`. This function's main purpose is to allow for very
  flexible perturbations. Outputs a `cloneSimsRestart` object.
* `rm_tibs`: Convert tibbles in a `cloneSims` or `cloneSimsRestart` object
  into normal data frames. 
* `sim_experiments`: This is the main function for simulations.
  It defaults to simulating experiments where we have two habitat patches, 
  one with and one without parasitoids.
  Outputs a `cloneSims` object.


# Usage

With all the defaults, you could simulate a resistant and susceptible clone
(with a cost for resistance) in two patches (only one with parasitoids)
for 250 days with the following:

```r
# set initial densities:
adults32 <- cbind(c(0,0,0,0,32), rep(0, 5))

line_s <- clonal_line("susceptible", adults32)
line_r <- clonal_line("resistant", adults32,
                      resistant = TRUE,
                      surv_juv_apterous = "low",
                      surv_adult_apterous = "low",
                      repro_apterous = "low")
sims <- sim_experiments(c(line_s, line_r))
print(sims)

```
