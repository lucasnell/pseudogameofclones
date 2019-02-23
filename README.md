
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

# clonewars <img src="logo.svg" align="right" />

This package contains code for the analysis and simulation of aphid
population-growth time series data. These data are used to help
understand the dynamics of clonal evolution in pea aphids.

## Project description

The aim for this project is to understand the mechanisms underlying
clonal evolution of pea aphid populations. Specifically, we’re
interested in 3 main questions regarding outcomes of evolution in this
system:

1.  How predictable are outcomes?
2.  Do traits that differ among clones help predict outcomes?
3.  How might complex patterns among multiple outcomes emerge?

To answer these questions, we’re combining detailed population assays
(growth, dispersal, etc.), computer simulations, and experiments.

## Assay description

The current data are from assays where we started two adults of the same
line on a fava bean plant that’s been cut down to 3 sets of leaves. We
then counted adult and juvenile aphids on the plant every day, keeping
track of where on the plant the aphids were and if they were not on the
plant at all. We stopped when the total number of aphids dropped for
three consecutive days or when aphid counts reached \< 80% of the
maximum number of aphids.

## Vignette descriptions

All vignettes are under the `Technical details` menu above.

  - `Model description`: Description of the model used to estimate
    growth rates and density dependences.
  - `Choosing priors`: Details of how priors were chosen for of the
    model used to estimate growth rates and density dependences.
  - `Estimating aphid population parameters`: The actual code deploying
    the model used to estimate growth rates and density dependences.
  - `Estimating dispersal`: Description of how I estimated dispersal.
  - `Plant death`: Description of how I estimated effects of plant
    death.
  - `Simulation description`: Description of how I simulated aphid
    clonal evolution.
  - `Simulating cages`: Code deploying my simulations of aphid clonal
    evolution.
