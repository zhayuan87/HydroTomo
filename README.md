
# HydroTomo

<!-- badges: start -->
[![R-CMD-check](https://github.com/zhayuan87/HydroTomo/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zhayuan87/HydroTomo/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of HydroTomo is to develop an R package for running Hydraulic tomography test. This integrates pumping test forward simulation and inverse simulation.

## Parallel computing

Both Monte Carlo (ensemble-based) and adjoint-based inversions are computationally demanding because they require solving the forward problem many times — once per ensemble member or once per observation point. HydroTomo addresses this through embarrassingly parallel execution using the `foreach` + `doParallel` framework.

Two parallelization strategies are employed:

1. **Ensemble parallelization** — Forward simulations for all ensemble members are distributed across available CPU cores. Each core independently runs the forward model for its assigned members, with no inter-process communication beyond the initial parameter distribution and final result collection.

2. **Adjoint parallelization** — Adjoint-state sensitivity computations for all observation points are executed in parallel across cores. Since each adjoint solve is independent of the others, this also achieves near-linear speedup.

The user controls the number of cores via the `ncore` parameter in applicable functions. The realized speedup depends on grid size, ensemble/observation count, solver choice, and hardware configuration.

## Installation

You can install the development version of HydroTomo like so:

``` r
library("devtools")
install_github("zhayuan87/HydroTomo")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(HydroTomo)
## basic example code
```

