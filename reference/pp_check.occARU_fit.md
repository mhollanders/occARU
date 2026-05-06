# Posterior predictive checks using rootograms for `occARU_fit` objects

Perform visual posterior predictive checks using rootograms using the
bayesplot package.

## Usage

``` r
# S3 method for class 'occARU_fit'
pp_check(object, level = c("Q", "y"), ndraws = NULL, ...)
```

## Arguments

- object:

  A fitted model object from
  [`occARU()`](https://mhollanders.github.io/occARU/reference/occARU.md).

- level:

  Character. One of `"Q"` or `"y"`. If `"Q"` (default), uses the
  aggregated counts per site and species as data. If `"y"`, uses the raw
  (survey-level) counts as input.

- ndraws:

  Positive integer. Number of draws to use. If `NULL` (default), uses
  all draws.

- ...:

  Additional arguments to be passed to
  [`bayesplot::ppc_rootogram_grouped()`](https://mc-stan.org/bayesplot/reference/PPC-discrete.html).

## Value

A ggplot object with posterior predictive rootograms faceted by species
which can be modified using the ggplot2 package.
