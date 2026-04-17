# Plot variance partitions

Plot variance partitions of the different model components for occupancy
and detection submodels.

## Usage

``` r
plot_partitions(fit, scales = FALSE, ...)
```

## Arguments

- fit:

  A fitted model object from
  [`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md).

- scales:

  `logical`. If `FALSE` (default), plots variance simplex partitions
  \\\boldsymbol{\phi}\\. If `TRUE`, produces component scales by
  plotting \\\sqrt{W \cdot \boldsymbol{\phi}}\\, where \\W\\ are
  variances of linear predictors. Useful for sparse simplexes, where few
  components account for most of the variance.

- ...:

  Additional arguments passed to
  [`ggdist::stat_pointinterval()`](https://mjskay.github.io/ggdist/reference/stat_pointinterval.html).

## Value

A `ggplot` object with occARU-specific attributes attached:

- `plot_data`:

  The tibble used to produce the plot.

## Details

The occARU model uses global-local shrinkage priors for the occupancy
and detection submodels, where half-Student-t priors are used for the
variances of both linear predictors which are simplex partitioned via
either Dirichlet or logistic-normal decomposition. Variance
decomposition only occurs when there is more than one model component.
Partitions exist for species-level intercepts, and species-level slopes,
sites effects, survey effects, and Poisson OLREs (with one mean
partition and one for species-level deviations). The species-level
scales for site and survey effects and OLREs are produced by additional
simplex decomposition of the species-level components.

## See also

[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md)
[`plot_intercepts()`](https://mhollanders.github.io/occARU/reference/plot_intercepts.md),
[`plot_coefficients()`](https://mhollanders.github.io/occARU/reference/plot_coefficients.md),
[`plot_sites()`](https://mhollanders.github.io/occARU/reference/plot_sites.md),
[`plot_surveys()`](https://mhollanders.github.io/occARU/reference/plot_surveys.md),
[`plot_correlations()`](https://mhollanders.github.io/occARU/reference/plot_correlations.md),
[`plot_realised()`](https://mhollanders.github.io/occARU/reference/plot_realised.md)
