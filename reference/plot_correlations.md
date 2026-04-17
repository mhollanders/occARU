# Plot interspecific correlations

Plot pairwise species correlations for different model components.For
occupancy, interspecific correlations are only estimated for responses
to site-level predictors. For detection, interspecific correlations are
potentially estimated for responses to site-level predictors,
survey-level predictors, random site effects, random survey effects, and
potentially Poisson OLREs.

## Usage

``` r
plot_correlations(
  fit,
  submodel = c("detection", "occupancy"),
  species = NULL,
  ...
)
```

## Arguments

- fit:

  A fitted model object from
  [`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md).

- submodel:

  `character`. Correlations of submodel to plot. One of `"detection"`
  (default) or `"occupancy"`.

- species:

  `character`. Vector of species to plot. If `NULL` (default), all
  species are plotted. Must be one of `attr(occARU_data, "species")`.

- ...:

  Additional arguments passed to
  [`ggdist::stat_pointinterval()`](https://mjskay.github.io/ggdist/reference/stat_pointinterval.html).

## Value

A `ggplot` object with occARU-specific attributes attached:

- `plot_data`:

  The tibble used to produce the plot.

## See also

[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md),
[`plot_intercepts()`](https://mhollanders.github.io/occARU/reference/plot_intercepts.md),
[`plot_coefficients()`](https://mhollanders.github.io/occARU/reference/plot_coefficients.md),
[`plot_sites()`](https://mhollanders.github.io/occARU/reference/plot_sites.md),
[`plot_surveys()`](https://mhollanders.github.io/occARU/reference/plot_surveys.md),
[`plot_partitions()`](https://mhollanders.github.io/occARU/reference/plot_partitions.md),
[`plot_realised()`](https://mhollanders.github.io/occARU/reference/plot_realised.md)
