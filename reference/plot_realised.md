# Plot realised occupancy proportions

Plots species-level proportions of occupied sites, \\\frac{\sum\_{i =
1}^I z\_{is}}{I}\\.

## Usage

``` r
plot_realised(fit, species = NULL, sites = NULL, ...)
```

## Arguments

- fit:

  A fitted model object from
  [`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md).

- species:

  `character`. Vector of species to plot. If `NULL` (default), all
  species are plotted. Must be one of `attr(occARU_data, "species")`.

- sites:

  `character`. Vector of sites to use. If `NULL` (default), all sites
  are used. Must be one of `attr(occARU_data, "sites")`.

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
[`plot_correlations()`](https://mhollanders.github.io/occARU/reference/plot_correlations.md),
[`plot_partitions()`](https://mhollanders.github.io/occARU/reference/plot_partitions.md)
