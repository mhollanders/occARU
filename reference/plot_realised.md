# Plot realised occupancy proportions

Plots species-level proportions of occupied sites by season,
\\\frac{\sum\_{i = 1}^I z\_{kis}}{I}\\, potentially by region.

## Usage

``` r
plot_realised(
  fit,
  by_region = FALSE,
  species = NULL,
  sites = NULL,
  regions = NULL,
  seasons = NULL,
  ...
)
```

## Arguments

- fit:

  A fitted model object from
  [`occARU()`](https://mhollanders.github.io/occARU/reference/occARU.md).

- by_region:

  `logical`. Whether to plot by region, if multiple regions were
  included.

- species:

  `character`. Vector of species to plot. If `NULL` (default), all
  species are plotted. Must be one of `attr(occARU_data, "species")`.

- sites:

  `character`. Vector of sites to use. If `NULL` (default), all sites
  are used. Must be one of `attr(occARU_data, "sites")`.

- regions:

  `character`. Vector of regions to use. If `NULL` (default), all
  regions are used. Must be one of `attr(occARU_data, "regions")`.

- seasons:

  `character`. Vector of seasons to use. If `NULL` (default), all
  seasons are used. Must be one of `attr(occARU_data, "seasons")`.

- ...:

  Additional arguments passed to
  [`ggdist::stat_pointinterval()`](https://mjskay.github.io/ggdist/reference/stat_pointinterval.html).

## Value

A `ggplot` object with occARU-specific attributes attached:

- `plot_data`:

  The tibble used to produce the plot.

## See also

[`occARU()`](https://mhollanders.github.io/occARU/reference/occARU.md),
[`plot_intercepts()`](https://mhollanders.github.io/occARU/reference/plot_intercepts.md),
[`plot_coefficients()`](https://mhollanders.github.io/occARU/reference/plot_coefficients.md),
[`plot_sites()`](https://mhollanders.github.io/occARU/reference/plot_sites.md),
[`plot_surveys()`](https://mhollanders.github.io/occARU/reference/plot_surveys.md),
[`plot_correlations()`](https://mhollanders.github.io/occARU/reference/plot_correlations.md),
[`plot_partitions()`](https://mhollanders.github.io/occARU/reference/plot_partitions.md)
