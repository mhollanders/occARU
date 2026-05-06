# Plot intercepts

Plot species-specific intercepts for occupancy and detection, by default
back-transformed to the original scale.

## Usage

``` r
plot_intercepts(
  fit,
  back_transform = TRUE,
  by_region = FALSE,
  species = NULL,
  regions = NULL,
  ...
)
```

## Arguments

- fit:

  A fitted model object from
  [`occARU()`](https://mhollanders.github.io/occARU/reference/occARU.md).

- back_transform:

  `logical`. If `TRUE` (default), intercepts are back-transformed to the
  natural scale via `inv_logit()` for occupancy and
  [`exp()`](https://rdrr.io/r/base/Log.html) for detection rates. If
  `FALSE`, values are left on the scale of the link functions (logit for
  occupancy and log for detection).

- by_region:

  `logical`. Whether to plot intercepts by region, if multiple regions
  were included.

- species:

  `character`. Vector of species to plot. If `NULL` (default), all
  species are plotted. Must be one of `attr(occARU_data, "species")`.

- regions:

  `character`. Vector of regions to plot. If `NULL` (default), all
  regions are plotted. Must be one of `attr(occARU_data, "regions")`.

- ...:

  Additional arguments passed to
  [`ggdist::stat_pointinterval()`](https://mjskay.github.io/ggdist/reference/stat_pointinterval.html).

## Value

A `ggplot` object with occARU-specific attributes attached:

- `plot_data`:

  The tibble used to produce the plot.

## See also

[`occARU()`](https://mhollanders.github.io/occARU/reference/occARU.md),
[`plot_coefficients()`](https://mhollanders.github.io/occARU/reference/plot_coefficients.md),
[`plot_sites()`](https://mhollanders.github.io/occARU/reference/plot_sites.md),
[`plot_surveys()`](https://mhollanders.github.io/occARU/reference/plot_surveys.md),
[`plot_correlations()`](https://mhollanders.github.io/occARU/reference/plot_correlations.md),
[`plot_partitions()`](https://mhollanders.github.io/occARU/reference/plot_partitions.md),
[`plot_realised()`](https://mhollanders.github.io/occARU/reference/plot_realised.md)
