# Plot intercepts

Plot species-specific intercepts for occupancy and detection, by default
back-transformed to the orginal scale.

## Usage

``` r
plot_intercepts(fit, species = NULL, back_transform = TRUE, ...)
```

## Arguments

- fit:

  A fitted model object from
  [`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md).

- species:

  Character vector of species to plot. If `NULL` (default), all species
  are plotted. Must be one of `attr(occARU_data, "species")`.

- back_transform:

  Logical. If `TRUE` (default), intercepts are back-transformed to the
  natural scale via `inv_logit()` for occupancy and
  [`exp()`](https://rdrr.io/r/base/Log.html) for detection rates. If
  `FALSE`, values are left on the scale of the link functions (logit for
  occupancy and log for detection).

- ...:

  Additional arguments passed to
  [`ggdist::stat_pointinterval()`](https://mjskay.github.io/ggdist/reference/stat_pointinterval.html).

## Value

A `ggplot` object with occARU-specific attributes attached:

- `plot_data`:

  The tibble used to produce the plot.

## See also

[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md),
[`plot_coefficients()`](https://mhollanders.github.io/occARU/reference/plot_coefficients.md),
[`plot_sites()`](https://mhollanders.github.io/occARU/reference/plot_sites.md),
[`plot_surveys()`](https://mhollanders.github.io/occARU/reference/plot_surveys.md),
[`plot_correlations()`](https://mhollanders.github.io/occARU/reference/plot_correlations.md),
[`plot_partitions()`](https://mhollanders.github.io/occARU/reference/plot_partitions.md)
