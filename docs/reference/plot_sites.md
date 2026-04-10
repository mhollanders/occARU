# Plot site occupancy and detection rates

Plots species-level occupancy and detection rates combining site
predictors (if included via
[`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md))
and site random effects (`iota`).

## Usage

``` r
plot_sites(
  fit,
  species = NULL,
  sites = NULL,
  map = TRUE,
  intercepts = TRUE,
  back_transform = TRUE,
  include_predictors = TRUE,
  restricted = TRUE,
  ndraws = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- fit:

  A fitted model object from
  [`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md).

- species:

  Character vector of species to plot. If `NULL` (default), all species
  are plotted. Must be one of `attr(occARU_data, "species")`.

- sites:

  Character vector of sites to plot. If `NULL` (default), all sites are
  plotted. Must be one of `attr(occARU_data, "sites")`.

- map:

  Logical. If `TRUE` (default), plot site effects summarised with
  posterior medians on a map using UTM coordinates. Requires site
  coordinates to have been supplied to
  [`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md).
  If `FALSE`, or if no coordinates are present, site effects are plotted
  with
  [`ggdist::stat_pointinterval()`](https://mjskay.github.io/ggdist/reference/stat_pointinterval.html).

- intercepts:

  Logical. If `TRUE` (default), species-level baseline log detection
  rates are added to the site effects. If `FALSE`, only the site
  deviations are plotted on the log scale.

- back_transform:

  Logical. Only used when `intercepts = TRUE`. If `TRUE` (default), log
  detection rates are back-transformed to the natural scale via
  [`exp()`](https://rdrr.io/r/base/Log.html). If `FALSE`, values are
  left on the log scale.

- include_predictors:

  Logical. If `TRUE` (default), includes predictors in the site effects,
  if included. If `FALSE`, only plots the random effects.

- restricted:

  Logical. If `TRUE` (default), when `include_predictors` is `FALSE`,
  plots random site effects with orthogonal projection, i.e.,
  \\\boldsymbol{\iota}(\boldsymbol{I} - \boldsymbol{P\_{X_2}})\\, where
  \\\boldsymbol{I} - \boldsymbol{P\_{X_2}}\\ is the orthogonal
  complement of the column space of the site design matrix. If `FALSE`,
  plots random effects without orthogonal projection, i.e.,
  \\\boldsymbol{\iota}\\ only. Has no effect when `include_predictors`
  is `TRUE` as the linear predictor is unaffected by orthogonal
  projection.

- ndraws:

  Number of draws to use for plotting, passed to
  [`tidybayes::spread_rvars()`](https://mjskay.github.io/tidybayes/reference/spread_rvars.html).
  Default: `NULL` (uses all draws).

- seed:

  Positive numeric. Seed to use when subsampling draws when `ndraws` is
  not `NULL`. Default: random integer.

- ...:

  Additional arguments passed to
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)
  when `map = TRUE` and
  [`ggdist::stat_pointinterval()`](https://mjskay.github.io/ggdist/reference/stat_pointinterval.html)
  when `map = FALSE`.

## Value

A `ggplot` object with occARU-specific attributes attached:

- `plot_data`:

  The tibble used to produce the plot.

## Details

If `map = TRUE` and site coordinates are present in the fitted object,
effects are displayed as points sized by the magnitude of the detection
rate on a map; otherwise site effects are plotted as point-intervals.
When `latent = TRUE` was set in
[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md),
sites with median posterior occupancy of 0 are shown as red crosses.
When `latent = FALSE`, detection rates are weighted by occupancy
probability (`inv_logit(logit_psi[s, i])`).

## See also

[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md),
[`plot_intercepts()`](https://mhollanders.github.io/occARU/reference/plot_intercepts.md),
[`plot_coefficients()`](https://mhollanders.github.io/occARU/reference/plot_coefficients.md),
[`plot_surveys()`](https://mhollanders.github.io/occARU/reference/plot_surveys.md),
[`plot_correlations()`](https://mhollanders.github.io/occARU/reference/plot_correlations.md),
[`plot_partitions()`](https://mhollanders.github.io/occARU/reference/plot_partitions.md)
