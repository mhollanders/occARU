# Plot temporal detection rates

Plots species-level temporal detection rates combining survey predictors
(if included via
[`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md))
and survey random effects (`kappa`).

## Usage

``` r
plot_surveys(
  fit,
  species = NULL,
  surveys = NULL,
  intercepts = TRUE,
  back_transform = TRUE,
  include_predictors = TRUE,
  restricted = TRUE,
  ndraws = NULL,
  seed = NULL,
  palette = "YlGn",
  ...
)
```

## Arguments

- fit:

  A fitted model object from
  [`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md).

- species:

  `character`. Vector of species to plot. If `NULL` (default), all
  species are plotted. Must be one of `attr(occARU_data, "species")`.

- surveys:

  `character`. Vector of survey dates to plot. If `NULL` (default), all
  surveys are plotted. Must be one of `attr(occARU_data, "surveys")`.

- intercepts:

  `logical`. If `TRUE` (default), species-level baseline log detection
  rates are added to the survey effects. If `FALSE`, only the temporal
  deviations are plotted on the log scale.

- back_transform:

  `logical`. Only used when `intercepts = TRUE`. If `TRUE` (default),
  log detection rates are back-transformed to the natural scale via
  [`exp()`](https://rdrr.io/r/base/Log.html). If `FALSE`, values are
  left on the log scale.

- include_predictors:

  `logical`. If `TRUE` (default), includes predictors in the survey
  effects, if included. If `FALSE`, only plots the random effects.

- restricted:

  `logical`. If `TRUE` (default), when `include_predictors` is `FALSE`
  and the model was fit with `project_kappa = TRUE`, plots random survey
  effects with orthogonal projection, i.e., \\(\boldsymbol{I} -
  \boldsymbol{P\_{X_3}}) \boldsymbol{\kappa}\\, where \\\boldsymbol{I} -
  \boldsymbol{P\_{X_3}}\\ is the orthogonal complement of the column
  space of the site-averaged survey design matrix. If `FALSE`, plots
  random effects without orthogonal projection, i.e.,
  \\\boldsymbol{\kappa}\\ only. Has no effect when `include_predictors`
  is `TRUE` as the linear predictor is unaffected by orthogonal
  projection.

- ndraws:

  Positive integer. Number of draws to use for plotting, passed to
  [`tidybayes::spread_rvars()`](https://mjskay.github.io/tidybayes/reference/spread_rvars.html).
  Default: `NULL` (uses all draws).

- seed:

  Positive numeric. Seed to use when subsampling draws when `ndraws` is
  not `NULL`. Default: random integer.

- palette:

  `character`. Colour palette to be passed to
  [`ggplot2::scale_fill_brewer()`](https://ggplot2.tidyverse.org/reference/scale_brewer.html).
  Default: `"YlGn"`.

- ...:

  Additional arguments passed to
  [`ggdist::stat_lineribbon()`](https://mjskay.github.io/ggdist/reference/stat_lineribbon.html)
  such as `.width` and `point_interval`.

## Value

A `ggplot` object with occARU-specific attributes attached:

- `plot_data`:

  The tibble used to produce the plot.

## See also

[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md),
[`plot_intercepts()`](https://mhollanders.github.io/occARU/reference/plot_intercepts.md),
[`plot_coefficients()`](https://mhollanders.github.io/occARU/reference/plot_coefficients.md),
[`plot_sites()`](https://mhollanders.github.io/occARU/reference/plot_sites.md),
[`plot_correlations()`](https://mhollanders.github.io/occARU/reference/plot_correlations.md),
[`plot_partitions()`](https://mhollanders.github.io/occARU/reference/plot_partitions.md),
[`plot_realised()`](https://mhollanders.github.io/occARU/reference/plot_realised.md)
