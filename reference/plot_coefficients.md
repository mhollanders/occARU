# Plot predictor coefficients

Plots coefficients of continuous, categorical, or ordinal predictors on
occupancy or detection submodels.

## Usage

``` r
plot_coefficients(
  fit,
  submodel = c("detection", "occupancy"),
  component = c("site", "survey"),
  type = c("continuous", "categorical", "ordinal"),
  level = c("species", "mean"),
  facet_by = c("predictor", "species"),
  species = NULL,
  restricted = TRUE,
  ordinal_categories = FALSE,
  ...
)
```

## Arguments

- fit:

  A fitted model object from
  [`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md).

- submodel:

  `character`. Predictors of submodel to plot. One of `"detection"`
  (default) or `"occupancy"`.

- component:

  `character`. Whether to plot `"site"` (default) or `"survey"`
  predictors. If `"survey"`, `submodel` must be `"detection"`.

- type:

  `character`. Type of predictors to plot. One of `"continuous"`
  (default), `"categorical"`, or `"ordinal"`.

- level:

  `character`. For multi-species models, whether to plot
  species-specific (`"species"`, default) or mean coefficients
  (`"mean"`).

- facet_by:

  `character`. Whether to use
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)
  or
  [`ggh4x::facet_grid2()`](https://teunbrand.github.io/ggh4x/reference/facet_grid2.html)
  to facet by `"predictor"` (default) or `"species"`. Only used if
  `level` is `"species"`.

- species:

  `character`. Vector of species to plot. If `NULL` (default), all
  species are plotted. Must be one of `attr(occARU_data, "species")`.

- restricted:

  `logical`. If `TRUE` (default), plots coefficients with orthogonal
  projection of the detection random site or survey effects, e.g.,
  \\\boldsymbol{\iota}(\boldsymbol{I} - \boldsymbol{P\_{X_2}})\\, where
  \\\boldsymbol{I} - \boldsymbol{P\_{X_2}}\\ is the orthogonal
  complement of the column space of the site or survey design matrix. If
  `FALSE`, recovers coefficients without orthogonal projection,
  \\\boldsymbol{\beta} - \boldsymbol{X_2}^+ \boldsymbol{\iota}\\, where
  \\\boldsymbol{X_2}^+\\ is the pseudo-inverse of the design matrix.
  Only used for site predictors if `submodel` is `"detection"`, or if
  survey random effects were also projected with `project_kappa = TRUE`
  in
  [`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md).

- ordinal_categories:

  `logical`. If `FALSE` (default), plots coefficients associated with
  maximum category (full effect). If `TRUE`, plots realised coefficient
  associated with each ordered category, where the first is used as the
  reference.

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
[`plot_sites()`](https://mhollanders.github.io/occARU/reference/plot_sites.md),
[`plot_surveys()`](https://mhollanders.github.io/occARU/reference/plot_surveys.md),
[`plot_correlations()`](https://mhollanders.github.io/occARU/reference/plot_correlations.md),
[`plot_partitions()`](https://mhollanders.github.io/occARU/reference/plot_partitions.md),
[`plot_realised()`](https://mhollanders.github.io/occARU/reference/plot_realised.md)
