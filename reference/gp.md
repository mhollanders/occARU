# Specify Gaussian process in occARU models

Specify Gaussian process (GP) random effects in occARU models, with a
choice of kernel and whether to estimate species-level length-scales.

## Usage

``` r
gp(
  kernel = c("exp_quad_cov", "matern32", "matern52"),
  species_length_scales = FALSE,
  periodic = FALSE,
  period = NULL
)
```

## Arguments

- kernel:

  `character`. The type of GP kernel to use. Must be one of `"exp_quad"`
  (exponentiatated quadratic), `"matern32"` (Matern 3/2), or
  `"matern52"` (Matern 5/2).

- species_length_scales:

  `logical`. If `TRUE`, species-specific length scales are estimated for
  that GP, each drawn independently from the shared length scale priors.
  Note that enabling this requires additional Cholesky decompositions
  per species and GP, which can substantially increase sampling time. If
  `FALSE` (default), only one length scale is estimated with one
  Cholesky decomposition performed per GP. Only used when multiple
  species are included.

- periodic:

  `logical`. If `TRUE`, a periodic kernel is added to the kernel for
  survey effects. Default: `FALSE`.

- period:

  Positive numeric. Period length in survey units (i.e. number of survey
  periods per cycle). Only used when `periodic = TRUE`. Defaults to
  `365 / survey_length`, corresponding to an annual cycle. For example,
  with `survey_length = 7` the default is `period = 52.1`. Override if
  your data span a different temporal cycle.

## Value

A `gp` object (a named list with class `"gp"`) for use in the `random`
argument of
[`occARU()`](https://mhollanders.github.io/occARU/reference/occARU.md).
Contains the following elements:

- `random`:

  Integer. Always `1L`, indicating a random effect.

- `kernel`:

  Integer. Kernel index (1 = exponentiated quadratic, 2 = Matern 3/2, 3
  = Matern 5/2).

- `species_length_scales`:

  Integer. Whether to estimate species-specific length scales (0 or 1).

- `periodic`:

  Integer. Whether a periodic kernel component is included (0 or 1).

- `period`:

  Numeric. Period length in survey units, or `0` if not periodic.

## See also

[`occARU()`](https://mhollanders.github.io/occARU/reference/occARU.md),
[`set_priors()`](https://mhollanders.github.io/occARU/reference/set_priors.md)
