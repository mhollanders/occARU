# Fit the occARU Stan model

Fits a Bayesian multispecies occupancy model with count observation
model to data prepared by
[`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md)
in Stan. Requires CmdStan \>= 2.36.0, which can be installed with
[`setup_occARU()`](https://mhollanders.github.io/occARU/reference/setup_occARU.md).

## Usage

``` r
fit_model(
  data,
  stan_file = NULL,
  spatial = c("gp", "mvn", "none"),
  temporal = c("gp", "mvn", "none"),
  periodic_gp = FALSE,
  period = NULL,
  species_length_scales = FALSE,
  project_kappa = TRUE,
  overdispersion = c("none", "nb", "olre"),
  variance_decomposition = c("dirichlet", "logistic-normal"),
  latent = TRUE,
  loo_draws = 100,
  ppc = c("Q", "y", "both", "none"),
  prior = set_priors(verbose = FALSE),
  init = "pathfinder",
  pathfinder_args = list(),
  threads = 1,
  grainsize = 1,
  ...
)
```

## Arguments

- data:

  A `occARU_data` object produced by
  [`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md).

- stan_file:

  Character. Path to a custom Stan file. If `NULL` (default), uses the
  built-in multispecies occARU model, or the single species version if
  only one species is included. Intended for advanced users who have
  modified the Stan program; note that custom models will likely require
  corresponding changes to the output of
  [`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md).

- spatial:

  Character. Structure of site-level random effects. `"gp"` (default)
  fits a hierarchical multi-species spatial Gaussian process with
  exponentiated quadratic kernel, which is the recommended option.
  `"mvn"` fits an unstructured multivariate normal, and `"none"` omits
  site-level random effects entirely.

- temporal:

  Character. Structure of survey-level random effects. `"gp"` (default)
  fits a hierarchical multi-species temporal Gaussian process with
  exponentiated quadratic kernel, which is the recommended option.
  `"mvn"` fits an unstructured multivariate normal, and `"none"` omits
  survey-level random effects entirely.

- periodic_gp:

  Logical. If `TRUE`, a periodic kernel is added to the temporal GP
  kernel. Only used when `temporal = "gp"`. Default: `FALSE`.

- period:

  Positive numeric. Period length in survey units (i.e. number of survey
  periods per cycle). Only used when `temporal = "gp"` and
  `periodic_gp = TRUE`. Defaults to `365 / survey_length`, corresponding
  to an annual cycle. For example, with `survey_length = 7` the default
  is `period = 52.1`. Override if your data span a different temporal
  cycle.

- species_length_scales:

  Logical. If `TRUE`, species-specific GP length scales are estimated
  for each kernel, each drawn independently from the shared length scale
  priors. Note that enabling this requires additional Cholesky
  decompositions per species and GP per iteration, which can
  substantially increase sampling time. If `FALSE` (default), only one
  Cholesky decomposition is performed per GP. Only used when multiple
  species are included, and `spatial = "gp"` or `temporal = "gp"`.

- project_kappa:

  Logical. If `TRUE` (default), uses orthogonal projection for random
  survey effects using the site-averaged survey predictor design matrix.
  Ignored when no survey predictors are provided.

- overdispersion:

  Character. Overdispersion model for the observation process. One of
  `"none"` (Poisson, default), `"nb"` (negative binomial), or `"olre"`
  (correlated observation-level random effects).

- variance_decomposition:

  Character. Prior for variance partitions. One of `"dirichlet"`
  (default) or `"logistic-normal"`.

- latent:

  Logical. If `TRUE` (default), latent occupancy states `z` are
  recovered for each species using the forward-backward sampling
  algorithm.

- loo_draws:

  Non-negative integer. Number of Monte Carlo draws for marginal
  log-likelihood estimation of site-level random effects and/or Poisson
  OLREs for PSI-LOO-CV via
  [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html). Default:
  `100`, which produces an additional `[S, I]` matrix `log_lik2` by
  marginalising over site effects (and OLRE residuals if applicable) via
  Monte Carlo integration. `log_lik2` is recommended over `log_lik` for
  PSIS-LOO-CV as it produces better Pareto-k diagnostics. Set to `0` to
  disable, returning only `log_lik`. Only used when `spatial` is not
  `"none"` or `overdispersion = "olre"`.

- ppc:

  Character. Posterior predictive checks to compute. One of `"Q"`
  (default), `"y"`, `"both"`, or `"none"`. `"y"` returns the full
  `[I, J, S]` prediction array (`yrep`); `"Q"` returns only aggregated
  counts `[I, S]` (`Qrep`). For large datasets, `"Q"` or `"none"` can
  substantially reduce memory usage and sampling time.

- prior:

  An `occARU_priors` object from
  [`set_priors()`](https://mhollanders.github.io/occARU/reference/set_priors.md).
  If omitted, default priors are used.

- init:

  Character, numeric, or list. Initialisation strategy passed to
  [cmdstanr::CmdStanModel](https://mc-stan.org/cmdstanr/reference/CmdStanModel.html)`$sample()`.
  One of:

  `"pathfinder"`

  :   Default. Use pathfinder to generate initial values (see
      [cmdstanr::CmdStanModel](https://mc-stan.org/cmdstanr/reference/CmdStanModel.html)`$pathfinder()`).
      Recommended for complex models as it can substantially reduce
      warmup time and improve convergence.

  A numeric scalar

  :   Initialise all parameters uniformly in \\\[-\\`init`\\,\\
      `init`\\\]\\.

  A list

  :   Custom initial values passed directly to
      [cmdstanr::CmdStanModel](https://mc-stan.org/cmdstanr/reference/CmdStanModel.html)`$sample()`.

- pathfinder_args:

  Named list of additional arguments passed to
  [cmdstanr::CmdStanModel](https://mc-stan.org/cmdstanr/reference/CmdStanModel.html)`$pathfinder()`
  when `init = "pathfinder"`. Overrides defaults (`refresh = 0`,
  `sig_figs = 14`, `init = 0.1`, `num_paths = chains`,
  `num_threads = chains`). Default:
  [`list()`](https://rdrr.io/r/base/list.html).

- threads:

  Positive integer. Number of threads for within-chain parallelisation
  via `reduce_sum()`. Default: `1` (no parallelisation). The total
  number of threads used is `threads * chains`, so for optimal
  performance set `threads = floor(available_cores / chains)`. For
  example, 8 cores with 4 chains gives `threads = 2`.

- grainsize:

  Positive integer. Chunk size (number of sites) for within-chain
  parallelisation via `reduce_sum()`. Only used when `threads > 1`.
  Default: `1`, which lets Stan automatically determine the optimal
  chunk size. Increase if you have many sites and want to reduce
  parallelisation overhead. See the [Stan User's
  Guide](https://mc-stan.org/docs/stan-users-guide/parallelization.html#reduce-sum-grainsize)
  for details on tuning grainsize.

- ...:

  Additional arguments passed to
  [cmdstanr::CmdStanModel](https://mc-stan.org/cmdstanr/reference/CmdStanModel.html)`$pathfinder()`.
  Uses parallel chains by default. All other sampling arguments use Stan
  defaults.

## Value

A [CmdStanFit](https://mc-stan.org/cmdstanr/reference/CmdStanFit.html)
object with occARU-specific attributes attached:

- `stan_data`:

  The full Stan data list passed to the model, including prior
  hyperparameters.

- `occARU_data`:

  The original `occARU_data` object from
  [`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md).

## See also

[`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md),
[`set_priors()`](https://mhollanders.github.io/occARU/reference/set_priors.md),
[`setup_occARU()`](https://mhollanders.github.io/occARU/reference/setup_occARU.md),
[cmdstanr::CmdStanMCMC](https://mc-stan.org/cmdstanr/reference/CmdStanMCMC.html)
for methods on the fitted model object. The statistical model is
described in
[`vignette("model", package = "occARU")`](https://mhollanders.github.io/occARU/articles/model.md).
