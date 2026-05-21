# Fit the occARU model

Fits the occARU model to data prepared by
[`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md)
in Stan. Requires CmdStan \>= 2.36.0, which can be installed with
[`setup_occARU()`](https://mhollanders.github.io/occARU/reference/setup_occARU.md).
Automatically determines number of species in the supplied data and
whether a single or multiseason model is required.

## Usage

``` r
occARU(
  data,
  prior = set_priors(verbose = FALSE),
  stan_file = NULL,
  random = list(site = gp(), survey = gp()),
  overdispersion = c("none", "nb", "olre"),
  variance_decomposition = c("dirichlet", "logistic-normal"),
  ppc = c("Q", "y", "both", "none"),
  latent = TRUE,
  loo_draws = 100L,
  init = 0.1,
  pathfinder_args = list(),
  threads = 1L,
  grainsize = 1L,
  ...
)
```

## Arguments

- data:

  A `occARU_data` object produced by
  [`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md).

- prior:

  An `occARU_priors` object produced by
  [`set_priors()`](https://mhollanders.github.io/occARU/reference/set_priors.md).
  If omitted, default priors are used.

- stan_file:

  `character`. Path to a custom Stan file. If `NULL` (default), uses the
  built-in occARU model. Intended for advanced users who have modified
  the Stan programs; note that custom models will likely require
  corresponding changes to the output of
  [`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md).

- random:

  A named `list` specifying random effect structures for site and survey
  effects on detection. For multiseason models, optionally accepts
  season effects on detection (`season`) and site (`site_occ`) and
  season (`season_occ`) effects on occupancy. Must be one of:

  - [`gp()`](https://mhollanders.github.io/occARU/reference/gp.md),
    which fits a (multispecies) Gaussian process. Default for site and
    survey effects. Also the default for season effects on both
    detection and occupancy and occupancy site effects in dynamic
    models.

  - `"mvn"`, which fits an unstructured (multispecies) normal.

  - `"none"`, which omits random effects entirely.

- overdispersion:

  `character`. Overdispersion model for the observation process. One of
  `"none"` (Poisson, default), `"nb"` (negative binomial), or `"olre"`
  (multispecies observation-level random effects).

- variance_decomposition:

  `character`. Prior for variance partitions. One of `"dirichlet"`
  (default) or `"logistic-normal"`.

- ppc:

  `character`. Posterior predictive checks to compute. One of `"Q"`
  (default), `"y"`, `"both"`, or `"none"`. `"y"` returns the full
  `[K, I, J, S]` prediction array (`yrep`); `"Q"` returns only
  aggregated counts `[K, I, S]` (`Qrep`). For large datasets, `"Q"` or
  `"none"` can substantially reduce memory usage and sampling time.

- latent:

  `logical`. If `TRUE` (default), latent occupancy states `z` are
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
  disable, returning only `log_lik`. Only used with site random effects
  or `overdispersion = "olre"`. Note: even high values (10K) seem
  inadequate for OLREs.

- init:

  `numeric`, `character`, or `list`. Initialisation strategy passed to
  [cmdstanr::CmdStanModel](https://mc-stan.org/cmdstanr/reference/CmdStanModel.html)`$sample()`.
  One of:

  - A numeric scalar (default). Initialises all parameters uniformly in
    `[-init, init]`. Default: 0.1.

  - `"pathfinder"`. Use Pathfinder to generate initial values (see
    [cmdstanr::CmdStanModel](https://mc-stan.org/cmdstanr/reference/CmdStanModel.html)`$pathfinder()`).

  - A named list with custom initial values passed directly to
    [cmdstanr::CmdStanModel](https://mc-stan.org/cmdstanr/reference/CmdStanModel.html)`$sample()`.

- pathfinder_args:

  Named list of additional arguments passed to
  [cmdstanr::CmdStanModel](https://mc-stan.org/cmdstanr/reference/CmdStanModel.html)`$pathfinder()`
  when `init = "pathfinder"`. Overrides occARU-specific defaults
  (`refresh = 0`, `init = 0.1`, `sig_figs = 14`, `num_threads = chains`,
  `num_paths = chains`, `max_lbfgs_iters = 200`,
  `psis_resample = FALSE`).

- threads:

  Positive integer. Number of threads for within-chain parallelisation
  via `reduce_sum()`. Default: `1` (no parallelisation). The total
  number of threads used is `threads * chains`, so for optimal
  performance set `threads = floor(available_cores / chains)`. For
  example, 8 cores with 4 chains gives `threads = 2`.

- grainsize:

  Positive integer. Chunk size for within-chain parallelisation when
  `threads > 1`. For data with multiple regions, chunks are number of
  regions; otherwise it is number of sites in each slice. Default: `1`,
  which lets Stan automatically determine the optimal chunk size.
  Increase if you have many sites or regions and want to reduce
  parallelisation overhead. See the [Stan User's
  Guide](https://mc-stan.org/docs/stan-users-guide/parallelization.html#reduce-sum-grainsize)
  for details on tuning grainsize.

- ...:

  Additional arguments passed to
  [cmdstanr::CmdStanModel](https://mc-stan.org/cmdstanr/reference/CmdStanModel.html)`$sample()`.
  Uses parallel chains by default. All other sampling arguments use Stan
  defaults.

## Value

An `occARU_fit` object, which extends the
[CmdStanFit](https://mc-stan.org/cmdstanr/reference/CmdStanFit.html)
class with occARU-specific attributes:

- `stan_data`:

  The full Stan data list passed to the model, including prior
  hyperparameters.

- `occARU_data`:

  The original `occARU_data` object from
  [`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md).

## See also

[`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md),
[`set_priors()`](https://mhollanders.github.io/occARU/reference/set_priors.md),
[`gp()`](https://mhollanders.github.io/occARU/reference/gp.md),
[`setup_occARU()`](https://mhollanders.github.io/occARU/reference/setup_occARU.md),
[`pp_check.occARU_fit()`](https://mhollanders.github.io/occARU/reference/pp_check.occARU_fit.md),
[cmdstanr::CmdStanMCMC](https://mc-stan.org/cmdstanr/reference/CmdStanMCMC.html)
for methods on the fitted model object. The statistical model is
described in
[`vignette("model")`](https://mhollanders.github.io/occARU/articles/model.md).
