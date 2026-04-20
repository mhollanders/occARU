#' Fit the occARU model
#'
#' Fits a Bayesian multispecies occupancy model with count observation model to
#' data prepared by [make_data()] in Stan. Requires CmdStan >= 2.36.0, which can
#' be installed with [setup_occARU()].
#'
#' @param data A `occARU_data` object produced by [make_data()].
#' @param stan_file `character`. Path to a custom Stan file. If `NULL` (default),
#'   uses the built-in multispecies occARU model, or the single species version
#'   if only one species is included. Intended for advanced users who have
#'   modified the Stan program; note that custom models will likely require
#'   corresponding changes to the output of [make_data()].
#' @param spatial `character`. Structure of site-level random effects. `"gp"`
#'   (default) fits a hierarchical multi-species spatial Gaussian process with
#'   exponentiated quadratic kernel, which is the recommended option. `"mvn"`
#'   fits an unstructured multivariate normal, and `"none"` omits site-level
#'   random effects entirely.
#' @param temporal `character`. Structure of survey-level random effects. `"gp"`
#'   (default) fits a hierarchical multi-species temporal Gaussian process with
#'   exponentiated quadratic kernel, which is the recommended option. `"mvn"`
#'   fits an unstructured multivariate normal, and `"none"` omits survey-level
#'   random effects entirely.
#' @param periodic_gp `logical`. If `TRUE`, a periodic kernel is added to the
#'   temporal GP kernel. Only used when `temporal = "gp"`. Default: `FALSE`.
#' @param period Positive numeric. Period length in survey units (i.e. number
#'   of survey periods per cycle). Only used when `temporal = "gp"` and
#'   `periodic_gp = TRUE`. Defaults to `365 / survey_length`, corresponding to
#'   an annual cycle. For example, with `survey_length = 7` the default is
#'   `period = 52.1`. Override if your data span a different temporal cycle.
#' @param species_length_scales `logical`. If `TRUE`, species-specific GP length
#'   scales are estimated for each kernel, each drawn independently from the
#'   shared length scale priors. Note that enabling this requires additional
#'   Cholesky decompositions per species and GP per iteration, which can
#'   substantially increase sampling time. If `FALSE` (default), only one
#'   Cholesky decomposition is performed per GP. Only used when multiple species
#'   are included, and `spatial = "gp"` or `temporal = "gp"`.
#' @param project_kappa `logical`. If `TRUE` (default), uses orthogonal
#'   projection for random survey effects using the site-averaged survey
#'   predictor design matrix. Ignored when no survey predictors are provided.
#' @param overdispersion `character`. Overdispersion model for the observation
#'   process. One of `"none"` (Poisson, default), `"nb"` (negative binomial),
#'   or `"olre"` (correlated observation-level random effects).
#' @param variance_decomposition `character`. Prior for variance partitions.
#'   One of `"dirichlet"` (default) or `"logistic-normal"`.
#' @param latent `logical`. If `TRUE` (default), latent occupancy states `z`
#'   are recovered for each species using the forward-backward sampling
#'   algorithm.
#' @param loo_draws Non-negative integer. Number of Monte Carlo draws for
#'   marginal log-likelihood estimation of site-level random effects and/or
#'   Poisson OLREs for PSI-LOO-CV via [loo::loo()]. Default: `100`, which
#'   produces an additional `[S, I]` matrix `log_lik2` by marginalising over
#'   site effects (and OLRE residuals if applicable) via Monte Carlo
#'   integration. `log_lik2` is recommended over `log_lik` for PSIS-LOO-CV as it
#'   produces better Pareto-k diagnostics. Set to `0` to disable, returning only
#'   `log_lik`. Only used when `spatial` is not `"none"` or `overdispersion =
#'   "olre"`.
#' @param ppc `character`. Posterior predictive checks to compute. One of `"Q"`
#'   (default), `"y"`, `"both"`, or `"none"`. `"y"` returns the full
#'   `[I, J, S]` prediction array (`yrep`); `"Q"` returns only aggregated
#'   counts `[I, S]` (`Qrep`). For large datasets, `"Q"` or `"none"` can
#'   substantially reduce memory usage and sampling time.
#' @param prior An `occARU_priors` object from [set_priors()]. If omitted,
#'   default priors are used.
#' @param init `character`, `numeric`, or `list`. Initialisation strategy passed
#'   to [cmdstanr::CmdStanModel]`$sample()`. One of:
#'   \describe{
#'     \item{`"pathfinder"`}{Default. Use pathfinder to generate initial values
#'       (see [cmdstanr::CmdStanModel]`$pathfinder()`). Recommended for complex
#'       models as it can substantially reduce warmup time and improve
#'       convergence.}
#'     \item{A numeric scalar}{Initialise all parameters uniformly in
#'       \eqn{[-}\code{init}\eqn{,} \code{init}\eqn{]}.}
#'     \item{A list}{Custom initial values passed directly to
#'       [cmdstanr::CmdStanModel]`$sample()`.}
#'   }
#' @param pathfinder_args Named list of additional arguments passed to
#'   [cmdstanr::CmdStanModel]`$pathfinder()` when `init = "pathfinder"`.
#'   Overrides defaults (`refresh = 0`, `sig_figs = 14`, `init = 0.1`,
#'   `num_paths = chains`, `num_threads = chains`). Default: `list()`.
#' @param threads Positive integer. Number of threads for within-chain
#'   parallelisation via `reduce_sum()`. Default: `1` (no parallelisation).
#'   The total number of threads used is `threads * chains`, so for optimal
#'   performance set `threads = floor(available_cores / chains)`. For example,
#'   8 cores with 4 chains gives `threads = 2`.
#' @param grainsize Positive integer. Chunk size (number of sites) for
#'   within-chain parallelisation via `reduce_sum()`. Only used when
#'   `threads > 1`. Default: `1`, which lets Stan automatically determine the
#'   optimal chunk size. Increase if you have many sites and want to reduce
#'   parallelisation overhead. See the
#'   \href{https://mc-stan.org/docs/stan-users-guide/parallelization.html#reduce-sum-grainsize}{Stan
#'   User's Guide} for details on tuning grainsize.
#' @param ... Additional arguments passed to
#'   [cmdstanr::CmdStanModel]`$pathfinder()`. Uses parallel chains by default.
#'   All other sampling arguments use Stan defaults.
#'
#' @return A \href{https://mc-stan.org/cmdstanr/reference/CmdStanFit.html}{CmdStanFit}
#'   object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`stan_data`}{The full Stan data list passed to the model,
#'       including prior hyperparameters.}
#'     \item{`occARU_data`}{The original `occARU_data` object from
#'       [make_data()].}
#'   }
#' @seealso [make_data()], [set_priors()], [setup_occARU()],
#'   [cmdstanr::CmdStanMCMC] for methods on the fitted model object.
#'   The statistical model is described in
#'   `vignette("model", package = "occARU")`.
#' @export
fit_model <- function(
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
  loo_draws = 100L,
  ppc = c("Q", "y", "both", "none"),
  prior = set_priors(verbose = FALSE),
  init = "pathfinder",
  pathfinder_args = list(),
  threads = 1L,
  grainsize = 1L,
  ...
) {
  # input checks
  if (!inherits(data, "occARU_data")) {
    cli::cli_abort(
      "{.arg data} must be a {.cls occARU_data} object from {.fun make_data}."
    )
  }
  if (!inherits(prior, "occARU_priors")) {
    cli::cli_abort(
      "{.arg prior} must be an {.cls occARU_priors} object from \\
      {.fun set_priors}."
    )
  }

  # check CmdStan version
  cmdstan_ver <- cmdstanr::cmdstan_version()
  if (package_version(cmdstan_ver) < package_version("2.36.0")) {
    cli::cli_abort(
      "occARU requires CmdStan >= 2.36.0. You have version {cmdstan_ver}. \\
       Run {.fun setup_occARU} to install the latest version."
    )
  }

  # match arguments
  spatial <- match.arg(spatial)
  temporal <- match.arg(temporal)
  overdispersion <- match.arg(overdispersion)
  variance_decomposition <- match.arg(variance_decomposition)
  ppc <- match.arg(ppc)

  # scalar argument checks
  if (!rlang::is_integerish(loo_draws) || loo_draws < 0) {
    cli::cli_abort("{.arg loo_draws} must be a non-negative integer.")
  } else if (!(spatial != "none" || overdispersion == "olre")) {
    cli::cli_warn(
      "{.arg loo_draws} only used when site effects or OLREs are included."
    )
  }
  if (!rlang::is_integerish(threads) || threads < 1) {
    cli::cli_abort("{.arg threads} must be a positive integer.")
  } else if (!rlang::is_integerish(grainsize) || grainsize < 0) {
    cli::cli_abort("{.arg grainsize} must be a positive integer.")
  } else if (grainsize >= data$I) {
    (cli::cli_abort(
      "{.arg grainsize} should be less than the number of sites."
    ))
  }

  # spatial GP check — warn and downgrade if no coordinates
  if (spatial == "gp" && all(data$XY == 0)) {
    cli::cli_abort(
      'No site coordinates found in {.arg data}. Set {.arg spatial = "mvn"} \\
      or supply {.arg latitude} and {.arg longitude} in {.fun make_data} to \\
      enable the spatial GP.'
    )
  }

  # periodic GP checks
  if (periodic_gp) {
    if (temporal != "gp") {
      cli::cli_abort(
        "{.arg periodic_gp = TRUE} requires {.arg temporal = 'gp'}."
      )
    }
    period_supplied <- period
    if (is.null(period)) {
      period <- 365 / attr(data, "survey_length")
    }
    if (!is.null(period) && (!is.numeric(period) || period <= 0)) {
      cli::cli_abort("{.arg period} must be a positive number.")
    }
  } else if (!is.null(period)) {
    cli::cli_abort("{.arg period} requires {.arg periodic_gp = TRUE}.")
  }

  # species-level GP checks
  if (data$S == 1 && species_length_scales) {
    cli::cli_abort(
      "{.arg species_length_scales} only applicable if multiple species are \\
      included."
    )
  }

  # orthogonal projection check
  # if (project_kappa) {
  #   if (temporal == "none") {
  #     cli::cli_abort(
  #       "Orthogonal projection only applicable if survey random effects are \\
  #       included. Set {.arg project_kappa = FALSE} or turn on survey random \\
  #       effects."
  #     )
  #   }
  #   if (!sum(c(data$P[3], data$P_cat[3], data$P_ord[3]))) {
  #     cli::cli_abort(
  #       "Orthogonal projection only applicable if survey predictors are \\
  #       included. Set {.arg project_kappa = FALSE} or include predictors."
  #     )
  #   }
  # }

  # build Stan data list
  stan_data <- c(
    unclass(data),
    list(
      period = if (!is.null(period)) period else 0,
      grainsize = if (threads > 1L) as.integer(grainsize) else 0L,
      D = as.integer(loo_draws),
      spatial = switch(spatial, "none" = 0L, "mvn" = 1L, "gp" = 2L),
      temporal = switch(temporal, "none" = 0L, "mvn" = 1L, "gp" = 2L),
      OD = switch(overdispersion, "none" = 0L, "olre" = 1L, "nb" = 2L),
      dirichlet = as.integer(variance_decomposition == "dirichlet"),
      SS = as.integer(species_length_scales),
      project_kappa = as.integer(project_kappa),
      latent = as.integer(latent),
      PPC_y = as.integer(ppc %in% c("both", "y")),
      PPC_Q = as.integer(ppc %in% c("both", "Q"))
    ),
    unclass(prior)
  )

  # strip non-Stan attributes that were attached by make_data()
  stan_data[c(
    "survey_length",
    "utm_crs",
    "thin_minutes",
    "reference_date",
    "site_levels",
    "species_levels",
    "surveys",
    "scaling"
  )] <- NULL

  # compile model
  if (is.null(stan_file)) {
    stan_file <- system.file(
      if (data$S > 1L) "stan/occARU.stan" else "stan/occARU-ss.stan",
      package = "occARU"
    )
  }
  mod <- cmdstanr::cmdstan_model(
    stan_file,
    cpp_options = list(stan_threads = TRUE)
  )

  # log model components
  cli::cli_inform(c(
    "Components: ",
    " " = "Spatial effects: {switch(spatial, 'gp' = 'GP', 'mvn' = 'MVN',
                                             'none' = 'None')}, \\
           Temporal effects: {switch(temporal, 'gp' = 'GP', 'mvn' = 'MVN',
                                               'none' = 'None')}\\
           {if (temporal == 'gp') paste0(', Periodic GP: ', {if (periodic_gp) \\
            'yes' else 'no'})}",
    if (temporal == "gp" && periodic_gp && is.null(period_supplied)) {
      c(
        " " = "\u00a0\u00a0Period not specified. Defaulting to period = \\
              {round(period, 1)} (annual cycle with \\
              {attr(data, 'survey_length')}-day surveys)"
      )
    },
    if ((spatial == "gp" || temporal == "gp") && data$S > 1) {
      c(
        " " = "\u00a0\u00a0Species-specific length scales: \\
              {if (species_length_scales) 'yes' else 'no'}"
      )
    },
    " " = "Overdispersion: {switch(overdispersion,
                                   'none' = 'None (Poisson)',
                                   'olre' = 'OLRE (Poisson)',
                                   'nb' = 'Negative Binomial')}",
    " " = "Variance decomposition: \\
          {switch(variance_decomposition,
                  'dirichlet' = 'Dirichlet',
                  'logistic-normal' = 'Logistic-normal')}",
    " " = "Latent occupancy: {if (latent) 'yes' else 'no'}",
    " " = "Posterior predictions: \\
           {switch(ppc, 'Q' = 'Q only', 'y' = 'y only',
                        'both' = 'y and Q', 'none' = 'none')}",
    if (spatial != "none" || overdispersion == "olre") {
      c(" " = "LOO Monte Carlo draws: {loo_draws}")
    }
  ))

  # initial values and sample
  dots <- list(...)
  chains <- dots[["chains"]] %||% 4L

  if (identical(init, "pathfinder")) {
    cli::cli_inform(c("i" = "Running pathfinder for initialisation..."))
    pathfinder_defaults <- list(
      data = stan_data,
      init = 0.1,
      refresh = 0,
      num_paths = chains,
      num_threads = chains,
      sig_figs = 14
    )
    init <- do.call(
      mod$pathfinder,
      utils::modifyList(pathfinder_defaults, pathfinder_args)
    )
  }

  fit <- mod$sample(
    data = stan_data,
    init = init,
    parallel_chains = chains,
    threads_per_chain = threads,
    ...
  )

  attr(fit, "stan_data") <- stan_data
  attr(fit, "occARU_data") <- data
  cli::cli_inform(c(
    "i" = 'Done. Stan data (including priors) and occARU data stored as \\
          attributes. Access with {.code attr(fit, "stan_data")} and \\
          {.code attr(fit, "occARU_data")}.'
  ))
  fit
}
