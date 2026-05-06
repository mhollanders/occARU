#' Fit the occARU model
#'
#' Fits the occARU model to data prepared by [make_data()] in Stan. Requires
#' CmdStan >= 2.36.0, which can be installed with [setup_occARU()].
#' Automatically determines number of species in the supplied data and whether a
#' single or multiseason model is required.
#'
#' @param data A `occARU_data` object produced by [make_data()].
#' @param prior An `occARU_priors` object produced by [set_priors()]. If
#'   omitted, default priors are used.
#' @param stan_file `character`. Path to a custom Stan file. If `NULL`
#'   (default), uses the built-in occARU models. Intended for advanced users who
#'   have modified the Stan programs; note that custom models will likely
#'   require corresponding changes to the output of [make_data()].
#' @param random A named `list` specifying random effect structures for site and
#'   survey effects on detection. For multiseason models, optionally accepts
#'   season effects on detection (`season`) and site (`site_occ`) and season
#'   (`season_occ`) effects on occupancy. Must be one of:
#'   * `gp()` (default for site, survey, and season effects on detection), which
#'     fits a (multispecies) Gaussian process. See [gp()] for details.
#'   * `"mvn"`, which fits an unstructured (multispecies) normal.
#'   * `"none"` (default for site and season effects on occupancy), which omits
#'     random effects entirely.
#' @param project A named `list` of logicals specifying whether to project
#'   random effects when predictors are included. Entries must be one of `site`
#'   or `survey`. If `TRUE` (default for both), orthogonally projects random
#'   effects. For survey effects, the design matrices are averaged across sites.
#'   For site effects in multiseason models, the design matrices are averaged
#'   across seasons.
#' @param overdispersion `character`. Overdispersion model for the observation
#'   process. One of `"none"` (Poisson, default), `"nb"` (negative binomial),
#'   or `"olre"` (multispecies observation-level random effects).
#' @param variance_decomposition `character`. Prior for variance partitions.
#'   One of `"dirichlet"` (default) or `"logistic-normal"`.
#' @param ppc `character`. Posterior predictive checks to compute. One of `"Q"`
#'   (default), `"y"`, `"both"`, or `"none"`. `"y"` returns the full
#'   `[I, J, S]` prediction array (`yrep`); `"Q"` returns only aggregated
#'   counts `[I, S]` (`Qrep`). For large datasets, `"Q"` or `"none"` can
#'   substantially reduce memory usage and sampling time.
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
#'   `log_lik`. Only used with site random effects or `overdispersion =
#'   "olre"`. Note: even high values (10K) seem inadequate for OLREs.
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
#' @param grainsize Positive integer. Chunk size for within-chain
#'   parallelisation via `reduce_sum()` when `threads > 1`. For sites nested in
#'   regions, chunks are number of regions; otherwise it is number of sites.
#'   Default: `1`, which lets Stan automatically determine the optimal chunk
#'   size. Increase if you have many sites or regions and want to reduce
#'   parallelisation overhead. See the
#'   \href{https://mc-stan.org/docs/stan-users-guide/parallelization.html#reduce-sum-grainsize}{Stan
#'   User's Guide} for details on tuning grainsize.
#' @param ... Additional arguments passed to
#'   [cmdstanr::CmdStanModel]`$sample()`. Uses parallel chains by default.
#'   All other sampling arguments use Stan defaults.
#'
#' @return An `occARU_fit` object, which extends the
#'   \href{https://mc-stan.org/cmdstanr/reference/CmdStanFit.html}{CmdStanFit}
#'   class with occARU-specific attributes:
#'   \describe{
#'     \item{`stan_data`}{The full Stan data list passed to the model,
#'       including prior hyperparameters.}
#'     \item{`occARU_data`}{The original `occARU_data` object from
#'       [make_data()].}
#'   }
#' @seealso [make_data()], [set_priors()], [gp()], [setup_occARU()],
#'   [pp_check.occARU_fit()], [cmdstanr::CmdStanMCMC] for methods on the fitted
#'   model object. The statistical model is described in `vignette("model")`.
#' @export
occARU <- function(
  data,
  prior = set_priors(verbose = FALSE),
  stan_file = NULL,
  random = list(site = gp(), survey = gp()),
  project = list(site = TRUE, survey = TRUE),
  overdispersion = c("none", "nb", "olre"),
  variance_decomposition = c("dirichlet", "logistic-normal"),
  ppc = c("Q", "y", "both", "none"),
  latent = TRUE,
  loo_draws = 100L,
  init = "pathfinder",
  pathfinder_args = list(),
  threads = 1L,
  grainsize = 1L,
  ...
) {
  # input checks
  if (!inherits(data, "occARU_data")) {
    cli::cli_warn(
      "{.arg data} is not a {.cls occARU_data} object from {.fun make_data}."
    )
  }
  if (!inherits(prior, "occARU_priors")) {
    cli::cli_warn(
      "{.arg prior} is not a {.cls occARU_priors} object from
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

  # prepare random effects
  random <- make_random(random, data$K)

  # spatial checks
  if (all(data$XY == 0)) {
    spatial <- intersect(c("site", "site_occ"), names(random))
    spatial_misuse <- names(purrr::map_lgl(
      random[spatial],
      ~ inherits(., "gp")
    ))
    if (length(spatial_misuse)) {
      cli::cli_abort(
        'No site coordinates found in {.arg data}. Set {.arg random = "mvn"}
        for element{?s} {.val {spatial_misuse}} or supply {.arg latitude} and
        {.arg longitude} in {.fun make_data} to enable the spatial GP{?s}.'
      )
    }
  }

  # add period
  default_period <- round(365 / attr(data, "survey_length"), 1)
  if (with(random$survey, kernel && periodic && period == 0)) {
    random$survey$period <- default_period
  }

  # site effects
  site_effects <- any(purrr::map_lgl(
    random[intersect(c("site", "site_occ"), names(random))],
    ~ .$random == 1
  ))
  site_gp <- any(purrr::map_lgl(
    random[intersect(c("site", "site_occ"), names(random))],
    ~ .$kernel > 0
  ))

  # check length scales
  if (any(purrr::map_lgl(random, ~ .$species_length_scales)) && data$S == 1) {
    cli::cli_abort(
      "{.arg species_length_scales = TRUE} is not applicable in single \\
      species models."
    )
  }

  # check projection and turn to 0 without predictors
  if (any(!is.logical(unlist(project)))) {
    cli::cli_abort("All entries in {.arg project} must be {.code logical}.")
  } else {
    project_defaults <- list(site = TRUE, survey = TRUE)
    project <- utils::modifyList(project_defaults, project)
    P_sum <- data$P + data$P_cat + data$P_ord
    if (!(random$site$random && P_sum[2]) && project$site) {
      project$site <- FALSE
    }
    if (!(random$survey$random && P_sum[3]) && project$survey) {
      project$survey <- FALSE
    }
  }

  # match arguments
  overdispersion <- match.arg(overdispersion)
  variance_decomposition <- match.arg(variance_decomposition)
  ppc <- match.arg(ppc)

  # scalar argument checks
  if (!rlang::is_integerish(loo_draws) || loo_draws < 0) {
    cli::cli_abort("{.arg loo_draws} must be a non-negative integer.")
  } else if (!(site_effects || overdispersion == "olre") && loo_draws != 100L) {
    cli::cli_warn(
      "{.arg loo_draws} only used when site effects or OLREs are included."
    )
  }

  if (!rlang::is_integerish(threads) || threads < 1) {
    cli::cli_abort("{.arg threads} must be a positive integer.")
  } else if (!rlang::is_integerish(grainsize) || grainsize < 0) {
    cli::cli_abort("{.arg grainsize} must be a positive integer.")
  } else if (data$R > 1 && random$site$kernel && grainsize > data$R) {
    cli::cli_abort(
      "{.arg grainsize} should be less than the number of regions."
    )
  } else if (grainsize > data$I) {
    cli::cli_abort(
      "{.arg grainsize} should be less than the number of sites."
    )
  }

  # build Stan data list
  stan_data <- c(
    unclass(data),
    list(
      random = purrr::map_int(random, ~ .$random),
      project = as.integer(unlist(project)),
      SS = purrr::map_int(random, ~ .$species_length_scales),
      kernel = purrr::map_int(random, ~ .$kernel),
      period = random$survey$period,
      OD = switch(overdispersion, "none" = 0L, "olre" = 1L, "nb" = 2L),
      dirichlet = as.integer(variance_decomposition == "dirichlet"),
      latent = as.integer(latent),
      PPC_y = as.integer(ppc %in% c("both", "y")),
      PPC_Q = as.integer(ppc %in% c("both", "Q")),
      grainsize = if (threads > 1L) as.integer(grainsize) else 0L,
      D = as.integer(loo_draws)
    ),
    unclass(prior)
  )

  # compile model
  mod <- cmdstanr::cmdstan_model(
    stan_file %||% system.file("stan/occARU.stan", package = "occARU"),
    cpp_options = list(stan_threads = TRUE)
  )

  # GP labels
  kernel_label <- \(k) c("exp. quad.", "Matern 3/2", "Matern 5/2")[k]
  fmt_effect <- \(x, role = NULL) {
    if (inherits(x, "gp")) {
      paste0(
        "GP (",
        paste(c(role, kernel_label(x$kernel)), collapse = ", "),
        ")"
      )
    } else if (x$random) {
      paste0("MVN", if (!is.null(role)) paste0(" (", role, ")"))
    } else {
      paste0("None", if (!is.null(role)) paste0(" (", role, ")"))
    }
  }

  # log model components
  cli::cli_inform(c(
    "Components: ",
    " " = "Site effects: {if (data$K == 1) {fmt_effect(random$site)}}\\
          {if (data$K > 1) paste0(', ', fmt_effect(random$site_occ, 'occupancy'))}",
    " " = "Survey effects: {fmt_effect(random$survey)}\\
          {if (inherits(random$survey, 'gp') && random$survey$periodic) ', periodic kernel: yes' else ''}",
    if (data$K > 1) {
      c(
        " " = "Seasonal effects: {fmt_effect(random$season, 'detection')}, \\
           {fmt_effect(random$season_occ, 'occupancy')}"
      )
    },
    if (with(random$survey, kernel && periodic && period == default_period)) {
      c(
        " " = "\u00a0\u00a0Period not specified. Defaulting to period = \
           {default_period} (annual cycle with \
           {attr(data, 'survey_length')}-day surveys)"
      )
    },

    if (site_effects && data$S > 1) {
      ss <- names(purrr::keep(
        random,
        ~ inherits(., "gp") && .$species_length_scales
      ))
      c(
        " " = "Species-specific length scales: {if (length(ss)) paste(ss, collapse = ', ') else 'no'}"
      )
    },
    " " = "Overdispersion: {switch(overdispersion,
                                 'none' = 'None (Poisson)',
                                 'olre' = 'OLRE (Poisson)',
                                 'nb' = 'Negative Binomial')}",
    " " = "Variance decomposition: {switch(variance_decomposition,
                                         'dirichlet' = 'Dirichlet',
                                         'logistic-normal' = 'Logistic-normal')}",
    " " = "Latent occupancy: {if (latent) 'yes' else 'no'}",
    " " = "Posterior predictions: {switch(ppc, 'Q' = 'Q only', 'y' = 'y only',
                                             'both' = 'y and Q', 'none' = 'none')}",
    if (site_effects || overdispersion == "olre") {
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
    init <- rlang::exec(
      mod$pathfinder,
      !!!utils::modifyList(pathfinder_defaults, pathfinder_args)
    )
  }

  fit <- mod$sample(
    data = stan_data,
    init = init,
    parallel_chains = chains,
    threads_per_chain = threads,
    ...
  )

  class(fit) <- c("occARU_fit", class(fit))
  attr(fit, "stan_data") <- stan_data
  attr(fit, "occARU_data") <- data
  cli::cli_inform(c(
    "i" = 'Done. Stan data (including priors) and occARU data stored as
          attributes. Access with {.code attr(fit, "stan_data")} and
          {.code attr(fit, "occARU_data")}.'
  ))
  fit
}
