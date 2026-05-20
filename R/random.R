#' Specify Gaussian process in occARU models
#'
#' Specify Gaussian process (GP) random effects in occARU models, with a choice
#' of kernel and whether to estimate species-level length-scales.
#'
#' @param kernel `character`. The type of GP kernel to use. Must be one of
#'   `"exp_quad"` (exponentiatated quadratic), `"matern32"` (Matern 3/2), or
#'   `"matern52"` (Matern 5/2).
#' @param species_length_scales `logical`. If `TRUE`, species-specific length
#'   scales are estimated for that GP, each drawn independently from the shared
#'   length scale priors. Note that enabling this requires additional Cholesky
#'   decompositions per species and GP, which can substantially increase
#'   sampling time. If `FALSE` (default), only one length scale is estimated
#'   with one Cholesky decomposition performed per GP. Only used when multiple
#'   species are included.
#' @param periodic `logical`. If `TRUE`, a periodic kernel is added to the
#'   kernel for survey effects. Default: `FALSE`.
#' @param period Positive numeric. Period length in survey units (i.e. number
#'   of survey periods per cycle). Only used when `periodic = TRUE`. Defaults to
#'   `365 / survey_length`, corresponding to an annual cycle. For example, with
#'   `survey_length = 7` the default is `period = 52.1`. Override if your data
#'   span a different temporal cycle.
#'
#' @return A `gp` object (a named list with class `"gp"`) for use in the
#'   `random` argument of [occARU()]. Contains the following elements:
#'   \describe{
#'     \item{`random`}{Integer. Always `1L`, indicating a random effect.}
#'     \item{`kernel`}{Integer. Kernel index (1 = exponentiated quadratic,
#'       2 = Matern 3/2, 3 = Matern 5/2).}
#'     \item{`species_length_scales`}{Integer. Whether to estimate
#'       species-specific length scales (0 or 1).}
#'     \item{`periodic`}{Integer. Whether a periodic kernel component is
#'       included (0 or 1).}
#'     \item{`period`}{Numeric. Period length in survey units, or `0` if not
#'       periodic.}
#'   }
#'
#' @seealso [occARU()], [set_priors()]
#' @export
#' @keywords internal
gp <- function(
  kernel = c("exp_quad_cov", "matern32", "matern52"),
  species_length_scales = FALSE,
  periodic = FALSE,
  period = NULL
) {
  kernel_choices <- c("exp_quad_cov", "matern32", "matern52")
  kernel <- match.arg(kernel)
  if (!is.logical(species_length_scales)) {
    cli::cli_abort("{.arg species_length_scales} must be TRUE or FALSE.")
  }
  if (!is.logical(periodic)) {
    cli::cli_abort("{.arg periodic} must be TRUE or FALSE.")
  }
  if (periodic) {
    if (!is.null(period) && (!is.numeric(period) || period <= 0)) {
      cli::cli_abort("{.arg period} must be a positive number.")
    }
  } else if (!is.null(period)) {
    cli::cli_abort("{.arg period} requires {.arg periodic = TRUE}.")
  }
  structure(
    list(
      random = 1L,
      kernel = which(kernel_choices == kernel),
      species_length_scales = as.integer(species_length_scales),
      periodic = as.integer(periodic),
      period = period %||% 0
    ),
    class = "gp"
  )
}

#' Organise random effects for use in [occARU()]
#'
#' Specify Gaussian process (GP) random effects in occARU models, with a choice
#' of kernel and whether to estimate species-level length-scales.
#'
#' @param random A named list of random effects.
#' @param K Number of seasons
#' @param dyn  Indicator of dynamic occupancy
#' @keywords internal
make_random <- function(random, K, dyn) {
  # check invalid list items
  invalid <- names(purrr::discard(
    random,
    ~ inherits(., "gp") || identical(., "mvn") || identical(., "none")
  ))
  if (length(invalid)) {
    cli::cli_abort(
      "Invalid {.arg random} element{?s} {.val {invalid}}: must be {.fn gp}, \\
      {.val mvn}, or {.val none}."
    )
  }

  # default list
  random_defaults <- list(
    site = gp(),
    survey = gp(),
    season = if (dyn) gp() else "none",
    site_occ = if (dyn) gp() else "none",
    season_occ = if (dyn) gp() else "none"
  )

  # validate multiseason
  if (K == 1) {
    wrong <- purrr::keep(
      names(random),
      ~ . %in% c("season", "site_occ", "season_occ")
    )
    if (length(wrong)) {
      cli::cli_abort(
        "{.arg random} element{?s} {.val {wrong}} not applicable in \\
        single season models."
      )
    }
  }

  # check names and modify input
  unknown <- setdiff(names(random), names(random_defaults))
  if (length(unknown)) {
    cli::cli_abort(
      "Unknown element{?s} in {.arg random}: {.val {unknown}}."
    )
  }
  random <- utils::modifyList(random_defaults, random)

  # check invalid periodic kernel
  periodic_misuse <- setdiff(
    names(purrr::keep(random, ~ inherits(., "gp") && .$periodic == 1)),
    "survey"
  )
  if (length(periodic_misuse)) {
    cli::cli_abort(
      "{.arg periodic} can only be used for survey effects, not \\
      {.val {periodic_misuse}}."
    )
  }

  # overwrite "mvn" and "none" option and return
  purrr::map_if(
    random,
    ~ !inherits(., "gp"),
    ~ list(
      random = ifelse(identical(., "mvn"), 1L, 0L),
      kernel = 0L,
      species_length_scales = 0L,
      period = 0L
    )
  )
}
