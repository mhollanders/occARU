# summary.occARU_fit <- function(fit) {
#   occARU_data <- attr(fit, "occARU_data")
#   stan_data <- attr(fit, "stan_data")
#   species_lvl <- attr(occARU_data, "species")
#
#   intercepts <- tidybayes::spread_rvars(fit, alpha[d, r, s]) |>
#     dplyr::summarise(alpha = posterior::rvar_mean(alpha), .by = c(d, s)) |>
#     dplyr::mutate(d = factor(d, labels = c("Occupancy Probability", "Detection Rate")),
#                   s = factor(s, labels = species_lvl)) |>
#     dplyr::rename(Submodel = d, Species = s)
#
#   fit$summary()
#
#   posterior::summarise_draws(intercepts, alpha)
#
#   gt::gt(intercepts)
#
# }

#' Posterior predictive checks using rootograms for `occARU_fit` objects
#'
#' Perform visual posterior predictive checks using rootograms using the
#' \pkg{bayesplot} package.
#'
#' @param object A fitted model object from [occARU()].
#' @param level Character. One of `"Q"` or `"y"`. If `"Q"` (default), uses the
#'   aggregated counts per site and species as data. If `"y"`, uses the raw
#'   (survey-level) counts as input.
#' @param group Character. Whether to group by `"species"` (default),
#'   `"region"`, or `"season"`.
#' @param ndraws Positive integer. Number of draws to use. If `NULL` (default),
#'   uses all draws.
#' @param ... Additional arguments to be passed to
#'   [bayesplot::ppc_rootogram_grouped()].
#' @return A ggplot object with posterior predictive rootograms faceted by
#'   species which can be modified using the \pkg{ggplot2} package.
#' @importFrom bayesplot pp_check
#' @aliases pp_check
#' @export pp_check
#' @export
pp_check.occARU_fit <- function(
  object,
  level = c("Q", "y"),
  group = c("species", "region", "season"),
  ndraws = NULL,
  ...
) {
  # checks and observed data
  fit <- object
  stan_data <- attr(fit, "stan_data")
  occARU_data <- attr(fit, "occARU_data")
  if (!(stan_data$PPC_Q || stan_data$PPC_y)) {
    cli::cli_abort(
      "Model was fit without posterior predictions."
    )
  }
  level <- match.arg(level, c("Q", "y"))
  group <- match.arg(group, c("species", "region", "season"))
  species_lvl <- attr(occARU_data, "species")
  region_lvl <- attr(occARU_data, "regions")
  season_lvl <- attr(occARU_data, "seasons")
  y <- stan_data$y

  # aggregated counts
  if (level == "Q") {
    draws <- if (stan_data$PPC_Q) {
      tidybayes::spread_rvars(fit, Qrep[k, i, s])
    } else {
      tidybayes::spread_rvars(fit, yrep[k, i, j, s]) |>
        dplyr::summarise(Qrep = posterior::rvar_sum(yrep), .by = c(k, i, s))
    }
    draws <- draws |>
      dplyr::mutate(y = apply(stan_data$y, c(1, 2, 4), sum)[cbind(k, i, s)]) |>
      dplyr::filter(apply(stan_data$Delta, c(1, 3), sum)[cbind(k, i)] > 0)
  } else if (!stan_data$PPC_y) {
    cli::cli_abort(
      'Model was fit without posterior predictions for `y`. Refit the model
      using with {.arg ppc = "y" or "both"}.'
    )
  } else {
    draws <- tidybayes::spread_rvars(fit, yrep[k, i, j, s]) |>
      dplyr::mutate(y = y[cbind(k, i, j, s)]) |>
      dplyr::filter(stan_data$Delta[cbind(k, j, i)] > 0)
  }

  # label and extract
  draws <- dplyr::mutate(
    draws,
    species = factor(species_lvl[s], species_lvl),
    region = factor(region_lvl[stan_data$region[i]], region_lvl),
    season = factor(season_lvl[k], season_lvl)
  ) |>
    dplyr::select(-dplyr::any_of(c("k", "i", "j", "s")))
  y <- draws$y
  yrep <- as_draws_matrix(if (level == "Q") draws$Qrep else draws$yrep)
  group <- draws[[group]]

  # subsample
  if (!is.null(ndraws)) {
    total_draws <- nrow(yrep)
    if (!rlang::is_integerish(ndraws) || ndraws > total_draws || ndraws <= 0) {
      cli::cli_abort(
        "{.arg ndraws} must be a positive integer less than the number of
        draws."
      )
    }
    yrep <- yrep[sample(1:total_draws, ndraws), ]
  }

  # plot
  bayesplot::ppc_rootogram_grouped(y, yrep, group = group, ...)
}
