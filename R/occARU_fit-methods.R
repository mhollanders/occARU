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
  ndraws = NULL,
  ...
) {
  # checks and observed data
  fit <- object
  stan_data <- attr(fit, "stan_data")
  if (!(stan_data$PPC_Q || stan_data$PPC_y)) {}
  level <- match.arg(level, c("Q", "y"))
  y <- stan_data$y

  # aggregated counts
  if (level == "Q") {
    obs <- c(apply(y, c(1, 3), sum))
    if (stan_data$PPC_Q) {
      rep <- fit$draws("Qrep", format = "draws_matrix")
    } else if (stan_data$PPC_y) {
      rep <- tidybayes::spread_rvars(fit, yrep[i, j, s]) |>
        dplyr::summarise(Qrep = posterior::rvar_sum(yrep), .by = c(i, s)) |>
        dplyr::mutate(name = paste0("Qrep[", i, ",", s, "]")) |>
        dplyr::select(name, Qrep) |>
        tibble::deframe() |>
        posterior::as_draws_matrix()
    } else {
      cli::cli_abort(
        'Model was fit without posterior predictions. Refit the model without
        {.arg ppc = "none"}.'
      )
    }

    # observed counts
  } else {
    if (!stan_data$PPC_y) {
      cli::cli_abort(
        'Model was fit without posterior predictions for `y`. Refit the model
        using with {.arg ppc = "y" or "both"}.'
      )
    }
    surveyed <- which(c(stan_data$Delta) > 0)
    obs <- c(y)[surveyed]
    rep <- fit$draws("yrep", format = "draws_matrix")[, surveyed]
  }

  # subsample
  if (!is.null(ndraws)) {
    total_draws <- nrow(rep)
    if (!rlang::is_integerish(ndraws) || ndraws > total_draws || ndraws <= 0) {
      cli::cli_abort(
        "{.arg ndraws} must be a positive integer less than the number of
        draws."
      )
    }
    rep <- rep[sample(1:total_draws, ndraws), ]
  }

  # plot
  occARU_data <- attr(fit, "occARU_data")
  bayesplot::ppc_rootogram_grouped(
    obs,
    rep,
    group = rep(
      attr(occARU_data, "species"),
      each = if (level == "y") length(surveyed) else stan_data$I
    ),
    ...
  )
}
