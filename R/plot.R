#' Plot predictor coefficients
#'
#' Plots coefficients of continuous, categorical, or ordinal predictors on
#' occupancy or detection submodels.
#'
#' occARU accommodates three types of predictors: continuous, categorical,
#' and ordinal (ordered categorical). Categorical predictors are modeled as
#' zero-sum Gaussians, and coefficients are plotted for each level in each
#' predictor. Ordinal predictors are modeled with a simplex decomposition of the
#' category differences. Ordinal coefficients can be plotted as the overall
#' effect (coefficient when the ordinal predictor has the maximum value) or
#' plotted for each level.
#'
#' @param fit A fitted model object from [occARU()].
#' @param submodel `character`. Predictors of submodel to plot. One of
#'   `"detection"` (default) or `"occupancy"`.
#' @param component `character`. Whether to plot `"site"` (default) or `"survey"`
#'   predictors. If `"survey"`, `submodel` must be `"detection"`.
#' @param type `character`. Type of predictors to plot. One of `"continuous"`
#'   (default), `"categorical"`, or `"ordinal"`.
#' @param level `character`. For multi-species models, whether to plot
#'   species-specific (`"species"`, default) or mean coefficients (`"mean"`).
#' @param facet_by `character`. Whether to use [ggplot2::facet_wrap()] or
#'   [ggh4x::facet_grid2()] to facet by `"predictor"` (default) or `"species"`.
#'   Only used if `level` is `"species"`.
#' @param species `character`. Vector of species to plot. If `NULL` (default),
#'   all species are plotted. Must be one of `attr(occARU_data, "species")`.
#' @param restricted `logical`. If `TRUE` (default), plots coefficients with
#'   orthogonal projection of the detection random site or survey effects, e.g.,
#'   \eqn{\boldsymbol{\iota}(\boldsymbol{I} - \boldsymbol{P_{X_2}})}, where
#'   \eqn{\boldsymbol{I} - \boldsymbol{P_{X_2}}} is the orthogonal complement of
#'   the column space of the design matrix. If `FALSE`, recovers coefficients
#'   without orthogonal projection, e.g., \eqn{\boldsymbol{\beta} -
#'   \boldsymbol{X_2}^+ \boldsymbol{\iota}}, where \eqn{\boldsymbol{X_2}^+} is
#'   the pseudo-inverse of the design matrix. Only used for site predictors if
#'   `submodel = "detection"`, or if survey random effects were projected with
#'   `project = list(survey = TRUE)` in `occARU()`.
#' @param ordinal_categories `logical`. If `FALSE` (default), plots coefficients
#'   associated with maximum category (full effect). If `TRUE`, plots realised
#'   coefficient associated with each ordered category, where the first
#'   is used as the reference.
#' @param ... Additional arguments passed to [ggdist::stat_pointinterval()].
#'
#' @return A `ggplot` object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`plot_data`}{The tibble used to produce the plot.}
#'   }
#'
#' @seealso [occARU()], [plot_intercepts()], [plot_sites()], [plot_surveys()],
#'   [plot_correlations()], [plot_partitions()], [plot_realised()]
#' @export
plot_coefficients <- function(
  fit,
  submodel = c("detection", "occupancy"),
  component = c("site", "survey"),
  type = c("continuous", "categorical", "ordinal"),
  level = c("species", "mean"),
  facet_by = c("species", "predictor"),
  species = NULL,
  restricted = TRUE,
  ordinal_categories = FALSE,
  ...
) {
  if (!inherits(fit, "CmdStanFit")) {
    cli::cli_abort(
      "{.arg fit} must be a {.cls CmdStanFit} object from {.fun occARU}."
    )
  }
  stan_data <- attr(fit, "stan_data")
  occARU_data <- attr(fit, "occARU_data")
  species_lvl <- attr(occARU_data, "species")
  MS <- length(species_lvl) > 1
  submodel <- match.arg(submodel, c("detection", "occupancy"))
  component <- match.arg(component, c("site", "survey"))
  type <- match.arg(type, c("continuous", "categorical", "ordinal"))
  level <- match.arg(level, c("species", "mean"))
  facet_by <- match.arg(facet_by, c("species", "predictor"))

  # get survey indicator
  survey <- FALSE
  if (component == "survey") {
    if (submodel == "occupancy") {
      cli::cli_abort(
        "Occupancy does not have survey-level predictors."
      )
    } else {
      survey <- TRUE
    }
  }

  # total number of predictors
  idx <- if (submodel == "occupancy") {
    1L
  } else if (survey) {
    3L
  } else {
    2L
  }
  P <- occARU_data$P[idx]
  P_cat <- occARU_data$P_cat[idx]
  P_ord <- occARU_data$P_ord[idx]
  where <- if (!survey) paste("for", submodel) else ""
  if (!sum(P, P_cat, P_ord)) {
    cli::cli_abort(
      "No {component} predictors were found {where} in {.arg fit}."
    )
  }

  # get restricted indicator
  res <- TRUE
  if (!restricted) {
    if (submodel == "occupancy") {
      cli::cli_abort(
        '{.arg restricted = FALSE} is only applicable when
        {.arg submodel = "detection"}.'
      )
    } else if (survey && !stan_data$project[2]) {
      cli::cli_abort(
        "{.arg restricted = FALSE} is only applicable when random survey
        effects were orthogonally projected. Requires
        {.arg project = list(survey = TRUE)} in {.fun occARU}."
      )
    } else {
      res <- FALSE
    }
    if (type == "categorical") {
      cli::cli_abort(
        "Recovery of unconditional categorical coefficients is currently not
        supported."
      )
    }
  }

  # get labels
  survey <- submodel == "detection" && component == "survey"
  prefix <- if (submodel == "occupancy") {
    "psi_"
  } else if (survey) {
    NULL
  } else {
    "mu_"
  }
  coef <- ifelse(survey, "gamma", "beta")
  suffix <- if (!res) "2" else NULL

  # get species indices
  if (level == "species") {
    species_idx <- indices(species, species_lvl)
  } else {
    if (!is.null(species)) {
      cli::cli_warn(
        "{.arg facet_by} and {.arg species} are ignored when
        {.arg level = {level}}."
      )
    }
  }

  # continuous
  if (type == "continuous") {
    if (P) {
      X <- stan_data[[
        if (submodel == "occupancy") {
          "X1"
        } else if (survey) {
          "X3"
        } else {
          "X2"
        }
      ]]
      X_lbl <- if (survey) dimnames(X)[[3]] else colnames(X)
      param <- rlang::sym(paste0(prefix, coef, suffix))
      draws <- tidybayes::spread_rvars(fit, (!!param)[p, s]) |>
        dplyr::filter(s %in% species_idx) |>
        dplyr::mutate(
          s = factor(species_lvl[s], levels = species_lvl[species_idx]),
          p = factor(p, labels = X_lbl)
        ) |>
        dplyr::rename(species = s, predictor = p)

      p <- ggplot2::ggplot(draws) +
        ggplot2::aes(
          xdist = !!param,
          y = if (facet_by == "predictor") {
            forcats::fct_rev(species)
          } else {
            predictor
          }
        ) +
        ggplot2::facet_wrap(
          ~ if (facet_by == "predictor") {
            predictor
          } else {
            species
          },
          scales = if (facet_by == "predictor") "free_x" else "fixed"
        ) +
        my_vline() +
        ggdist::stat_pointinterval(...) +
        ggplot2::labs(
          x = "Coefficient",
          y = if (facet_by == "predictor") "Species" else "Predictor"
        )
    } else {
      cli::cli_abort(
        "No {type} {component} predictors were found {where} in {.arg fit}."
      )
    }

    # factors
  } else {
    lvls <- attr(occARU_data, "levels")[[
      if (submodel == "occupancy") {
        "X1"
      } else if (survey) {
        "X3"
      } else {
        "X2"
      }
    ]]
    if (type == "categorical") {
      if (P_cat) {
        lvls_df <- tibble::enframe(
          lvls$cat,
          name = "predictor",
          value = "levels"
        ) |>
          dplyr::mutate(p = dplyr::row_number()) |>
          tidyr::unnest_longer(levels, values_to = "category") |>
          dplyr::mutate(c = dplyr::row_number(), .by = predictor)

        param <- rlang::sym(paste0(prefix, coef, "_cat"))
        draws <- tidybayes::spread_rvars(fit, (!!param)[p, c, s]) |>
          dplyr::filter(s %in% species_idx) |>
          dplyr::right_join(lvls_df, by = dplyr::join_by(p, c)) |>
          dplyr::mutate(
            species = factor(
              species_lvl[s],
              levels = species_lvl[species_idx]
            ),
            predictor = factor(predictor, levels = names(lvls$cat))
          ) |>
          dplyr::mutate(
            category = factor(category, levels = lvls$ord[[predictor[1]]]),
            .by = predictor
          ) |>
          dplyr::select(-c(p, s, c)) |>
          dplyr::relocate(
            species,
            predictor,
            category,
            .before = !!param
          )

        p <- ggplot2::ggplot(draws) +
          ggplot2::aes(xdist = !!param, y = forcats::fct_rev(category)) +
          ggh4x::facet_grid2(
            if (facet_by == "predictor") {
              predictor ~ species
            } else {
              species ~ predictor
            },
            scales = "free_y",
            independent = "y"
          ) +
          my_vline() +
          ggdist::stat_pointinterval(...) +
          ggplot2::labs(x = "Coefficient", y = "Category")
      } else {
        cli::cli_abort(
          "No {type} {component} predictors were found {where} in {.arg fit}."
        )
      }
    } else {
      if (P_ord) {
        lvls_df <- tibble::enframe(
          lvls$ord,
          name = "predictor",
          value = "levels"
        ) |>
          dplyr::mutate(p = dplyr::row_number()) |>
          tidyr::unnest_longer(levels, values_to = "category") |>
          dplyr::slice(-1, .by = predictor) |>
          dplyr::mutate(o = dplyr::row_number(), .by = predictor)

        cs <- rlang::sym(paste0(prefix, coef, "_ord_cs"))
        param <- rlang::sym(paste0(prefix, coef, "_ord", suffix))
        draws <- tidybayes::spread_rvars(
          fit,
          (!!param)[p, s],
          (!!cs)[o, p]
        ) |>
          dplyr::filter(s %in% species_idx) |>
          dplyr::right_join(lvls_df, by = dplyr::join_by(o, p)) |>
          dplyr::mutate(
            species = factor(
              species_lvl[s],
              levels = species_lvl[species_idx]
            ),
            predictor = factor(predictor, levels = names(lvls$ord))
          ) |>
          dplyr::mutate(
            category = factor(category, levels = lvls$ord[[predictor[1]]]),
            .by = predictor
          ) |>
          dplyr::select(-c(p, s, o)) |>
          dplyr::relocate(
            species,
            predictor,
            category,
            .before = !!param
          )

        if (ordinal_categories) {
          p <- ggplot2::ggplot(draws) +
            ggplot2::aes(
              xdist = !!param * !!cs,
              y = forcats::fct_rev(category)
            ) +
            ggh4x::facet_grid2(
              if (facet_by == "predictor") {
                predictor ~ species
              } else {
                species ~ predictor
              },
              scales = "free_y",
              independent = "y"
            ) +
            my_vline() +
            ggdist::stat_pointinterval(...) +
            ggplot2::labs(x = "Coefficient", y = "Category")
        } else {
          p <- ggplot2::ggplot(draws) +
            ggplot2::aes(
              xdist = !!param,
              y = if (facet_by == "predictor") {
                forcats::fct_rev(species)
              } else {
                forcats::fct_rev(predictor)
              }
            ) +
            ggplot2::facet_wrap(
              ~ if (facet_by == "predictor") {
                predictor
              } else {
                species
              },
              scales = if (facet_by == "predictor") "free_x" else "fixed"
            ) +
            my_vline() +
            ggdist::stat_pointinterval(...) +
            ggplot2::labs(
              x = "Coefficient",
              y = if (facet_by == "predictor") "Species" else "Predictor"
            )
        }
      } else {
        cli::cli_abort(
          "No {type} {component} predictors were found {where} in {.arg fit}."
        )
      }
    }
  }
  attr(p, "plot_data") <- draws
  p
}

#' Plot interspecific correlations
#'
#' Plot pairwise species correlations for different model components. For
#' occupancy, interspecific correlations are only estimated for responses
#' to site-level predictors. For detection, interspecific correlations are
#' potentially estimated for responses to site-level predictors, survey-level
#' predictors, random site effects, random survey effects, and potentially
#' Poisson OLREs.
#'
#' @param fit A fitted model object from [occARU()].
#' @param submodel `character`. Correlations of submodel to plot. One of
#'   `"detection"` (default) or `"occupancy"`.
#' @param species `character`. Vector of species to plot. If `NULL` (default),
#'   all species are plotted. Must be one of `attr(occARU_data, "species")`.
#' @param ... Additional arguments passed to [ggdist::stat_pointinterval()].
#'
#' @return A `ggplot` object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`plot_data`}{The tibble used to produce the plot.}
#'   }
#'
#' @seealso [occARU()], [plot_intercepts()], [plot_coefficients()],
#'   [plot_sites()], [plot_surveys()], [plot_partitions()],
#'   [plot_realised()]
#' @export
plot_correlations <- function(
  fit,
  submodel = c("detection", "occupancy"),
  species = NULL,
  ...
) {
  if (!inherits(fit, "CmdStanFit")) {
    cli::cli_abort(
      "{.arg fit} must be a {.cls CmdStanFit} object from {.fun occARU}."
    )
  }
  stan_data <- attr(fit, "stan_data")
  if (stan_data$S == 1) {
    cli::cli_abort(
      "No interspecific correlations exist in single species models."
    )
  }
  occARU_data <- attr(fit, "occARU_data")
  species_lvl <- attr(occARU_data, "species")
  species_idx <- indices(species, species_lvl)
  submodel <- match.arg(submodel, c("detection", "occupancy"))

  # extract correlations
  draws <- tidyr::expand_grid(s = species_idx, ss = species_idx)
  if (submodel == "occupancy") {
    P <- occARU_data$P[1] + occARU_data$P_cat[1] + occARU_data$P_ord[1]
    if (P) {
      draws <- dplyr::left_join(
        draws,
        tidybayes::spread_rvars(fit, psi_beta_O[s, ss])
      )
    } else {
      cli::cli_abort("No predictors were found for occupancy in {.arg fit}.")
    }
  } else {
    P <- occARU_data$P[2:3] + occARU_data$P_cat[2:3] + occARU_data$P_ord[2:3]
    random <- stan_data$random
    OLRE <- stan_data$OD == 1
    if (sum(P) + sum(random) + OLRE) {
      if (P[1]) {
        draws <- dplyr::left_join(
          draws,
          tidybayes::spread_rvars(fit, mu_beta_O[s, ss]),
          by = dplyr::join_by(s, ss)
        )
      }
      if (P[2]) {
        draws <- dplyr::left_join(
          draws,
          tidybayes::spread_rvars(fit, gamma_O[s, ss]),
          by = dplyr::join_by(s, ss)
        )
      }
      if (random[1]) {
        draws <- dplyr::left_join(
          draws,
          tidybayes::spread_rvars(fit, iota_O[s, ss]),
          by = dplyr::join_by(s, ss)
        )
      }
      if (random[1]) {
        draws <- dplyr::left_join(
          draws,
          tidybayes::spread_rvars(fit, kappa_O[s, ss]),
          by = dplyr::join_by(s, ss)
        )
      }
      if (OLRE) {
        draws <- dplyr::left_join(
          draws,
          tidybayes::spread_rvars(fit, epsilon_O[s, ss]),
          by = dplyr::join_by(s, ss)
        )
      }
    } else {
      cli::cli_abort(
        "No predictors or random effects were found for detection in {.arg fit}."
      )
    }
  }

  # pivot and label
  draws <- tidyr::pivot_longer(
    draws,
    -c(s, ss),
    names_to = ".variable",
    values_to = ".value"
  ) |>
    dplyr::mutate(
      dplyr::across(c(s, ss), ~ factor(., labels = species_lvl[species_idx])),
      .variable = factor(
        .variable,
        levels = paste0(
          c("psi_beta", "mu_beta", "gamma", "iota", "kappa", "epsilon"),
          "_O"
        ),
        labels = c(
          "Site Predictors",
          "Site Predictors",
          "Survey Predictors",
          "Site Effects",
          "Survey Effects",
          "OLRE"
        )
      )
    ) |>
    dplyr::filter(s != ss) |>
    dplyr::rename(species = s, species2 = ss)

  # plot
  p <- ggplot2::ggplot(draws) +
    ggplot2::aes(xdist = .value, y = forcats::fct_rev(species)) +
    ggplot2::facet_grid(species2 ~ .variable) +
    my_vline() +
    ggdist::stat_pointinterval(...) +
    ggplot2::scale_x_continuous(
      breaks = seq(-0.5, 0.5, 0.5),
      limits = c(-1, 1),
      expand = c(0, 0)
    ) +
    ggplot2::labs(x = "Correlation", y = "Species")
  attr(p, "plot_data") <- draws
  p
}

#' Plot intercepts
#'
#' Plot species-specific intercepts for occupancy and detection, by default
#' back-transformed to the original scale.
#'
#' @param fit A fitted model object from [occARU()].
#' @param back_transform `logical`. If `TRUE` (default), intercepts are
#'   back-transformed to the natural scale via `inv_logit()` for occupancy and
#'   `exp()` for detection rates. If `FALSE`, values are left on the scale of
#'   the link functions (logit for occupancy and log for detection).
#' @param by_region `logical`. Whether to plot intercepts by region, if multiple
#'   regions were included.
#' @param species `character`. Vector of species to plot. If `NULL` (default),
#'   all species are plotted. Must be one of `attr(occARU_data, "species")`.
#' @param regions `character`. Vector of regions to plot. If `NULL` (default),
#'   all regions are plotted. Must be one of `attr(occARU_data, "regions")`.
#'
#' @param ... Additional arguments passed to [ggdist::stat_pointinterval()].
#'
#' @return A `ggplot` object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`plot_data`}{The tibble used to produce the plot.}
#'   }
#'
#' @seealso [occARU()], [plot_coefficients()], [plot_sites()], [plot_surveys()],
#'   [plot_correlations()], [plot_partitions()], [plot_realised()]
#' @export
plot_intercepts <- function(
  fit,
  back_transform = TRUE,
  by_region = FALSE,
  species = NULL,
  regions = NULL,
  ...
) {
  if (!inherits(fit, "CmdStanFit")) {
    cli::cli_abort(
      "{.arg fit} must be a {.cls CmdStanFit} object from {.fun occARU}."
    )
  }
  stan_data <- attr(fit, "stan_data")
  occARU_data <- attr(fit, "occARU_data")
  species_lvl <- attr(occARU_data, "species")
  MS <- length(species_lvl) > 1
  region_lvl <- attr(occARU_data, "regions")
  MR <- length(region_lvl) > 1
  if (!MR && by_region) {
    cli::cli_abort(
      "{.arg by_region = TRUE} only applicable with multiple regions."
    )
  }

  # get indices
  species_idx <- indices(species, species_lvl)
  region_idx <- indices(regions, region_lvl)

  # extract intercepts, transform, and label
  draws <- tidybayes::spread_rvars(fit, alpha[d, r, s]) |>
    dplyr::filter(r %in% region_idx, s %in% species_idx) |>
    dplyr::mutate(
      intercept = if (back_transform) {
        dplyr::if_else(d == 1, inv_logit(alpha), exp(alpha))
      } else {
        alpha
      }
    ) |>
    dplyr::mutate(
      d = factor(
        d,
        labels = if (back_transform) {
          c("Occupancy Probability", "Detection Rate")
        } else {
          c("Occupancy Log Odds", "Log Detection Rate")
        }
      ),
      r = factor(r, labels = region_lvl[region_idx]),
      s = factor(s, labels = species_lvl[species_idx])
    ) |>
    dplyr::rename(submodel = d, region = r, species = s)
  if (!by_region) {
    draws <- dplyr::summarise(
      draws |> dplyr::select(-region),
      intercept = posterior::rvar_mean(intercept),
      .by = c(submodel, species)
    )
  }

  # plot
  p <- ggplot2::ggplot(draws) +
    ggplot2::aes(xdist = intercept) +
    ggplot2::labs(x = "Estimate")
  if (MS) {
    p <- p +
      ggdist::stat_pointinterval(
        ggplot2::aes(y = forcats::fct_rev(species)),
        ...
      ) +
      ggplot2::labs(y = "Species")
    p <- if (by_region) {
      p + ggplot2::facet_grid(region ~ submodel, scales = "free_x")
    } else {
      p + ggplot2::facet_wrap(~submodel, nrow = 1, scales = "free_x")
    }
  } else {
    p <- if (by_region) {
      p +
        ggplot2::facet_wrap(~submodel, nrow = 1, scales = "free_x") +
        ggdist::stat_pointinterval(
          ggplot2::aes(y = forcats::fct_rev(region))
        ) +
        ggplot2::labs(y = "Region")
    } else {
      p +
        ggdist::stat_pointinterval(
          ggplot2::aes(y = forcats::fct_rev(submodel)),
          ...
        ) +
        ggplot2::labs(y = "Submodel")
    }
  }
  attr(p, "plot_data") <- draws
  p
}


#' Plot variance partitions
#'
#' Plot variance partitions of the different model components for occupancy
#' and detection submodels.
#'
#' The occARU model uses global-local shrinkage priors for the occupancy and
#' detection submodels, where half-Student-t priors are used for the variances
#' of both linear predictors which are simplex partitioned via either Dirichlet
#' or logistic-normal decomposition. Variance decomposition only occurs when
#' there is more than one model component. Partitions exist for species-level
#' intercepts, and species-level slopes, sites effects, survey effects, and
#' Poisson OLREs (with one mean partition and one for species-level deviations).
#' The species-level scales for site and survey effects and OLREs are produced
#' by additional simplex decomposition of the species-level components.
#'
#' @param fit A fitted model object from [occARU()].
#' @param scales `logical`. If `FALSE` (default), plots variance simplex
#'   partitions \eqn{\boldsymbol{\phi}}. If `TRUE`, produces component scales by
#'   plotting \eqn{\sqrt{W \cdot \boldsymbol{\phi}}}, where \eqn{W} are
#'   variances of linear predictors. Useful for sparse simplexes, where few
#'   components account for most of the variance.
#' @param ... Additional arguments passed to [ggdist::stat_pointinterval()].
#'
#' @return A `ggplot` object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`plot_data`}{The tibble used to produce the plot.}
#'   }
#'
#' @seealso [occARU()] [plot_intercepts()], [plot_coefficients()],
#'   [plot_sites()], [plot_surveys()], [plot_correlations()],
#'   [plot_realised()]
#' @export
plot_partitions <- function(
  fit,
  scales = FALSE,
  ...
) {
  if (!inherits(fit, "CmdStanFit")) {
    cli::cli_abort(
      "{.arg fit} must be a {.cls CmdStanFit} object from {.fun occARU}."
    )
  }
  stan_data <- attr(fit, "stan_data")
  occARU_data <- attr(fit, "occARU_data")
  S <- stan_data$S
  MS <- S > 1
  R <- stan_data$R
  MR <- R > 1
  random <- stan_data$random
  OLRE <- stan_data$OD == 1
  P_sum <- stan_data$P + stan_data$P_cat + stan_data$P_ord

  # partitions and labels
  suffix <- if (MS) paste0(" (", c("mean", "species"), ")")
  psi_V <- MS + (1 + MS) * (MR + P_sum[1])
  mu_V <- MS + (1 + MS) * (MR + sum(P_sum[2:3]) + sum(random) + OLRE)
  if (psi_V == 1 && mu_V == 1) {
    cli::cli_abort(
      "The model was fit without any predictors or random effects and contains
      only random intercepts, so the variance partitions are 1."
    )
  } else {
    if (psi_V > 1) {
      V_lbl <- if (MS) "intercept: species"
      if (MR) {
        V_lbl <- c(V_lbl, paste0("intercept: region", suffix))
      }
      if (P_sum[1]) {
        X_lbl <- c(
          colnames(stan_data$X1),
          colnames(stan_data$X1_cat),
          colnames(stan_data$X1_ord)
        )
        V_lbl <- c(
          V_lbl,
          paste0(
            rep(X_lbl, ifelse(MS, 2, 1)),
            rep(suffix, each = length(X_lbl))
          )
        )
      }
      psi_phi <- tidybayes::spread_rvars(fit, psi_phi[v], psi_W) |>
        dplyr::mutate(
          v = factor(v, labels = V_lbl),
          submodel = "Occupancy",
          .before = 1
        ) |>
        dplyr::rename(partition = v, phi = psi_phi, W = psi_W)
    } else if (psi_V == 1) {
      cli::cli_warn(
        "The model was fit without any predictors for occupancy and contains
        only random intercepts, so the variance partition is 1 and not plotted."
      )
      psi_phi <- NULL
    }
    if (mu_V > 1) {
      V_lbl <- if (MS) "intercept: species"
      if (MR) {
        V_lbl <- c(V_lbl, paste0("intercept: region", suffix))
      }
      if (P_sum[2]) {
        X_lbl <- list(
          colnames(stan_data$X2),
          colnames(stan_data$X_cat2),
          colnames(stan_data$X_ord2)
        ) |>
          sapply(\(x) {
            if (!is.null(x)) {
              paste0(
                "site:",
                rep(x, ifelse(MS, 2, 1)),
                rep(suffix, each = length(x))
              )
            }
          }) |>
          unlist()
        V_lbl <- c(V_lbl, X_lbl)
      }
      if (P_sum[3]) {
        X_lbl <- list(
          dimnames(stan_data$X3)[[3]],
          dimnames(stan_data$X_cat3)[[3]],
          dimnames(stan_data$X_ord3)[[3]]
        ) |>
          sapply(\(x) {
            if (!is.null(x)) {
              paste0(
                "survey:",
                rep(x, ifelse(MS, 2, 1)),
                rep(suffix, each = length(x))
              )
            }
          }) |>
          unlist()
        V_lbl <- c(V_lbl, X_lbl)
      }
      if (random[1]) {
        V_lbl <- c(V_lbl, paste0("random: site", suffix))
      }
      if (random[2]) {
        V_lbl <- c(V_lbl, paste0("random: survey", suffix))
      }
      if (OLRE) {
        V_lbl <- c(V_lbl, paste0("random: OLRE", suffix))
      }
      mu_phi <- tidybayes::spread_rvars(fit, mu_phi[v], mu_W) |>
        dplyr::mutate(
          v = factor(v, labels = V_lbl),
          submodel = "Detection",
          .before = 1
        ) |>
        dplyr::rename(partition = v, phi = mu_phi, W = mu_W)
    } else if (mu_V == 1) {
      cli::cli_warn(
        "The model was fit without any predictors for detection and \\
          contains only random intercepts, so the variance partition is 1 and \\
          not plotted."
      )
      mu_phi <- NULL
    }
  }

  # join psi and mu and produce scales
  phi <- dplyr::bind_rows(psi_phi, mu_phi) |>
    dplyr::mutate(
      submodel = factor(submodel, levels = c("Occupancy", "Detection")),
      tau = sqrt(W * phi)
    )

  # plot
  p <- ggplot2::ggplot(phi) +
    ggplot2::aes(
      xdist = if (scales) tau else phi,
      y = forcats::fct_rev(partition)
    ) +
    ggplot2::facet_wrap(
      ~submodel,
      nrow = 1,
      scales = if (scales) "free" else "free_y"
    ) +
    ggdist::stat_pointinterval(...) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(
      x = "Estimate",
      y = if (scales) "Scale" else "Variance Partition"
    )
  attr(p, "plot_data") <- phi
  p
}


#' Plot realised occupancy proportions
#'
#' Plots species-level proportions of occupied sites,
#' \eqn{\frac{\sum_{i = 1}^I z_{is}}{I}}, potentially by region.
#'
#' @param fit A fitted model object from [occARU()].
#' @param by_region `logical`. Whether to plot by region, if multiple regions
#' were included.
#' @param species `character`. Vector of species to plot. If `NULL` (default),
#'   all species are plotted. Must be one of `attr(occARU_data, "species")`.
#' @param sites `character`. Vector of sites to use. If `NULL` (default), all
#'   sites are used. Must be one of `attr(occARU_data, "sites")`.
#' @param regions `character`. Vector of regions to use. If `NULL` (default),
#'   all regions are used. Must be one of `attr(occARU_data, "regions")`.
#' @param ... Additional arguments passed to [ggdist::stat_pointinterval()].
#'
#' @return A `ggplot` object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`plot_data`}{The tibble used to produce the plot.}
#'   }
#'
#' @seealso [occARU()], [plot_intercepts()], [plot_coefficients()],
#'   [plot_sites()], [plot_surveys()], [plot_correlations()],
#'   [plot_partitions()]
#'
#' @export
plot_realised <- function(
  fit,
  by_region = FALSE,
  species = NULL,
  sites = NULL,
  regions = NULL,
  ...
) {
  if (!inherits(fit, "CmdStanFit")) {
    cli::cli_abort(
      "{.arg fit} must be a {.cls CmdStanFit} object from {.fun occARU}."
    )
  }
  stan_data <- attr(fit, "stan_data")
  if (!stan_data$latent) {
    cli::cli_abort(
      "{.arg fit} does not have latent occupancy states because {.fun occARU}
      was run with {.arg latent = FALSE}."
    )
  }
  occARU_data <- attr(fit, "occARU_data")
  species_lvl <- attr(occARU_data, "species")
  species_idx <- indices(species, species_lvl)
  MS <- length(species_idx) > 1
  site_lvl <- attr(occARU_data, "sites")
  site_idx <- indices(sites, site_lvl)
  region_lvl <- attr(occARU_data, "regions")
  MR <- length(region_lvl) > 1
  if (!MR && by_region) {
    cli::cli_abort(
      "{.arg by_region = TRUE} is only applicable with multiple regions."
    )
  }
  region_idx <- indices(regions, region_lvl)

  # extract draws, label, and summarise
  draws <- tidybayes::spread_rvars(fit, z[i, s]) |>
    dplyr::mutate(r = stan_data$region[i]) |>
    dplyr::filter(i %in% site_idx, s %in% species_idx, r %in% region_idx) |>
    dplyr::mutate(
      region = factor(region_lvl[r], labels = region_lvl[region_idx]),
      species = factor(species_lvl[s], levels = species_lvl[species_idx])
    ) |>
    dplyr::summarise(
      realised = posterior::rvar_sum(z) / n(),
      .by = c("species", if (by_region) "region")
    )

  # plot
  p <- ggplot2::ggplot(draws) +
    ggplot2::aes(
      xdist = realised,
      y = forcats::fct_rev(if (by_region && !MS) region else species)
    ) +
    ggdist::stat_pointinterval(...) +
    ggplot2::scale_x_continuous(
      breaks = seq(0.2, 0.8, 0.2),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x = "Realised Occupancy",
      y = if (by_region && !MS) "Region" else "Species"
    )
  attr(p, "plot_data") <- draws
  if (by_region && MS) {
    p + facet_wrap(~region)
  } else {
    p
  }
}


#' Plot site occupancy and detection rates
#'
#' Plots species-level occupancy and detection rates combining site predictors
#' (if included via [make_data()]) and site random effects (`iota`).
#'
#' If `map = TRUE` and site coordinates are present in the fitted object,
#' effects are displayed as points sized by the magnitude of the detection rate
#' on a map; otherwise site effects are plotted as point-intervals. When
#' `latent = TRUE` was set in [occARU()], sites with median posterior
#' occupancy of 0 are shown as red crosses. When `latent = FALSE`, detection
#' rates are weighted by occupancy probability (`inv_logit(logit_psi[i, s])`).
#'
#' @param fit A fitted model object from [occARU()].
#' @param map `logical`. If `TRUE` (default), plot site effects summarised with
#'   posterior medians on a map using UTM coordinates. Requires site coordinates
#'   to have been supplied to [make_data()]. If `FALSE`, or if no coordinates
#'   are present, site effects are plotted with
#'   [ggdist::stat_pointinterval()].
#' @param intercepts `logical`. If `TRUE` (default), species-level baseline log
#'   detection rates are added to the site effects. If `FALSE`, only the site
#'   deviations are plotted on the log scale.
#' @param back_transform `logical`. Only used when `intercepts = TRUE`. If `TRUE`
#'   (default), log detection rates are back-transformed to the natural scale
#'   via `exp()`. If `FALSE`, values are left on the log scale.
#' @param include_predictors `logical`. If `TRUE` (default), includes predictors
#'   in the site effects, if included. If `FALSE`, only plots the random
#'   effects.
#' @param restricted `logical`. If `TRUE` (default), when `include_predictors` is
#'   `FALSE`, plots random site effects with orthogonal projection, i.e.,
#'   \eqn{(\boldsymbol{I} - \boldsymbol{P_{X_2}}) \boldsymbol{\iota}}, where
#'   \eqn{\boldsymbol{I} - \boldsymbol{P_{X_2}}} is the orthogonal complement of
#'   the column space of the site design matrix. If `FALSE`, plots random
#'   effects without orthogonal projection, i.e., \eqn{\boldsymbol{\iota}} only.
#'   Has no effect when `include_predictors` is `TRUE` as the linear predictor
#'   is unaffected by orthogonal projection.
#' @param species `character`. Vector of species to plot. If `NULL` (default),
#'   all species are plotted. Must be one of `attr(occARU_data, "species")`.
#' @param sites `character`. Vector of sites to plot. If `NULL` (default), all
#'   sites are plotted. Must be one of `attr(occARU_data, "sites")`.
#' @param regions `character`. Vector of regions to plot. If `NULL` (default),
#'   all regions are plotted. Must be one of `attr(occARU_data, "regions")`.
#' @param ndraws Positive integerish. Number of draws to use for plotting,
#'   passed to [tidybayes::spread_rvars()]. Default: `NULL` (uses all draws).
#' @param seed Positive numeric. Seed to use when subsampling draws when
#'   `ndraws` is not `NULL`. Default: random integer.
#' @param ... Additional arguments passed to [ggplot2::geom_point()] when
#'   `map = TRUE` and  [ggdist::stat_pointinterval()] when `map = FALSE`.
#'
#' @return A `ggplot` object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`plot_data`}{The tibble used to produce the plot.}
#'   }
#'
#' @seealso [occARU()], [plot_intercepts()], [plot_coefficients()],
#'   [plot_surveys()], [plot_correlations()], [plot_partitions()],
#'   [plot_realised()]
#'
#' @export
plot_sites <- function(
  fit,
  map = TRUE,
  intercepts = TRUE,
  back_transform = TRUE,
  include_predictors = TRUE,
  restricted = TRUE,
  species = NULL,
  sites = NULL,
  regions = NULL,
  ndraws = NULL,
  seed = NULL,
  ...
) {
  if (!inherits(fit, "CmdStanFit")) {
    cli::cli_abort(
      "{.arg fit} must be a {.cls CmdStanFit} object from {.fun occARU}."
    )
  }
  stan_data <- attr(fit, "stan_data")
  occARU_data <- attr(fit, "occARU_data")
  site_lvl <- attr(occARU_data, "sites")
  species_lvl <- attr(occARU_data, "species")
  region_lvl <- attr(occARU_data, "regions")
  survey_length <- attr(occARU_data, "survey_length")
  P <- occARU_data$P[2]
  P_cat <- occARU_data$P_cat[2]
  P_ord <- occARU_data$P_ord[2]
  P_sum <- sum(c(P, P_cat, P_ord))
  random <- stan_data$random[1]
  has_latent <- stan_data$latent == 1L
  transform <- back_transform && intercepts
  has_coords <- any(occARU_data$XY != 0)
  use_map <- map && has_coords
  if (!random && !P_sum) {
    cli::cli_abort(
      "Model was fit without site predictors or random effects."
    )
  }
  if (map && !has_coords) {
    cli::cli_warn(
      "No site coordinates found in {.arg data}. Falling back to
      point-intervals."
    )
  }
  if (!is.null(ndraws)) {
    if (!rlang::is_integerish(ndraws) || ndraws < 0) {
      cli::cli_abort(
        "{.arg ndraws} must be a positive integer."
      )
    } else if (is.null(seed)) {
      seed <- sample.int(.Machine$integer.max, 1)
    }
  }

  # get indices
  species_idx <- indices(species, species_lvl)
  site_idx <- indices(sites, site_lvl)
  region_idx <- indices(regions, region_lvl)

  # determine projection
  res <- TRUE
  if (!restricted) {
    if (P_sum && include_predictors) {
      cli::cli_abort(
        "{.arg restricted = FALSE} is only applicable when predictors are in
        the model but excluded from the plot. Set
        {.arg include_predictors = FALSE}."
      )
    } else if (!P_sum) {
      cli::cli_abort(
        "{.arg restricted = FALSE} is only applicable when predictors are
        included."
      )
    } else if (!stan_data$project[1]) {
      cli::cli_abort(
        "{.arg restricted = FALSE} is ignored when random site effects were
        not orthogonally projected. Requires
        {.arg project = list(site = TRUE)} in {.fun occARU}."
      )
    } else {
      res <- FALSE
    }
  }

  # initialise log_mu
  draws <- tidyr::expand_grid(i = site_idx, s = species_idx) |>
    dplyr::mutate(r = stan_data$region[i], .before = i) |>
    dplyr::filter(r %in% region_idx) |>
    dplyr::mutate(log_mu = 0)

  # increment predictor effects
  if (P_sum && include_predictors) {
    if (P) {
      betaX <- dplyr::left_join(
        tidybayes::spread_rvars(
          fit,
          mu_beta[p, s],
          ndraws = ndraws,
          seed = seed
        ) |>
          dplyr::filter(s %in% species_idx),
        predictor_matrix_to_tibble(stan_data$X2),
        by = dplyr::join_by(p),
        relationship = "many-to-many"
      ) |>
        dplyr::summarise(
          betaX = posterior::rvar_sum(mu_beta * x),
          .by = c(i, s)
        )

      draws <- dplyr::left_join(draws, betaX, by = dplyr::join_by(i, s)) |>
        dplyr::mutate(log_mu = log_mu + betaX) |>
        dplyr::relocate(betaX, .before = log_mu)
    }
    if (P_cat) {
      betaX_cat <- dplyr::left_join(
        tidybayes::spread_rvars(
          fit,
          mu_beta_cat[p, x, s],
          ndraws = ndraws,
          seed = seed
        ) |>
          dplyr::filter(s %in% species_idx, !is.nan(median(mu_beta_cat))),
        predictor_matrix_to_tibble(stan_data$X_cat2),
        by = dplyr::join_by(x, p),
        relationship = "many-to-many"
      ) |>
        dplyr::summarise(
          betaX_cat = posterior::rvar_sum(mu_beta_cat),
          .by = c(i, s)
        )

      draws <- dplyr::left_join(
        draws,
        betaX_cat,
        by = dplyr::join_by(i, s)
      ) |>
        dplyr::mutate(log_mu = log_mu + betaX_cat) |>
        dplyr::relocate(betaX_cat, .before = log_mu)
    }
    if (P_ord) {
      betaX_ord <- dplyr::left_join(
        tidybayes::spread_rvars(
          fit,
          mu_beta_ord_realised[p, x, s],
          ndraws = ndraws,
          seed = seed
        ) |>
          dplyr::filter(
            s %in% species_idx,
            !is.na(median(mu_beta_ord_realised))
          ),
        predictor_matrix_to_tibble(stan_data$X_ord2 + 1),
        by = dplyr::join_by(x, p),
        relationship = "many-to-many"
      ) |>
        dplyr::summarise(
          betaX_ord = posterior::rvar_sum(mu_beta_ord_realised),
          .by = c(i, s)
        )

      draws <- dplyr::left_join(
        draws,
        betaX_ord,
        by = dplyr::join_by(i, s)
      ) |>
        dplyr::mutate(log_mu = log_mu + betaX_ord) |>
        dplyr::relocate(betaX_ord, .before = log_mu)
    }
  }

  # increment random site effects
  if (random) {
    param <- rlang::sym(paste0("iota", if (!res) "2"))
    draws <- dplyr::left_join(
      draws,
      tidybayes::spread_rvars(
        fit,
        (!!param)[i, s],
        ndraws = ndraws,
        seed = seed
      ),
      by = dplyr::join_by(i, s)
    )

    draws <- draws |>
      dplyr::mutate(log_mu = log_mu + !!param) |>
      dplyr::relocate(!!param, .before = log_mu)
  }

  # join intercepts
  if (intercepts) {
    draws <- dplyr::left_join(
      draws,
      tidybayes::spread_rvars(
        fit,
        alpha[d, r, s],
        ndraws = ndraws,
        seed = seed
      ) |>
        dplyr::filter(d == 2) |>
        dplyr::select(-d),
      by = dplyr::join_by(r, s)
    )

    draws <- draws |>
      dplyr::mutate(log_mu = alpha + log_mu) |>
      dplyr::relocate(alpha, .after = s)
  } else if (back_transform) {
    cli::cli_warn(
      "{.arg back_transform = TRUE} has no effect when
      {.arg intercepts = FALSE}. Plotting log detection rates."
    )
  }

  # join with occupancy states if recovered or occupancy log odds
  if (has_latent) {
    draws <- dplyr::left_join(
      draws,
      tidybayes::spread_rvars(fit, z[i, s], ndraws = ndraws, seed = seed),
      by = dplyr::join_by(i, s)
    )

    draws <- draws |>
      dplyr::mutate(
        occupied = factor(median(z) == 1, levels = c(TRUE, FALSE)),
        pred = log(z * exp(log_mu))
      ) |>
      dplyr::relocate(z, occupied, .after = s)
  } else {
    draws <- dplyr::left_join(
      draws,
      tidybayes::spread_rvars(
        fit,
        logit_psi[i, s],
        ndraws = ndraws,
        seed = seed
      ),
      by = dplyr::join_by(i, s)
    )

    draws <- draws |>
      dplyr::relocate(logit_psi, .after = s) |>
      dplyr::mutate(pred = log(inv_logit(logit_psi) * exp(log_mu)))
  }

  # add labels
  draws <- draws |>
    dplyr::mutate(
      s = factor(species_lvl[s], levels = species_lvl[species_idx]),
      i = factor(site_lvl[i], levels = site_lvl[site_idx]),
      r = factor(region_lvl[r], labels = region_lvl[unique(r)]),
    ) |>
    dplyr::rename(species = s, site = i, region = r)

  # transform
  draws <- dplyr::mutate(draws, pred = if (transform) exp(log_mu) else log_mu)
  mu_lab <- paste0(
    "Detection rate per ",
    if (use_map) "\n",
    survey_length,
    "-day survey",
    if (!transform) " (log)"
  )

  # plot
  if (use_map) {
    draws <- dplyr::left_join(
      draws,
      as.data.frame(occARU_data$XY) |>
        tibble::rownames_to_column(var = "site"),
      by = "site"
    )
    p <- ggplot2::ggplot(draws) +
      ggplot2::aes(X, Y) +
      ggplot2::labs(x = "Easting (km)", y = "Northing (km)") +
      ggplot2::coord_fixed()
    if (has_latent) {
      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(
            colour = occupied,
            shape = occupied,
            size = dplyr::if_else(occupied == TRUE, median(pred), 0)
          ),
          ...
        ) +
        ggplot2::scale_colour_manual(values = c("#333333", "red4")) +
        ggplot2::scale_shape_manual(values = c(16, 4)) +
        ggplot2::labs(colour = "Occupied", shape = "Occupied", size = mu_lab)
    } else {
      p <- p +
        ggplot2::geom_point(ggplot2::aes(size = median(pred)), ...) +
        ggplot2::labs(size = mu_lab)
    }
  } else {
    p <- ggplot2::ggplot(draws) +
      ggplot2::aes(y = forcats::fct_rev(site)) +
      ggplot2::labs(x = mu_lab, y = "Site")
    if (has_latent) {
      p <- p +
        ggdist::stat_pointinterval(
          ggplot2::aes(
            xdist = dplyr::if_else(occupied == TRUE, pred, 0),
            colour = occupied,
            shape = occupied
          ),
          ...
        ) +
        ggplot2::scale_colour_manual(values = c("#333333", "red4")) +
        ggplot2::scale_shape_manual(values = c(16, 4)) +
        ggplot2::labs(colour = "Occupied", shape = "Occupied")
    } else {
      p <- p +
        ggdist::stat_pointinterval(
          ggplot2::aes(xdist = pred),
          ...
        )
    }
  }
  p <- p +
    ggplot2::facet_wrap(
      ~species,
      scales = if (!use_map && transform) "free_x"
    )
  attr(p, "plot_data") <- draws
  p
}

#' Plot temporal detection rates
#'
#' Plots species-level temporal detection rates combining survey predictors
#' averaged over sites (if included via [make_data()]) and survey random effects
#' (`kappa`).
#'
#' @param fit A fitted model object from [occARU()].
#' @param species `character`. Vector of species to plot. If `NULL` (default),
#'   all species are plotted. Must be one of `attr(occARU_data, "species")`.
#' @param surveys `character`. Vector of survey dates to plot. If `NULL`
#'   (default), all surveys are plotted. Must be one of
#'   `attr(occARU_data, "surveys")`.
#' @param intercepts `logical`. If `TRUE` (default), species-level baseline log
#'   detection rates are added to the survey effects. If `FALSE`, only the
#'   temporal deviations are plotted on the log scale.
#' @param back_transform `logical`. Only used when `intercepts = TRUE`. If `TRUE`
#'   (default), log detection rates are back-transformed to the natural scale
#'   via `exp()`. If `FALSE`, values are left on the log scale.
#' @param include_predictors `logical`. If `TRUE` (default), includes predictors
#'   in the survey effects, if included. If `FALSE`, only plots the random
#'   effects.
#' @param restricted `logical`. If `TRUE` (default), when `include_predictors`
#'   is `FALSE` and the model was fit `project = list(survey = TRUE)`, plots
#'   random survey effects with orthogonal projection, i.e.,
#'   \eqn{(\boldsymbol{I} - \boldsymbol{P_{X_3}}) \boldsymbol{\kappa}}, where
#'   \eqn{\boldsymbol{I} - \boldsymbol{P_{X_3}}} is the orthogonal complement of
#'   the column space of the site-averaged survey design matrix. If `FALSE`,
#'   plots random effects without orthogonal projection, i.e.,
#'   \eqn{\boldsymbol{\kappa}} only. Has no effect when `include_predictors`
#'   is `TRUE` as the linear predictor is unaffected by orthogonal projection.
#' @param ndraws Positive integer. Number of draws to use for plotting, passed
#'   to [tidybayes::spread_rvars()]. Default: `NULL` (uses all draws).
#' @param seed Positive numeric. Seed to use when subsampling draws when
#'   `ndraws` is not `NULL`. Default: random integer.
#' @param palette `character`. Colour palette to be passed to
#'   [ggplot2::scale_fill_brewer()]. Default: `"YlGn"`.
#' @param ... Additional arguments passed to [ggdist::stat_lineribbon()]
#'   such as `.width` and `point_interval`.
#'
#' @return A `ggplot` object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`plot_data`}{The tibble used to produce the plot.}
#'   }
#'
#' @seealso [occARU()], [plot_intercepts()], [plot_coefficients()],
#'   [plot_sites()], [plot_correlations()], [plot_partitions()],
#'   [plot_realised()]
#' @export
plot_surveys <- function(
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
) {
  if (!inherits(fit, "CmdStanFit")) {
    cli::cli_abort(
      "{.arg fit} must be a {.cls CmdStanFit} object from {.fun occARU}."
    )
  }
  stan_data <- attr(fit, "stan_data")
  occARU_data <- attr(fit, "occARU_data")
  survey_lvl <- attr(occARU_data, "surveys")$.survey
  species_lvl <- attr(occARU_data, "species")
  survey_length <- attr(occARU_data, "survey_length")
  transform <- back_transform && intercepts
  P <- occARU_data$P[3]
  P_cat <- occARU_data$P_cat[3]
  P_ord <- occARU_data$P_ord[3]
  P_sum <- sum(c(P, P_cat, P_ord))
  random <- stan_data$random[2]
  if (!random && !P_sum) {
    cli::cli_abort(
      "Model was fit without survey predictors or random effects."
    )
  }
  transform <- back_transform && intercepts
  if (!is.null(ndraws)) {
    if (!rlang::is_integerish(ndraws) || ndraws < 0) {
      cli::cli_abort(
        "{.arg ndraws} must be a positive integer."
      )
    } else if (is.null(seed)) {
      seed <- sample.int(.Machine$integer.max, 1)
    }
  }

  # get species and survey indices
  species_idx <- indices(species, species_lvl)
  survey_idx <- indices(surveys, survey_lvl)

  # determine projection
  res <- TRUE
  if (!restricted) {
    if (P_sum && include_predictors) {
      cli::cli_abort(
        "{.arg restricted = FALSE} is only applicable when predictors are in the
        model but excluded from the plot."
      )
    } else if (!P_sum) {
      cli::cli_abort(
        "{.arg restricted = FALSE} is only applicable when predictors are
        included."
      )
    } else if (!stan_data$project[2]) {
      cli::cli_abort(
        "{.arg restricted = FALSE} is ignored when random survey effects were
        not orthogonally projected. Requires
        {.arg project = list(survey = TRUE)} in {.fun occARU}."
      )
    } else {
      res <- FALSE
    }
  }

  # initialise log_mu
  draws <- tidyr::expand_grid(s = species_idx, j = survey_idx) |>
    dplyr::mutate(log_mu = 0)

  # increment predictor effects
  if (P_sum && include_predictors) {
    if (P) {
      gammaX <- dplyr::left_join(
        tidybayes::spread_rvars(
          fit,
          gamma[p, s],
          ndraws = ndraws,
          seed = seed
        ) |>
          dplyr::filter(s %in% species_idx),
        predictor_matrix_to_tibble(apply(stan_data$X3, 2:3, mean)) |>
          dplyr::rename(j = i) |>
          dplyr::filter(j %in% survey_idx),
        by = dplyr::join_by(p),
        relationship = "many-to-many"
      ) |>
        dplyr::summarise(
          gammaX = posterior::rvar_sum(gamma * x),
          .by = c(j, s)
        )

      draws <- dplyr::left_join(draws, gammaX, by = dplyr::join_by(j, s)) |>
        dplyr::mutate(log_mu = log_mu + gammaX) |>
        dplyr::relocate(gammaX, .before = log_mu)
    }
    if (P_cat) {
      gammaX_cat <- dplyr::left_join(
        tidybayes::spread_rvars(fit, gamma_cat[p, x, s]) |>
          dplyr::filter(s %in% species_idx),
        predictor_matrix_to_tibble(apply(stan_data$X_cat3, 2:3, int_mode)) |>
          dplyr::rename(j = i) |>
          dplyr::filter(j %in% survey_idx),
        by = dplyr::join_by(x, p),
        relationship = "many-to-many"
      ) |>
        dplyr::summarise(
          gammaX_cat = posterior::rvar_sum(gamma_cat),
          .by = c(j, s)
        )

      draws <- dplyr::left_join(
        draws,
        gammaX_cat,
        by = dplyr::join_by(j, s)
      ) |>
        dplyr::mutate(log_mu = log_mu + gammaX_cat) |>
        dplyr::relocate(gammaX_cat, .before = log_mu)
    }
    if (P_ord) {
      gammaX_ord <- dplyr::left_join(
        tidybayes::spread_rvars(
          fit,
          gamma_ord_realised[p, x, s],
          ndraws = ndraws,
          seed = seed
        ) |>
          dplyr::filter(s %in% species_idx),
        predictor_matrix_to_tibble(apply(
          stan_data$X_ord3 + 1,
          2:3,
          int_mode
        )) |>
          dplyr::rename(j = i) |>
          dplyr::filter(j %in% survey_idx),
        by = dplyr::join_by(x, p),
        relationship = "many-to-many"
      ) |>
        dplyr::summarise(
          gammaX_ord = posterior::rvar_sum(gamma_ord_realised),
          .by = c(j, s)
        )

      draws <- dplyr::left_join(
        draws,
        gammaX_ord,
        by = dplyr::join_by(j, s)
      ) |>
        dplyr::mutate(log_mu = log_mu + gammaX_ord) |>
        dplyr::relocate(gammaX_ord, .before = log_mu)
    }
  }

  # increment random survey effects
  if (random) {
    param <- rlang::sym(paste0("kappa", if (!res) "2"))
    draws <- dplyr::left_join(
      draws,
      tidybayes::spread_rvars(
        fit,
        (!!param)[j, s],
        ndraws = ndraws,
        seed = seed
      ),
      by = dplyr::join_by(j, s)
    ) |>
      dplyr::mutate(log_mu = log_mu + !!param) |>
      dplyr::relocate(!!param, .before = log_mu)
  }

  # join intercepts
  if (intercepts) {
    draws <- dplyr::left_join(
      draws,
      tidybayes::spread_rvars(
        fit,
        alpha[d, r, s],
        ndraws = ndraws,
        seed = seed
      ) |>
        dplyr::filter(d == 2) |>
        dplyr::select(-d) |>
        dplyr::summarise(alpha = posterior::rvar_mean(alpha), .by = s),
      by = dplyr::join_by(s)
    ) |>
      dplyr::mutate(log_mu = alpha + log_mu) |>
      dplyr::relocate(alpha, .after = j)
  } else if (back_transform) {
    cli::cli_warn(
      "{.arg back_transform = TRUE} has no effect when \\
      {.arg intercepts = FALSE}. Plotting log detection rates."
    )
  }

  # add labels
  draws <- draws |>
    dplyr::mutate(
      s = factor(species_lvl[s], levels = species_lvl[species_idx]),
      j = lubridate::ymd(survey_lvl[j])
    ) |>
    dplyr::rename(species = s, survey = j)

  # transform and plot
  draws <- dplyr::mutate(draws, pred = if (transform) exp(log_mu) else log_mu)
  p <- ggplot2::ggplot(draws) +
    ggplot2::aes(x = survey, ydist = pred) +
    ggdist::stat_lineribbon(...) +
    ggplot2::facet_wrap(~species, scales = if (transform) "free_y") +
    ggplot2::scale_x_date(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::labs(
      x = "Date",
      y = paste0(
        "Detection rate per ",
        survey_length,
        "-day survey",
        if (!transform) " (log)"
      )
    )
  attr(p, "plot_data") <- draws
  p
}

#' Convert `[P, I]` design matrix in `stan_data` object to long format
#'
#' @param X Design matrix
#' @return A `tibble` with continuous or integer valued predictors
#' @noRd
predictor_matrix_to_tibble <- function(X) {
  categorical <- all(X == round(X))
  P <- ncol(X)
  tidyr::as_tibble(X) |>
    purrr::set_names(1:P) |>
    dplyr::mutate(i = dplyr::row_number()) |>
    tidyr::pivot_longer(-i, names_to = "p", values_to = "x") |>
    dplyr::mutate(p = as.integer(p), x = if (categorical) as.integer(x) else x)
}


#' Produce indices based on supplied character vector of names
#'
#' @param names Character vector of names. Default: `NULL`
#' @param levels Levels of orginal data
#' @return A vector of integers corresponding to names to retain
#' @noRd
indices <- function(names = NULL, levels) {
  arg_name <- deparse(substitute(names))
  if (!is.null(names)) {
    unknown <- setdiff(names, levels)
    if (length(unknown)) {
      cli::cli_abort(
        "The following {arg_name} were not found in the model: \\
        {.val {unknown}}."
      )
    }
    idx <- which(levels %in% names)
  } else {
    idx <- 1:length(levels)
  }
  idx
}

#' occARU ggplot2 theme
#'
#' A clean ggplot2 theme used by occARU plotting functions. Can also be applied
#' to custom plots with `+ theme_occARU()`.
#'
#' @param base_size Base font size in pts. Default: `11`.
#' @param base_family Base font family. Default: `""`.
#'
#' @return A `ggplot2` theme object.
#' @importFrom ggplot2 %+replace%
#' @noRd
theme_occARU <- function(base_size = 11, base_family = "") {
  my_black <- "#333333"
  half_line <- base_size / 2
  base_line_size <- base_size / 22

  ggplot2::theme_bw(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      # axes
      axis.text = ggplot2::element_text(
        colour = my_black,
        size = ggplot2::rel(0.9)
      ),
      axis.title = ggplot2::element_text(colour = my_black),
      axis.ticks = ggplot2::element_line(
        colour = my_black,
        linewidth = base_line_size
      ),

      # legend
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.9)),

      # panel
      panel.border = ggplot2::element_rect(
        fill = NA,
        colour = my_black,
        linewidth = base_line_size
      ),
      panel.grid = ggplot2::element_blank(),

      # facet strips
      strip.background = ggplot2::element_rect(
        fill = my_black,
        colour = my_black
      ),
      strip.text = ggplot2::element_text(
        colour = "white",
        size = ggplot2::rel(1),
        margin = ggplot2::margin(
          0.8 * half_line,
          0.8 * half_line,
          0.8 * half_line,
          0.8 * half_line
        )
      ),

      # text
      text = ggplot2::element_text(colour = my_black),
      complete = TRUE
    )
}

#' Inverse logit function
#'
#' @param x Logit
#' @return Probability
#' @noRd
inv_logit <- function(x) 1 / (1 + exp(-x))

#' Plot red dashed verticle 0-line
#'
#' @param p Plot
#' @return Plot with ggplot2::geom_vline()
#' @noRd
my_vline <- function() {
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "red4")
}
