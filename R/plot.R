#' Plot predictor coefficients
#'
#' Plots coefficients of continuous, categorical, or ordinal predictors on
#' occupancy or detection submodels.
#'
#' @param fit A fitted model object from [fit_model()].
#' @param submodel Character. Predictors of submodel to plot. One of
#'   `"detection"` (default) or `"occupancy"`.
#' @param component Character. Whether to plot `"site"` (default) or `"survey"`
#'   predictors. If `"survey"`, `submodel` must be `"detection"`.
#' @param type Character. Type of predictors to plot. One of `"continuous"`
#'   (default), `"categorical"`, or `"ordinal"`.
#' @param level Character. For multi-species models, whether to plot
#'   species-specific (`"species"`, default) or mean coefficients (`"mean"`).
#' @param facet_by Character. Whether to use [ggplot2::facet_wrap()] or
#'   [ggh4x::facet_grid2()] to facet by `"predictor"` (default) or `"species"`.
#'   Only used if `level` is `"species"`.
#' @param species Character vector of species to plot. If `NULL` (default), all
#'   species are plotted. Must be one of `attr(occARU_data, "species")`.
#' @param restricted Logical. If `TRUE` (default), plots coefficients with
#'   orthogonal projection of the detection random site or survey effects, i.e.,
#'   \eqn{\boldsymbol{\iota}(\boldsymbol{I} - \boldsymbol{P}_X)}, where
#'   \eqn{\boldsymbol{I} - \boldsymbol{P}_X} is the orthogonal complement of the
#'   column space of the site or survey design matrix. If `FALSE`, recovers
#'   coefficients without orthogonal projection, \eqn{\boldsymbol{\beta} -
#'   \boldsymbol{\iota} \boldsymbol{X}^+}, where  \eqn{\boldsymbol{X}^+} is the
#'   pseudo-inverse of the design matrix. Only used if `submodel` is
#'   `"detection"`.
#' @param ordinal_categories Logical. If `FALSE` (default), plots coefficients
#'   associated with maximum category (full effect). If `TRUE`, plots realised
#'   coefficient associated with each ordered category, where the first
#'   is used as the reference.
#' @param ... Additional arguments passed to [tidybayes::stat_pointinterval()].
#'
#' @return A `ggplot` object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`plot_data`}{The data structure used to produce the plot.}
#'   }
#'
#' @seealso [fit_model()], [plot_intercepts()], [plot_sites()],
#'   [plot_surveys()], [plot_correlations()], [plot_partitions()]
#' @export
plot_coefficients <- function(
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
) {
  if (!inherits(fit, "CmdStanFit")) {
    cli::cli_abort(
      "{.arg fit} must be a {.cls CmdStanFit} object from {.fun fit_model}."
    )
  }
  rlang::check_installed(c("tidybayes", "ggplot2"))
  stan_data <- attr(fit, "stan_data")
  occARU_data <- attr(fit, "occARU_data")
  species_lvl <- attr(occARU_data, "species")
  MS <- length(species_lvl) > 1
  submodel <- match.arg(submodel, c("detection", "occupancy"))
  component <- match.arg(component, c("site", "survey"))
  type <- match.arg(type, c("continuous", "categorical", "ordinal"))
  level <- match.arg(level, c("species", "mean"))
  facet_by <- match.arg(facet_by, c("predictor", "species"))

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

  # get restricted indicator
  res <- TRUE
  if (!restricted) {
    if (submodel == "occupancy") {
      cli::cli_warn(
        "{.arg restricted} is ignored when {.arg submodel} is \\
        {.val {submodel}}."
      )
    } else {
      res <- FALSE
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
  idx <- if (submodel == "occupancy") {
    1L
  } else if (survey) {
    3L
  } else {
    2L
  }
  suffix <- if (!res) "_unc" else NULL

  # total number of predictors
  P <- occARU_data$P[idx]
  P_cat <- occARU_data$P_cat[idx]
  P_ord <- occARU_data$P_ord[idx]
  where <- if (!survey) paste("for", submodel) else ""
  if (!sum(P, P_cat, P_ord)) {
    cli::cli_abort(
      "No {component} predictors were found {where} in {.arg fit}."
    )
  }

  # get species indices
  if (level == "species") {
    species_idx <- species_indices(species, species_lvl)
  } else {
    if (!is.null(species)) {
      cli::cli_warn(
        "{.arg facet_by} and {.arg species} are ignored when {.arg level} is \\
        {.val {level}}."
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
      X_lbl <- if (survey) colnames(X) else rownames(X)
      if (MS && level == "species") {
        param <- rlang::sym(paste0(prefix, coef, suffix))
        draws <- tidybayes::spread_rvars(fit, (!!param)[s, p]) |>
          dplyr::filter(s %in% species_idx) |>
          dplyr::mutate(
            s = factor(species_lvl[s], levels = species_lvl[species_idx]),
            p = factor(p, labels = X_lbl)
          ) |>
          dplyr::rename(species = s, predictor = p)
        p <- ggplot2::ggplot(draws) +
          ggplot2::aes(xdist = !!param) +
          ggplot2::facet_wrap(
            ~ if (facet_by == "species") {
              species
            } else {
              predictor
            },
            scales = if (facet_by == "species") "fixed" else "free_x"
          ) +
          my_vline() +
          tidybayes::stat_pointinterval(
            ggplot2::aes(
              y = if (facet_by == "species") {
                forcats::fct_rev(predictor)
              } else {
                forcats::fct_rev(species)
              }
            ),
            ...
          ) +
          ggplot2::labs(
            x = "Coefficient",
            y = if (facet_by == "species") "Predictor" else "Species"
          )
      } else {
        param <- rlang::sym(paste0(prefix, coef, if (MS) "_bar", suffix))
        draws <- tidybayes::spread_rvars(fit, (!!param)[p]) |>
          dplyr::mutate(
            p = factor(p, labels = X_lbl)
          ) |>
          dplyr::rename(predictor = p)
        p <- ggplot2::ggplot(draws) +
          ggplot2::aes(xdist = !!param, y = forcats::fct_rev(predictor)) +
          my_vline() +
          tidybayes::stat_pointinterval(...) +
          ggplot2::labs(
            x = "Coefficient",
            y = "Predictor"
          )
      }
    } else {
      cli::cli_abort(
        "No {type} {component} predictors were found {where} in {.arg fit}."
      )
    }

    # factors
  } else {
    rlang::check_installed("ggh4x")
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

        if (MS && level == "species") {
          param <- rlang::sym(paste0(prefix, coef, "_cat", suffix))
          draws <- tidybayes::spread_rvars(fit, (!!param)[p, s, c]) |>
            dplyr::filter(s %in% species_idx) |>
            dplyr::right_join(lvls_df, by = dplyr::join_by(p, c)) |>
            dplyr::mutate(
              species = factor(
                species_lvl[s],
                levels = species_lvl[species_idx]
              ),
              predictor = factor(predictor, levels = names(lvls$cat)),
              category = factor(category, levels = unlist(lvls$cat)),
            ) |>
            dplyr::select(-c(p, s, c)) |>
            dplyr::relocate(
              species,
              predictor,
              category,
              .before = !!param
            )

          p <- ggplot2::ggplot(draws) +
            ggplot2::aes(xdist = !!param, y = category) +
            ggh4x::facet_grid2(
              if (facet_by == "species") {
                species ~ predictor
              } else {
                predictor ~ species
              },
              scales = "free_y",
              independent = "y"
            ) +
            my_vline() +
            tidybayes::stat_pointinterval(...) +
            ggplot2::labs(x = "Coefficient", y = "Category")
        } else {
          param <- rlang::sym(paste0(
            prefix,
            coef,
            "_cat",
            if (MS) "_bar",
            suffix
          ))
          draws <- tidybayes::spread_rvars(fit, (!!param)[p, c]) |>
            dplyr::filter(!is.na(median(!!param))) |>
            dplyr::left_join(lvls_df, by = dplyr::join_by(p, c)) |>
            dplyr::mutate(
              predictor = factor(predictor, levels = names(lvls$cat)),
              category = factor(category, levels = unlist(lvls$cat))
            ) |>
            dplyr::select(-c(p, c)) |>
            dplyr::relocate(predictor, category, .before = !!param)

          p <- ggplot2::ggplot(draws) +
            ggplot2::aes(xdist = !!param, y = forcats::fct_rev(category)) +
            ggplot2::facet_wrap(~predictor, scales = "free_y") +
            my_vline() +
            tidybayes::stat_pointinterval(...) +
            ggplot2::labs(x = "Coefficient", y = "Category")
        }
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
        if (MS && level == "species") {
          param <- rlang::sym(paste0(prefix, coef, "_ord", suffix))
          draws <- tidybayes::spread_rvars(
            fit,
            (!!param)[s, p],
            (!!cs)[p, o]
          ) |>
            dplyr::filter(s %in% species_idx) |>
            dplyr::right_join(lvls_df, by = dplyr::join_by(p, o)) |>
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
                if (facet_by == "species") {
                  species ~ predictor
                } else {
                  predictor ~ species
                },
                scales = "free_y",
                independent = "y"
              ) +
              my_vline() +
              tidybayes::stat_pointinterval(...) +
              ggplot2::labs(x = "Coefficient", y = "Category")
          } else {
            p <- ggplot2::ggplot(draws) +
              ggplot2::aes(
                xdist = !!param,
                y = if (facet_by == "species") {
                  forcats::fct_rev(predictor)
                } else {
                  forcats::fct_rev(species)
                }
              ) +
              ggplot2::facet_wrap(
                ~ if (facet_by == "species") {
                  species
                } else {
                  predictor
                },
                scales = if (facet_by == "species") "fixed" else "free_x"
              ) +
              my_vline() +
              tidybayes::stat_pointinterval(...) +
              ggplot2::labs(
                x = "Coefficient",
                y = if (facet_by == "species") "Predictor" else "Species"
              )
          }
        } else {
          param <- rlang::sym(paste0(
            prefix,
            coef,
            "_ord",
            if (MS) "_bar",
            suffix
          ))
          draws <- tidybayes::spread_rvars(fit, (!!param)[p], (!!cs)[p, o]) |>
            dplyr::filter(!is.na(median(!!param))) |>
            dplyr::right_join(lvls_df, by = dplyr::join_by(p, o)) |>
            dplyr::mutate(
              predictor = factor(predictor, levels = names(lvls$ord))
            ) |>
            dplyr::mutate(
              category = factor(category, levels = lvls$ord[[predictor[1]]]),
              .by = predictor
            ) |>
            dplyr::select(-c(p, o)) |>
            dplyr::relocate(predictor, category, .before = !!param)

          if (ordinal_categories) {
            p <- ggplot2::ggplot(draws) +
              ggplot2::aes(
                xdist = !!param * !!cs,
                y = forcats::fct_rev(category)
              ) +
              ggplot2::facet_wrap(~predictor, scales = "free_y") +
              my_vline() +
              tidybayes::stat_pointinterval(...) +
              ggplot2::labs(x = "Coefficient", y = "Category")
          } else {
            p <- ggplot2::ggplot(draws) +
              ggplot2::aes(xdist = !!param, y = forcats::fct_rev(predictor)) +
              my_vline() +
              tidybayes::stat_pointinterval(...) +
              ggplot2::labs(x = "Coefficient", y = "Predictor")
          }
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
#' Plot pairwise species correlations for different model components.For
#' occupancy, interspecific correlations are only estimated for responses
#' to site-level predictors. For detection, interspecific correlations are
#' potentially estimated for responses to site-level predictors, survey-level
#' predictors, random site effects, random survey effects, and potentially
#' Poisson OLREs.
#'
#' @param fit A fitted model object from [fit_model()].
#' @param submodel Character. Correlations of submodel to plot. One of
#'   `"detection"` (default) or `"occupancy"`.
#' @param species Character vector of species to plot. If `NULL` (default), all
#'   species are plotted. Must be one of `attr(occARU_data, "species")`.
#' @param ... Additional arguments passed to [tidybayes::stat_pointinterval()].
#'
#' @return A `ggplot` object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`plot_data`}{The data structure used to produce the plot.}
#'   }
#'
#' @seealso [fit_model()], [plot_intercepts()], [plot_coefficients()],
#'    [plot_sites()], [plot_surveys()], [plot_partitions()]
#' @export
plot_correlations <- function(
  fit,
  submodel = c("detection", "occupancy"),
  species = NULL,
  ...
) {
  if (!inherits(fit, "CmdStanFit")) {
    cli::cli_abort(
      "{.arg fit} must be a {.cls CmdStanFit} object from {.fun fit_model}."
    )
  }
  rlang::check_installed(c("tidybayes", "ggplot2"))
  stan_data <- attr(fit, "stan_data")
  if (stan_data$S == 1) {
    cli::cli_abort(
      "No interspecific correlations exist in single species models."
    )
  }
  occARU_data <- attr(fit, "occARU_data")
  species_lvl <- attr(occARU_data, "species")
  species_idx <- species_indices(species, species_lvl)
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
    SP <- stan_data$spatial > 0
    TE <- stan_data$temporal > 0
    OLRE <- stan_data$OD == 1
    if (sum(P) + SP + TE + OLRE) {
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
      if (SP) {
        draws <- dplyr::left_join(
          draws,
          tidybayes::spread_rvars(fit, iota_O[s, ss]),
          by = dplyr::join_by(s, ss)
        )
      }
      if (TE) {
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
        "No predictors or random effects were found for \\
                     detection in {.arg fit}."
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
      across(c(s, ss), ~ factor(., labels = species_lvl)),
      .variable = factor(
        .variable,
        levels = c(str_c(
          c("psi_beta", "mu_beta", "gamma", "iota", "kappa", "epsilon"),
          "_O"
        )),
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
    tidybayes::stat_pointinterval(...) +
    ggplot2::scale_x_continuous(
      breaks = seq(-0.5, 0.5, 0.5),
      limits = c(-1, 1),
      expand = c(0, 0)
    ) +
    labs(x = "Correlation", y = "Species")
  attr(p, "plot_data") <- draws
  p
}

#' Plot intercepts
#'
#' Plot species-specific intercepts for occupancy and detection, by default
#' back-transformed to the orginal scale.
#'
#' @param fit A fitted model object from [fit_model()].
#' @param species Character vector of species to plot. If `NULL` (default), all
#'   species are plotted. Must be one of `attr(occARU_data, "species")`.
#' @param back_transform Logical. If `TRUE` (default), intercepts are
#'   back-transformed to the natural scale via `inv_logit()` for occupancy and
#'   `exp()` for detection rates. If `FALSE`, values are left on the scale of
#'   the link functions (logit for occupancy and log for detection).
#'
#' @param ... Additional arguments passed to [tidybayes::stat_pointinterval()].
#'
#' @return A `ggplot` object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`plot_data`}{The data structure used to produce the plot.}
#'   }
#'
#' @seealso [fit_model()], [plot_coefficients()], [plot_sites()],
#'    [plot_surveys()], [plot_correlations()], [plot_partitions()]
#' @export
plot_intercepts <- function(
  fit,
  species = NULL,
  back_transform = TRUE,
  ...
) {
  if (!inherits(fit, "CmdStanFit")) {
    cli::cli_abort(
      "{.arg fit} must be a {.cls CmdStanFit} object from {.fun fit_model}."
    )
  }
  rlang::check_installed(c("tidybayes", "ggplot2"))
  stan_data <- attr(fit, "stan_data")
  occARU_data <- attr(fit, "occARU_data")
  MS <- stan_data$S > 1

  # extract intercepts and label
  draws <- if (MS) {
    species_lvl <- attr(occARU_data, "species")
    species_idx <- species_indices(species, species_lvl)
    tidybayes::spread_rvars(fit, alpha[s, d]) |>
      dplyr::filter(s %in% species_idx) |>
      dplyr::mutate(
        s = factor(species_lvl[s], levels = species_lvl[species_idx])
      ) |>
      dplyr::rename(species = s)
  } else {
    if (!is.null(species)) {
      cli::cli_warn(
        "{.arg species} is ignored when single species model was fit."
      )
    }
    tidybayes::spread_rvars(fit, alpha[d])
  }

  # transform
  if (back_transform) {
    inv_logit <- function(x) 1 / (1 + exp(-x))
    draws <- dplyr::mutate(
      draws,
      alpha = dplyr::if_else(d == 1, inv_logit(alpha), exp(alpha)),
      d = factor(d, labels = c("Occupancy Probability", "Detection Rate"))
    )
  } else {
    draws <- dplyr::mutate(
      draws,
      d = factor(d, labels = c("Occupancy Log Odds", "Log Detection Rate"))
    )
  }
  draws <- dplyr::rename(draws, submodel = d)

  # plot
  p <- ggplot2::ggplot(draws) +
    ggplot2::aes(xdist = alpha) +
    ggplot2::labs(x = "Estimate")
  if (MS) {
    p <- p +
      ggplot2::facet_wrap(~submodel, nrow = 1, scales = "free_x") +
      tidybayes::stat_pointinterval(
        ggplot2::aes(y = forcats::fct_rev(species)),
        ...
      ) +
      ggplot2::labs(y = "Species")
  } else {
    p <- p +
      tidybayes::stat_pointinterval(
        ggplot2::aes(y = forcats::fct_rev(submodel)),
        ...
      ) +
      ggplot2::labs(y = "Submodel")
  }
  attr(p, "plot_data") <- draws
  p
}


#' Plot variance partitions
#'
#' The occARU model uses global-local shrinkage priors for the occupancy and
#' detection submodels, where half Student-t priors are used for the variances
#' of both linear predictors which are simplex partitioned via either Dirichlet
#' or logistic-normal decomposition. Variance decomposition only occurs when
#' there is more than one model component. Partitions exist for species-level
#' intercepts, species-level slopes (with one mean partition and one partition
#' for species-level deviations), species-level site effects (with mean and
#' species-specific components), species-level survey effects (with mean and
#' species-specific components), and Poisson OLREs (with mean and
#' species-specific components).
#'
#' @param fit A fitted model object from [fit_model()].
#' @param scales Logical. If `FALSE` (default), plots variance simplex
#'   partitions \eqn{\phi}. If `TRUE`, produces component scales by plotting
#'   \eqn{\sqrt{W \cdot \phi}}, where \eqn{W} are variances of linear
#'   predictors.
#' @param ... Additional arguments passed to [tidybayes::stat_pointinterval()].
#'
#' @return A `ggplot` object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`plot_data`}{The data structure used to produce the plot.}
#'   }
#'
#' @seealso [fit_model()] [plot_intercepts()], [plot_coefficients()],
#'    [plot_sites()], [plot_surveys()], [plot_correlations()]
#' @export
plot_partitions <- function(
  fit,
  scales = FALSE,
  ...
) {
  if (!inherits(fit, "CmdStanFit")) {
    cli::cli_abort(
      "{.arg fit} must be a {.cls CmdStanFit} object from {.fun fit_model}."
    )
  }
  rlang::check_installed(c("tidybayes", "ggplot2"))
  stan_data <- attr(fit, "stan_data")
  occARU_data <- attr(fit, "occARU_data")
  species_lvl <- attr(occARU_data, "species")
  S <- stan_data$S
  MS <- S > 1
  SP <- stan_data$spatial > 0
  TE <- stan_data$temporal > 0
  OLRE <- stan_data$OD == 1
  P_sum <- stan_data$P + stan_data$P_cat + stan_data$P_ord

  # labels
  species_lbl <- paste0("(", species_lvl, ")")
  suffix <- paste0("(", c("Mean", "Species"), ")")

  # multispecies
  if (MS) {
    psi_V <- 1 + 2 * P_sum[1]
    mu_V <- 1 + 2 * sum(P_sum[2:3]) + (SP + TE + OLRE) * (S + 1)
    if (psi_V == 1 && mu_V == 1) {
      cli::cli_abort(
        "The model was fit without any predictors or random effects and \\
        contains only random intercepts, so the variance partitions are 1."
      )
    } else {
      if (psi_V > 1) {
        V_lbl <- "Intercept"
        if (P_sum[1]) {
          X_lbl <- c(
            rownames(stan_data$X1),
            rownames(stan_data$X1_cat),
            rownames(stan_data$X1_ord)
          )
          V_lbl <- c(
            V_lbl,
            paste(rep(X_lbl, 2), rep(suffix, each = length(X_lbl)))
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
          "The model was fit without any predictors for occupancy and \\
          contains only random intercepts, so the variance partition is 1 and \\
          not plotted."
        )
        psi_phi <- NULL
      }
      if (mu_V > 1) {
        V_lbl <- "Intercept"
        if (P_sum[2]) {
          X_lbl <- list(
            rownames(stan_data$X2),
            rownames(stan_data$X_cat2),
            rownames(stan_data$X_ord2)
          ) |>
            sapply(\(x) {
              if (!is.null(x)) {
                paste("site:", rep(x, 2), rep(suffix, each = length(x)))
              }
            }) |>
            unlist()
          V_lbl <- c(V_lbl, X_lbl)
        }
        if (P_sum[3]) {
          X_lbl <- list(
            colnames(stan_data$X3),
            colnames(stan_data$X_cat3),
            colnames(stan_data$X_ord3)
          ) |>
            sapply(\(x) {
              if (!is.null(x)) {
                paste("survey:", rep(x, 2), rep(suffix, each = length(x)))
              }
            }) |>
            unlist()
          V_lbl <- c(V_lbl, X_lbl)
        }
        if (SP) {
          V_lbl <- c(V_lbl, paste("Site Effects", c("(Global)", species_lbl)))
        }
        if (TE) {
          V_lbl <- c(V_lbl, paste("Survey Effects", c("(Global)", species_lbl)))
        }
        if (OLRE) {
          V_lbl <- c(V_lbl, paste("OLRE", c("(Global)", species_lbl)))
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

    # single species
  } else {
    psi_V <- P_sum[1]
    mu_V <- sum(P_sum[2:3]) + SP + TE + OLRE
    if (psi_V == 0 && mu_V == 0) {
      cli::cli_abort(
        "The model was fit without any predictors or random effects, \\
        so there are no variance partitions."
      )
    }
    psi_phi <- NULL
    if (psi_V > 1) {
      V_lbl <- NULL
      if (P_sum[1]) {
        X_lbl <- c(
          rownames(stan_data$X1),
          rownames(stan_data$X1_cat),
          rownames(stan_data$X1_ord)
        )
        V_lbl <- c(V_lbl, X_lbl)
      }
      psi_phi <- tidybayes::spread_rvars(fit, psi_phi[v], psi_W[d]) |>
        dplyr::mutate(
          v = factor(v, labels = V_lbl),
          submodel = "Occupancy",
          .before = 1
        ) |>
        dplyr::rename(partition = v, phi = psi_phi, W = psi_W)
    } else if (psi_V == 1) {
      cli::cli_warn(
        "The model was fit with one predictor for occupancy, so the \\
          variance partition is 1 and not plotted."
      )
    }
    mu_phi <- NULL
    if (mu_V > 1) {
      V_lbl <- NULL
      if (P_sum[2]) {
        X_lbl <- c(
          rownames(stan_data$X2),
          rownames(stan_data$X_cat2),
          rownames(stan_data$X_ord2)
        )
        V_lbl <- c(V_lbl, paste("site:", X_lbl))
      }
      if (P_sum[3]) {
        X_lbl <- c(
          colnames(stan_data$X3),
          colnames(stan_data$X_cat3),
          colnames(stan_data$X_ord3)
        )
        V_lbl <- c(V_lbl, paste("survey:", X_lbl))
      }
      if (SP) {
        V_lbl <- c(V_lbl, "Site Effects")
      }
      if (TE) {
        V_lbl <- c(V_lbl, "Survey Effects")
      }
      if (OLRE) {
        V_lbl <- c(V_lbl, "OLRE")
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
        "The model was fit with only one component for detection, \\
          so the variance partition is 1 and not plotted."
      )
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
    ggplot2::facet_wrap(~submodel, nrow = 1, scales = "free_y") +
    tidybayes::stat_pointinterval(...) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(
      x = "Estimate",
      y = if (scales) "Scale" else "Variance Partition"
    )
  attr(p, "plot_data") <- phi
  p
}


#' Plot site occupancy and detection rates
#'
#' Plots species-level detection rates combining site predictors (if included
#' via [make_data()]) and site random effects (`iota`). If `map = TRUE` and
#' site coordinates are present in the fitted object, effects are displayed as
#' points sized by the magnitude of the detection rate on a map; otherwise site
#' effects are plotted as point-intervals. When `latent = TRUE` was set in
#' [fit_model()], sites with median posterior occupancy of 0 are shown as red
#' crosses. When `latent = FALSE`, detection rates are weighted by occupancy
#' probability (`inv_logit(logit_psi[s, i])`).
#'
#' @param fit A fitted model object from [fit_model()].
#' @param species Character vector of species to plot. If `NULL` (default), all
#'   species are plotted. Must be one of `attr(occARU_data, "species")`.
#' @param sites Character vector of sites to plot. If `NULL` (default), all
#'   sites are plotted. Must be one of `attr(occARU_data, "sites")`.
#' @param map Logical. If `TRUE` (default), plot site effects summarised with
#'   posterior medians on a map using UTM coordinates. Requires site coordinates
#'   to have been supplied to [make_data()]. If `FALSE`, or if no coordinates
#'   are present, site effects are plotted with
#'   [tidybayes::stat_pointinterval()].
#' @param intercepts Logical. If `TRUE` (default), species-level baseline log
#'   detection rates are added to the site effects. If `FALSE`, only the site
#'   deviations are plotted on the log scale.
#' @param back_transform Logical. Only used when `intercepts = TRUE`. If `TRUE`
#'   (default), log detection rates are back-transformed to the natural scale
#'   via `exp()`. If `FALSE`, values are left on the log scale.
#' @param include_predictors Logical. If `TRUE` (default), includes predictors
#'   in the site effects, if included. If `FALSE`, only plots the random
#'   effects.
#' @param restricted Logical. If `TRUE` (default), when `include_predictors` is
#'   `FALSE`, plots random site effects with orthogonal projection, i.e.,
#'   \eqn{\boldsymbol{\iota}(\boldsymbol{I} - \boldsymbol{P}_X)}, where
#'   \eqn{\boldsymbol{I} - \boldsymbol{P}_X} is the orthogonal complement of the
#'   column space of the site design matrix. If `FALSE`, plots random effects
#'   without orthogonal projection, i.e., \eqn{\boldsymbol{\iota}} only. Has
#'   no effect when `include_predictors` is `TRUE` as the linear predictor
#'   is unaffected by orthogonal projection.
#' @param ndraws  Number of draws to use for plotting, passed to
#'   [tidybayes::spread_rvars()]. Default: `NULL` (uses all draws).
#' @param seed Positive numeric. Seed to use when subsampling draws when
#'   `ndraws` is not `NULL`. Default: random integer.
#' @param ... Additional arguments passed to [ggplot2::geom_point()] when
#'   `map = TRUE` and  [tidybayes::stat_pointinterval()] when `map = FALSE`.
#'
#' @return A `ggplot` object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`plot_data`}{The data structure used to produce the plot.}
#'   }
#'
#' @seealso [fit_model()], [plot_intercepts()], [plot_coefficients()],
#'   [plot_surveys()], [plot_correlations()], [plot_partitions()]
#'
#' @export
plot_sites <- function(
  fit,
  species = NULL,
  sites = NULL,
  map = TRUE,
  intercepts = TRUE,
  back_transform = TRUE,
  include_predictors = TRUE,
  restricted = TRUE,
  ndraws = NULL,
  seed = NULL,
  ...
) {
  if (!inherits(fit, "CmdStanFit")) {
    cli::cli_abort(
      "{.arg fit} must be a {.cls CmdStanFit} object from {.fun fit_model}."
    )
  }
  rlang::check_installed(c("tidybayes", "ggplot2"))
  stan_data <- attr(fit, "stan_data")
  occARU_data <- attr(fit, "occARU_data")
  site_lvl <- attr(occARU_data, "sites")
  species_lvl <- attr(occARU_data, "species")
  MS <- length(species_lvl) > 1
  survey_length <- attr(occARU_data, "survey_length")
  P <- occARU_data$P[2]
  P_cat <- occARU_data$P_cat[2]
  P_ord <- occARU_data$P_ord[2]
  P_sum <- sum(c(P, P_cat, P_ord))
  SP <- stan_data$spatial > 0
  has_latent <- stan_data$latent == 1L
  transform <- back_transform && intercepts
  has_coords <- any(occARU_data$XY != 0)
  use_map <- map && has_coords
  if (!SP && !P_sum) {
    cli::cli_abort(
      "Model was fit without site predictors or random effects."
    )
  }
  if (map && !has_coords) {
    cli::cli_warn(
      "No site coordinates found in {.arg data}. Falling back to \\
      point-intervals."
    )
  }
  if (!is.null(ndraws) && is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1)
  }

  # get species and site indices
  species_idx <- species_indices(species, species_lvl)
  if (!is.null(sites)) {
    unknown <- setdiff(sites, site_lvl)
    if (length(unknown)) {
      cli::cli_abort(
        "The following sites were not found in the model: {.val {unknown}}."
      )
    }
    site_idx <- which(site_lvl %in% sites)
  } else {
    site_idx <- 1:length(site_lvl)
  }

  # determine projection
  if (!restricted) {
    if (P_sum && include_predictors) {
      cli::cli_abort(
        "{.arg restricted} = {restricted} is only applicable when predictors \\
        are in the model but excluded from the plot."
      )
    } else if (!P_sum) {
      cli::cli_warn(
        "{.arg restricted} is only applicable when predictors are included."
      )
    }
  }
  res <- !(P_sum && !restricted && !include_predictors)

  # initialise log_mu
  draws <- tidyr::expand_grid(s = species_idx, i = site_idx) |>
    dplyr::mutate(log_mu = 0)

  # increment predictor effects
  if (P_sum && include_predictors) {
    if (P) {
      if (MS) {
        betaX <- dplyr::left_join(
          tidybayes::spread_rvars(
            fit,
            mu_beta[s, p],
            ndraws = ndraws,
            seed = seed
          ) |>
            dplyr::filter(s %in% species_idx),
          predictor_matrix_to_tibble(stan_data$X2) |>
            dplyr::filter(i %in% site_idx),
          by = dplyr::join_by(p),
          relationship = "many-to-many"
        ) |>
          dplyr::summarise(
            betaX = posterior::rvar_sum(mu_beta * x),
            .by = c(s, i)
          )

        draws <- dplyr::left_join(draws, betaX, by = dplyr::join_by(s, i)) |>
          dplyr::mutate(log_mu = log_mu + betaX) |>
          dplyr::relocate(betaX, .before = log_mu)
      } else {
        betaX <- dplyr::left_join(
          tidybayes::spread_rvars(
            fit,
            mu_beta[p],
            ndraws = ndraws,
            seed = seed
          ),
          predictor_matrix_to_tibble(stan_data$X2) |>
            dplyr::filter(i %in% site_idx),
          by = dplyr::join_by(p),
          relationship = "many-to-many"
        ) |>
          dplyr::summarise(
            betaX = posterior::rvar_sum(mu_beta * x),
            .by = i
          )

        draws <- dplyr::left_join(draws, betaX, by = dplyr::join_by(i)) |>
          dplyr::mutate(log_mu = log_mu + betaX) |>
          dplyr::relocate(betaX, .before = log_mu)
      }
    }
    if (P_cat) {
      if (MS) {
        betaX_cat <- dplyr::left_join(
          tidybayes::spread_rvars(
            fit,
            mu_beta_cat[p, s, x],
            ndraws = ndraws,
            seed = seed
          ) |>
            dplyr::filter(s %in% species_idx, !is.nan(median(mu_beta_cat))),
          predictor_matrix_to_tibble(stan_data$X_cat2) |>
            dplyr::filter(i %in% site_idx),
          by = dplyr::join_by(p, x),
          relationship = "many-to-many"
        ) |>
          dplyr::summarise(
            betaX_cat = posterior::rvar_sum(mu_beta_cat),
            .by = c(s, i)
          )

        draws <- dplyr::left_join(
          draws,
          betaX_cat,
          by = dplyr::join_by(s, i)
        ) |>
          dplyr::mutate(log_mu = log_mu + betaX_cat) |>
          dplyr::relocate(betaX_cat, .before = log_mu)
      } else {
        betaX_cat <- dplyr::left_join(
          tidybayes::spread_rvars(
            fit,
            mu_beta_cat[p, x],
            ndraws = ndraws,
            seed = seed
          ) |>
            dplyr::filter(!is.na(median(mu_beta_cat))),
          predictor_matrix_to_tibble(stan_data$X_cat2) |>
            dplyr::filter(i %in% site_idx),
          by = dplyr::join_by(p, x),
          relationship = "many-to-many"
        ) |>
          dplyr::summarise(
            betaX_cat = posterior::rvar_sum(mu_beta_cat),
            .by = i
          )

        draws <- dplyr::left_join(draws, betaX_cat, by = dplyr::join_by(i)) |>
          dplyr::mutate(log_mu = log_mu + betaX_cat) |>
          dplyr::relocate(betaX_cat, .before = log_mu)
      }
    }
    if (P_ord) {
      if (MS) {
        betaX_ord <- dplyr::left_join(
          tidybayes::spread_rvars(
            fit,
            mu_beta_ord[s, p],
            mu_beta_ord_cs[p, x],
            ndraws = ndraws,
            seed = seed
          ) |>
            dplyr::filter(s %in% species_idx, !is.na(median(mu_beta_ord_cs))),
          predictor_matrix_to_tibble(stan_data$X_ord2) |>
            dplyr::filter(i %in% site_idx),
          by = dplyr::join_by(p, x),
          relationship = "many-to-many"
        ) |>
          dplyr::summarise(
            betaX_ord = posterior::rvar_sum(mu_beta_ord * mu_beta_ord_cs),
            .by = c(s, i)
          )

        draws <- dplyr::left_join(
          draws,
          betaX_ord,
          by = dplyr::join_by(s, i)
        ) |>
          dplyr::mutate(log_mu = log_mu + betaX_ord) |>
          dplyr::relocate(betaX_ord, .before = log_mu)
      } else {
        betaX_ord <- dplyr::left_join(
          tidybayes::spread_rvars(
            fit,
            mu_beta_ord[p],
            mu_beta_ord_cs[p, x],
            ndraws = ndraws,
            seed = seed
          ) |>
            dplyr::filter(!is.na(median(mu_beta_ord_cs))),
          predictor_matrix_to_tibble(stan_data$X_ord2) |>
            dplyr::filter(i %in% site_idx),
          by = dplyr::join_by(p, x),
          relationship = "many-to-many"
        ) |>
          dplyr::summarise(
            betaX_ord = posterior::rvar_sum(mu_beta_ord * mu_beta_ord_cs),
            .by = i
          )

        draws <- dplyr::left_join(draws, betaX_ord, by = dplyr::join_by(i)) |>
          dplyr::mutate(
            log_mu = dplyr::if_else(
              is.na(median(betaX_ord)),
              log_mu,
              log_mu + betaX_ord
            )
          ) |>
          dplyr::relocate(betaX_ord, .before = log_mu)
      }
    }
  }

  # increment random site effects
  if (SP) {
    param <- rlang::sym(paste0("iota", if (!res) "_unc"))
    if (MS) {
      draws <- dplyr::left_join(
        draws,
        tidybayes::spread_rvars(
          fit,
          (!!param)[s, i],
          ndraws = ndraws,
          seed = seed
        ),
        by = dplyr::join_by(s, i)
      )
    } else {
      draws <- dplyr::left_join(
        draws,
        tidybayes::spread_rvars(
          fit,
          (!!param)[i],
          ndraws = ndraws,
          seed = seed
        ),
        by = dplyr::join_by(i)
      )
    }

    draws <- draws |>
      dplyr::mutate(log_mu = log_mu + !!param) |>
      dplyr::relocate(!!param, .before = log_mu)
  }

  # join intercepts
  if (intercepts) {
    if (MS) {
      draws <- dplyr::left_join(
        draws,
        tidybayes::spread_rvars(
          fit,
          alpha[s, d],
          ndraws = ndraws,
          seed = seed
        ) |>
          dplyr::filter(s %in% species_idx, d == 2) |>
          dplyr::select(-d),
        by = dplyr::join_by(s)
      )
    } else {
      draws <- dplyr::cross_join(
        draws,
        tidybayes::spread_rvars(fit, alpha[d], ndraws = ndraws, seed = seed) |>
          dplyr::filter(d == 2) |>
          dplyr::select(-d)
      )
    }

    draws <- draws |>
      dplyr::mutate(log_mu = alpha + log_mu) |>
      dplyr::relocate(alpha, .after = i)
  } else if (back_transform) {
    cli::cli_warn(
      "{.arg back_transform = TRUE} has no effect when{.arg intercepts = FALSE}. \\
      Plotting log detection rates."
    )
  }

  # join with occupancy states if recovered or occupancy log odds
  if (has_latent) {
    if (MS) {
      draws <- dplyr::left_join(
        draws,
        tidybayes::spread_rvars(fit, z[i, s], ndraws = ndraws, seed = seed) |>
          dplyr::filter(s %in% species_idx, i %in% site_idx),
        by = dplyr::join_by(i, s)
      )
    } else {
      draws <- dplyr::left_join(
        draws,
        tidybayes::spread_rvars(fit, z[i], ndraws = ndraws, seed = seed) |>
          dplyr::filter(i %in% site_idx),
        by = dplyr::join_by(i)
      )
    }
    draws <- draws |>
      dplyr::mutate(
        occupied = factor(median(z) == 1, levels = c(TRUE, FALSE)),
        pred = log(z * exp(log_mu))
      ) |>
      dplyr::relocate(z, occupied, .after = i)
  } else {
    inv_logit <- function(x) 1 / (1 + exp(-x))
    if (MS) {
      draws <- dplyr::left_join(
        draws,
        tidybayes::spread_rvars(
          fit,
          logit_psi[s, i],
          ndraws = ndraws,
          seed = seed
        ) |>
          dplyr::filter(s %in% species_idx, i %in% site_idx),
        by = dplyr::join_by(i, s)
      )
    } else {
      draws <- dplyr::left_join(
        draws,
        tidybayes::spread_rvars(
          fit,
          logit_psi[i],
          ndraws = ndraws,
          seed = seed
        ) |>
          dplyr::filter(i %in% site_idx),
        by = dplyr::join_by(i)
      )
    }

    draws <- draws |>
      dplyr::relocate(logit_psi, .after = i) |>
      dplyr::mutate(pred = log(inv_logit(logit_psi) * exp(log_mu)))
  }

  # add labels
  draws <- draws |>
    dplyr::mutate(
      s = factor(species_lvl[s], levels = species_lvl[species_idx]),
      i = factor(site_lvl[i], levels = site_lvl[site_idx])
    ) |>
    dplyr::rename(species = s, site = i)

  # transform
  draws <- dplyr::mutate(draws, pred = if (transform) exp(log_mu) else log_mu)
  mu_lab <- paste0(
    "Detection rate per ",
    survey_length,
    "-day survey",
    if (!transform) " (log)"
  )

  # plot
  if (use_map) {
    draws <- dplyr::left_join(
      draws,
      as.data.frame(occARU_data$XY[site_idx, ]) |>
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
        tidybayes::stat_pointinterval(
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
        tidybayes::stat_pointinterval(
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
#' Plots species-level temporal detection rates combining survey predictors (if
#' included via [make_data()]) and survey random effects (`kappa`).
#'
#' @param fit A fitted model object from [fit_model()].
#' @param species Character vector of species to plot. If `NULL` (default), all
#'   species are plotted. Must be one of `attr(occARU_data, "species")`.
#' @param surveys Character vector of survey dates to plot. If `NULL` (default),
#'   all surveys are plotted. Must be one of `attr(occARU_data, "surveys")`.
#' @param intercepts Logical. If `TRUE` (default), species-level baseline log
#'   detection rates are added to the survey effects. If `FALSE`, only the
#'   temporal deviations are plotted on the log scale.
#' @param back_transform Logical. Only used when `intercepts = TRUE`. If `TRUE`
#'   (default), log detection rates are back-transformed to the natural scale
#'   via `exp()`. If `FALSE`, values are left on the log scale.
#' @param include_predictors Logical. If `TRUE` (default), includes predictors
#'   in the survey effects, if included. If `FALSE`, only plots the random
#'   effects.
#' @param restricted Logical. If `TRUE` (default), when `include_predictors` is
#'   `FALSE`, plots random survey effects with orthogonal projection, i.e.,
#'   \eqn{\boldsymbol{\kappa}(\boldsymbol{I} - \boldsymbol{P_X})}, where
#'   \eqn{\boldsymbol{I} - \boldsymbol{P_X}} is the orthogonal complement of the
#'   column space of the survey design matrix. If `FALSE`, plots random effects
#'   without orthogonal projection, i.e., \eqn{\boldsymbol{\kappa}} only. Has
#'   no effect when `include_predictors` is `TRUE` as the linear predictor
#'   is unaffected by orthogonal projection. The orthogonal complement of the
#'   design matrix is computed using the site-averaged survey predictor values.
#' @param ndraws  Number of draws to use for plotting, passed to
#'   [tidybayes::spread_rvars()]. Default: `NULL` (uses all draws).
#' @param seed Positive numeric. Seed to use when subsampling draws when
#'   `ndraws` is not `NULL`. Default: random integer.
#' @param palette Colour palette to be passed to [ggplot2::scale_fill_brewer()].
#'   Default: `"YlGn"`.
#' @param ... Additional arguments passed to [tidybayes::stat_lineribbon()]
#'   such as `.width` and `point_interval`.
#'
#' @return A `ggplot` object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`plot_data`}{The data structure used to produce the plot.}
#'   }
#'
#' @seealso [fit_model()], [plot_intercepts()], [plot_coefficients()],
#'   [plot_sites()], [plot_correlations()], [plot_partitions()]
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
      "{.arg fit} must be a {.cls CmdStanFit} object from {.fun fit_model}."
    )
  }
  rlang::check_installed(c("tidybayes", "ggplot2"))
  stan_data <- attr(fit, "stan_data")
  occARU_data <- attr(fit, "occARU_data")
  survey_lvl <- attr(occARU_data, "surveys")
  species_lvl <- attr(occARU_data, "species")
  MS <- length(species_lvl) > 1
  survey_length <- attr(occARU_data, "survey_length")
  transform <- back_transform && intercepts
  P <- occARU_data$P[3]
  P_cat <- occARU_data$P_cat[3]
  P_ord <- occARU_data$P_ord[3]
  P_sum <- sum(c(P, P_cat, P_ord))
  TE <- stan_data$temporal > 0
  if (!TE && !P_sum) {
    cli::cli_abort(
      "Model was fit without survey predictors or random effects."
    )
  }
  transform <- back_transform && intercepts
  if (!is.null(ndraws) && is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1)
  }

  # get species and survey indices
  species_idx <- species_indices(species, species_lvl)
  if (!is.null(surveys)) {
    unknown <- setdiff(surveys, survey_lvl)
    if (length(unknown)) {
      cli::cli_abort(
        "The following surveys were not found in the model: {.val {unknown}}."
      )
    }
    survey_idx <- which(survey_lvl %in% surveys)
  } else {
    survey_idx <- 1:length(survey_lvl)
  }

  # determine projection
  if (!restricted) {
    if (P_sum && include_predictors) {
      cli::cli_abort(
        "{.arg restricted} = {restricted} is only applicable when predictors \\
      are in the model but excluded from the plot."
      )
    } else if (!P_sum) {
      cli::cli_warn(
        "{.arg restricted} is only applicable when predictors are included."
      )
    }
  }
  res <- !(P_sum && !restricted && !include_predictors)

  # initialise log_mu
  draws <- tidyr::expand_grid(s = species_idx, j = survey_idx) |>
    dplyr::mutate(log_mu = 0)

  # increment predictor effects
  if (P_sum && include_predictors) {
    if (P) {
      if (MS) {
        gammaX <- dplyr::left_join(
          tidybayes::spread_rvars(
            fit,
            gamma[s, p],
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
            .by = c(s, j)
          )

        draws <- dplyr::left_join(draws, gammaX, by = dplyr::join_by(s, j)) |>
          dplyr::mutate(log_mu = log_mu + gammaX) |>
          dplyr::relocate(gammaX, .before = log_mu)
      } else {
        gammaX <- dplyr::left_join(
          tidybayes::spread_rvars(
            fit,
            gamma[p],
            ndraws = ndraws,
            seed = seed
          ),
          predictor_matrix_to_tibble(apply(stan_data$X3, 2:3, mean)) |>
            dplyr::rename(j = i) |>
            dplyr::filter(j %in% survey_idx),
          by = dplyr::join_by(p),
          relationship = "many-to-many"
        ) |>
          dplyr::summarise(
            gammaX = posterior::rvar_sum(gamma * x),
            .by = j
          )

        draws <- dplyr::left_join(draws, gammaX, by = dplyr::join_by(j)) |>
          dplyr::mutate(log_mu = log_mu + gammaX) |>
          dplyr::relocate(gammaX, .before = log_mu)
      }
    }
    if (P_cat) {
      if (MS) {
        gammaX_cat <- dplyr::left_join(
          tidybayes::spread_rvars(fit, gamma_cat[p, s, x]) |>
            dplyr::filter(s %in% species_idx),
          predictor_matrix_to_tibble(apply(stan_data$X_cat3, 2:3, int_mode)) |>
            dplyr::rename(j = i) |>
            dplyr::filter(j %in% survey_idx),
          by = dplyr::join_by(p, x),
          relationship = "many-to-many"
        ) |>
          dplyr::summarise(
            gammaX_cat = posterior::rvar_sum(gamma_cat),
            .by = c(s, j)
          )

        draws <- dplyr::left_join(
          draws,
          gammaX_cat,
          by = dplyr::join_by(s, j)
        ) |>
          dplyr::mutate(log_mu = log_mu + gammaX_cat) |>
          dplyr::relocate(gammaX_cat, .before = log_mu)
      } else {
        gammaX_cat <- dplyr::left_join(
          tidybayes::spread_rvars(fit, gamma_cat[p, x]),
          predictor_matrix_to_tibble(apply(stan_data$X_cat3, 2:3, int_mode)) |>
            dplyr::rename(j = i) |>
            dplyr::filter(j %in% survey_idx),
          by = dplyr::join_by(p, x),
          relationship = "many-to-many"
        ) |>
          dplyr::summarise(
            gammaX_cat = posterior::rvar_sum(gamma_cat),
            .by = j
          )

        draws <- dplyr::left_join(draws, gammaX_cat, by = dplyr::join_by(j)) |>
          dplyr::mutate(log_mu = log_mu + gammaX_cat) |>
          dplyr::relocate(gammaX_cat, .before = log_mu)
      }
    }
    if (P_ord) {
      if (MS) {
        gammaX_ord <- dplyr::left_join(
          tidybayes::spread_rvars(
            fit,
            gamma_ord[s, p],
            gamma_ord_cs[p, x],
            ndraws = ndraws,
            seed = seed
          ) |>
            dplyr::filter(s %in% species_idx),
          predictor_matrix_to_tibble(apply(stan_data$X_ord3, 2:3, int_mode)) |>
            dplyr::rename(j = i) |>
            dplyr::filter(j %in% survey_idx),
          by = dplyr::join_by(p, x),
          relationship = "many-to-many"
        ) |>
          dplyr::summarise(
            gammaX_ord = posterior::rvar_sum(gamma_ord * gamma_ord_cs),
            .by = c(s, j)
          )

        draws <- dplyr::left_join(
          draws,
          gammaX_ord,
          by = dplyr::join_by(s, j)
        ) |>
          dplyr::mutate(
            log_mu = dplyr::if_else(
              is.na(median(betaX_ord)),
              log_mu,
              log_mu + betaX_ord
            )
          ) |>
          dplyr::relocate(gammaX_ord, .before = log_mu)
      } else {
        gammaX_ord <- dplyr::left_join(
          tidybayes::spread_rvars(
            fit,
            gamma_ord[p],
            gamma_ord_cs[p, x],
            ndraws = ndraws,
            seed = seed
          ),
          predictor_matrix_to_tibble(apply(stan_data$X_ord3, 2:3, int_mode)) |>
            dplyr::rename(j = i) |>
            dplyr::filter(j %in% survey_idx),
          by = dplyr::join_by(p, x),
          relationship = "many-to-many"
        ) |>
          dplyr::summarise(
            gammaX_ord = posterior::rvar_sum(gamma_ord * gamma_ord_cs),
            .by = j
          )
        draws <- dplyr::left_join(draws, gammaX_ord, by = dplyr::join_by(j)) |>
          dplyr::mutate(
            log_mu = dplyr::if_else(
              is.na(median(betaX_ord)),
              log_mu,
              log_mu + betaX_ord
            )
          ) |>
          dplyr::relocate(gammaX_ord, .before = log_mu)
      }
    }
  }

  # increment random survey effects
  if (TE) {
    param <- rlang::sym(paste0("kappa", if (!res) "_unc"))
    if (MS) {
      draws <- dplyr::left_join(
        draws,
        tidybayes::spread_rvars(
          fit,
          (!!param)[s, j],
          ndraws = ndraws,
          seed = seed
        ),
        by = dplyr::join_by(s, j)
      ) |>
        dplyr::mutate(log_mu = log_mu + !!param) |>
        dplyr::relocate(!!param, .before = log_mu)
    } else {
      draws <- dplyr::left_join(
        draws,
        tidybayes::spread_rvars(
          fit,
          (!!param)[j],
          ndraws = ndraws,
          seed = seed
        ),
        by = dplyr::join_by(j)
      ) |>
        dplyr::mutate(log_mu = log_mu + !!param) |>
        dplyr::relocate(!!param, .before = log_mu)
    }
  }

  # join intercepts
  if (intercepts) {
    if (MS) {
      draws <- dplyr::left_join(
        draws,
        tidybayes::spread_rvars(
          fit,
          alpha[s, d],
          ndraws = ndraws,
          seed = seed
        ) |>
          dplyr::filter(s %in% species_idx, d == 2) |>
          dplyr::select(-d),
        by = dplyr::join_by(s)
      ) |>
        dplyr::mutate(log_mu = alpha + log_mu) |>
        dplyr::relocate(alpha, .after = j)
    } else {
      draws <- dplyr::cross_join(
        draws,
        tidybayes::spread_rvars(fit, alpha[d], ndraws = ndraws, seed = seed) |>
          dplyr::filter(d == 2) |>
          dplyr::select(-d)
      ) |>
        dplyr::mutate(log_mu = alpha + log_mu) |>
        dplyr::relocate(alpha, .after = j)
    }
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
    tidybayes::stat_lineribbon(...) +
    ggplot2::facet_wrap(~species, scales = if (transform) "free_y", ncol = 1) +
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
  P <- nrow(X)
  tidyr::as_tibble(t(X)) |>
    setNames(1:P) |>
    dplyr::mutate(i = dplyr::row_number()) |>
    tidyr::pivot_longer(-i, names_to = "p", values_to = "x") |>
    dplyr::mutate(p = as.integer(p), x = if (categorical) as.integer(x) else x)
}

#' Produce species indices based on supplied character vector of species
#'
#' @param species Character vector species names. Default: `NULL`
#' @param species_lvl Species levels of orginal data
#' @return A vector of integers corresponding to species to retain
#' @noRd
species_indices <- function(species = NULL, species_lvl) {
  if (!is.null(species)) {
    unknown <- setdiff(species, species_lvl)
    if (length(unknown)) {
      cli::cli_abort(
        "The following species were not found in the model: {.val {unknown}}."
      )
    }
    species_idx <- which(species_lvl %in% species)
  } else {
    species_idx <- 1:length(species_lvl)
  }
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


#' Plot red dashed verticle 0-line
#'
#' @param p Plot
#' @return Plot with ggplot2::geom_vline()
#' @noRd
my_vline <- function() {
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "red4")
}
