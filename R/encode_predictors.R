#' Encode a predictor dataframe for the occARU Stan model
#'
#' Splits a predictor dataframe into separate matrices/arrays for continuous,
#' categorical, and ordinal predictors, optionally scaling continuous
#' predictors, and extracts factor levels for use in the Stan model and for
#' post-processing.
#'
#' @param df A dataframe of predictors including the deployment ID column,
#'   and for survey mode a date column. If `NULL`, returns `NULL` immediately.
#' @param mode Character. Either `"site"` (returns `[I, P]` matrices) or
#'   `"survey"` (returns `[I, P, J]` arrays).
#' @param locationID <[`data-masking`][rlang::args_data_masking]> Name of the
#'   site column. Default: `locationID`.
#' @param season <[`data-masking`][rlang::args_data_masking]> Name of the season
#'   column. Default: `season`.
#' @param date `Date`. <[`data-masking`][rlang::args_data_masking]> Name of the
#'   date column. Default: `NULL`.
#' @param site_lvl `character`. Vector of site factor levels.
#' @param season_lvl `character`. Vector of season factor levels.
#' @param reference_dates A dataframe of dates defining the start of the first
#'   survey period per season. Default: `NULL`.
#' @param surveys A dataframe of survey dates and indexes per season. Default:
#'   `NULL`
#' @param summary_functions An optional named list mapping continuous survey
#'   predictor column names to summary functions, used when aggregating survey
#'   predictors over `survey_length`-length periods. Each value can be a
#'   function name as a string (e.g. `"sum"`) or a function object (e.g. `sum`).
#'   Numeric predictors not named in `summary_functions` are summarised with
#'   `mean`; categorical and ordinal
#' @param survey_length Positive integer. Number of days per survey period. Only
#'   used when `mode = "survey"`.
#' @param scale_predictors Logical. If `TRUE` (default), continuous predictors
#'   are scaled to zero mean and unit variance. When `mode = "survey"`, scaling
#'   parameters are derived from the site-averaged values per survey period
#'   rather than the raw site-by-survey values.
#'
#' @return A named list or `NULL` if `df` is `NULL`:
#'   \describe{
#'     \item{`X`}{Continuous predictor matrix `[I, P]` or array `[I, P, J]`.}
#'     \item{`X_cat`}{Categorical integer matrfx `[I, P_cat]` or array
#'       `[I, P_cat, J]`.}
#'     \item{`X_ord`}{Ordinal integer matrix `[I, P_ord]` or array
#'       `[I, P_ord, J]`.}
#'     \item{`C`}{Named list of factor levels per categorical predictor.}
#'     \item{`O`}{Named list of factor levels per ordinal predictor.}
#'     \item{`scale_params`}{Named list of `list(mean, sd)` per continuous
#'       predictor, or `NULL`.}
#'   }
#' @noRd
encode_predictors <- function(
  df,
  mode = c("site", "survey"),
  locationID = locationID,
  .season = .season,
  date = date,
  site_lvl = site_lvl,
  season_lvl = season_lvl,
  surveys = NULL,
  reference_dates = NULL,
  summary_functions = NULL,
  survey_length = 1,
  scale_predictors = TRUE
) {
  mode <- match.arg(mode)

  # produce empty matrix and array
  I <- length(site_lvl)
  K <- length(season_lvl)
  if (!is.null(surveys)) {
    J_max <- max(as.integer(surveys$.survey_idx))
  }

  empty_array <- if (mode == "site") {
    array(0L, c(K, I, 0L), dimnames = list(season_lvl, site_lvl, NULL))
  } else {
    array(
      0L,
      c(K, I, J_max, 0),
      dimnames = list(
        season_lvl,
        site_lvl,
        if (K == 1) as.character(surveys$.survey),
        NULL
      )
    )
  }

  # return empties
  if (is.null(df)) {
    list(
      X = empty_array,
      X_cat = empty_array,
      X_ord = empty_array
    )
  } else {
    # classify
    df_pred <- df |>
      dplyr::select(-c({{ locationID }}, {{ .season }}))
    if (mode == "survey") {
      df_pred <- df_pred |>
        dplyr::select(-{{ date }})
    }
    pred_cols <- colnames(df_pred)
    types <- classify_predictors(df_pred)
    num_cols <- types$numeric
    cat_cols <- types$categorical
    ord_cols <- types$ordinal
    fct_cols <- c(cat_cols, ord_cols)

    P <- length(num_cols)
    P_cat <- length(cat_cols)
    C <- lapply(df[cat_cols], levels)
    P_ord <- length(ord_cols)
    O <- lapply(df[ord_cols], levels)

    # check factors and convert to integers
    single_cat <- names(C)[sapply(C, length) < 2]
    if (length(single_cat)) {
      cli::cli_abort(
        "The following categorical predictor{?s} {?has/have} fewer than 2 \\
        levels: {.val {single_cat}}. Remove {?it/them} before passing to \\
        {.fun make_data}."
      )
    }

    double_ord <- names(O)[sapply(O, length) < 3]
    if (length(double_ord)) {
      cli::cli_abort(
        "The following ordinal predictor{?s} {?has/have} fewer than 3 levels: \\
        {.val {double_ord}}. Should be supplied as categorical instead?"
      )
    }

    check_empty_levels(df, !!!rlang::syms(fct_cols))

    df <- df |>
      dplyr::mutate(
        dplyr::across(dplyr::all_of(fct_cols), as.integer)
      )

    # survey aggregation
    if (mode == "survey") {
      # check summary functions
      if (!is.null(summary_functions)) {
        unknown <- setdiff(names(summary_functions), pred_cols)
        if (length(unknown)) {
          cli::cli_abort(
            "The following predictors in {.arg summary_functions} were not \\
            found in {.arg survey_predictors}: {.val {unknown}}."
          )
        }

        invalid_fns <- names(purrr::discard(
          summary_functions,
          ~ is.function(.) || exists(as.character(.), mode = "function")
        ))
        if (length(invalid_fns)) {
          cli::cli_abort(
            "The following entries in {.arg summary_functions} have invalid \\
            functions: {.val {invalid_fns}}."
          )
        }

        overrides <- intersect(names(summary_functions), fct_cols)
        if (length(overrides)) {
          cli::cli_warn(
            "Ignoring {.arg summary_functions} entries for non-numeric \\
            predictor{?s} {.val {overrides}}. Factor and ordinal predictors \\
            always use modal aggregation."
          )
        }
      }

      # produce final summary functions
      summary_fns <- purrr::set_names(pred_cols) |>
        purrr::map(
          ~ {
            if (
              . %in%
                num_cols &&
                !is.null(summary_functions) &&
                . %in% names(summary_functions)
            ) {
              match.fun(summary_functions[[.]])
            } else if (. %in% num_cols) {
              mean
            } else {
              int_mode
            }
          }
        )

      # aggregate and summarise
      df <- dplyr::left_join(df, reference_dates, by = ".season") |>
        dplyr::mutate(
          .survey = aggregate_by_days(
            {{ date }},
            reference,
            survey_length
          ),
          .by = ".season"
        ) |>
        # dplyr::right_join(surveys, by = c(".season", ".survey")) |>
        dplyr::summarise(
          dplyr::across(
            dplyr::all_of(pred_cols),
            ~ summary_fns[[dplyr::cur_column()]](.)
          ),
          .by = c({{ locationID }}, .survey)
        )

      # join survey indices
      df <- dplyr::left_join(df, surveys, by = ".survey")
    }

    # scale continuous predictors
    scale_params <- NULL
    if (scale_predictors && P) {
      df_num <- df |>
        dplyr::select(-dplyr::all_of(fct_cols))
      if (mode == "site") {
        # average over seasons for scaling
        site_means <- df_num |>
          dplyr::summarise(
            dplyr::across(dplyr::all_of(num_cols), mean),
            .by = {{ locationID }}
          )
        scale_params <- scaling_parameters(site_means, num_cols)
        df <- df |>
          dplyr::mutate(
            dplyr::across(
              dplyr::all_of(num_cols),
              ~ (. - scale_params[[dplyr::cur_column()]][1]) /
                scale_params[[dplyr::cur_column()]][2]
            )
          )
      } else {
        # average over sites for scaling
        df_num <- df |>
          dplyr::select(-dplyr::all_of(fct_cols))
        survey_means <- df_num |>
          dplyr::summarise(
            dplyr::across(dplyr::all_of(num_cols), mean),
            .by = .survey
          )
        scale_params <- scaling_parameters(survey_means, num_cols)
        df <- df |>
          dplyr::mutate(
            dplyr::across(
              dplyr::all_of(num_cols),
              ~ (. - scale_params[[dplyr::cur_column()]][1]) /
                scale_params[[dplyr::cur_column()]][2]
            )
          )
      }
    }

    # build return arrays
    to_array <- if (mode == "site") {
      function(df, cols) {
        P <- length(cols)
        if (P) {
          df |>
            dplyr::select(
              {{ locationID }},
              {{ .season }},
              dplyr::all_of(cols)
            ) |>
            tidyr::complete(
              {{ .season }} := factor(season_lvl, season_lvl),
              {{ locationID }} := factor(site_lvl, site_lvl),
              fill = as.list(purrr::set_names(
                ifelse(cols %in% num_cols, 0L, 1L),
                cols
              ))
            ) |>
            tidyr::pivot_longer(
              dplyr::all_of(cols),
              names_to = "p",
              values_to = "x"
            ) |>
            dplyr::mutate(p = factor(p, levels = cols)) |>
            dplyr::arrange(p, {{ locationID }}, {{ .season }}) |>
            dplyr::pull(x) |>
            array(c(K, I, P), dimnames = list(season_lvl, site_lvl, cols))
        } else {
          empty_array
        }
      }
    } else {
      function(df, cols) {
        P <- length(cols)
        if (P) {
          df |>
            dplyr::select(
              {{ locationID }},
              {{ .season }},
              .survey_idx,
              dplyr::all_of(cols)
            ) |>
            tidyr::complete(
              {{ locationID }},
              {{ .season }},
              .survey_idx,
              fill = as.list(purrr::set_names(
                ifelse(cols %in% num_cols, 0L, 1L),
                cols
              ))
            ) |>
            tidyr::pivot_longer(
              dplyr::all_of(cols),
              names_to = "p",
              values_to = "x"
            ) |>
            dplyr::mutate(p = factor(p, levels = cols)) |>
            dplyr::arrange(p, .survey_idx, {{ locationID }}, {{ .season }}) |>
            dplyr::pull(x) |>
            array(
              c(K, I, J_max, P),
              dimnames = list(
                season_lvl,
                site_lvl,
                if (K == 1) as.character(surveys$.survey),
                cols
              )
            )
        } else {
          empty_array
        }
      }
    }

    X <- to_array(df, num_cols)
    X_cat <- to_array(df, cat_cols)
    X_ord <- to_array(df, ord_cols)

    list(
      X = X,
      X_cat = X_cat,
      X_ord = X_ord - 1, # 0 is reference
      C = C,
      O = O,
      scale_params = scale_params
    )
  }
}

#' Extract factor levels from an encoder output
#'
#' @param enc Output of `encode_predictors()` or `NULL`.
#' @return A named list with elements `cat` and `ord`, each a named list of
#'   character vectors.
#' @noRd
enc_levels <- function(enc) {
  if (is.null(enc)) {
    NULL
  }
  list(
    cat = if (!rlang::is_empty(enc$C)) enc$C,
    ord = if (!rlang::is_empty(enc$O)) enc$O
  )
}
