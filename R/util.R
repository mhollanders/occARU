#' Aggregate dates to survey periods
#'
#' @param dates A vector of dates.
#' @param reference Reference date for aggregation.
#' @param survey_length Integer number of days per survey period.
#' @return A vector of dates floored to the nearest survey period.
#' @noRd
aggregate_by_days <- function(dates, reference, survey_length = 1) {
  reference +
    as.integer(as.Date(dates) - reference) %/% survey_length * survey_length
}

#' Align a factor column to reference levels, warning or erroring on missing
#' values
#'
#' @param df A dataframe.
#' @param col <[`data-masking`][rlang::args_data_masking]> Column name of
#'   factor.
#' @param ref_levels Reference levels.
#' @param strict If `TRUE`, missing values trigger an error. If `FALSE`,
#'   missing values trigger a warning and affected rows are removed.
#' @return The dataframe with `col` realigned and missing values removed.
#' @noRd
align_factor <- function(df, col, ref_levels, strict = FALSE) {
  df_name <- deparse(substitute(df))
  check_cols_exist(df, {{ col }})
  df <- dplyr::mutate(df, {{ col }} := factor({{ col }}, levels = ref_levels))
  n_missing <- dplyr::filter(df, is.na({{ col }})) |> nrow()
  if (n_missing) {
    if (strict) {
      cli::cli_abort(
        "{.arg {df_name}} contains {n_missing} record{?s} with {.val {col}} \\
        values not found."
      )
    } else {
      cli::cli_warn(
        "{.arg {df_name}} contains {n_missing} record{?s} with {.val {col}} \\
        values not found. {?This record/These records} will be ignored."
      )
      df <- tidyr::drop_na(df, col)
    }
  }
  df
}

#' Assign detections to temporal clusters
#'
#' Assigns each detection to a cluster based on a minimum time gap between
#' consecutive detections. A new cluster begins whenever the time since the
#' previous detection exceeds `gap` minutes.
#'
#' @param times A vector of POSIXt timestamps, assumed to be sorted.
#' @param gap Numeric. Minimum gap in minutes to start a new cluster.
#' @return An integer vector of cluster assignments.
#' @noRd
assign_clusters <- function(times, gap) {
  gaps <- c(
    0,
    as.double(difftime(times[-1], times[-length(times)], units = "mins"))
  )
  cumsum(gaps >= gap) + 1L
}

#' Check that columns inherit a given class
#'
#' @param df A dataframe.
#' @param class Expected class as a character string.
#' @param ... <[`data-masking`][rlang::args_data_masking]> Column names.
#' @noRd
check_cols_class <- function(df, class, ...) {
  df_name <- deparse(substitute(df))
  col_names <- purrr::map_chr(rlang::enquos(...), rlang::as_name)
  wrong_cols <- purrr::discard(col_names, ~ inherits(df[[.]], class))
  if (length(wrong_cols)) {
    cli::cli_abort(
      "The following columns in {.arg {df_name}} must be {.cls {class}}: \\
       {.val {wrong_cols}}."
    )
  }
}

#' Check if a dataframe has duplicate values in columns
#'
#' @param df A dataframe.
#' @param ... <[`data-masking`][rlang::args_data_masking]> Column names.
#' @noRd
check_cols_duplicates <- function(df, ...) {
  df_name <- deparse(substitute(df))
  col_names <- purrr::map_chr(rlang::enquos(...), rlang::as_name)
  purrr::walk(
    col_names,
    ~ {
      dupes <- df |>
        dplyr::count(.data[[.]]) |>
        dplyr::filter(n > 1) |>
        dplyr::pull(.data[[.]])
      if (length(dupes)) {
        cli::cli_abort(
          "{.arg {df_name}} contains duplicate {.val {.}} values: \\
          {.val {dupes}}."
        )
      }
    }
  )
}

#' Check that columns exist in a dataframe
#'
#' @param df A dataframe.
#' @param ... <[`data-masking`][rlang::args_data_masking]> Column names.
#' @noRd
check_cols_exist <- function(df, ...) {
  df_name <- deparse(substitute(df))
  col_names <- purrr::map_chr(rlang::enquos(...), rlang::as_name)
  missing_cols <- purrr::discard(col_names, ~ rlang::has_name(df, .))
  if (length(missing_cols)) {
    cli::cli_abort(
      "The following columns were not found in {.arg {df_name}}: \\
       {.val {missing_cols}}."
    )
  }
}

#' Check that all factor levels are present in the data
#'
#' Aborts if any factor column contains levels with no corresponding rows.
#'
#' @param df A dataframe.
#' @param ... <[`data-masking`][rlang::args_data_masking]> Factor column names.
#' @param strict If `TRUE`, missing values trigger an error. If `FALSE`,
#'   missing values trigger a warning.
#' @noRd
check_empty_levels <- function(df, ..., strict = TRUE) {
  df_name <- deparse(substitute(df))
  col_names <- purrr::map_chr(rlang::enquos(...), rlang::as_name)
  purrr::walk(col_names, \(col_name) {
    col_vec <- df[[col_name]]
    empty <- purrr::discard(levels(col_vec), ~ . %in% as.character(col_vec))
    if (length(empty)) {
      if (strict) {
        cli::cli_abort(
          "{.arg {df_name}} has missing {.arg {col_name}} factor level{?s}: \\
          {.val {empty}}."
        )
      } else {
        cli::cli_warn(
          "{.arg {df_name}} has missing {.arg {col_name}} factor level{?s}: \\
          {.val {empty}}. This is not necessarily an error."
        )
      }
    }
  })
}

#' Check that all predictor columns are numeric, factor, or ordered factor
#'
#' @param df A dataframe of predictors (ID column already removed).
#' @param ... <[`data-masking`][rlang::args_data_masking]> Column names to
#'  ignore.
#' @noRd
check_mixed_predictors <- function(df, ...) {
  df_name <- deparse(substitute(df))
  bad <- df |>
    dplyr::select(-c(...)) |>
    dplyr::select(dplyr::where(\(x) !is.numeric(x) && !is.factor(x))) |>
    colnames()
  if (length(bad)) {
    cli::cli_abort(
      "{.arg {df_name}} contains columns that are not {.cls numeric}, \\
       {.cls factor}, or {.cls ordered}: {.val {bad}}. \\
       Convert them before passing to {.fun make_data}."
    )
  }
}

#' Check that survey predictors cover the full deployment period for each site
#'
#' @param survey_predictors A dataframe of survey-level predictors. Must
#'    contain columns `deploymentID` and `date`.
#' @param deployments A dataframe of deployment information, one row per
#'   site. Must contain columns `deploymentID`, `deploymentStart`, and
#'   `deploymentEnd`.
#' @param deploymentID <[`data-masking`][rlang::args_data_masking]> Column
#'   name for sites (ARUs). Must be a factor with identical levels in
#'   `deployments` and `survey_predictors`. Default: `deploymentID`.
#' @param deploymentStart <[`data-masking`][rlang::args_data_masking]> Column
#'   name for deployment start dates in `deployments`. Must be a `Date`.
#'   Default: `deploymentStart`.
#' @param deploymentEnd <[`data-masking`][rlang::args_data_masking]> Column
#'   name for deployment end dates in `deployments`. Must be a `Date`.
#'   Default: `deploymentEnd`.
#' @param date <[`data-masking`][rlang::args_data_masking]> Column name for
#'   dates in `survey_predictors`. Default: `date`.
#' @noRd
check_survey_predictors_coverage <- function(
  survey_predictors,
  deployments,
  deploymentID = deploymentID,
  deploymentStart = deploymentStart,
  deploymentEnd = deploymentEnd,
  date = date,
  verbose = TRUE
) {
  first_last <- survey_predictors |>
    dplyr::summarise(
      first = min({{ date }}),
      last = max({{ date }}),
      .by = {{ deploymentID }}
    ) |>
    dplyr::left_join(
      deployments |>
        dplyr::select(
          {{ deploymentID }},
          {{ deploymentStart }},
          {{ deploymentEnd }}
        ),
      by = dplyr::join_by({{ deploymentID }})
    )

  incomplete <- first_last |>
    dplyr::filter(
      first > {{ deploymentStart }} |
        last < {{ deploymentEnd }}
    ) |>
    dplyr::pull({{ deploymentID }})

  if (length(incomplete)) {
    cli::cli_abort(
      "{.arg survey_predictors} does not cover the full deployment period \\
       for the following site{?s}: {.val {incomplete}}."
    )
  } else if (verbose) {
    cli::cli_alert_success(
      "{.arg survey_predictors} covers the full deployment period."
    )
  }
}

#' Classify predictor columns by type
#'
#' @param df A dataframe of predictor columns only.
#' @return A named list with elements `numeric`, `categorical`, `ordinal`.
#' @noRd
classify_predictors <- function(df) {
  list(
    numeric = names(purrr::keep(df, is.numeric)),
    categorical = names(purrr::keep(df, \(x) is.factor(x) && !is.ordered(x))),
    ordinal = names(purrr::keep(df, is.ordered))
  )
}

#' Project site coordinates from WGS84 to UTM
#'
#' Converts longitude/latitude columns in a dataframe to UTM coordinates
#' (km), with the zone auto-detected from the mean longitude. Returns both
#' the projected coordinate matrix and the CRS string used.
#'
#' @param deployments A dataframe containing longitude and latitude columns.
#' @param deploymentID <[`data-masking`][rlang::args_data_masking]> Column
#'   name for sites (ARUs). Must be a factor with identical levels in
#'   `deployments` and `observations`. Default: `deploymentID`.
#' @param latitude <[`data-masking`][rlang::args_data_masking]> Column name
#'   for WGS84 latitude in `deployments`. Default: `latitude`.
#' @param longitude <[`data-masking`][rlang::args_data_masking]> Column name
#'   for WGS84 longitude in `deployments`. Default: `longitude`.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{`XY`}{Numeric matrix of UTM coordinates in km, with columns
#'       \code{X} and \code{Y}.}
#'     \item{`utm_crs`}{Character string of the PROJ CRS used for the
#'       transformation.}
#'   }
#' @noRd
coords_to_utm <- function(deployments, deploymentID, latitude, longitude) {
  lats <- dplyr::pull(deployments, {{ latitude }})
  lons <- dplyr::pull(deployments, {{ longitude }})

  mean_lon <- mean(lons)
  mean_lat <- mean(lats)
  zone <- floor((mean_lon + 180) / 6) + 1
  hemisphere <- ifelse(mean_lat >= 0, "north", "south")
  utm_crs <- paste0(
    "+proj=utm +zone=",
    zone,
    " +",
    hemisphere,
    " +datum=WGS84 +units=km"
  )

  zones <- floor((lons + 180) / 6) + 1
  if (length(unique(zones)) > 1) {
    cli::cli_warn(
      "Sites span multiple UTM zones ({.val {unique(zones)}}). \\
       Projecting all sites to zone {zone}."
    )
  }

  XY <- deployments |>
    sf::st_as_sf(
      coords = c(
        rlang::as_name(rlang::enquo(longitude)),
        rlang::as_name(rlang::enquo(latitude))
      ),
      crs = 4326
    ) |>
    sf::st_transform(utm_crs) |>
    sf::st_coordinates()
  row.names(XY) <- dplyr::pull(deployments, {{ deploymentID }}) |> levels()

  list(XY = XY, utm_crs = utm_crs)
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
#' @param deploymentID <[`data-masking`][rlang::args_data_masking]> Name of the
#'   site column.
#' @param date <[`data-masking`][rlang::args_data_masking]> Name of the date
#'   column.
#' @param reference_date A `Date` defining the start of the first survey
#'   period.
#' @param scale_predictors Logical. If `TRUE` (default), continuous predictors
#'   are scaled to zero mean and unit variance. When `mode = "survey"`, scaling
#'   parameters are derived from the site-averaged values per survey period
#'   rather than the raw site-by-survey values.
#' @param summary_functions An optional named list mapping continuous survey
#'   predictor column names to summary functions, used when aggregating survey
#'   predictors over `survey_length`-length periods. Each value can be a
#'   function name as a string (e.g. `"sum"`) or a function object (e.g. `sum`).
#'   Numeric predictors not named in `summary_functions` are summarised with
#'   `mean`; categorical and ordinal
#' @param survey_length Positive integer. Number of days per survey period. Only
#'   used when `mode = "survey"`.
#'
#' @return A named list or `NULL` if `df` is `NULL`:
#'   \describe{
#'     \item{`X`}{Continuous predictor matrix `[I, P]` or array `[I, P, J]`.}
#'     \item{`X_cat`}{Categorical integer matrix `[I, P_cat]` or array
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
  deploymentID,
  date = NULL,
  site_lvl,
  surveys = NULL,
  reference_date,
  scale_predictors = TRUE,
  summary_functions = NULL,
  survey_length = 1
) {
  mode <- match.arg(mode)

  # produce empty matrix and array
  I <- length(site_lvl)
  if (mode == "site") {
    empty_matrix <- matrix(
      0L,
      nrow = 0,
      ncol = I,
      dimnames = list(NULL, site_lvl)
    )
  } else {
    J <- length(surveys)
    empty_array <- array(
      0L,
      dim = c(I, 0, J),
      dimnames = list(site_lvl, NULL, as.character(surveys))
    )
  }

  # return empties
  if (is.null(df)) {
    if (mode == "site") {
      list(
        X = empty_matrix,
        X_cat = empty_matrix,
        X_ord = empty_matrix
      )
    } else {
      list(
        X = empty_array,
        X_cat = empty_array,
        X_ord = empty_array
      )
    }
  } else {
    # classify
    df_pred <- df |>
      dplyr::select(-{{ deploymentID }})
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

      if (survey_length == 1L) {
        df <- dplyr::mutate(df, survey = {{ date }})
      } else {
        df <- df |>
          dplyr::mutate(
            survey = aggregate_by_days(
              {{ date }},
              reference_date,
              survey_length
            )
          ) |>
          dplyr::summarise(
            dplyr::across(
              dplyr::all_of(pred_cols),
              ~ summary_fns[[dplyr::cur_column()]](.)
            ),
            .by = c({{ deploymentID }}, survey)
          )
      }
      surveys <- unique(df$survey)
      J <- length(surveys)
    }

    # scale continuous predictors
    scale_params <- NULL
    if (scale_predictors && P) {
      if (mode == "site") {
        scale_params <- scaling_parameters(df, num_cols)
        df <- df |>
          dplyr::mutate(
            dplyr::across(
              all_of(num_cols),
              ~ (. - scale_params[[dplyr::cur_column()]][1]) /
                scale_params[[dplyr::cur_column()]][2]
            )
          )
      } else {
        # use survey values averaged across sites for scaling
        df_num <- df |>
          dplyr::select(-dplyr::all_of(fct_cols))
        survey_means <- df_num |>
          dplyr::summarise(
            dplyr::across(dplyr::all_of(num_cols), mean),
            .by = survey
          )
        scale_params <- scaling_parameters(survey_means, num_cols)
        df <- df |>
          dplyr::mutate(
            dplyr::across(
              all_of(num_cols),
              ~ (. - scale_params[[dplyr::cur_column()]][1]) /
                scale_params[[dplyr::cur_column()]][2]
            )
          )
      }
    }

    # build return matrices/arrays
    if (mode == "site") {
      to_matrix <- function(df, cols) {
        if (length(cols)) {
          df |>
            dplyr::select({{ deploymentID }}, dplyr::all_of(cols)) |>
            tibble::column_to_rownames(rlang::as_name(rlang::enquo(
              deploymentID
            ))) |>
            t()
        } else {
          empty_matrix
        }
      }

      X <- to_matrix(df, num_cols)
      X_cat <- to_matrix(df, cat_cols)
      X_ord <- to_matrix(df, ord_cols)
    } else {
      to_array <- function(df, cols) {
        P <- length(cols)
        if (P) {
          df |>
            dplyr::select({{ deploymentID }}, survey, dplyr::all_of(cols)) |>
            tidyr::complete(
              {{ deploymentID }},
              survey,
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
            dplyr::mutate(p = factor(p, levels = num_cols)) |>
            dplyr::arrange(survey, p, {{ deploymentID }}) |>
            dplyr::pull(x) |>
            array(
              c(I, P, J),
              dimnames = list(site_lvl, cols, as.character(surveys))
            )
        } else {
          empty_array
        }
      }

      X <- to_array(df, num_cols)
      X_cat <- to_array(df, cat_cols)
      X_ord <- to_array(df, ord_cols)
    }

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

#' Remove observations outside deployment window
#'
#' @param deployments A dataframe of deployment information, one row per
#'   site. Must contain columns `deploymentID`, `deploymentStart`, and
#'   `deploymentEnd` (or equivalents specified via the corresponding
#'   arguments).
#' @param observations A dataframe of observation records. Must contain
#'   columns `deploymentID`, `eventStart`, `scientificName`, and `count`
#'   (or equivalents specified via the corresponding arguments).
#' @param deploymentID <[`data-masking`][rlang::args_data_masking]> Column
#'   name for sites (ARUs). Must be a factor with identical levels in
#'   `deployments` and `observations`. Default: `deploymentID`.
#' @param deploymentStart <[`data-masking`][rlang::args_data_masking]> Column
#'   name for deployment start dates in `deployments`. Must be a `Date`.
#'   Default: `deploymentStart`.
#' @param deploymentEnd <[`data-masking`][rlang::args_data_masking]> Column
#'   name for deployment end dates in `deployments`. Must be a `Date`.
#'   Default: `deploymentEnd`.
#' @param eventStart <[`data-masking`][rlang::args_data_masking]> Column name
#'   for observation timestamps in `observations`. Must be `POSIXt`. Default:
#'   `eventStart`.
#' @return The observations dataframe with out-of-window records removed.
#' @noRd
filter_observations_window <- function(
  deployments,
  observations,
  deploymentID = deploymentID,
  deploymentStart = deploymentStart,
  deploymentEnd = deploymentEnd,
  eventStart = eventStart
) {
  n_before <- nrow(observations)
  cli::cli_inform("Starting with {n_before} events...")
  observations <- observations |>
    dplyr::left_join(
      deployments |>
        dplyr::select(
          {{ deploymentID }},
          {{ deploymentStart }},
          {{ deploymentEnd }}
        ),
      by = dplyr::join_by({{ deploymentID }})
    ) |>
    dplyr::filter(dplyr::between(
      {{ eventStart }},
      {{ deploymentStart }},
      {{ deploymentEnd }}
    )) |>
    dplyr::select(-c({{ deploymentStart }}, {{ deploymentEnd }}))
  n_removed <- n_before - nrow(observations)
  if (n_removed) {
    cli::cli_inform(
      "{n_removed} event{?s} outside the deployment window \\
       {?was/were} removed ({nrow(observations)} remaining)..."
    )
  } else {
    cli::cli_inform("No events outside the deployment window...")
  }
  observations
}

#' Find the mode of an integer vector
#'
#' @param x An integer or numeric vector.
#' @return The most frequent value in `x`. If there are ties, the smallest
#'   value is returned.
#' @noRd
int_mode <- function(x) {
  tab <- tabulate(x + 1)
  which.max(tab) - 1
}

#' Return mean and SD of numeric columns
#'
#' @param df Dataframe of predictors and site identifiers.
#' @return Mean and SD of each predictor
#' @noRd
scaling_parameters <- function(df, cols) {
  sds <- sapply(df[cols], sd)
  zero_var <- names(sds[sds == 0])
  if (length(zero_var)) {
    cli::cli_abort(
      "The following continuous predictor{?s} {?has/have} zero variance: \\
      {.val {zero_var}}. Remove {?it/them} before passing to {.fun make_data}."
    )
  }
  dplyr::bind_rows(colMeans(df[cols]), sds) |>
    dplyr::mutate(component = c("mean", "sd"), .before = 1)
}
