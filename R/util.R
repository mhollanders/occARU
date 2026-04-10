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

#' Align a factor column to reference levels, warning or erroring on missing values
#'
#' @param df A data frame.
#' @param col Character. Name of the column to align.
#' @param levels Reference levels.
#' @param strict If `TRUE`, missing values trigger an error. If `FALSE`,
#'   missing values trigger a warning and affected rows are removed.
#' @return The data frame with `col` realigned and missing values removed.
#' @noRd
align_factor <- function(df, col, levels, strict = FALSE) {
  df_name <- deparse(substitute(df))
  df[[col]] <- factor(df[[col]], levels = levels)
  n_missing <- sum(is.na(df[[col]]))
  if (n_missing) {
    if (strict) {
      cli::cli_abort(
        "{.arg {df_name}} contains {n_missing} record{?s} with \\
        {.val {col}} values not found in {.arg deployments}."
      )
    } else {
      cli::cli_warn(
        "{.arg {df_name}} contains {n_missing} record{?s} with \\
        {.val {col}} values not found in {.arg deployments}. \\
        {?This record/These records} will be ignored."
      )
      df <- df[!is.na(df[[col]]), ]
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

#' Check that columns exist in a data frame
#'
#' @param df A data frame.
#' @param ... Column names as strings.
#' @noRd
check_cols_exist <- function(df, ...) {
  df_name <- deparse(substitute(df))
  cols <- c(...)
  missing_cols <- cols[!rlang::has_name(df, cols)]
  if (length(missing_cols)) {
    cli::cli_abort(
      "The following columns were not found in {.arg {df_name}}: \\
       {.val {missing_cols}}."
    )
  }
}

#' Check that columns inherit a given class
#'
#' @param df A data frame.
#' @param class Expected class as a character string.
#' @param ... Column names as strings.
#' @noRd
check_cols_class <- function(df, class, ...) {
  df_name <- deparse(substitute(df))
  cols <- c(...)
  wrong_cols <- cols[!sapply(cols, \(col) inherits(df[[col]], class))]
  if (length(wrong_cols)) {
    cli::cli_abort(
      "The following columns in {.arg {df_name}} must be {.cls {class}}: \\
       {.val {wrong_cols}}."
    )
  }
}

#' Check that all predictor columns are numeric, factor, or ordered factor
#'
#' @param df A data frame of predictors (ID column already removed).
#' @param arg_name Name of the dataframe to return.
#' @noRd
check_mixed_predictors <- function(df, arg_name) {
  bad <- df |>
    dplyr::select(dplyr::where(\(x) !is.numeric(x) && !is.factor(x))) |>
    colnames()
  if (length(bad)) {
    cli::cli_abort(
      "{.arg {arg_name}} contains columns that are not {.cls numeric}, \\
       {.cls factor}, or {.cls ordered}: {.val {bad}}. \\
       Convert them before passing to {.fun make_data}."
    )
  }
}

#' Check that a data frame has no duplicate values in a column
#'
#' @param df A data frame.
#' @param col Character. Name of the column to check.
#' @noRd
check_no_duplicates <- function(df, col) {
  df_name <- deparse(substitute(df))
  dupes <- df |>
    dplyr::count(.data[[col]]) |>
    dplyr::filter(n > 1) |>
    dplyr::pull(.data[[col]])
  if (length(dupes)) {
    cli::cli_abort(
      "{.arg {df_name}} contains duplicate {.val {col}} values: {.val {dupes}}."
    )
  }
}

#' Check that all factor levels are present in the data
#'
#' Aborts if a factor column contains levels with no corresponding rows,
#' which typically occurs when a data frame has been filtered after the
#' factor was defined.
#'
#' @param df A data frame.
#' @param col_chr Column name as a string.
#'
#' @return `df`, invisibly, if all levels are present.
#' @keywords internal
check_no_empty_levels <- function(df, col_chr) {
  df_name <- deparse(substitute(df))
  col <- df[[col_chr]]
  if (!is.factor(col)) {
    return(invisible(df))
  }
  empty <- setdiff(levels(col), as.character(col))
  if (length(empty)) {
    cli::cli_abort(
      c(
        "{.arg {df_name}} has unused {.code {col_chr}} factor level{?s}:",
        "x" = "{.val {empty}}",
        "i" = "Drop empty levels with {.code droplevels()} before calling
               {.fn make_data}."
      )
    )
  }
  invisible(df)
}

#' Check that survey predictors cover the full deployment period for each site
#'
#' @param survey_predictors A data frame of survey-level predictors.
#' @param deployments A data frame of deployment information.
#' @param dep_id_chr Character. Name of the deployment ID column.
#' @param dep_start_chr Character. Name of the deployment start date column.
#' @param dep_end_chr Character. Name of the deployment end date column.
#' @param date_chr Character. Name of the date column in survey predictors.
#' @noRd
check_survey_predictors_coverage <- function(
  survey_predictors,
  deployments,
  dep_id_chr,
  dep_start_chr,
  dep_end_chr,
  date_chr
) {
  first_last <- survey_predictors |>
    dplyr::summarise(
      first = min(.data[[date_chr]]),
      last = max(.data[[date_chr]]),
      .by = dplyr::all_of(dep_id_chr)
    ) |>
    dplyr::left_join(
      deployments |>
        dplyr::select(dplyr::all_of(c(dep_id_chr, dep_start_chr, dep_end_chr))),
      by = dep_id_chr
    )

  incomplete <- first_last |>
    dplyr::filter(
      first > .data[[dep_start_chr]] |
        last < .data[[dep_end_chr]]
    ) |>
    dplyr::pull(.data[[dep_id_chr]])

  if (length(incomplete)) {
    cli::cli_abort(
      "{.arg survey_predictors} does not cover the full deployment period \\
       for the following site{?s}: {.val {incomplete}}."
    )
  }
}

#' Classify predictor columns by type
#'
#' @param df A data frame of predictor columns only.
#' @return A named list with elements `numeric`, `categorical`, `ordinal`.
#' @noRd
classify_predictors <- function(df) {
  nms <- colnames(df)
  list(
    numeric = nms[sapply(nms, \(n) is.numeric(df[[n]]))],
    categorical = nms[sapply(nms, \(n) {
      is.factor(df[[n]]) &&
        !is.ordered(df[[n]])
    })],
    ordinal = nms[sapply(nms, \(n) is.ordered(df[[n]]))]
  )
}

#' Project site coordinates from WGS84 to UTM
#'
#' Converts longitude/latitude columns in a data frame to UTM coordinates
#' (km), with the zone auto-detected from the mean longitude. Returns both
#' the projected coordinate matrix and the CRS string used.
#'
#' @param deployments A data frame containing longitude and latitude columns.
#' @param lon_chr String. Name of the longitude column (WGS84 decimal degrees).
#' @param lat_chr String. Name of the latitude column (WGS84 decimal degrees).
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{`XY`}{Numeric matrix of UTM coordinates in km, with columns
#'       \code{X} and \code{Y}.}
#'     \item{`utm_crs`}{Character string of the PROJ CRS used for the
#'       transformation.}
#'   }
#' @noRd
coords_to_utm <- function(deployments, dep_id_chr, lon_chr, lat_chr) {
  lons <- deployments[[lon_chr]]
  lats <- deployments[[lat_chr]]

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
    sf::st_as_sf(coords = c(lon_chr, lat_chr), crs = 4326) |>
    sf::st_transform(utm_crs) |>
    sf::st_coordinates()
  row.names(XY) <- levels(deployments[[dep_id_chr]])

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

#' Encode a predictor data frame for the occARU Stan model
#'
#' Splits a predictor data frame into separate matrices/arrays for continuous,
#' categorical, and ordinal predictors, optionally scaling continuous
#' predictors, and extracts factor levels for use in the Stan model and for
#' post-processing.
#'
#' @param df A data frame of predictors including the deployment ID column,
#'   and for survey mode a date column. If `NULL`, returns `NULL` immediately.
#' @param mode Character. Either `"site"` (returns `[I, P]` matrices) or
#'   `"survey"` (returns `[I, P, J]` arrays).
#' @param dep_id_chr Character. Name of the deployment ID column.
#' @param date_chr Character. Name of the date column. Only used when
#'   `mode = "survey"`.
#' @param date_chr Character. Name of the date column. Only used when
#'   `mode = "survey"`.
#' @param reference_date A `Date` defining the start of the first survey
#'   period.
#' @param scale_predictors Logical. If `TRUE` (default), continuous predictors
#'   are scaled to zero mean and unit variance. When `mode = "survey"`, scaling
#'   parameters are derived from the site-averaged values per survey period
#'   rather than the raw site-by-survey values.
#' @param survey_summary Optional named list mapping continuous survey predictor
#'   column names to aggregation functions, used when aggregating survey
#'   predictors over `survey_length`-length periods. Each value can be a
#'   function name as a string (e.g. `"sum"`) or a function object (e.g. `sum`).
#'   Numeric predictors not named in `survey_summary` default to `mean`; factor
#'   and ordered factor predictors default to the modal value. Default: `NULL`.
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
  dep_id_chr,
  date_chr = NULL,
  site_lvl,
  surveys = NULL,
  reference_date,
  scale_predictors = TRUE,
  survey_summary = NULL,
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
    # --- classify -------------------------------------------------------------
    ignore <- c(dep_id_chr, if (mode == "survey") date_chr)
    pred_cols <- df |>
      dplyr::select(-dplyr::all_of(ignore)) |>
      colnames()
    types <- classify_predictors(df |> dplyr::select(dplyr::all_of(pred_cols)))
    num_cols <- types$numeric
    cat_cols <- types$categorical
    ord_cols <- types$ordinal

    P <- length(num_cols)
    P_cat <- length(cat_cols)
    C <- lapply(df[cat_cols], levels)
    P_ord <- length(ord_cols)
    O <- lapply(df[ord_cols], levels)

    # --- checks ---------------------------------------------------------------
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
        {.val {double_ord}}. Remove {?it/them} before passing to \\
        {.fun make_data}."
      )
    }

    all_factor_cols <- c(cat_cols, ord_cols)
    empty_levels <- Filter(
      length,
      lapply(
        setNames(all_factor_cols, all_factor_cols),
        \(col) levels(df[[col]])[!levels(df[[col]]) %in% df[[col]]]
      )
    )
    if (length(empty_levels)) {
      detail <- mapply(
        \(col, lvls) cli::format_inline("{.val {col}}: {.val {lvls}}"),
        names(empty_levels),
        empty_levels
      )
      cli::cli_warn(c(
        "Some factor predictors have levels not present in the data:",
        setNames(detail, rep("*", length(detail)))
      ))
    }

    # --- convert factors to integers ------------------------------------------
    df <- df |>
      dplyr::mutate(
        dplyr::across(dplyr::all_of(all_factor_cols), as.integer)
      )

    # --- survey aggregation ---------------------------------------------------
    if (mode == "survey") {
      if (!is.null(survey_summary)) {
        unknown <- setdiff(names(survey_summary), pred_cols)
        if (length(unknown)) {
          cli::cli_abort(
            "The following predictors in {.arg survey_summary} were not found \\
           in {.arg survey_predictors}: {.val {unknown}}."
          )
        }
        invalid_fns <- names(survey_summary)[
          !sapply(survey_summary, \(f) {
            tryCatch(
              {
                match.fun(f)
                TRUE
              },
              error = \(e) FALSE
            )
          })
        ]
        if (length(invalid_fns)) {
          cli::cli_abort(
            "The following functions in {.arg survey_summary} are not valid: \\
            {.val {invalid_fns}}."
          )
        }
        non_num_overrides <- intersect(
          names(survey_summary),
          c(cat_cols, ord_cols)
        )
        if (length(non_num_overrides)) {
          cli::cli_warn(
            "Ignoring {.arg survey_summary} entries for non-numeric \\
           predictor{?s} {.val {non_num_overrides}}. Factor and ordinal \\
           predictors always use modal aggregation."
          )
        }
      }

      summary_fns <- setNames(
        lapply(pred_cols, \(col) {
          if (
            col %in%
              num_cols &&
              !is.null(survey_summary) &&
              col %in% names(survey_summary)
          ) {
            match.fun(survey_summary[[col]])
          } else if (col %in% num_cols) {
            mean
          } else {
            int_mode
          }
        }),
        pred_cols
      )

      df <- df |>
        dplyr::mutate(
          survey = aggregate_by_days(
            .data[[date_chr]],
            reference_date,
            survey_length
          )
        ) |>
        dplyr::summarise(
          dplyr::across(dplyr::all_of(pred_cols), \(x) {
            summary_fns[[dplyr::cur_column()]](x)
          }),
          .by = c(dplyr::all_of(dep_id_chr), "survey")
        )

      surveys <- unique(df$survey)
      J <- length(surveys)
    }

    # --- scale ------------------------------------------------------------------
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
        df_num <- df |>
          dplyr::select(-dplyr::all_of(c(cat_cols, ord_cols)))
        survey_means <- df_num |>
          dplyr::summarise(
            dplyr::across(dplyr::all_of(num_cols), mean),
            .by = "survey"
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

    # --- build matrices / arrays ------------------------------------------------
    if (mode == "site") {
      to_matrix <- function(df, cols) {
        if (length(cols)) {
          df |>
            dplyr::select(dplyr::all_of(c(dep_id_chr, cols))) |>
            column_to_rownames(dep_id_chr) |>
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
            dplyr::select(all_of(c(dep_id_chr, "survey", cols))) |>
            tidyr::complete(
              !!rlang::sym(dep_id_chr),
              survey,
              fill = setNames(
                lapply(cols, \(c) if (c %in% num_cols) 0 else 1L),
                cols
              )
            ) |>
            tidyr::pivot_longer(dplyr::all_of(cols), names_to = "p") |>
            dplyr::mutate(p = factor(p, levels = num_cols)) |>
            dplyr::arrange(survey, p, {{ dep_id_chr }}) |>
            dplyr::pull(value) |>
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
      X_ord = X_ord - 1,
      C = C,
      O = O,
      scale_params = scale_params
    )
  }
}

#' Remove observations outside deployment window
#'
#' @param deployments A data frame of deployment information.
#' @param observations A data frame of observation records.
#' @param dep_id_chr Character. Name of the deployment ID column.
#' @param dep_start_chr Character. Name of the deployment start date column.
#' @param dep_end_chr Character. Name of the deployment end date column.
#' @param event_chr Character. Name of the observation timestamp column.
#' @return The observations data frame with out-of-window records removed.
#' @noRd
filter_observations_window <- function(
  deployments,
  observations,
  dep_id_chr,
  dep_start_chr,
  dep_end_chr,
  event_chr
) {
  n_before <- nrow(observations)
  cli::cli_inform("Starting with {n_before} events...")
  observations <- observations |>
    dplyr::left_join(
      deployments |>
        dplyr::select(dplyr::all_of(c(dep_id_chr, dep_start_chr, dep_end_chr))),
      by = dep_id_chr
    ) |>
    dplyr::filter(dplyr::between(
      lubridate::as_date(.data[[event_chr]]),
      .data[[dep_start_chr]],
      .data[[dep_end_chr]]
    )) |>
    dplyr::select(-dplyr::all_of(c(dep_start_chr, dep_end_chr)))
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

#' Compare user-supplied priors against defaults and return names of
#' non-default priors
#'
#' @param prior A named list of user-supplied prior hyperparameters.
#' @param defaults A named list of default prior hyperparameters.
#' @return A character vector of prior names that differ from defaults.
#' @noRd
non_default_priors <- function(prior, defaults) {
  names(prior)[
    sapply(names(prior), \(nm) !isTRUE(all.equal(prior[[nm]], defaults[[nm]])))
  ]
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

#' List sites in deployments with no observations
#'
#' @param deployments A data frame of deployment information.
#' @param observations A data frame of observation records.
#' @param dep_id_chr Character. Name of the deployment ID column.
#' @return A character vector of site IDs with no observations, invisibly.
#' @noRd
sites_without_observations <- function(deployments, observations, dep_id_chr) {
  deployment_sites <- deployments[[dep_id_chr]]
  observation_sites <- unique(observations[[dep_id_chr]])
  empty_sites <- setdiff(deployment_sites, observation_sites)
  if (length(empty_sites)) {
    cli::cli_inform(c(
      "{length(empty_sites)} site{?s} in {.arg deployments} have no observations:",
      "*" = "{.val {empty_sites}}"
    ))
  }
  invisible(empty_sites)
}
