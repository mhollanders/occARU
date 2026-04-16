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
  col_name <- rlang::as_name(rlang::enquo(col))
  df <- dplyr::mutate(df, {{ col }} := factor({{ col }}, levels = ref_levels))
  n_missing <- dplyr::filter(df, is.na({{ col }})) |> nrow()
  if (n_missing) {
    if (strict) {
      cli::cli_abort(
        "{.arg {df_name}} contains {n_missing} record{?s} with \\
        {.val {col_name}} values not found."
      )
    } else {
      cli::cli_warn(
        "{.arg {df_name}} contains {n_missing} record{?s} with \\
        {.val {col_name}} values not found. \\
        {cli::qty(n_missing)}{?This record/These records} will be ignored."
      )
      df <- tidyr::drop_na(df, {{ col }})
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

#' Check that unique deployments are aligned with the right season
#'
#' Aborts if deployments in one season start before the last deployment ended
#' in the previous season.
#'
#' @param deployments A dataframe of deployment information, one row per
#'   site and season. Must contain columns `deploymentID`, `deploymentStart`,
#'   `deploymentEnd`, and `.season`.
#' @param deploymentID <[`data-masking`][rlang::args_data_masking]> Column
#'   name for sites (ARUs). Default: `deploymentID`.
#' @param deploymentStart <[`data-masking`][rlang::args_data_masking]> Column
#'   name for deployment start dates. Must be a `Date`. Default:
#'   `deploymentStart`.
#' @param deploymentEnd <[`data-masking`][rlang::args_data_masking]> Column
#'   name for deployment end dates. Must be a `Date`. Default: `deploymentEnd`.
#' @param .season <[`data-masking`][rlang::args_data_masking]> Column
#'   name for season. Default: `.season`.
#' @noRd
check_deployment_times <- function(
  deployments,
  deploymentID = deploymentID,
  deploymentStart = deploymentStart,
  deploymentEnd = deploymentEnd,
  .season = .season
) {
  deployments <- dplyr::arrange(deployments, {{ deploymentID }}, {{ .season }})

  last_first <- deployments |>
    dplyr::summarise(
      last = max({{ deploymentEnd }}),
      first = min({{ deploymentStart }}),
      .by = c({{ deploymentID }}, {{ .season }})
    ) |>
    dplyr::mutate(
      prev_last = dplyr::lag(last),
      .by = {{ deploymentID }}
    )

  wrong_deps <- dplyr::left_join(
    deployments,
    last_first,
    by = dplyr::join_by({{ deploymentID }}, {{ .season }})
  ) |>
    dplyr::filter(
      {{ deploymentStart }} < prev_last,
      .by = {{ deploymentID }}
    ) |>
    dplyr::pull({{ deploymentID }})

  if (length(wrong_deps)) {
    cli::cli_abort(
      "The following {.arg deploymentID}{?s} in {.arg deployments} contain \\
      {.arg deploymentStart} that occur before the last {.arg deploymentEnd} \\
      of the previous season: {.val {wrong_deps}}. All deployment periods in \\
      each season must occur before the start of the next season."
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
#' @param site_predictors A dataframe of site-level predictors. Must contain
#'   columns `deploymentID` and `.season`.
#' @param deployments A dataframe of deployment information, one row per
#'   site. Must contain columns `deploymentID` and `season`.
#' @param deploymentID <[`data-masking`][rlang::args_data_masking]> Column
#'   name for sites (ARUs). Must be a factor with identical levels in
#'   `deployments` and `survey_predictors`. Default: `deploymentID`.
#' @param .season <[`data-masking`][rlang::args_data_masking]> Column
#'   name for season. Must be a factor with identical levels in `deployments`
#'   and `survey_predictors`. Default: `.season`.
#' @noRd
check_site_predictors_coverage <- function(
  site_predictors,
  deployments,
  deploymentID = deploymentID,
  .season = .season,
  verbose = TRUE
) {
  df_name <- deparse(substitute(site_predictors))
  site_lvl <- dplyr::pull(deployments, {{ deploymentID }}) |> levels()
  complete <- dplyr::right_join(
    site_predictors,
    dplyr::distinct(deployments, {{ deploymentID }}, {{ .season }}),
    by = dplyr::join_by({{ deploymentID }}, {{ .season }})
  ) |>
    tidyr::drop_na() |>
    dplyr::pull({{ deploymentID }})
  missing <- setdiff(site_lvl, complete)
  if (length(missing)) {
    cli::cli_abort(
      "{.arg {df_name}} has incomplete predictor values for the following \\
      {.arg deploymentID}: {.vals {missing}}. Predictor values must be \\
      supplied for each season a {.arg deploymentID} was active."
    )
  } else if (verbose) {
    cli::cli_alert_success(
      "{.arg {df_name}} covers the full deployment period."
    )
  }
}

#' Check that survey predictors cover the full deployment period for each site
#'
#' @param survey_predictors A dataframe of survey-level predictors. Must
#'    contain columns `deploymentID`, `.season`, and `date`.
#' @param deployments A dataframe of deployment information, one row per
#'   site. Must contain columns `deploymentID`, `deploymentStart`,
#'   `deploymentEnd`, and `season`.
#' @param deploymentID <[`data-masking`][rlang::args_data_masking]> Column
#'   name for sites (ARUs). Must be a factor with identical levels in
#'   `deployments` and `survey_predictors`. Default: `deploymentID`.
#' @param deploymentStart `Date`. <[`data-masking`][rlang::args_data_masking]>
#'   Column name for deployment start dates in `deployments`. Default:
#'   `deploymentStart`.
#' @param deploymentEnd `Date`. <[`data-masking`][rlang::args_data_masking]>
#'   Column name for deployment end dates in `deployments`. Default:
#'   `deploymentEnd`.
#' @param .season <[`data-masking`][rlang::args_data_masking]> Column
#'   name for season. Must be a factor with identical levels in `deployments`
#'   and `survey_predictors`. Default: `.season`.
#' @param date `Date`. <[`data-masking`][rlang::args_data_masking]> Column name
#'   for dates in `survey_predictors`. Default: `date`.
#' @noRd
check_survey_predictors_coverage <- function(
  survey_predictors,
  deployments,
  deploymentID = deploymentID,
  deploymentStart = deploymentStart,
  deploymentEnd = deploymentEnd,
  .season = .season,
  date = date,
  verbose = TRUE
) {
  first_last <- survey_predictors |>
    dplyr::summarise(
      first = min({{ date }}),
      last = max({{ date }}),
      .by = c({{ deploymentID }}, {{ .season }})
    ) |>
    dplyr::left_join(
      deployments |>
        dplyr::select(
          {{ deploymentID }},
          {{ deploymentStart }},
          {{ deploymentEnd }},
          {{ .season }}
        ),
      by = dplyr::join_by({{ deploymentID }}, {{ .season }})
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

#' Check that end dates are after start dates
#'
#' @param df A dataframe.
#' @param deploymentID <[`data-masking`][rlang::args_data_masking]> Column
#'   name for sites (ARUs). Default: `deploymentID`.
#' @param start <[`data-masking`][rlang::args_data_masking]> `Date`. Column name
#'   for start dates to check. Default: `deploymentStart`.
#' @param end <[`data-masking`][rlang::args_data_masking]> `Date`. Column name
#'   for end dates to check. Default: `deploymentEnd`.
#' @noRd
check_dates <- function(
  df,
  deploymentID = deploymentID,
  start = deploymentStart,
  end = deploymentEnd
) {
  df_name <- deparse(substitute(df))
  wrong_sites <- df |>
    dplyr::filter({{ end }} <= {{ start }}) |>
    dplyr::distinct({{ deploymentID }}) |>
    dplyr::pull({{ deploymentID }})
  if (length(wrong_sites)) {
    cli::cli_abort(
      "The following {.arg deploymentID} in {.arg {df_name}} have end \\
      times before start times: {.val {wrong_sites}}."
    )
  }
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
  eventStart = eventStart,
  .season = .season
) {
  n_before <- nrow(observations)
  cli::cli_inform("Starting with {n_before} events...")
  observations <- dplyr::left_join(
    observations,
    deployments |>
      dplyr::select(
        {{ deploymentID }},
        {{ deploymentStart }},
        {{ deploymentEnd }},
        {{ .season }}
      ),
    by = dplyr::join_by({{ deploymentID }}, {{ .season }})
  ) |>
    dplyr::filter(
      dplyr::between(
        lubridate::as_date({{ eventStart }}),
        {{ deploymentStart }},
        {{ deploymentEnd }},
      ),
    ) |>
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
