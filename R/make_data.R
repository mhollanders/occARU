#' Prepare data for the occARU Stan model
#'
#' Transforms raw deployment and observation data into a named list suitable
#' for passing to [fit_model()]. Follows the
#' [camtrapDP](https://tdwg.github.io/camtrap-dp/) data format by default.
#' Site coordinates are automatically projected from WGS84 latitude/longitude
#' to UTM (km), with the zone auto-detected from the mean longitude.
#'
#' @param deployments A data frame of deployment information, one row per
#'   site. Must contain columns `deploymentID`, `deploymentStart`, and
#'   `deploymentEnd` (or equivalents specified via the corresponding
#'   arguments). Optionally, `latitude` and `longitude` columns enable
#'   the spatial Gaussian process.
#' @param observations A data frame of observation records. Must contain
#'   columns `deploymentID`, `eventStart`, `scientificName`, and `count`
#'   (or equivalents specified via the corresponding arguments).
#' @param failures Optional data frame of recorder failure periods. If
#'   supplied, must have the same `deploymentID` factor levels as
#'   `deployments`, with each row corresponding to one failure period at a
#'   `deploymentID` from `failureStart` to `failureEnd`.
#' @param deploymentID <[`data-masking`][rlang::args_data_masking]> Column
#'   name for site IDs. Must be a factor with identical levels in `deployments`
#'   and `observations`. Default: `deploymentID`.
#' @param deploymentStart <[`data-masking`][rlang::args_data_masking]> Column
#'   name for deployment start dates in `deployments`. Must be a `Date`.
#'   Default: `deploymentStart`.
#' @param deploymentEnd <[`data-masking`][rlang::args_data_masking]> Column
#'   name for deployment end dates in `deployments`. Must be a `Date`.
#'   Default: `deploymentEnd`.
#' @param eventStart <[`data-masking`][rlang::args_data_masking]> Column name
#'   for observation timestamps in `observations`. Must be `POSIXt`. Default:
#'   `eventStart`.
#' @param scientificName <[`data-masking`][rlang::args_data_masking]> Column
#'   name for species names in `observations`. Must be a factor. Default:
#'   `scientificName`.
#' @param count <[`data-masking`][rlang::args_data_masking]> Column name for
#'   number of individuals per observation record. Default: `count`.
#' @param latitude <[`data-masking`][rlang::args_data_masking]> Column name
#'   for WGS84 latitude in `deployments`. If omitted alongside `longitude`, no
#'   spatial Gaussian process is fitted. Default: `latitude`.
#' @param longitude <[`data-masking`][rlang::args_data_masking]> Column name
#'   for WGS84 longitude in `deployments`. Default: `longitude`.
#' @param survey_length Positive integer defining the length of each survey
#'   period in days. Observations are aggregated within each survey period by
#'   summing `count`, and recording effort (`Delta`) is computed as the fraction
#'   of days within each period that the recorder was active. For example,
#'   `days = 7` aggregates to weekly survey periods, with `Delta` ranging from 0
#'   (recorder failed all week) to 1 (recorder active all week). Longer periods
#'   reduce the number of surveys `J` but increase the counts per survey,
#'   trading off temporal resolution against model complexity and the assumption
#'   of closure within a survey period. Default: `7`.
#' @param thin_minutes Positive numeric. If supplied, observations within
#'   `thin_minutes` minutes of each other (per site and species) are thinned to
#'   a single observation, retaining the record with the highest `count`.
#'   Default: `30`.
#' @param day_start Character. Whether survey days start at `"midnight"` or
#'   `"midday"`. Default: `"midday"`.
#' @param reference_date A `Date` defining the start of the first survey
#'   period. Defaults to the earliest deployment start date.
#' @param failureStart <[`data-masking`][rlang::args_data_masking]> Column
#'   name for failure start dates in `failures`. Must be a `Date`. Default:
#'   `failureStart`.
#' @param failureEnd <[`data-masking`][rlang::args_data_masking]> Column
#'   name for failure end dates (inclusive) in `failures`. Must be a `Date`.
#'   Default: `failureEnd`.
#' @param occupancy_site_predictors Optional data frame of site-level
#'   covariates for the occupancy submodel. Must contain a `deploymentID`
#'   column with the same factor levels as `deployments`. Predictor columns
#'   must be `numeric` (continuous), `factor` (unordered categorical), or
#'   `ordered factor` (ordinal).
#' @param detection_site_predictors Optional data frame of site-level
#'   covariates for the detection submodel. Same column-type rules as
#'   `occupancy_site_predictors`. If identical to
#'   `occupancy_site_predictors`, the same matrices are reused.
#' @param survey_predictors Optional data frame of site-by-survey level
#'   covariates, with one row per site and date. Must contain `deploymentID`
#'   and `date` columns. Predictor columns follow the same type rules as the
#'   site-level predictor data frames. Must cover the full deployment period.
#' @param date <[`data-masking`][rlang::args_data_masking]> Column name for
#'   dates in `survey_predictors`. Default: `date`.
#' @param survey_summary An optional named list mapping continuous survey
#'   predictor column names to summary functions, used when aggregating survey
#'   predictors over `survey_length`-length periods. Each value can be a
#'   function name as a string (e.g. `"sum"`) or a function object (e.g. `sum`).
#'   Numeric predictors not named in `survey_summary` are summarised with
#'   `mean`; categorical and ordinal predictors are summarised with the
#'   modal value. Default: `NULL`.
#' @param scale_predictors Logical. If `TRUE`, continuous predictors are are
#'   scaled to zero mean and unit variance. Scaling parameters are stored as an
#'   attribute. Default: `TRUE`.
#'
#' @return A named list of class `"occARU_data"` containing all inputs
#'   required by the occARU Stan model, except for model specification
#'   arguments which are added by [fit_model()]. The list contains:
#'   \describe{
#'     \item{`I`}{Number of sites.}
#'     \item{`J`}{Number of survey periods.}
#'     \item{`S`}{Number of species.}
#'     \item{`Delta`}{`[I, J]` matrix of recording effort (0-1).}
#'     \item{`y`}{`[I, J, S]` array of detection counts.}
#'     \item{`XY`}{`[I, 2]` matrix of UTM coordinates in km, or zeros if
#'       coordinates not supplied.}
#'     \item{`P`}{Integer vector of length 3: number of continuous predictors
#'       for occupancy, and site and survey detection.}
#'     \item{`P_cat`}{Integer vector of length 3: number of categorical
#'       predictors for each component.}
#'     \item{`P_ord`}{Integer vector of length 3: number of ordinal predictors
#'       for each component.}
#'     \item{`X1`}{`[P[1], I]` occupancy continuous design matrix.}
#'     \item{`X_cat1`}{`[P_cat[1], I]` occupancy categorical integer matrix.}
#'     \item{`X_ord1`}{`[P_ord[1], I]` occupancy ordinal integer matrix.}
#'     \item{`X2`}{`[P[2], I]` site-level detection continuous design matrix.}
#'     \item{`X_cat2`}{`[P_cat[2], I]` site-level detection categorical integer
#'        matrix.}
#'     \item{`X_ord2`}{`[P_ord[2], I]` site-level detection ordinal integer
#'       matrix.}
#'     \item{`X3`}{`[I, P[3], J]` site-by-survey level detection continuous
#'       array.}
#'     \item{`X_cat3`}{`[I, P_cat[3], J]` site-by-survey categorical integer
#'       array.}
#'     \item{`X_ord3`}{`[I, P_ord[3], J]` site-by-survey survey ordinal integer
#'       array.}
#'   }
#'
#' @seealso [fit_model()]
#' @export
make_data <- function(
  deployments,
  observations,
  failures = NULL,
  deploymentID = deploymentID,
  latitude = latitude,
  longitude = longitude,
  deploymentStart = deploymentStart,
  deploymentEnd = deploymentEnd,
  eventStart = eventStart,
  scientificName = scientificName,
  count = count,
  failureStart = failureStart,
  failureEnd = failureEnd,
  survey_length = 7,
  thin_minutes = 30,
  day_start = "midday",
  reference_date,
  occupancy_site_predictors = NULL,
  detection_site_predictors = NULL,
  survey_predictors = NULL,
  date = date,
  survey_summary = NULL,
  scale_predictors = TRUE
) {
  # --- capture column names as strings ----------------------------------------
  dep_id_chr <- rlang::as_name(rlang::enquo(deploymentID))
  dep_start_chr <- rlang::as_name(rlang::enquo(deploymentStart))
  dep_end_chr <- rlang::as_name(rlang::enquo(deploymentEnd))
  event_chr <- rlang::as_name(rlang::enquo(eventStart))
  species_chr <- rlang::as_name(rlang::enquo(scientificName))
  count_chr <- rlang::as_name(rlang::enquo(count))
  date_chr <- rlang::as_name(rlang::enquo(date))
  fail_start_chr <- rlang::as_name(rlang::enquo(failureStart))
  fail_end_chr <- rlang::as_name(rlang::enquo(failureEnd))
  lat_chr <- tryCatch(
    rlang::as_name(rlang::enquo(latitude)),
    error = function(e) NULL
  )
  lon_chr <- tryCatch(
    rlang::as_name(rlang::enquo(longitude)),
    error = function(e) NULL
  )

  # --- deployments checks -----------------------------------------------------
  check_no_duplicates(deployments, dep_id_chr)
  check_cols_exist(deployments, dep_id_chr, dep_start_chr, dep_end_chr)
  check_cols_class(deployments, "Date", dep_start_chr, dep_end_chr)
  if (!is.factor(deployments |> dplyr::pull({{ deploymentID }}))) {
    deployments <- deployments |>
      dplyr::mutate({{ deploymentID }} := factor({{ deploymentID }}))
  } else {
    check_no_empty_levels(deployments, dep_id_chr)
  }
  site_lvl <- deployments |> dplyr::pull({{ deploymentID }}) |> levels()
  I <- length(site_lvl)
  if (I == 1L) {
    cli::cli_abort("occARU requires more than one site.")
  }

  # --- observations -----------------------------------------------------------
  check_cols_exist(observations, dep_id_chr, event_chr, species_chr, count_chr)
  check_no_empty_levels(observations, dep_id_chr)
  check_cols_class(observations, "POSIXt", event_chr)
  observations <- align_factor(observations, dep_id_chr, site_lvl)
  if (!is.factor(observations |> dplyr::pull({{ scientificName }}))) {
    observations <- observations |>
      dplyr::mutate({{ scientificName }} := factor({{ scientificName }}))
  } else {
    check_no_empty_levels(observations, species_chr)
  }
  species_lvl <- observations |> dplyr::pull({{ scientificName }}) |> levels()
  S <- length(species_lvl)

  # --- failures ---------------------------------------------------------------
  if (!is.null(failures)) {
    check_cols_exist(failures, dep_id_chr, fail_start_chr, fail_end_chr)
    check_cols_class(failures, "Date", fail_start_chr, fail_end_chr)
    failures <- align_factor(failures, dep_id_chr, site_lvl)
  }

  # --- scalar arguments -------------------------------------------------------
  day_start <- match.arg(day_start, c("midnight", "midday"))
  if (survey_length %% 1 != 0 || survey_length < 1) {
    cli::cli_abort("{.arg survey_length} must be a positive integer.")
  }
  if (!is.null(thin_minutes) && thin_minutes <= 0) {
    cli::cli_abort("{.arg thin_minutes} must be a positive number of minutes.")
  }

  # --- reference date ---------------------------------------------------------
  if (missing(reference_date)) {
    reference_date <- deployments |>
      dplyr::pull({{ deploymentStart }}) |>
      min()
  } else if (!inherits(reference_date, "Date")) {
    cli::cli_abort("{.arg reference_date} must be a {.cls Date}.")
  }

  # --- predictors -------------------------------------------------------------
  if (!is.null(occupancy_site_predictors)) {
    check_cols_exist(occupancy_site_predictors, dep_id_chr)
    check_no_duplicates(occupancy_site_predictors, dep_id_chr)
    check_mixed_predictors(
      occupancy_site_predictors |> dplyr::select(-{{ deploymentID }})
    )
    occupancy_site_predictors <- align_factor(
      occupancy_site_predictors,
      dep_id_chr,
      site_lvl,
      strict = TRUE
    )
  }
  if (
    !is.null(detection_site_predictors) &&
      !isTRUE(all.equal(occupancy_site_predictors, detection_site_predictors))
  ) {
    check_cols_exist(detection_site_predictors, dep_id_chr)
    check_no_duplicates(detection_site_predictors, dep_id_chr)
    check_mixed_predictors(
      detection_site_predictors |> dplyr::select(-{{ deploymentID }})
    )
    detection_site_predictors <- align_factor(
      detection_site_predictors,
      dep_id_chr,
      site_lvl,
      strict = TRUE
    )
  }
  if (!is.null(survey_predictors)) {
    check_cols_exist(survey_predictors, dep_id_chr, date_chr)
    check_mixed_predictors(
      survey_predictors |> dplyr::select(-c({{ deploymentID }}, {{ date }}))
    )
    survey_predictors <- align_factor(
      survey_predictors,
      dep_id_chr,
      site_lvl,
      strict = TRUE
    )
    check_survey_predictors_coverage(
      survey_predictors,
      deployments,
      dep_id_chr,
      dep_start_chr,
      dep_end_chr,
      date_chr
    )
  }

  # --- deployment grid and Delta ----------------------------------------------
  daily_grid <- deployments |>
    dplyr::select(
      {{ deploymentID }},
      {{ deploymentStart }},
      {{ deploymentEnd }}
    ) |>
    dplyr::mutate(
      date = list(seq.Date(
        {{ deploymentStart }},
        {{ deploymentEnd }},
        by = "day"
      )),
      .by = {{ deploymentID }}
    ) |>
    tidyr::unnest(c(date))

  if (is.null(failures)) {
    daily_grid <- daily_grid |> dplyr::mutate(Delta = 1L)
  } else {
    failure_dates <- failures |>
      dplyr::select({{ deploymentID }}, {{ failureStart }}, {{ failureEnd }}) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        date = list(seq.Date({{ failureStart }}, {{ failureEnd }}, by = "day"))
      ) |>
      tidyr::unnest(c(date)) |>
      dplyr::mutate(failure = 1L) |>
      dplyr::select({{ deploymentID }}, date, failure)

    daily_grid <- daily_grid |>
      dplyr::left_join(
        failure_dates,
        by = dplyr::join_by({{ deploymentID }}, date)
      ) |>
      dplyr::mutate(Delta = dplyr::if_else(is.na(failure), 1L, 0L)) |>
      dplyr::select(-failure)
  }

  deployments_aggregated <- daily_grid |>
    dplyr::mutate(
      survey = aggregate_by_days(date, reference_date, survey_length)
    ) |>
    dplyr::summarise(
      Delta = sum(Delta) / survey_length,
      .by = c({{ deploymentID }}, survey)
    )

  Delta <- deployments_aggregated |>
    tidyr::pivot_wider(
      names_from = {{ deploymentID }},
      values_from = Delta,
      values_fill = 0,
      names_sort = TRUE
    ) |>
    dplyr::arrange(survey) |>
    tibble::column_to_rownames("survey") |>
    t()
  surveys <- colnames(Delta)
  J <- length(surveys)
  if (J == 1L) {
    cli::cli_abort(
      "occARU requires more than one survey. Did you make the survey length \\
      too long? Reduce {.arg survey_length} and try again."
    )
  }

  # --- detection history ------------------------------------------------------
  if (day_start == "midday") {
    observations <- observations |>
      dplyr::mutate(
        {{ eventStart }} := {{ eventStart }} -
          lubridate::hours(12)
      )
  }

  sites_without_observations(deployments, observations, dep_id_chr)
  observations <- filter_observations_window(
    deployments,
    observations,
    dep_id_chr,
    dep_start_chr,
    dep_end_chr,
    event_chr
  )

  if (!is.null(thin_minutes)) {
    n_before <- nrow(observations)
    observations <- observations |>
      dplyr::mutate(
        cluster = assign_clusters({{ eventStart }}, thin_minutes),
        .by = c({{ deploymentID }}, {{ scientificName }})
      ) |>
      dplyr::slice_max(
        {{ count }},
        with_ties = FALSE,
        by = c({{ deploymentID }}, {{ scientificName }}, cluster)
      )
    n_removed <- n_before - nrow(observations)
    cli::cli_inform(
      "{n_removed} event{?s} removed by {thin_minutes}-minute thinning window \\
      ({nrow(observations)} remaining)..."
    )
  }

  observations_aggregated <- observations |>
    dplyr::mutate(
      date = lubridate::as_date({{ eventStart }}),
      survey = aggregate_by_days(date, reference_date, survey_length)
    ) |>
    dplyr::filter(survey %in% lubridate::ymd(surveys)) |>
    dplyr::summarise(
      n = sum({{ count }}),
      .by = c({{ deploymentID }}, survey, {{ scientificName }})
    )

  y <- observations_aggregated |>
    tidyr::complete(
      {{ scientificName }} := factor(species_lvl, species_lvl),
      survey := lubridate::ymd(surveys),
      {{ deploymentID }} := factor(site_lvl, site_lvl),
      fill = list(n = 0)
    ) |>
    dplyr::arrange({{ scientificName }}, survey, {{ deploymentID }}) |>
    dplyr::pull(n) |>
    array(c(I, J, S), dimnames = list(site_lvl, surveys, species_lvl))

  # --- coordinates ------------------------------------------------------------
  if (!is.null(lat_chr) && !is.null(lon_chr)) {
    utm <- coords_to_utm(deployments, dep_id_chr, lon_chr, lat_chr)
    XY <- utm$XY
    utm_crs <- utm$utm_crs
  } else {
    XY <- matrix(0, I, 2, dimnames = list(site_lvl, c("X", "Y")))
    utm_crs <- NA_character_
  }

  # --- encode predictors ------------------------------------------------------
  enc1 <- encode_predictors(
    occupancy_site_predictors,
    "site",
    dep_id_chr,
    scale_predictors = scale_predictors
  )
  enc2 <- if (
    isTRUE(all.equal(occupancy_site_predictors, detection_site_predictors))
  ) {
    enc1
  } else {
    encode_predictors(
      detection_site_predictors,
      "site",
      dep_id_chr,
      scale_predictors = scale_predictors
    )
  }
  enc3 <- encode_predictors(
    survey_predictors,
    "survey",
    dep_id_chr,
    date_chr = date_chr,
    scale_predictors = scale_predictors,
    survey_summary = survey_summary,
    survey_length = survey_length,
    reference_date = reference_date
  )

  empty_site_matrix <- matrix(
    0L,
    nrow = 0,
    ncol = I,
    dimnames = list(NULL, site_lvl)
  )
  empty_survey_array <- array(
    0L,
    dim = c(I, 0, J),
    dimnames = list(site_lvl, NULL, as.character(surveys))
  )

  # --- predictor dimensions ---------------------------------------------------
  P <- c(
    if (!is.null(enc1)) nrow(enc1$X) else 0L,
    if (!is.null(enc2)) nrow(enc2$X) else 0L,
    if (!is.null(enc3)) dim(enc3$X)[2] else 0L
  )
  P_cat <- c(
    if (!is.null(enc1)) nrow(enc1$X_cat) else 0L,
    if (!is.null(enc2)) nrow(enc2$X_cat) else 0L,
    if (!is.null(enc3)) dim(enc3$X_cat)[2] else 0L
  )
  P_ord <- c(
    if (!is.null(enc1)) nrow(enc1$X_ord) else 0L,
    if (!is.null(enc2)) nrow(enc2$X_ord) else 0L,
    if (!is.null(enc3)) dim(enc3$X_ord)[2] else 0L
  )

  # --- return -----------------------------------------------------------------
  structure(
    list(
      I = I,
      J = J,
      S = S,
      Delta = Delta,
      XY = XY,
      y = y[,, 1:S],
      P = P,
      P_cat = P_cat,
      P_ord = P_ord,
      X1 = enc1$X %||% empty_site_matrix,
      X_cat1 = enc1$X_cat %||% empty_site_matrix,
      X_ord1 = enc1$X_ord %||% empty_site_matrix,
      X2 = enc2$X %||% empty_site_matrix,
      X_cat2 = enc2$X_cat %||% empty_site_matrix,
      X_ord2 = enc2$X_ord %||% empty_site_matrix,
      X3 = enc3$X %||% empty_survey_array,
      X_cat3 = enc3$X_cat %||% empty_survey_array,
      X_ord3 = enc3$X_ord %||% empty_survey_array
    ),
    class = "occARU_data",
    sites = site_lvl,
    surveys = surveys,
    species = species_lvl,
    survey_length = survey_length,
    thin_minutes = thin_minutes,
    reference_date = reference_date,
    utm_crs = utm_crs,
    scale = if (scale_predictors) {
      list(
        X1 = enc1$scale_params,
        X2 = enc2$scale_params,
        X3 = enc3$scale_params
      )
    } else {
      NULL
    },
    levels = list(
      X1 = enc_levels(enc1),
      X2 = enc_levels(enc2),
      X3 = enc_levels(enc3)
    )
  ) |>
    print()
}

#' Print method for occARU_data objects
#'
#' @param x A `occARU_data` object.
#' @param ... Ignored.
#' @noRd
print.occARU_data <- function(x, ...) {
  cli::cli_h1("occARU data")
  P <- sapply(1:3, \(p) x$P[p] + x$P_cat[p] + x$P_ord[p])
  cli::cli_dl(c(
    "Sites (I)" = "{x$I}",
    "Surveys (J)" = "{x$J}",
    "Species (S)" = "{x$S}",
    "Detections" = "{sum(x$y)}",
    "Site coordinates" = "{ifelse(any(x$XY != 0), 'yes', 'no')}",
    "Deployment span" = "{attr(x, 'surveys')[1]} to \\
                        {attr(x, 'surveys')[x$J]}",
    "Survey length" = "{attr(x, 'survey_length')} day{?s}",
    "Thinning" = "{ifelse(!is.null(attr(x, 'thin_minutes')),
                            paste(attr(x, 'thin_minutes'), 'minutes'), 'none')}",
    "Occupancy predictors" = P[1]
  ))
  if (P[1]) {
    cli::cli_bullets(c(
      " " = "Continuous: {x$P[1]}",
      " " = "Categorical: {x$P_cat[1]}",
      " " = "Ordinal: {x$P_ord[1]}"
    ))
  }
  cli::cli_dl(c(
    "Detection site predictors" = P[2]
  ))
  if (P[2]) {
    cli::cli_bullets(c(
      " " = "Continuous: {x$P[2]}",
      " " = "Categorical: {x$P_cat[2]}",
      " " = "Ordinal: {x$P_ord[2]}"
    ))
  }
  cli::cli_dl(c(
    "Survey predictors" = P[3]
  ))
  if (P[3]) {
    cli::cli_bullets(c(
      " " = "Continuous: {x$P[3]}",
      " " = "Categorical: {x$P_cat[3]}",
      " " = "Ordinal: {x$P_ord[3]}"
    ))
  }
  invisible(x)
}
