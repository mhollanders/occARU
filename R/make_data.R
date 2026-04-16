#' Prepare data for the occARU model
#'
#' Transforms raw deployment and observation data into a named list suitable
#' for passing to [fit_model()]. Follows the
#' [camtrapDP](https://tdwg.github.io/camtrap-dp/) data format by default.
#' Site coordinates are automatically projected from WGS84 latitude/longitude
#' to UTM (km), with the zone auto-detected from the mean longitude.
#'
#' @param deployments A dataframe of deployment information, one row per
#'   site (and potentially season). Must contain columns `deploymentID`,
#'   `deploymentStart`, and `deploymentEnd` (or equivalents specified via the
#'   corresponding arguments). Optionally, `latitude` and `longitude` columns
#'   enable the spatial Gaussian process. If multiple seasons, must also contain
#'   column `season`.
#' @param observations A dataframe of observation records. Must contain
#'   columns `deploymentID`, `eventStart`, `scientificName`, and `count`
#'   (or equivalents specified via the corresponding arguments).
#' @param failures Optional dataframe of ARU failure periods. Must contain
#'   columns `deploymentID`, `failureStart`, and `failureEnd`, with each row
#'   corresponding to one failure period at a `deploymentID` from `failureStart`
#'   to `failureEnd` (inclusive). See [find_failures()].
#' @param deploymentID <[`data-masking`][rlang::args_data_masking]> Column
#'   name for sites (ARUs). Retains levels if supplied as factor.  Default:
#'   `deploymentID`.
#' @param deploymentStart <[`data-masking`][rlang::args_data_masking]> `Date`.
#'   Column name for deployment start dates in `deployments`. Default:
#'   `deploymentStart`.
#' @param deploymentEnd <[`data-masking`][rlang::args_data_masking]> `Date.`
#'   Column name for deployment end dates in `deployments`. Default:
#'   `deploymentEnd`.
#' @param latitude <[`data-masking`][rlang::args_data_masking]> `numeric`.
#'   Column name for WGS84 latitude in `deployments`. If omitted alongside
#'   `longitude`, no spatial Gaussian process is fitted. Default: `latitude`.
#' @param longitude <[`data-masking`][rlang::args_data_masking]> `numeric`.
#'   Column name for WGS84 longitude in `deployments`. Default: `longitude`.
#' @param season <[`data-masking`][rlang::args_data_masking]> Optional column
#'   specifying season. The column must be a factor to ensure correct ordering.
#'   If the column is not present in `deployments`, all observations are treated
#'   as a single season. Default: `season`.
#' @param eventStart <[`data-masking`][rlang::args_data_masking]> `POSIXt`.
#'   Column name for observation timestamps in `observations`. Default:
#'   `eventStart`.
#' @param scientificName <[`data-masking`][rlang::args_data_masking]> Column
#'   name for species names in `observations`. Retains levels if supplied as
#'   factor. Default: `scientificName`.
#' @param count <[`data-masking`][rlang::args_data_masking]> `integerish`.
#'   Column name for number of individuals per observation record. Default:
#'   `count`.
#' @param failureStart <[`data-masking`][rlang::args_data_masking]> `Date`.
#'   Column name for failure start dates in `failures`. Default:
#'   `failureStart`.
#' @param failureEnd <[`data-masking`][rlang::args_data_masking]> `Date`. Column
#'   name for failure end dates (inclusive) in `failures`. Default:
#'   `failureEnd`.
#' @param survey_length Positive integer. Defines the length of each survey
#'   period in days. Observations are aggregated within each survey period by
#'   summing `count`, and recording effort (`Delta`) is computed as the fraction
#'   of the survey length the ARU was active. For example, `survey_length = 7`
#'   aggregates to weekly survey periods, with `Delta` ranging from 0
#'   (ARU failed all week) to 1 (ARU active all week). Longer periods reduce the
#'   number of surveys `J` but increase the counts per survey, trading off
#'   temporal resolution against model complexity and the closure assumption
#'   within a survey period. Default: `1L`.
#' @param thin_minutes Non-negative numeric. If supplied, observations within
#'   `thin_minutes` minutes of each other (per site and species) are thinned to
#'   a single observation, retaining the record with the highest `count`.
#'   Thinning is performed via [thin_observations()]. Default: `30`.
#' @param day_start Whether survey days start at `"midnight"` or `"midday"`.
#'    Default: `"midday"`.
#' @param occupancy_site_predictors Optional dataframe of site-level
#'   covariates for the occupancy submodel. Must contain a `deploymentID` column
#'   with the same entries as `deployments`. Predictor columns must be `numeric`
#'   (continuous), `factor` (unordered categorical), or `ordered factor`
#'   (ordinal). If multiple seasons, each `deploymentID` requires a value for
#'   each `season` it was deployed.
#' @param detection_site_predictors Optional dataframe of site-level
#'   covariates for the detection submodel. Same column-type rules as
#'   `occupancy_site_predictors`. If identical to
#'   `occupancy_site_predictors`, the same matrices are reused.
#' @param survey_predictors Optional dataframe of site-by-survey level
#'   covariates, with one row per site and date. Must contain `deploymentID`
#'   and `date` columns. Predictor columns follow the same type rules as the
#'   site-level predictor dataframes. Must cover the full deployment period for
#'   each `deploymentID`.
#' @param date <[`data-masking`][rlang::args_data_masking]> Column name for
#'   dates in `survey_predictors`. Default: `date`.
#' @param summary_functions An optional named list mapping continuous survey
#'   predictor column names to summary functions, used when aggregating survey
#'   predictors over `survey_length`-length periods. Each value can be a
#'   function name as a string (e.g. `"sum"`) or a function object (e.g. `sum`).
#'   Numeric predictors not named in `summary_functions` are summarised with
#'   `mean`; categorical and ordinal predictors are summarised with the
#'   modal value. Default: `NULL`.
#' @param scale_predictors Logical. If `TRUE`, continuous predictors are scaled
#'   to zero mean and unit variance. Survey predictors are scaled using
#'   parameters derived from site-averaged values per survey period (a `[P, J]`
#'   matrix) rather than the raw `[I, P, J]` array, so that spatial variation
#'   across sites does not inflate the scaling. Scaling parameters (means and
#'   SDs) are stored as an attribute. Default: `TRUE`.
#'
#' @return A named list of class `"occARU_data"` containing all inputs
#'   required by the occARU Stan model, except for model specification
#'   arguments which are added by [fit_model()]. The list contains:
#'   \describe{
#'     \item{`I`}{Number of sites (ARUs).}
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
#'   The object also carries the following attributes, accessible via
#'   [attr()]:
#'   \describe{
#'     \item{`sites`}{Character vector of site identifiers.}
#'     \item{`surveys`}{Character vector of start dates for each survey
#'       period).}
#'     \item{`species`}{Character vector species names.}
#'     \item{`utm_crs`}{Character. PROJ string of the UTM coordinate reference
#'       system used to transform site coordinates, or `NULL` if no coordinates
#'       were supplied.}
#'     \item{`scaling`}{tibble of means and standard deviations used to
#'       standardise continuous predictors, or `NULL` if
#'       `scale_predictors = FALSE`.}
#'     \item{`levels`}{Named list of category levels for categorical and ordinal
#'       predictors}
#'     \item{`survey_length`}{}
#'     \item{`thin_minutes`}{}
#'     \item{`reference_date`}{}
#'   }
#' @importFrom rlang :=
#'
#' @seealso [fit_model()], [thin_observations()], [find_failures()]
#'   The model is described in detail in `vignette("model", package = "occARU")`.
#' @export
make_data <- function(
  deployments,
  observations,
  failures = NULL,
  deploymentID = deploymentID,
  deploymentStart = deploymentStart,
  deploymentEnd = deploymentEnd,
  latitude = latitude,
  longitude = longitude,
  season = season,
  eventStart = eventStart,
  scientificName = scientificName,
  count = count,
  failureStart = failureStart,
  failureEnd = failureEnd,
  survey_length = 1L,
  thin_minutes = 30,
  day_start = c("midday", "midnight"),
  occupancy_site_predictors = NULL,
  detection_site_predictors = NULL,
  survey_predictors = NULL,
  date = date,
  summary_functions = NULL,
  scale_predictors = TRUE
) {
  # DEPLOYMENTS
  if (anyNA(deployments)) {
    cli::cli_abort("{.arg deployments} cannot have missing values.")
  }
  check_cols_exist(
    deployments,
    {{ deploymentID }},
    {{ deploymentStart }},
    {{ deploymentEnd }}
  )
  check_cols_class(
    deployments,
    "Date",
    {{ deploymentStart }},
    {{ deploymentEnd }}
  )
  check_dates(
    deployments,
    {{ deploymentID }},
    {{ deploymentStart }},
    {{ deploymentEnd }}
  )
  if (is.factor(deployments |> dplyr::pull({{ deploymentID }}))) {
    check_empty_levels(deployments, {{ deploymentID }})
  } else {
    deployments <- deployments |>
      dplyr::mutate({{ deploymentID }} := factor({{ deploymentID }}))
  }
  site_lvl <- deployments |> dplyr::pull({{ deploymentID }}) |> levels()
  I <- length(site_lvl)
  if (I == 1L) {
    cli::cli_abort("occARU requires more than one site.")
  }

  # seasons
  season_chr <- rlang::as_name(rlang::enquo(season))
  if (!rlang::has_name(deployments, season_chr)) {
    deployments <- dplyr::mutate(deployments, .season = factor(1L))
    check_cols_duplicates(deployments, {{ deploymentID }})
  } else if (!is.factor(dplyr::pull(deployments, {{ season }}))) {
    cli::cli_abort(
      "{.arg season} must be a factor to ensure correct ordering."
    )
  } else {
    deployments <- dplyr::rename(deployments, .season = {{ season }})
    check_empty_levels(deployments, .season)
    if (nlevels(deployments$.season) == 1) {
      cli::cli_warn(
        "{.arg season} has one level. Proceeding with single season model."
      )
    } else {
      check_deployment_times(deployments)
    }
  }
  season_lvl <- deployments |> dplyr::pull(.season) |> levels()
  K <- length(season_lvl)

  # expand deployments to individual dates
  daily_grid <- deployments |>
    dplyr::select(
      {{ deploymentID }},
      {{ deploymentStart }},
      {{ deploymentEnd }},
      .season
    ) |>
    dplyr::mutate(
      .date = purrr::map2(
        {{ deploymentStart }},
        {{ deploymentEnd }},
        ~ seq.Date(.x, .y, by = "day")
      ),
      .by = c({{ deploymentID }}, .season)
    ) |>
    tidyr::unnest(.date)

  # identify failure dates
  if (is.null(failures)) {
    daily_grid <- daily_grid |> dplyr::mutate(.Delta = 1L)
  } else {
    # checks
    if (anyNA(failures)) {
      cli::cli_abort("{.arg failures} cannot have missing values.")
    }
    check_cols_exist(
      failures,
      {{ deploymentID }},
      {{ failureStart }},
      {{ failureEnd }}
    )
    failures <- align_factor(failures, {{ deploymentID }}, site_lvl)
    check_cols_class(failures, "Date", {{ failureStart }}, {{ failureEnd }})

    failure_dates <- failures |>
      dplyr::select({{ deploymentID }}, {{ failureStart }}, {{ failureEnd }}) |>
      dplyr::mutate(
        .date = purrr::map2(
          {{ failureStart }},
          {{ failureEnd }},
          ~ seq.Date(.x, .y, by = "day")
        ),
        .failure = 1L
      ) |>
      tidyr::unnest(.date) |>
      dplyr::select({{ deploymentID }}, .date, .failure)

    daily_grid <- daily_grid |>
      dplyr::left_join(
        failure_dates,
        by = dplyr::join_by({{ deploymentID }}, .date)
      ) |>
      dplyr::mutate(.Delta = dplyr::if_else(is.na(.failure), 1L, 0L)) |>
      dplyr::select(-.failure)
  }

  # reference date
  reference_dates <- deployments |>
    dplyr::summarise(reference = min({{ deploymentStart }}), .by = .season)

  # aggregate surveys
  if (!rlang::is_integerish(survey_length) || survey_length < 1) {
    cli::cli_abort("{.arg survey_length} must be a positive integer.")
  }
  deployments_aggregated <- dplyr::left_join(
    daily_grid,
    reference_dates,
    by = ".season"
  ) |>
    dplyr::mutate(
      .survey = aggregate_by_days(.date, reference, survey_length),
      .by = .season
    ) |>
    dplyr::summarise(
      .Delta = sum(.Delta) / survey_length,
      .by = c({{ deploymentID }}, .survey, .season)
    )

  surveys <- deployments_aggregated |>
    dplyr::distinct(.season, .survey) |>
    dplyr::arrange(.season, .survey) |>
    dplyr::mutate(.survey_idx = dplyr::dense_rank(.survey), .by = .season)
  J <- surveys |>
    dplyr::count(.season) |>
    dplyr::pull(n)
  J_max <- max(J)
  if (any(J == 1L)) {
    cli::cli_abort(
      "occARU requires more than one survey. Is {.arg survey_length} too long?"
    )
  }

  # produce Delta
  Delta <- if (K > 1) {
    deployments_aggregated |>
      dplyr::left_join(
        surveys,
        by = c(".season", ".survey")
      ) |>
      dplyr::arrange(.season, .survey, {{ deploymentID }}) |>
      tidyr::complete(
        {{ deploymentID }} := factor(site_lvl, site_lvl),
        .survey_idx = 1:J_max,
        .season = factor(season_lvl, levels = season_lvl),
        fill = list(.Delta = 0)
      ) |>
      dplyr::pull(.Delta) |>
      array(c(I, J_max, K), dimnames = list(site_lvl, NULL, season_lvl))
  } else {
    deployments_aggregated |>
      dplyr::select(-.season) |>
      tidyr::pivot_wider(
        names_from = {{ deploymentID }},
        values_from = .Delta,
        values_fill = 0,
        names_sort = TRUE
      ) |>
      dplyr::arrange(.survey) |>
      tibble::column_to_rownames(".survey") |>
      t()
  }

  # OBSERVATIONS
  if (anyNA(observations)) {
    cli::cli_abort("{.arg observations} cannot have missing values.")
  }
  check_cols_exist(
    observations,
    {{ deploymentID }},
    {{ eventStart }},
    {{ scientificName }},
    {{ count }}
  )
  observations <- align_factor(observations, {{ deploymentID }}, site_lvl)
  if (!rlang::has_name(observations, season_chr)) {
    observations <- dplyr::mutate(observations, .season = factor(1L))
  } else {
    observations <- dplyr::rename(observations, .season = {{ season }})
    observations <- align_factor(observations, .season, season_lvl)
  }
  check_cols_class(observations, "POSIXt", {{ eventStart }})
  check_empty_levels(observations, {{ deploymentID }}, .season, strict = FALSE)
  if (is.factor(observations |> dplyr::pull({{ scientificName }}))) {
    check_empty_levels(observations, {{ scientificName }})
  } else {
    observations <- dplyr::mutate(
      observations,
      {{ scientificName }} := factor({{ scientificName }})
    )
  }
  species_lvl <- dplyr::pull(observations, {{ scientificName }}) |> levels()
  S <- length(species_lvl)
  observations <- dplyr::arrange(
    observations,
    {{ deploymentID }},
    {{ eventStart }}
  )

  # modify eventStart for midday
  day_start <- match.arg(day_start, c("midday", "midnight"))
  if (day_start == "midday") {
    observations <- dplyr::mutate(
      observations,
      {{ eventStart }} := {{ eventStart }} -
        lubridate::hours(12)
    )
  }

  # remove observations outside of supplied deployments
  observations <- filter_observations_window(
    deployments,
    observations,
    {{ deploymentID }},
    {{ deploymentStart }},
    {{ deploymentEnd }},
    {{ eventStart }},
    .season
  )

  # thin observations
  if (!is.null(thin_minutes)) {
    observations <- thin_observations(
      observations,
      {{ deploymentID }},
      {{ scientificName }},
      {{ eventStart }},
      {{ count }},
      thin_minutes
    )
  }

  # aggregate and fill array
  observations_aggregated <- dplyr::left_join(
    observations,
    reference_dates,
    by = ".season"
  ) |>
    dplyr::mutate(
      .date = lubridate::as_date({{ eventStart }}),
      .survey = aggregate_by_days(.date, reference, survey_length),
      .by = ".season"
    ) |>
    dplyr::semi_join(surveys, by = c(".season", ".survey")) |>
    dplyr::summarise(
      .y = sum({{ count }}),
      .by = c({{ deploymentID }}, .survey, .season, {{ scientificName }})
    )

  y <- dplyr::left_join(
    observations_aggregated,
    surveys,
    by = c(".season", ".survey")
  ) |>
    tidyr::complete(
      {{ deploymentID }} := factor(site_lvl, site_lvl),
      .survey_idx,
      .season = factor(season_lvl, season_lvl),
      {{ scientificName }} := factor(species_lvl, species_lvl),
      fill = list(.y = 0)
    ) |>
    dplyr::arrange(
      {{ scientificName }},
      .season,
      .survey_idx,
      {{ deploymentID }}
    ) |>
    dplyr::pull(.y) |>
    array(
      c(I, J_max, K, S),
      dimnames = list(site_lvl, NULL, season_lvl, species_lvl)
    )

  # SITE COORDINATES
  has_coords <- rlang::has_name(
    deployments,
    rlang::as_name(rlang::enquo(latitude))
  ) &&
    rlang::has_name(deployments, rlang::as_name(rlang::enquo(latitude)))
  if (has_coords) {
    wrong_coords <- dplyr::distinct(
      deployments,
      {{ deploymentID }},
      {{ latitude }},
      {{ longitude }}
    ) |>
      dplyr::add_count({{ deploymentID }}) |>
      dplyr::filter(n > 1) |>
      pull({{ deploymentID }})
    if (length(wrong_coords)) {
      cli::cli_abort(
        "The following {.arg deploymentID}{?s} have different {.arg latitude} \\
        or {.arg longitude} across seasons."
      )
    }
    utm <- coords_to_utm(
      dplyr::distinct(deployments, {{ deploymentID }}, .keep_all = TRUE),
      {{ deploymentID }},
      {{ latitude }},
      {{ longitude }}
    )
    XY <- utm$XY
    utm_crs <- utm$utm_crs
  } else {
    XY <- matrix(0, I, 2, dimnames = list(site_lvl, c("X", "Y")))
    utm_crs <- NA_character_
  }

  # PREDICTORS

  # occupancy site predictors
  if (!is.null(occupancy_site_predictors)) {
    if (anyNA(occupancy_site_predictors)) {
      cli::cli_abort(
        "{.arg occupancy_site_predictors} cannot have missing values."
      )
    }
    check_cols_exist(occupancy_site_predictors, {{ deploymentID }})
    occupancy_site_predictors <- align_factor(
      occupancy_site_predictors,
      {{ deploymentID }},
      site_lvl,
      strict = TRUE
    )
    if (!rlang::has_name(occupancy_site_predictors, season_chr)) {
      occupancy_site_predictors <- dplyr::mutate(
        occupancy_site_predictors,
        .season = factor(1L)
      )
      check_cols_duplicates(occupancy_site_predictors, {{ deploymentID }})
    } else {
      occupancy_site_predictors <- dplyr::rename(
        occupancy_site_predictors,
        .season = {{ season }}
      )
      occupancy_site_predictors <- align_factor(
        occupancy_site_predictors,
        .season,
        season_lvl
      )
    }
    check_site_predictors_coverage(
      occupancy_site_predictors,
      deployments,
      {{ deploymentID }},
      .season,
      verbose = FALSE
    )
    check_mixed_predictors(occupancy_site_predictors, {{ deploymentID }})
  }

  # detection site predictors
  if (
    !is.null(detection_site_predictors) &&
      !isTRUE(all.equal(occupancy_site_predictors, detection_site_predictors))
  ) {
    if (anyNA(detection_site_predictors)) {
      cli::cli_abort(
        "{.arg detection_site_predictors} cannot have missing values."
      )
    }
    check_cols_exist(detection_site_predictors, {{ deploymentID }})
    detection_site_predictors <- align_factor(
      detection_site_predictors,
      {{ deploymentID }},
      site_lvl,
      strict = TRUE
    )
    if (!rlang::has_name(detection_site_predictors, season_chr)) {
      detection_site_predictors <- dplyr::mutate(
        detection_site_predictors,
        .season = factor(1L)
      )
      check_cols_duplicates(detection_site_predictors, {{ deploymentID }})
    } else {
      detection_site_predictors <- dplyr::rename(
        detection_site_predictors,
        .season = {{ season }}
      )
      detection_site_predictors <- align_factor(
        detection_site_predictors,
        .season,
        season_lvl
      )
    }
    check_site_predictors_coverage(
      detection_site_predictors,
      deployments,
      {{ deploymentID }},
      .season,
      verbose = FALSE
    )
    check_mixed_predictors(detection_site_predictors, {{ deploymentID }})
  }

  # survey predictors
  if (!is.null(survey_predictors)) {
    if (anyNA(survey_predictors)) {
      cli::cli_abort("{.arg survey_predictors} cannot have missing values.")
    }
    check_cols_exist(survey_predictors, {{ deploymentID }}, {{ date }})
    survey_predictors <- align_factor(
      survey_predictors,
      {{ deploymentID }},
      site_lvl,
      strict = TRUE
    )
    if (!rlang::has_name(survey_predictors, season_chr)) {
      survey_predictors <- dplyr::mutate(
        survey_predictors,
        .season = factor(1L)
      )
    } else {
      survey_predictors <- dplyr::rename(
        survey_predictors,
        .season = {{ season }}
      )
      survey_predictors <- align_factor(
        survey_predictors,
        .season,
        season_lvl
      )
    }
    check_survey_predictors_coverage(
      survey_predictors,
      deployments,
      {{ deploymentID }},
      {{ deploymentStart }},
      {{ deploymentEnd }},
      .season,
      {{ date }},
      verbose = FALSE
    )
    check_mixed_predictors(survey_predictors, {{ deploymentID }}, {{ date }})
  }

  # encode predictors

  enc1 <- encode_predictors(
    occupancy_site_predictors,
    "site",
    {{ deploymentID }},
    .season,
    site_lvl = site_lvl,
    season_lvl = season_lvl,
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
      {{ deploymentID }},
      .season,
      site_lvl = site_lvl,
      season_lvl = season_lvl,
      scale_predictors = scale_predictors
    )
  }

  if (!is.null(survey_predictors)) {
    survey_predictors <- dplyr::semi_join(
      survey_predictors,
      daily_grid |> dplyr::filter(.Delta == 1L),
      by = dplyr::join_by({{ deploymentID }}, {{ date }} == .date)
    )
  }
  enc3 <- encode_predictors(
    survey_predictors,
    "survey",
    {{ deploymentID }},
    .season,
    {{ date }},
    site_lvl = site_lvl,
    season_lvl = season_lvl,
    surveys,
    reference_dates,
    summary_functions = summary_functions,
    survey_length = survey_length,
    scale_predictors = scale_predictors
  )

  # predictor dimensions
  idx <- ifelse(K == 1, 2, 3)
  P <- c(nrow(enc1$X), nrow(enc2$X), dim(enc3$X)[idx])
  P_cat <- c(nrow(enc1$X_cat), nrow(enc2$X_cat), dim(enc3$X_cat)[idx])
  P_ord <- c(nrow(enc1$X_ord), nrow(enc2$X_ord), dim(enc3$X_ord)[idx])

  # return
  structure(
    c(
      list(I = I),
      if (K == 1) {
        list(J = J_max)
      } else {
        list(J_max = J_max, K = K, J = J)
      },
      list(
        S = S,
        Delta = Delta,
        XY = XY,
        y = y[,, 1:K, 1:S],
        P = P,
        P_cat = P_cat,
        P_ord = P_ord,
        X1 = enc1$X,
        X_cat1 = if (K == 1) enc1$X_cat,
        X_ord1 = enc1$X_ord,
        X2 = enc2$X,
        X_cat2 = enc2$X_cat,
        X_ord2 = enc2$X_ord,
        X3 = enc3$X,
        X_cat3 = enc3$X_cat,
        X_ord3 = enc3$X_ord
      )
    ),
    class = "occARU_data",
    sites = site_lvl,
    surveys = surveys,
    seasons = season_lvl,
    species = species_lvl,
    utm_crs = utm_crs,
    scaling = if (scale_predictors) {
      list(
        X1 = enc1$scale_params,
        X2 = enc2$scale_params,
        X3 = enc3$scale_params
      )
    },
    levels = list(
      X1 = if (P_cat[1] || P_ord[1]) enc_levels(enc1),
      X2 = if (P_cat[2] || P_ord[2]) enc_levels(enc2),
      X3 = if (P_cat[3] || P_ord[3]) enc_levels(enc3)
    ),
    survey_length = survey_length,
    thin_minutes = thin_minutes
  ) |>
    print()
}

#' Print method for occARU_data objects
#'
#' @param x A `occARU_data` object.
#' @param ... Ignored.
#' @keywords internal
#' @export
print.occARU_data <- function(x, ...) {
  cli::cli_h1("occARU data")
  P <- sapply(1:3, \(p) x$P[p] + x$P_cat[p] + x$P_ord[p])
  K <- !is.null(x$K)
  surveys <- attr(x, "surveys")

  cli::cli_dl(c(
    "Sites (I)" = "{x$I}",
    if (K) {
      c("Seasons (K)" = "{x$K}")
    },
    "Surveys (J)" = "{paste(x$J, collapse = ', ')}"
  ))

  if (K) {
    cli::cli_text("Deployment spans:")
    cli::cli_div(theme = list(dl = list("margin-left" = 2)))
    cli::cli_dl(
      purrr::map_chr(
        levels(surveys$.season),
        ~ {
          s <- dplyr::filter(surveys, .season == .x)
          paste(
            as.character(min(s$.survey)),
            "to",
            as.character(max(s$.survey))
          )
        }
      ) |>
        purrr::set_names(paste0("  ", levels(surveys$.season)))
    )
    cli::cli_end()
  } else {
    cli::cli_text(
      "Deployment span: {as.character(min(surveys$.survey))} to \\
      {as.character(max(surveys$.survey))}"
    )
  }

  cli::cli_dl(c(
    "Species (S)" = "{x$S}",
    "Detections" = "{sum(x$y)}",
    "Site coordinates" = "{ifelse(any(x$XY != 0), 'yes', 'no')}",
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
