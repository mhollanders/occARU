#' Prepare data for the occARU model
#'
#' Transforms raw deployment and observation data into a named list suitable
#' for passing to [occARU()]. Follows the
#' [camtrapDP](https://tdwg.github.io/camtrap-dp/) data format by default.
#' Site coordinates are automatically projected from WGS84 latitude/longitude
#' to UTM (km), with the zone auto-detected from the mean longitude.
#'
#' @param deployments A dataframe of deployment information, one row per
#'   site (and potentially season). Must contain columns `locationID`,
#'   `deploymentStart`, and `deploymentEnd` (or equivalents specified via the
#'   corresponding arguments). Optionally, `latitude` and `longitude` columns
#'   enable the spatial Gaussian process. If multiple seasons, must also contain
#'   column `season`.
#' @param observations A dataframe of observation records. Must contain
#'   columns `locationID`, `eventStart`, `scientificName`, and `count`
#'   (or equivalents specified via the corresponding arguments). If multiple
#'   seasons, must also contain column `season`.
#' @param failures Optional dataframe of ARU failure periods. Must contain
#'   columns `locationID`, `failureStart`, and `failureEnd`, with each row
#'   corresponding to one failure period at a `locationID` from `failureStart`
#'   to `failureEnd` (inclusive). See [find_failures()].
#' @param locationID <[`data-masking`][rlang::args_data_masking]> Column
#'   name for sites (ARUs). Retains levels if supplied as factor.  Default:
#'   `locationID`.
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
#' @param region <[`data-masking`][rlang::args_data_masking]> Optional column in
#'   `deployments` specifying region, defined as a cluster of ARUs. Leads to
#'   faster model fits when spatial site effects are included in [occARU()].
#'   If the column is not present in `deployments`, all observations are treated
#'   as a single region. Default: `region`.
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
#'   covariates for the occupancy submodel. Must contain a `locationID` column
#'   with the same entries as `deployments`. Predictor columns must be `numeric`
#'   (continuous), `factor` (unordered categorical), or `ordered factor`
#'   (ordinal). If multiple seasons, each `locationID` requires a value for
#'   each `season` it was deployed.
#' @param detection_site_predictors Optional dataframe of site-level
#'   covariates for the detection submodel. Same column-type rules as
#'   `occupancy_site_predictors`. If identical to
#'   `occupancy_site_predictors`, the same matrices are reused.
#' @param survey_predictors Optional dataframe of site-by-survey level
#'   covariates, with one row per site and date. Must contain `locationID`
#'   and `date` columns. Predictor columns follow the same type rules as the
#'   site-level predictor dataframes. Must cover the full deployment period for
#'   each `locationID`.
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
#' @param verbose Logical. If `TRUE` (default), prints data.
#'
#' @return A named list of class `"occARU_data"` containing all inputs
#'   required by the occARU Stan model, except for model specification
#'   arguments which are added by [occARU()]. The list contains:
#'   \describe{
#'     \item{`I`}{Number of sites (ARUs).}
#'     \item{`R`}{Number of regions (groups of sites).}
#'     \item{`J`}{Number of survey periods.}
#'     \item{`K`}{Number of seasons (if multiseason).}
#'     \item{`tau`}{Interval length in years between end of previous deploymment
#'       and start of current deployment (if multiseason).}
#'     \item{`dyn`}{Indicator for dynamic occupancy, when at least one site was
#'       deployed over multiple seasons.}
#'     \item{`S`}{Number of species.}
#'     \item{`Delta`}{`[I, J(, K)]` array of recording effort (0-1).}
#'     \item{`y`}{`[I, J(, K), S]` array of detection counts.}
#'     \item{`XY`}{`[I, 2]` matrix of UTM coordinates in km, or zeros if
#'       coordinates not supplied.}
#'     \item{`P`}{Integer vector of length 3: number of continuous predictors
#'       for occupancy, and site and survey detection.}
#'     \item{`P_cat`}{Integer vector of length 3: number of categorical
#'       predictors for each component.}
#'     \item{`P_ord`}{Integer vector of length 3: number of ordinal predictors
#'       for each component.}
#'     \item{`X1`}{`[I(, K), P[1]]` occupancy continuous design array.}
#'     \item{`X_cat1`}{`[I(, K), P_cat[1]]` occupancy categorical integer
#'       array.}
#'     \item{`X_ord1`}{`[I(, K), P_ord[1]]` occupancy ordinal integer array.}
#'     \item{`X2`}{`[I(, K), P[2]]` site-level detection continuous design
#'       array.}
#'     \item{`X_cat2`}{`[I(, K), P_cat[2]]` site-level detection categorical
#'       integer array.}
#'     \item{`X_ord2`}{`[I(, K), P_ord[2]]` site-level detection ordinal integer
#'       array.}
#'     \item{`X3`}{`[I(, K), J, P[3]]` site-by-survey level detection continuous
#'       array.}
#'     \item{`X_cat3`}{`[I(, K), J, P_cat[3]]` site-by-survey categorical
#'       integer array.}
#'     \item{`X_ord3`}{`[I(, K), J, P_ord[3]]` site-by-survey survey ordinal
#'       integer array.}
#'   }
#'   The object also carries the following attributes, accessible via
#'   [attr()]:
#'   \describe{
#'     \item{`sites`}{Character vector of site identifiers.}
#'     \item{`surveys`}{tibble of start dates and indices for each survey
#'       period per season.}
#'     \item{`seasons`}{Character vector of season identifiers.}
#'     \item{`regions`}{Character vector of region identifiers.}
#'     \item{`species`}{Character vector species names.}
#'     \item{`utm_crs`}{Character. PROJ string of the UTM coordinate reference
#'       system used to transform site coordinates, or `NULL` if no coordinates
#'       were supplied.}
#'     \item{`scaling`}{tibble of means and standard deviations used to
#'       standardise continuous predictors, or `NULL` if
#'       `scale_predictors = FALSE`.}
#'     \item{`levels`}{Named list of category levels for categorical and ordinal
#'       predictors.}
#'     \item{`survey_length`}{}
#'     \item{`thin_minutes`}{}
#'     \item{`reference_dates`}{}
#'     \item{`day_start`}{}
#'   }
#' @importFrom rlang :=
#'
#' @seealso [occARU()], [thin_observations()], [find_failures()]
#'   The model is described in detail in
#'   `vignette("model", package = "occARU")`.
#' @export
make_data <- function(
  deployments,
  observations,
  failures = NULL,
  locationID = locationID,
  deploymentStart = deploymentStart,
  deploymentEnd = deploymentEnd,
  latitude = latitude,
  longitude = longitude,
  region = region,
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
  scale_predictors = TRUE,
  verbose = TRUE
) {
  # deployments
  deployments <- make_deployments(
    deployments,
    {{ locationID }},
    {{ deploymentStart }},
    {{ deploymentEnd }},
    {{ region }},
    {{ season }}
  )

  # levels and dimensions
  site_lvl <- deployments |> dplyr::pull({{ locationID }}) |> levels()
  I <- length(site_lvl)
  if (I == 1L) {
    cli::cli_abort("occARU requires more than one site (ARU).")
  }
  region_lvl <- deployments |> dplyr::pull(.region) |> levels()
  R <- length(region_lvl)
  season_lvl <- deployments |> dplyr::pull(.season) |> levels()
  K <- length(season_lvl)

  # reference dates
  reference_dates <- deployments |>
    dplyr::summarise(reference = min({{ deploymentStart }}), .by = .season)

  # expand deployments to dates and incorporate failures
  daily_grid <- make_daily_grid(
    deployments,
    failures,
    {{ locationID }},
    {{ deploymentStart }},
    {{ deploymentEnd }},
    .season,
    {{ failureStart }},
    {{ failureEnd }}
  )

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
      .by = c({{ locationID }}, .survey, .season)
    )

  # surveys
  surveys <- deployments_aggregated |>
    dplyr::distinct(.season, .survey) |>
    dplyr::arrange(.season, .survey) |>
    dplyr::mutate(
      .survey_idx = factor(dplyr::dense_rank(.survey)),
      .by = .season
    )
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
  Delta <- deployments_aggregated |>
    dplyr::left_join(
      surveys,
      by = c(".season", ".survey")
    ) |>
    tidyr::complete(
      .season,
      .survey_idx,
      {{ locationID }},
      fill = list(.Delta = 0)
    ) |>
    dplyr::pull(.Delta) |>
    array(
      c(I, J_max, K),
      dimnames = list(
        site_lvl,
        if (K == 1) as.character(surveys$.survey) else NULL,
        season_lvl
      )
    )

  # produce tau
  if (K > 1) {
    tau <- deployments |>
      dplyr::arrange({{ locationID }}, .season) |>
      dplyr::mutate(
        tau = difftime(
          {{ deploymentStart }},
          dplyr::lag({{ deploymentEnd }}),
          units = "weeks"
        ) |>
          as.numeric(),
        .by = {{ locationID }}
      ) |>
      tidyr::drop_na() |>
      tidyr::complete(
        {{ locationID }},
        .season = season_lvl[-1],
        fill = list(tau = 0)
      ) |>
      dplyr::pull(tau) |>
      matrix(K - 1, I, dimnames = list(season_lvl[-1], site_lvl))
    dyn <- any(colSums(tau) > 0)
  }

  # observations
  observations <- make_observations(
    observations,
    daily_grid,
    {{ locationID }},
    {{ eventStart }},
    {{ scientificName }},
    {{ count }},
    {{ season }}
  )
  species_lvl <- dplyr::pull(observations, {{ scientificName }}) |> levels()
  S <- length(species_lvl)

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
    {{ locationID }},
    {{ deploymentStart }},
    {{ deploymentEnd }},
    {{ eventStart }},
    .season
  )

  # thin observations
  observations <- thin_observations(
    observations,
    {{ locationID }},
    {{ scientificName }},
    {{ eventStart }},
    {{ count }},
    thin_minutes
  )
  species_lvl <- dplyr::pull(observations, {{ scientificName }}) |> levels()
  S <- length(species_lvl)

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
    dplyr::summarise(
      .y = sum({{ count }}),
      .by = c({{ locationID }}, .survey, .season, {{ scientificName }})
    )

  # make detection history
  y <- dplyr::left_join(
    observations_aggregated,
    surveys,
    by = c(".season", ".survey")
  ) |>
    tidyr::complete(
      {{ scientificName }},
      .season,
      .survey_idx,
      {{ locationID }},
      fill = list(.y = 0)
    ) |>
    dplyr::pull(.y) |>
    array(
      c(I, J_max, K, S),
      dimnames = list(
        site_lvl,
        if (K == 1) as.character(surveys$.survey) else NULL,
        season_lvl,
        species_lvl
      )
    )

  # site_coordinates
  XY <- make_coordinates(
    deployments,
    {{ locationID }},
    {{ latitude }},
    {{ longitude }}
  )

  # predictors
  occupancy_site_predictors <- make_site_predictors(
    occupancy_site_predictors,
    deployments,
    {{ locationID }},
    .season
  )
  enc1 <- encode_predictors(
    occupancy_site_predictors,
    "site",
    {{ locationID }},
    .season,
    site_lvl = site_lvl,
    season_lvl = season_lvl,
    scale_predictors = scale_predictors
  )

  if (isTRUE(all.equal(occupancy_site_predictors, detection_site_predictors))) {
    enc2 <- enc1
  } else {
    detection_site_predictors <- make_site_predictors(
      detection_site_predictors,
      deployments,
      {{ locationID }},
      .season
    )
    enc2 <- encode_predictors(
      detection_site_predictors,
      "site",
      {{ locationID }},
      .season,
      site_lvl = site_lvl,
      season_lvl = season_lvl,
      scale_predictors = scale_predictors
    )
  }

  survey_predictors <- make_survey_predictors(
    survey_predictors,
    deployments,
    {{ locationID }},
    {{ season }},
    {{ date }},
    {{ deploymentStart }},
    {{ deploymentEnd }}
  )
  if (!is.null(survey_predictors)) {
    survey_predictors <- dplyr::semi_join(
      survey_predictors,
      daily_grid |> dplyr::filter(.Delta == 1L),
      by = dplyr::join_by({{ locationID }}, {{ date }} == .date)
    )
  }
  enc3 <- encode_predictors(
    survey_predictors,
    "survey",
    {{ locationID }},
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
  P <- c(dim(enc1$X)[idx], dim(enc2$X)[idx], dim(enc3$X)[idx + 1])
  P_cat <- c(
    dim(enc1$X_cat)[idx],
    dim(enc2$X_cat)[idx],
    dim(enc3$X_cat)[idx + 1]
  )
  P_ord <- c(
    dim(enc1$X_ord)[idx],
    dim(enc2$X_ord)[idx],
    dim(enc3$X_ord)[idx + 1]
  )

  # return
  data <- structure(
    c(
      list(I = I, J = J_max, R = R, K = K),
      if (K > 1) {
        list(tau = tau / 52, dyn = dyn, J_i = J)
      },
      list(
        S = S,
        Delta = Delta[,, 1:K],
        region = as.integer(deployments$.region),
        XY = XY,
        y = y, # y[,, 1:K, ],
        P = P,
        P_cat = P_cat,
        P_ord = P_ord,
        X1 = enc1$X,
        X_cat1 = enc1$X_cat,
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
    regions = region_lvl,
    species = species_lvl,
    utm_crs = attr(XY, "utm_crs"),
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
    thin_minutes = thin_minutes,
    day_start = day_start
  )
  if (verbose) {
    print(data)
  }
  invisible(data)
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
  surveys <- attr(x, "surveys")

  cli::cli_dl(c(
    "Sites (I)" = "{x$I}",
    if (x$K > 1) {
      c(
        "Seasons (K)" = "{x$K}",
        "Surveys (J)" = "{paste(x$J_i, collapse = ', ')}"
      )
    } else {
      c("Surveys (J)" = "{x$J}")
    }
  ))

  if (x$K > 1) {
    cli::cli_text("Deployment spans:")
    cli::cli_div(theme = list(dl = list("margin-left" = 2)))
    cli::cli_dl(
      purrr::map_chr(
        levels(surveys$.season),
        ~ {
          s <- dplyr::filter(surveys, .season == .)
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
      "Deployment span: {as.character(min(surveys$.survey))} to
      {as.character(max(surveys$.survey))}"
    )
  }

  cli::cli_dl(c(
    if (x$R > 1) {
      c("Regions (R)" = "{x$R}")
    },
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
}


#' Prepare deployments dataframe
#' @noRd
make_deployments <- function(
  deployments,
  locationID,
  deploymentStart,
  deploymentEnd,
  region,
  season
) {
  # checks
  check_cols_exist(
    deployments,
    {{ locationID }},
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
    {{ locationID }},
    {{ deploymentStart }},
    {{ deploymentEnd }}
  )

  # regions
  if (!rlang::has_name(deployments, rlang::as_name(rlang::enquo(region)))) {
    deployments <- dplyr::mutate(deployments, .region = factor(1L))
  } else {
    deployments <- dplyr::rename(deployments, .region = {{ region }})
    if (!is.factor(dplyr::pull(deployments, .region))) {
      deployments <- dplyr::mutate(deployments, .region = factor(.region))
    }
  }

  # seasons
  if (!rlang::has_name(deployments, rlang::as_name(rlang::enquo(season)))) {
    deployments <- dplyr::mutate(deployments, .season = factor(1L))
    check_cols_duplicates(deployments, {{ locationID }})
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

  # sites
  if (is.factor(deployments |> dplyr::pull({{ locationID }}))) {
    check_empty_levels(deployments, {{ locationID }})
  } else {
    deployments <- dplyr::mutate(
      deployments,
      {{ locationID }} := factor({{ locationID }})
    )
  }
  deployments <- dplyr::mutate(
    deployments,
    {{ locationID }} := forcats::fct_reorder(
      {{ locationID }},
      as.integer(.region)
    )
  )
  check_missing(deployments)
  deployments |>
    dplyr::arrange({{ locationID }}, .season)
}

#' Produce daily grid with failure dates
#' @noRd
make_daily_grid <- function(
  deployments,
  failures = NULL,
  locationID,
  deploymentStart,
  deploymentEnd,
  season,
  failureStart,
  failureEnd
) {
  # expand deployment dates
  daily_grid <- deployments |>
    dplyr::select(
      {{ locationID }},
      {{ deploymentStart }},
      {{ deploymentEnd }},
      {{ season }}
    ) |>
    dplyr::mutate(
      .date = purrr::map2(
        {{ deploymentStart }},
        {{ deploymentEnd }},
        ~ seq.Date(.x, .y, by = "day")
      ),
      .by = c({{ locationID }}, {{ season }})
    ) |>
    tidyr::unnest(.date) |>
    dplyr::mutate(.Delta = 1L)

  # join failures
  if (!is.null(failures)) {
    check_missing(failures)
    check_cols_exist(
      failures,
      {{ locationID }},
      {{ failureStart }},
      {{ failureEnd }}
    )
    failures <- align_factor(failures, {{ locationID }}, site_lvl)
    check_cols_class(failures, "Date", {{ failureStart }}, {{ failureEnd }})
    check_dates(
      failures,
      {{ locationID }},
      {{ failureStart }},
      {{ failureEnd }}
    )

    failure_dates <- failures |>
      dplyr::select({{ locationID }}, {{ failureStart }}, {{ failureEnd }}) |>
      dplyr::mutate(
        .date = purrr::map2(
          {{ failureStart }},
          {{ failureEnd }},
          ~ seq.Date(.x, .y, by = "day")
        ),
        .failure = 1L
      ) |>
      tidyr::unnest(.date) |>
      dplyr::select({{ locationID }}, .date, .failure)

    daily_grid <- dplyr::left_join(
      daily_grid,
      failure_dates,
      by = dplyr::join_by({{ locationID }}, .date)
    ) |>
      dplyr::mutate(
        .Delta = dplyr::replace_when(.Delta, .failure == 1L ~ 0L)
      ) |>
      dplyr::select(-.failure)
  }
  daily_grid
}

#' Prepare observations dataframe
#' @noRd
make_observations <- function(
  observations,
  daily_grid,
  locationID,
  eventStart,
  scientificName,
  count,
  season
) {
  # checks
  check_cols_exist(
    observations,
    {{ locationID }},
    {{ eventStart }},
    {{ scientificName }},
    {{ count }}
  )
  site_lvl <- dplyr::pull(daily_grid, {{ locationID }}) |> levels()
  observations <- align_factor(observations, {{ locationID }}, site_lvl)
  season_lvl <- dplyr::pull(daily_grid, .season) |> levels()
  if (length(season_lvl) == 1) {
    observations <- dplyr::mutate(
      observations,
      .season = factor(1L, labels = season_lvl)
    )
  } else {
    check_cols_exist(observations, {{ season }})
    observations <- dplyr::rename(observations, .season = {{ season }})
    observations <- align_factor(observations, .season, season_lvl)
  }
  check_cols_class(observations, "POSIXt", {{ eventStart }})
  check_empty_levels(observations, {{ locationID }}, .season, strict = FALSE)
  if (is.factor(observations |> dplyr::pull({{ scientificName }}))) {
    check_empty_levels(observations, {{ scientificName }})
  } else {
    observations <- dplyr::mutate(
      observations,
      {{ scientificName }} := factor({{ scientificName }})
    )
  }
  observations <- dplyr::select(
    observations,
    {{ locationID }},
    {{ eventStart }},
    {{ scientificName }},
    {{ count }},
    .season
  ) |>
    dplyr::arrange(
      {{ locationID }},
      {{ eventStart }}
    )
  check_missing(observations)
  observations
}

#' Make coordinates matrix
#' @noRd
make_coordinates <- function(
  deployments,
  locationID,
  latitude,
  longitude
) {
  has_coords <- rlang::has_name(
    deployments,
    rlang::as_name(rlang::enquo(latitude))
  ) &&
    rlang::has_name(deployments, rlang::as_name(rlang::enquo(longitude)))
  if (has_coords) {
    wrong_coords <- deployments |>
      dplyr::distinct({{ locationID }}, {{ latitude }}, {{ longitude }}) |>
      dplyr::add_count({{ locationID }}) |>
      dplyr::filter(n > 1) |>
      dplyr::pull({{ locationID }})
    n_wrong <- length(wrong_coords)
    if (n_wrong) {
      cli::cli_abort(
        "The following {n_wrong} {.arg locationID} {?has/have} different
        {.arg latitude} or {.arg longitude} across seasons:
        {.val {wrong_coords}}."
      )
    }
    deployments <- deployments |>
      dplyr::distinct({{ locationID }}, {{ latitude }}, {{ longitude }}) |>
      coords_to_utm({{ latitude }}, {{ longitude }})

    XY <- tibble::column_to_rownames(
      deployments |>
        dplyr::select(-c({{ latitude }}, {{ longitude }})) |>
        dplyr::arrange({{ locationID }}),
      rlang::as_name(rlang::enquo(locationID))
    )
    utm_crs <- attr(deployments, "utm_crs")
  } else {
    XY <- matrix(0, I, 2, dimnames = list(site_lvl, c("X", "Y")))
    utm_crs <- NA_character_
  }
  attr(XY, "utm_crs") <- utm_crs
  XY
}

#' Prepare site predictors dataframe
#' @noRd
make_site_predictors <- function(predictors, deployments, locationID, season) {
  # checks
  if (!is.null(predictors)) {
    check_missing(predictors)
    check_cols_exist(predictors, {{ locationID }})
    site_lvl <- dplyr::pull(deployments, {{ locationID }}) |> levels()
    predictors <- align_factor(
      predictors |> dplyr::filter({{ locationID }} %in% site_lvl),
      {{ locationID }},
      site_lvl,
      strict = TRUE
    )
    season_lvl <- dplyr::pull(deployments, .season) |> levels()
    if (length(season_lvl == 1)) {
      predictors <- dplyr::mutate(
        predictors,
        .season = factor(1L, labels = season_lvl)
      )
      check_cols_duplicates(predictors, {{ locationID }})
    } else {
      check_cols_exist(predictors, {{ season }})
      predictors <- dplyr::rename(
        predictors,
        .season = {{ season }}
      )
      predictors <- align_factor(
        predictors,
        .season,
        season_lvl
      )
    }
    check_site_predictors_coverage(
      predictors,
      deployments,
      {{ locationID }},
      .season,
      verbose = FALSE
    )
    check_mixed_predictors(predictors, {{ locationID }})
  }
  predictors
}

#' Prepare survey predictors dataframe
#' @noRd
make_survey_predictors <- function(
  predictors,
  deployments,
  locationID,
  season,
  date,
  deploymentStart,
  deploymentEnd
) {
  if (!is.null(predictors)) {
    check_missing(predictors)
    check_cols_exist(predictors, {{ locationID }}, {{ date }})
    site_lvl <- dplyr::pull(deployments, {{ locationID }}) |> levels()
    predictors <- align_factor(
      predictors |> dplyr::filter({{ locationID }} %in% site_lvl),
      {{ locationID }},
      site_lvl,
      strict = TRUE
    )
    season_lvl <- dplyr::pull(deployments, .season) |> levels()
    if (length(season_lvl)) {
      predictors <- dplyr::mutate(
        predictors,
        .season = factor(1L, labels = season_lvl)
      )
    } else {
      check_cols_exist(predictors, {{ season }})
      predictors <- dplyr::rename(predictors, .season = {{ season }})
      predictors <- align_factor(
        predictors,
        .season,
        season_lvl
      )
    }
    check_survey_predictors_coverage(
      predictors,
      deployments,
      {{ locationID }},
      {{ deploymentStart }},
      {{ deploymentEnd }},
      .season,
      {{ date }},
      verbose = FALSE
    )
    check_mixed_predictors(predictors, {{ locationID }}, {{ date }})
  }
  predictors
}
