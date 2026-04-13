#' Find potential ARU failure periods
#'
#' Identifies potential failure periods by searching for unusually long gaps
#' between consecutive detections at each site. Best used on a dataframe of
#' all records including null records, as a suitable gap threshold requires
#' knowledge of expected detection rates.
#'
#' @param df A dataframe of records. Must contain columns `deploymentID` and
#'   `eventStart` (or equivalents).
#' @param deploymentID <[`data-masking`][rlang::args_data_masking]> Column
#'   name for site identifiers. Default: `deploymentID`.
#' @param eventStart <[`data-masking`][rlang::args_data_masking]> `POSIXt`.
#'   Observation timestamps. Default: `eventStart`.
#' @param buffer_days Positive integer. Number of days after the last detection
#'   and before the next detection to exclude from the inferred failure period.
#'   For example, `buffer_days = 2` means the failure period starts 2 days
#'   after the last detection and ends 2 days before the next detection.
#'
#' @return A dataframe with columns `deploymentID`, `failureStart`, and
#'   `failureEnd`, one row per inferred failure period.
#' @export
find_failures <- function(
  df,
  deploymentID = deploymentID,
  eventStart = eventStart,
  buffer_days
) {
  # checks
  check_cols_exist(
    observations,
    {{ deploymentID }},
    {{ eventStart }}
  )
  if (!rlang::is_integerish(buffer_days) || buffer_days < 0) {
    cli::cli_abort(
      "{.arg buffer_days} must be a positive integer."
    )
  }

  # produce failures
  df |>
    dplyr::arrange({{ deploymentID }}, {{ eventStart }}) |>
    dplyr::mutate(
      gap = difftime(
        dplyr::lead({{ eventStart }}),
        {{ eventStart }},
        units = "days"
      ),
      failureStart = as.Date({{ eventStart }}) + buffer_days,
      failureEnd = as.Date(dplyr::lead({{ eventStart }})) - buffer_days,
      .by = {{ deploymentID }}
    ) |>
    dplyr::filter(gap > 2 * buffer_days) |>
    dplyr::select({{ deploymentID }}, failureStart, failureEnd)
}
