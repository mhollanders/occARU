#' Thin observations
#'
#' Removes redundant detections by grouping consecutive records within
#' `thin_minutes` minutes of each other into clusters (per site and species),
#' retaining only the record with the highest count from each cluster.
#'
#' @param observations A dataframe of observation records. Must contain
#'   columns `deploymentID`, `eventStart`, `scientificName`, and `count`
#'   (or equivalents specified via the corresponding arguments).
#' @param deploymentID <[`data-masking`][rlang::args_data_masking]> Column
#'   name for sites (ARUs). Default: `deploymentID`.
#' @param eventStart <[`data-masking`][rlang::args_data_masking]> `POSIXt`.
#'   Column name for observation timestamps. Default: `eventStart`.
#' @param scientificName <[`data-masking`][rlang::args_data_masking]> Column
#'   name for species names. Default: `scientificName`.
#' @param count <[`data-masking`][rlang::args_data_masking]> `integerish`.Column
#'   name for number of individuals per observation record. Default: `count`.
#' @param thin_minutes Non-negative numeric. Observations within `thin_minutes`
#'   minutes of each other (per site and species) are thinned to a single
#'   observation, retaining the record with the highest `count`. Default: `30`.
#'
#' @return A dataframe of thinned observation records, sorted by `deploymentID`,
#'   `scientificName`, and `eventStart`, or the original dataframe if
#'   `thin_minutes = 0`, with `thin_minutes` stored as attribute.
#' @seealso [make_data()]
#' @export
thin_observations <- function(
  observations,
  deploymentID = deploymentID,
  scientificName = scientificName,
  eventStart = eventStart,
  count = count,
  thin_minutes = 30
) {
  if (thin_minutes < 0) {
    cli::cli_abort("{.arg thin_minutes} must be a positive scalar.")
  } else if (thin_minutes > 0) {
    check_cols_exist(
      observations,
      {{ deploymentID }},
      {{ scientificName }},
      {{ eventStart }}
    )
    check_cols_class(observations, "POSIXt", {{ eventStart }})
    n_before <- nrow(observations)
    if (thin_minutes > 0) {
      observations <- observations |>
        dplyr::arrange(
          {{ deploymentID }},
          {{ scientificName }},
          {{ eventStart }}
        ) |>
        dplyr::mutate(
          cluster = assign_clusters({{ eventStart }}, thin_minutes),
          .by = c({{ deploymentID }}, {{ scientificName }})
        ) |>
        dplyr::slice_max(
          {{ count }},
          with_ties = FALSE,
          by = c({{ deploymentID }}, {{ scientificName }}, cluster)
        ) |>
        dplyr::select(-cluster)
    } else if (thin_minutes < 0) {
      cli::cli_abort(
        "{.arg thin_minutes} must be a positive number of minutes, \\
      or 0 if no thinning is required."
      )
    }
    n_removed <- n_before - nrow(observations)
    cli::cli_inform(
      "{n_removed} event{?s} removed by {thin_minutes}-minute thinning window \\
    ({nrow(observations)} remaining)..."
    )
  }
  attr(observations, "thin_minutes") <- thin_minutes
  observations
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
