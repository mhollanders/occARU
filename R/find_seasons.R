#' Find seasons of deployments
#'
#' Assigns seasons such that the last `deploymentEnd` within any season strictly
#' precedes the first `deploymentStart` of the next.
#'
#' @param deployments A dataframe of deployment information. Must contain
#'   columns `locationID`, `deploymentStart`, and `deploymentEnd` (or
#'   equivalents specified via the corresponding arguments).
#' @param locationID <[`data-masking`][rlang::args_data_masking]> Column
#'   name for site identifiers. Default: `locationID`.
#' @param deploymentStart <[`data-masking`][rlang::args_data_masking]> `Date`.
#'   Column name for deployment start dates. Default: `deploymentStart`.
#' @param deploymentEnd <[`data-masking`][rlang::args_data_masking]> `Date.`
#'   Column name for deployment end dates. Default: `deploymentEnd`.
#'
#' @return The `deployments` dataframe with additional factor column `season`,
#'   labeled with `year_quarter` of the midpoint of each season, arranged by
#'   `deploymentStart`.
#' @seealso [make_data()], [find_failures()]
#' @export
find_seasons <- function(
  deployments,
  locationID = locationID,
  deploymentStart = deploymentStart,
  deploymentEnd = deploymentEnd
) {
  deployments |>
    dplyr::arrange(deploymentStart) |>
    dplyr::mutate(
      season = purrr::accumulate2(
        deploymentStart[-1],
        deploymentEnd[-n()],
        .init = 1L,
        \(season, start, end) {
          if (start > end) season + 1L else season
        }
      )
    ) |>
    dplyr::mutate(
      .mid = min(deploymentStart) +
        (max(deploymentEnd) - min(deploymentStart)) / 2,
      season = factor(paste0(
        lubridate::year(.mid),
        "_Q",
        lubridate::quarter(.mid)
      )),
      .mid = NULL,
      .by = season
    )
}
