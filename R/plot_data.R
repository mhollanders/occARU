#' Plot deployment dates
#'
#' Plot deployment dates for each site across multiple seasons.
#'
#' @param data An `occARU_data` object from [make_data()].
#' @param ... Additional arguments passed to [ggplot2::geom_rect()].
#' @return A `ggplot` object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`plot_data`}{The tibble used to produce the plot.}
#'   }
#' @seealso [make_data()] [plot_observations()]
#' @export
plot_deployments <- function(data, ...) {
  if (!inherits(data, "occARU_data")) {
    cli::cli_abort(
      "{.arg data} must be an {.cls occARU_data} object from {.fun make_data}."
    )
  }
  failures <- attr(data, "failures")
  deployments <- attr(data, "deployments")
  p <- deployments |>
    ggplot2::ggplot() +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = deploymentStart,
        xmax = deploymentEnd,
        y = forcats::fct_rev(locationID),
        height = 0.8
      ),
      ...
    ) +
    ggplot2::scale_x_date() +
    ggplot2::labs(x = "Deployment", y = "Site")
  p <- if (!is.null(failures)) {
    p +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = failureStart,
          xmax = failureEnd,
          y = forcats::fct_rev(locationID),
          height = 1
        ),
        data = failures,
        fill = "white"
      )
  } else {
    p
  }
  attr(p, "plot_data") <- deployments
  p
}

#' Plot observations
#'
#' Plot observations of the data used in the occARU model, aggregated to the
#' chosen survey length and summed across sites.
#'
#' @param data An `occARU_data` object from [make_data()].
#' @param by_region `logical.` Whether to summarise observations by region, if
#'   multiple regions were included. Default: `FALSE`.
#' @param species `character`. Vector of species to plot. If `NULL` (default),
#'   all species are plotted. Must be one of `attr(occARU_data, "species")`.
#' @param regions `character`. Vector of regions to use. If `NULL` (default),
#'   all regions are used. Must be one of `attr(occARU_data, "regions")`.
#' @param ... Additional arguments passed to [ggplot2::geom_point()].
#' @return A `ggplot` object with occARU-specific attributes attached:
#'   \describe{
#'     \item{`plot_data`}{The tibble used to produce the plot.}
#'   }
#' @seealso [make_data()] [plot_deployments()]
#' @export
plot_observations <- function(
  data,
  by_region = FALSE,
  species = NULL,
  regions = NULL,
  ...
) {
  if (!inherits(data, "occARU_data")) {
    cli::cli_abort(
      "{.arg data} must be an {.cls occARU_data} object from {.fun make_data}."
    )
  }
  species_lvl <- attr(data, "species")
  species_idx <- indices(species, species_lvl)
  region_lvl <- attr(data, "regions")
  region_idx <- indices(regions, region_lvl)
  MR <- data$R > 1
  if (by_region && !MR) {
    cli::cli_warn(
      "{.arg by_region = TRUE} ignored with only one region."
    )
    by_region <- FALSE
  }
  K <- data$K
  S <- data$S

  df <- attr(data, "observations") |>
    dplyr::filter(
      as.integer(scientificName) %in% species_idx,
      as.integer(.region) %in% region_idx
    ) |>
    dplyr::summarise(
      .Y = sum(.y),
      .by = c(scientificName, if (by_region) ".region", .survey, .season)
    )

  p <- df |>
    ggplot2::ggplot(ggplot2::aes(.survey, .Y)) +
    ggplot2::geom_point(...) +
    ggplot2::labs(x = "Survey", y = "Aggregated Count")
  p <- if (by_region) {
    if (S > 1) {
      p + ggplot2::facet_grid(scientificName ~ .region, scales = "free")
    } else {
      p + ggplot2::facet_wrap(~.region, scales = "free_x")
    }
  } else if (S > 1) {
    p + ggplot2::facet_wrap(~scientificName, ncol = 1, scales = "free_y")
  } else {
    p
  }
  attr(p, "plot_data") <- df
  p
}
