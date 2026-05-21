# Plot observations

Plot observations of the data used in the occARU model, aggregated to
the chosen survey length and summed across sites.

## Usage

``` r
plot_observations(data, by_region = FALSE, species = NULL, regions = NULL, ...)
```

## Arguments

- data:

  An `occARU_data` object from
  [`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md).

- by_region:

  `logical.` Whether to summarise observations by region, if multiple
  regions were included. Default: `FALSE`.

- species:

  `character`. Vector of species to plot. If `NULL` (default), all
  species are plotted. Must be one of `attr(occARU_data, "species")`.

- regions:

  `character`. Vector of regions to use. If `NULL` (default), all
  regions are used. Must be one of `attr(occARU_data, "regions")`.

- ...:

  Additional arguments passed to
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html).

## Value

A `ggplot` object with occARU-specific attributes attached:

- `plot_data`:

  The tibble used to produce the plot.

## See also

[`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md)
[`plot_deployments()`](https://mhollanders.github.io/occARU/reference/plot_deployments.md)
