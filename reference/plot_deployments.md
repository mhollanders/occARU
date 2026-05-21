# Plot deployment dates

Plot deployment dates for each site across multiple seasons.

## Usage

``` r
plot_deployments(data, ...)
```

## Arguments

- data:

  An `occARU_data` object from
  [`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md).

- ...:

  Additional arguments passed to
  [`ggplot2::geom_rect()`](https://ggplot2.tidyverse.org/reference/geom_tile.html).

## Value

A `ggplot` object with occARU-specific attributes attached:

- `plot_data`:

  The tibble used to produce the plot.

## See also

[`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md)
[`plot_observations()`](https://mhollanders.github.io/occARU/reference/plot_observations.md)
