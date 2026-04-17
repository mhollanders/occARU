# Thin observations

Removes redundant detections by grouping consecutive records within
`thin_minutes` minutes of each other into clusters (per site and
species), retaining only the record with the highest count from each
cluster.

## Usage

``` r
thin_observations(
  observations,
  deploymentID = deploymentID,
  scientificName = scientificName,
  eventStart = eventStart,
  count = count,
  thin_minutes = 30
)
```

## Arguments

- observations:

  A dataframe of observation records. Must contain columns
  `deploymentID`, `eventStart`, `scientificName`, and `count` (or
  equivalents specified via the corresponding arguments).

- deploymentID:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column name for sites (ARUs). Default: `deploymentID`.

- scientificName:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column name for species names. Default: `scientificName`.

- eventStart:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  `POSIXt`. Column name for observation timestamps. Default:
  `eventStart`.

- count:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  `integerish`.Column name for number of individuals per observation
  record. Default: `count`.

- thin_minutes:

  Non-negative numeric. Observations within `thin_minutes` minutes of
  each other (per site and species) are thinned to a single observation,
  retaining the record with the highest `count`. Default: `30`.

## Value

A dataframe of thinned observation records, sorted by `deploymentID`,
`scientificName`, and `eventStart`, or the original dataframe if
`thin_minutes = 0`, with `thin_minutes` stored as attribute.

## See also

[`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md)
