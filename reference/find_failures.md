# Find potential ARU failure periods

Identifies potential failure periods by searching for unusually long
gaps between consecutive detections at each site. Best used on a
dataframe of all records including null records, as a suitable gap
threshold requires knowledge of expected detection rates.

## Usage

``` r
find_failures(
  df,
  locationID = locationID,
  eventStart = eventStart,
  buffer_days
)
```

## Arguments

- df:

  A dataframe of records. Must contain columns `locationID` and
  `eventStart` (or equivalents).

- locationID:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column name for site identifiers. Default: `locationID`.

- eventStart:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  `POSIXt`. Observation timestamps. Default: `eventStart`.

- buffer_days:

  Positive integer. Number of days after the last detection and before
  the next detection to exclude from the inferred failure period. For
  example, `buffer_days = 2` means the failure period starts 2 days
  after the last detection and ends 2 days before the next detection.

## Value

A dataframe with columns `locationID`, `failureStart`, and `failureEnd`,
one row per inferred failure period, with `buffer_days` stored as
attribute.

## See also

[`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md),
[`find_seasons()`](https://mhollanders.github.io/occARU/reference/find_seasons.md)
