# Find potential ARU failure periods

Identifies potential failure periods by searching for unusually long
gaps between consecutive detections at each site. Best used on a
dataframe of all records including null records, as a suitable gap
threshold requires knowledge of expected detection rates.

## Usage

``` r
find_failures(
  df,
  deploymentID = deploymentID,
  eventStart = eventStart,
  buffer_days
)
```

## Arguments

- df:

  A dataframe of records. Must contain columns `deploymentID` and
  `eventStart` (or equivalents).

- deploymentID:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column name for site identifiers. Default: `deploymentID`.

- eventStart:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  `POSIXt`. Observation timestamps. Default: `eventStart`.

- buffer_days:

  Positive integer. Number of days after the last detection and before
  the next detection to exclude from the inferred failure period. For
  example, `buffer_days = 2` means the failure period starts 2 days
  after the last detection and ends 2 days before the next detection.

## Value

A dataframe with columns `deploymentID`, `failureStart`, and
`failureEnd`, one row per inferred failure period.
