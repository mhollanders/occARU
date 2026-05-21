# Find seasons of deployments

Assigns seasons such that the last `deploymentEnd` within any season
strictly precedes the first `deploymentStart` of the next.

## Usage

``` r
find_seasons(
  deployments,
  locationID = locationID,
  deploymentStart = deploymentStart,
  deploymentEnd = deploymentEnd
)
```

## Arguments

- deployments:

  A dataframe of deployment information. Must contain columns
  `locationID`, `deploymentStart`, and `deploymentEnd` (or equivalents
  specified via the corresponding arguments).

- locationID:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column name for site identifiers. Default: `locationID`.

- deploymentStart:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  `Date`. Column name for deployment start dates. Default:
  `deploymentStart`.

- deploymentEnd:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  `Date.` Column name for deployment end dates. Default:
  `deploymentEnd`.

## Value

The `deployments` dataframe with additional factor column `season`,
labeled with `year_quarter` of the midpoint of each season, arranged by
`deploymentStart`.

## See also

[`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md),
[`find_failures()`](https://mhollanders.github.io/occARU/reference/find_failures.md)
