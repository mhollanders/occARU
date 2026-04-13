# Prepare data for the occARU model

Transforms raw deployment and observation data into a named list
suitable for passing to
[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md).
Follows the [camtrapDP](https://tdwg.github.io/camtrap-dp/) data format
by default. Site coordinates are automatically projected from WGS84
latitude/longitude to UTM (km), with the zone auto-detected from the
mean longitude.

## Usage

``` r
make_data(
  deployments,
  observations,
  failures = NULL,
  deploymentID = deploymentID,
  deploymentStart = deploymentStart,
  deploymentEnd = deploymentEnd,
  latitude = latitude,
  longitude = longitude,
  eventStart = eventStart,
  scientificName = scientificName,
  count = count,
  failureStart = failureStart,
  failureEnd = failureEnd,
  survey_length = 1L,
  thin_minutes = 30,
  day_start = c("midday", "midnight"),
  reference_date = NULL,
  occupancy_site_predictors = NULL,
  detection_site_predictors = NULL,
  survey_predictors = NULL,
  date = date,
  summary_functions = NULL,
  scale_predictors = TRUE
)
```

## Arguments

- deployments:

  A dataframe of deployment information, one row per site. Must contain
  columns `deploymentID`, `deploymentStart`, and `deploymentEnd` (or
  equivalents specified via the corresponding arguments). Optionally,
  `latitude` and `longitude` columns enable the spatial Gaussian
  process.

- observations:

  A dataframe of observation records. Must contain columns
  `deploymentID`, `eventStart`, `scientificName`, and `count` (or
  equivalents specified via the corresponding arguments).

- failures:

  Optional dataframe of ARU failure periods. Must contain columns
  `deploymentID`, `failureStart`, and `failureEnd`, with each row
  corresponding to one failure period at a `deploymentID` from
  `failureStart` to `failureEnd` (inclusive). See
  [`find_failures()`](https://mhollanders.github.io/occARU/reference/find_failures.md).

- deploymentID:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column name for sites (ARUs). Retains levels if supplied as factor.
  Default: `deploymentID`.

- deploymentStart:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  `Date`. Column name for deployment start dates in `deployments`.
  Default: `deploymentStart`.

- deploymentEnd:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  `Date.` Column name for deployment end dates in `deployments`.
  Default: `deploymentEnd`.

- latitude:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  `numeric`. Column name for WGS84 latitude in `deployments`. If omitted
  alongside `longitude`, no spatial Gaussian process is fitted. Default:
  `latitude`.

- longitude:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  `numeric`. Column name for WGS84 longitude in `deployments`. Default:
  `longitude`.

- eventStart:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  `POSIXt`. Column name for observation timestamps in `observations`.
  Default: `eventStart`.

- scientificName:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column name for species names in `observations`. Retains levels if
  supplied as factor. Default: `scientificName`.

- count:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  `integerish`. Column name for number of individuals per observation
  record. Default: `count`.

- failureStart:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  `Date`. Column name for failure start dates in `failures`. Default:
  `failureStart`.

- failureEnd:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  `Date`. Column name for failure end dates (inclusive) in `failures`.
  Default: `failureEnd`.

- survey_length:

  Positive integer. Defines the length of each survey period in days.
  Observations are aggregated within each survey period by summing
  `count`, and recording effort (`Delta`) is computed as the fraction of
  the survey length the ARU was active. For example, `survey_length = 7`
  aggregates to weekly survey periods, with `Delta` ranging from 0 (ARU
  failed all week) to 1 (ARU active all week). Longer periods reduce the
  number of surveys `J` but increase the counts per survey, trading off
  temporal resolution against model complexity and the closure
  assumption within a survey period. Default: `1L`.

- thin_minutes:

  Non-negative numeric. If supplied, observations within `thin_minutes`
  minutes of each other (per site and species) are thinned to a single
  observation, retaining the record with the highest `count`. Thinning
  is performed via
  [`thin_observations()`](https://mhollanders.github.io/occARU/reference/thin_observations.md).
  Default: `30`.

- day_start:

  Whether survey days start at `"midnight"` or `"midday"`. Default:
  `"midday"`.

- reference_date:

  `Date`. Defines the start of the first survey period. If `NULL`
  (default), uses the earliest `deploymentStart`.

- occupancy_site_predictors:

  Optional dataframe of site-level covariates for the occupancy
  submodel. Must contain a `deploymentID` column with the same entries
  as `deployments`. Predictor columns must be `numeric` (continuous),
  `factor` (unordered categorical), or `ordered factor` (ordinal).

- detection_site_predictors:

  Optional dataframe of site-level covariates for the detection
  submodel. Same column-type rules as `occupancy_site_predictors`. If
  identical to `occupancy_site_predictors`, the same matrices are
  reused.

- survey_predictors:

  Optional dataframe of site-by-survey level covariates, with one row
  per site and date. Must contain `deploymentID` and `date` columns.
  Predictor columns follow the same type rules as the site-level
  predictor dataframes. Must cover the full deployment period for each
  ARU.

- date:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column name for dates in `survey_predictors`. Default: `date`.

- summary_functions:

  An optional named list mapping continuous survey predictor column
  names to summary functions, used when aggregating survey predictors
  over `survey_length`-length periods. Each value can be a function name
  as a string (e.g. `"sum"`) or a function object (e.g. `sum`). Numeric
  predictors not named in `summary_functions` are summarised with
  `mean`; categorical and ordinal predictors are summarised with the
  modal value. Default: `NULL`.

- scale_predictors:

  Logical. If `TRUE`, continuous predictors are scaled to zero mean and
  unit variance. Survey predictors are scaled using parameters derived
  from site-averaged values per survey period (a `[P, J]` matrix) rather
  than the raw `[I, P, J]` array, so that spatial variation across sites
  does not inflate the scaling. Scaling parameters (means and SDs) are
  stored as an attribute. Default: `TRUE`.

## Value

A named list of class `"occARU_data"` containing all inputs required by
the occARU Stan model, except for model specification arguments which
are added by
[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md).
The list contains:

- `I`:

  Number of sites (ARUs).

- `J`:

  Number of survey periods.

- `S`:

  Number of species.

- `Delta`:

  `[I, J]` matrix of recording effort (0-1).

- `y`:

  `[I, J, S]` array of detection counts.

- `XY`:

  `[I, 2]` matrix of UTM coordinates in km, or zeros if coordinates not
  supplied.

- `P`:

  Integer vector of length 3: number of continuous predictors for
  occupancy, and site and survey detection.

- `P_cat`:

  Integer vector of length 3: number of categorical predictors for each
  component.

- `P_ord`:

  Integer vector of length 3: number of ordinal predictors for each
  component.

- `X1`:

  `[P[1], I]` occupancy continuous design matrix.

- `X_cat1`:

  `[P_cat[1], I]` occupancy categorical integer matrix.

- `X_ord1`:

  `[P_ord[1], I]` occupancy ordinal integer matrix.

- `X2`:

  `[P[2], I]` site-level detection continuous design matrix.

- `X_cat2`:

  `[P_cat[2], I]` site-level detection categorical integer matrix.

- `X_ord2`:

  `[P_ord[2], I]` site-level detection ordinal integer matrix.

- `X3`:

  `[I, P[3], J]` site-by-survey level detection continuous array.

- `X_cat3`:

  `[I, P_cat[3], J]` site-by-survey categorical integer array.

- `X_ord3`:

  `[I, P_ord[3], J]` site-by-survey survey ordinal integer array.

The object also carries the following attributes, accessible via
[`attr()`](https://rdrr.io/r/base/attr.html):

- `sites`:

  Character vector of site identifiers.

- `surveys`:

  Character vector of start dates for each survey period).

- `species`:

  Character vector species names.

- `utm_crs`:

  Character. PROJ string of the UTM coordinate reference system used to
  transform site coordinates, or `NULL` if no coordinates were supplied.

- `scaling`:

  tibble of means and standard deviations used to standardise continuous
  predictors, or `NULL` if `scale_predictors = FALSE`.

- `levels`:

  Named list of category levels for categorical and ordinal predictors

- `survey_length`:

- `thin_minutes`:

- `reference_date`:

## See also

[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md),
[`thin_observations()`](https://mhollanders.github.io/occARU/reference/thin_observations.md),
[`find_failures()`](https://mhollanders.github.io/occARU/reference/find_failures.md)
The model is described in detail in
[`vignette("model", package = "occARU")`](https://mhollanders.github.io/occARU/articles/model.md).
