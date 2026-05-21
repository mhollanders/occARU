# Prepare data for the occARU model

Transforms raw deployment and observation data into a named list
suitable for passing to
[`occARU()`](https://mhollanders.github.io/occARU/reference/occARU.md).
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
  locationID = locationID,
  deploymentStart = deploymentStart,
  deploymentEnd = deploymentEnd,
  latitude = latitude,
  longitude = longitude,
  region = region,
  season = season,
  eventStart = eventStart,
  scientificName = scientificName,
  count = count,
  failureStart = failureStart,
  failureEnd = failureEnd,
  survey_length = 1L,
  thin_minutes = 30,
  day_start = c("midday", "midnight"),
  occupancy_site_predictors = NULL,
  detection_site_predictors = NULL,
  survey_predictors = NULL,
  date = date,
  summary_functions = NULL,
  scale_predictors = TRUE,
  verbose = TRUE
)
```

## Arguments

- deployments:

  A dataframe of deployment information, one row per site (and
  potentially season). Must contain columns `locationID`,
  `deploymentStart`, and `deploymentEnd` (or equivalents specified via
  the corresponding arguments). Optionally, `latitude` and `longitude`
  columns enable the spatial Gaussian process. If multiple seasons, must
  also contain column `season`.

- observations:

  A dataframe of observation records. Must contain columns `locationID`,
  `eventStart`, `scientificName`, and `count` (or equivalents specified
  via the corresponding arguments). If multiple seasons, must also
  contain column `season`.

- failures:

  Optional dataframe of ARU failure periods. Must contain columns
  `locationID`, `failureStart`, and `failureEnd`, with each row
  corresponding to one failure period at a `locationID` from
  `failureStart` to `failureEnd` (inclusive). See
  [`find_failures()`](https://mhollanders.github.io/occARU/reference/find_failures.md).

- locationID:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column name for sites (ARUs). Retains levels if supplied as factor.
  Default: `locationID`.

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

- region:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Optional column in `deployments` specifying region, defined as a
  cluster of ARUs. Leads to faster model fits when spatial site effects
  are included in
  [`occARU()`](https://mhollanders.github.io/occARU/reference/occARU.md).
  If the column is not present in `deployments`, all observations are
  treated as a single region. Default: `region`.

- season:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Optional column specifying season in `deployments`. The column must be
  a factor to ensure correct ordering. If the column is not present, a
  single season is assumed. Default: `season`. See
  [`find_seasons()`](https://mhollanders.github.io/occARU/reference/find_seasons.md).

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

- occupancy_site_predictors:

  Optional dataframe of site-level covariates for the occupancy
  submodel. Must contain a `locationID` column with the same entries as
  `deployments`. Predictor columns must be `numeric` (continuous),
  `factor` (unordered categorical), or `ordered factor` (ordinal). If
  multiple seasons, each `locationID` requires a value for each `season`
  it was deployed.

- detection_site_predictors:

  Optional dataframe of site-level covariates for the detection
  submodel. Same column-type rules as `occupancy_site_predictors`. If
  identical to `occupancy_site_predictors`, the same matrices are
  reused.

- survey_predictors:

  Optional dataframe of site-by-survey level covariates, with one row
  per site and date. Must contain `locationID` and `date` columns.
  Predictor columns follow the same type rules as the site-level
  predictor dataframes. Must cover the full deployment period for each
  `locationID`.

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

- verbose:

  Logical. If `TRUE` (default), prints data.

## Value

A named list of class `"occARU_data"` containing all inputs required by
the occARU Stan model, except for model specification arguments which
are added by
[`occARU()`](https://mhollanders.github.io/occARU/reference/occARU.md).
The list contains:

- `I`:

  Number of sites (ARUs).

- `R`:

  Number of regions (groups of sites).

- `J`:

  Number of survey periods (maximum).

- `K`:

  Number of seasons.

- `S`:

  Number of species.

- `tau`:

  Interval length in years between end of previous deploymment and start
  of current deployment (if multiseason).

- `dyn`:

  Indicator for dynamic occupancy, when at least one site was deployed
  over multiple seasons.

- `Delta`:

  `[K, J, I]` array of recording effort (0-1).

- `y`:

  `[K, I, J, S]` array of detection counts.

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

  `[K, I, P[1]]` occupancy continuous design array.

- `X_cat1`:

  `[K, I, P_cat[1]]` occupancy categorical integer array.

- `X_ord1`:

  `[K, I, P_ord[1]]` occupancy ordinal integer array.

- `X2`:

  `[K, I, P[2]]` site-level detection continuous design array.

- `X_cat2`:

  `[K, I, P_cat[2]]` site-level detection categorical integer array.

- `X_ord2`:

  `[K, I, P_ord[2]]` site-level detection ordinal integer array.

- `X3`:

  `[K, I, J, P[3]]` site-by-survey level detection continuous array.

- `X_cat3`:

  `[K, I, J, P_cat[3]]` site-by-survey categorical integer array.

- `X_ord3`:

  `[K, I, J, P_ord[3]]` site-by-survey survey ordinal integer array.

The object also carries the following attributes, accessible via
[`attr()`](https://rdrr.io/r/base/attr.html):

- `deployments`:

  The processed `deployments`.

- `observations`:

  The processed and thinned `observations` aggregated to the chosen
  survey length.

- `sites`:

  Character vector of site identifiers.

- `regions`:

  Character vector of region identifiers.

- `surveys`:

  tibble of start dates and indices for each survey period per season.

- `seasons`:

  Character vector of season identifiers.

- `species`:

  Character vector species names.

- `utm_crs`:

  Character. PROJ string of the UTM coordinate reference system used to
  transform site coordinates, or `NULL` if no coordinates were supplied.

- `scaling`:

  tibble of means and standard deviations used to standardise continuous
  predictors, or `NULL` if `scale_predictors = FALSE`.

- `levels`:

  Named list of category levels for categorical and ordinal predictors.

- `reference_dates`:

  First deploymentStart per season.

- `survey_length`:

- `thin_minutes`:

- `day_start`:

## See also

[`occARU()`](https://mhollanders.github.io/occARU/reference/occARU.md),
[`plot_deployments()`](https://mhollanders.github.io/occARU/reference/plot_deployments.md),
[`plot_observations()`](https://mhollanders.github.io/occARU/reference/plot_observations.md),
[`thin_observations()`](https://mhollanders.github.io/occARU/reference/thin_observations.md),
[`find_failures()`](https://mhollanders.github.io/occARU/reference/find_failures.md),
[`find_seasons()`](https://mhollanders.github.io/occARU/reference/find_seasons.md)
The model is described in detail in
[`vignette("model")`](https://mhollanders.github.io/occARU/articles/model.md).
