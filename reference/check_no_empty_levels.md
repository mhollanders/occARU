# Check that all factor levels are present in the data

Aborts if a factor column contains levels with no corresponding rows,
which typically occurs when a data frame has been filtered after the
factor was defined.

## Usage

``` r
check_no_empty_levels(df, col_chr)
```

## Arguments

- df:

  A data frame.

- col_chr:

  Column name as a string.

## Value

`df`, invisibly, if all levels are present.
