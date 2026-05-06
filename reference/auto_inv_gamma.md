# Automatically determine inverse gamma hyperparameters, used for Gaussian process length scales in occARU models.

Finds inverse gamma shape (alpha) and rate (beta) parameters such that
`tail_prob` prior probability mass falls below `bounds[1]` and
`tail_prob` falls above `bounds[2]`, placing `1 - 2 * tail_prob` of the
mass within the bounds.

## Usage

``` r
auto_inv_gamma(bounds, tail_prob)
```

## Arguments

- bounds:

  Numeric vector of length 2. Lower and upper bounds on the length
  scale. For a spatial GP, the bounds could be the minimum and maximum
  distance between sites

  For a periodic kernel with a fixed annual cycle, the lower bound
  should be the minimum temporal spacing between observations (e.g.
  `1/52` for weekly data) and the upper bound should be the period (e.g.
  `1` for an annual cycle).

- tail_prob:

  Numeric scalar in (0, 0.5). The probability mass to place in each
  tail. For example, `0.05` places 5% of the prior below `bounds[1]` and
  5% above `bounds[2]`.

## Value

Numeric vector of length 2, `c(alpha, beta)`, giving the shape and rate
of the fitted inverse gamma prior.
