# Set priors for the occARU model

Constructs and validates a named list of prior hyperparameters for use
with
[`occARU()`](https://mhollanders.github.io/occARU/reference/occARU.md).
Any unspecified priors use the defaults listed below.

## Usage

``` r
set_priors(
  psi_bar = c(1, 1),
  q_bar = c(1, 3),
  mu_bar = c(1, 1),
  psi_W = c(3, 0, 1),
  mu_W = c(3, 0, 2.5),
  psi_theta = c(1, 1),
  mu_theta = c(1, 1),
  iota_ell = c(1, 1),
  kappa_ell = c(1, 1),
  kappa_ell_periodic = c(1, 1),
  K_phi = c(1, 1),
  nu_ell = c(1, 1),
  phi = c(0.4, 0.3),
  alpha_O_L = 1,
  psi_beta_O_L = 1,
  mu_beta_O_L = 1,
  gamma_O_L = 1,
  iota_O_L = 1,
  kappa_O_L = 1,
  epsilon_O_L = 1,
  verbose = TRUE
)
```

## Arguments

- psi_bar:

  Numeric vector of length 2. `c(a, b)` for a Beta(a, b) prior on mean
  occupancy probability. Default: `c(1, 1)`.

- q_bar:

  Numeric vector of length 2. `c(shape, rate)` for Gamma(shape, rate)
  priors on the mean annual colonisation and emigration rates for
  multiseason models. Default: `c(1, 3)`.

- mu_bar:

  Numeric vector of length 2. `c(shape, rate)` for a Gamma(shape, rate)
  prior on mean detection rate. Default: `c(1, 1)`.

- psi_W:

  Numeric vector of length 3. `c(df, mu, sigma)` for a Student-t+(df,
  mu, sigma) prior on total occupancy log odds variance. Default:
  `c(3, 0, 1)`.

- mu_W:

  Numeric vector of length 3. `c(df, mu, sigma)` for a Student-t+(df,
  mu, sigma) prior on total log detection variance. Default:
  `c(3, 0, 2.5)`.

- psi_theta:

  Numeric vector of length 2. `c(shape, rate)` for a Gamma(shape, rate)
  prior on the occupancy variance partition sparsity parameter. Default:
  `c(1, 1)`.

- mu_theta:

  Numeric vector of length 2. `c(shape, rate)` for a Gamma(shape, rate)
  prior on the detection variance partition sparsity parameter. Default:
  `c(1, 1)`.

- iota_ell:

  Numeric vector of length 2. `c(alpha, beta)` for an InvGamma(alpha,
  beta) prior on the spatial GP length scale(s). Default: `c(1, 1)`.

- kappa_ell:

  Numeric vector of length 2. `c(alpha, beta)` for an InvGamma(alpha,
  beta) prior on the temporal GP length scale. Default: `c(1, 1)`.

- kappa_ell_periodic:

  Numeric vector of length 2. `c(alpha, beta)` for an InvGamma(alpha,
  beta) prior on the periodic temporal GP length scale. Default:
  `c(1, 1)`.

- K_phi:

  Numeric vector of length 2. `c(alpha[1], alpha[2])` for a
  Dirichlet(alpha) prior on the temporal GP variance partitions of the
  Matern and periodic kernels. Default: `c(1, 1)`.

- nu_ell:

  Numeric vector of length 2. `c(alpha, beta)` for an InvGamma(alpha,
  beta) prior on the seasonal GP length scale. Default: `c(1, 1)`.

- phi:

  Numeric vector of length 2. `c(alpha, beta)` for an InvGamma(alpha,
  beta) prior on species-specific negative binomial overdispersion
  parameters. Default: `c(0.4, 0.3)`.

- alpha_O_L:

  Positive scalar. LKJ prior on the \\\[2, 2\]\\ occupancy log odds/log
  detection rate correlation matrix. Default: `1`.

- psi_beta_O_L:

  Positive scalar. LKJ prior on the \\\[S, S\]\\ correlation matrix of
  species-specific occupancy site coefficients. Only used when site
  predictors are supplied for occupancy in
  [`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md).
  Default: `1`.

- mu_beta_O_L:

  Positive scalar. LKJ prior on the \\\[S, S\]\\ correlation matrix of
  species-specific detection site coefficients. Only used when site
  predictors are supplied for detection in
  [`make_data()`](https://mhollanders.github.io/occARU/reference/make_data.md).
  Default: `1`.

- gamma_O_L:

  Positive scalar. LKJ prior on the \\\[S, S\]\\ correlation matrix of
  species-specific survey coefficients. Only used when survey predictors
  are supplied. Default: `1`.

- iota_O_L:

  Positive scalar. LKJ prior on the \\\[S, S\]\\ correlation matrix of
  species-specific site effects. Default: `1`.

- kappa_O_L:

  Positive scalar. LKJ prior on the \\\[S, S\]\\ correlation matrix of
  species-specific survey effects. Default: `1`.

- epsilon_O_L:

  Positive scalar. LKJ prior on the \\\[S, S\]\\ correlation matrix of
  species-specific OLRE residuals. Default: `1`.

- verbose:

  Logical. If `TRUE` (default), prints list of priors.

## Value

An `occARU_priors` object (a named list) for use with
[`occARU()`](https://mhollanders.github.io/occARU/reference/occARU.md).

## See also

[`occARU()`](https://mhollanders.github.io/occARU/reference/occARU.md)
