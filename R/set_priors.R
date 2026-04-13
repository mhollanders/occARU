#' Set priors for an occARU model
#'
#' Constructs and validates a named list of prior hyperparameters for use with
#' [fit_model()]. Any unspecified priors use the defaults listed below.
#'
#' @param psi_bar Numeric vector of length 2. `c(a, b)` for a Beta(a, b) prior
#'   on mean occupancy probability. Default: `c(1, 1)`.
#' @param mu_bar Numeric vector of length 2. `c(shape, rate)` for a
#'   Gamma(shape, rate) prior on mean detection rate. Default: `c(1, 1)`.
#' @param psi_W Numeric vector of length 3. `c(df, mu, sigma)` for a
#'   Student-t+(df, mu, sigma) prior on total occupancy log odds variance.
#'   Default: `c(3, 0, 1)`.
#' @param mu_W Numeric vector of length 3. `c(df, mu, sigma)` for a
#'   Student-t+(df, mu, sigma) prior on total log detection variance.
#'   Default: `c(3, 0, 2.5)`.
#' @param psi_theta Numeric vector of length 2. `c(shape, rate)` for a
#'   Gamma(shape, rate) prior on the occupancy variance partition sparsity
#'   parameter. Default: `c(1, 1)`.
#' @param mu_theta Numeric vector of length 2. `c(shape, rate)` for a
#'   Gamma(shape, rate) prior on the detection variance partition sparsity
#'   parameter. Default: `c(1, 1)`.
#' @param iota_ell Numeric vector of length 2. `c(alpha, beta)` for an
#'   InvGamma(alpha, beta) prior on the spatial GP length scale(s). Only used
#'   when `spatial = "gp"` in [fit_model()]. Default: `c(1, 1)`.
#' @param kappa_ell Numeric vector of length 2. `c(alpha, beta)` for an
#'   InvGamma(alpha, beta) prior on the exp. quad. temporal GP length scale.
#'   Only used when `temporal = "gp"` in [fit_model()]. Default: `c(1, 1)`.
#' @param kappa_ell_periodic Numeric vector of length 2. `c(alpha, beta)` for an
#'   InvGamma(alpha, beta) prior on the periodic temporal GP length scale. Only
#'   used when `temporal = "gp"` and `periodic = TRUE` in [fit_model()].
#'   Default: `c(1, 1)`.
#' @param K_phi Numeric vector of length 2. `c(alpha[1], alpha[2])` for a
#'   Dirichlet(alpha) prior on the temporal GP variance partitions of the
#'   exp. quad. and periodic kernels. Only used when `temporal_gp = TRUE` and
#'   `periodic_GP = TRUE`. Default: `c(1, 1)`.
#' @param phi Numeric vector of length 2. `c(alpha, beta)` for an
#'   InvGamma(alpha, beta) prior on species-specific negative binomial
#'   overdispersion parameters. Only used when `overdispersion = "nb"` in
#'   [fit_model()]. Default: `c(0.4, 0.3)`.
#' @param alpha_O_L Positive scalar. LKJ prior on the \eqn{[2, 2]} occupancy
#'   log odds/log detection rate correlation matrix. Default: `1`.
#' @param psi_beta_O_L Positive scalar. LKJ prior on the \eqn{[S, S]}
#'   correlation matrix of species-specific occupancy site coefficients. Only
#'   used when site predictors are supplied for occupancy in [make_data()].
#'   Default: `1`.
#' @param mu_beta_O_L Positive scalar. LKJ prior on the \eqn{[S, S]}
#'   correlation matrix of species-specific detection site coefficients. Only
#'   used when site predictors are supplied for detection in [make_data()].
#'   Default: `1`.
#' @param gamma_O_L Positive scalar. LKJ prior on the \eqn{[S, S]} correlation
#'   matrix of species-specific survey coefficients. Only used when survey
#'   predictors are supplied. Default: `1`.
#' @param iota_O_L Positive scalar. LKJ prior on the \eqn{[S, S]} correlation
#'   matrix of species-specific site effects. Default: `1`.
#' @param kappa_O_L Positive scalar. LKJ prior on the \eqn{[S, S]} correlation
#'   matrix of species-specific survey effects. Default: `1`.
#' @param epsilon_O_L Positive scalar. LKJ prior on the \eqn{[S, S]}
#'   correlation matrix of species-specific OLRE residuals. Only used when
#'   `overdispersion = "olre"` in [fit_model()]. Default: `1`.
#' @param verbose Logical. If `TRUE` (default), prints list of priors.
#'
#' @return An `occARU_priors` object (a named list) for use with [fit_model()].
#'
#' @seealso [fit_model()]
#' @export
set_priors <- function(
  psi_bar = c(1, 1),
  mu_bar = c(1, 1),
  psi_W = c(3, 0, 1),
  mu_W = c(3, 0, 2.5),
  psi_theta = c(1, 1),
  mu_theta = c(1, 1),
  iota_ell = c(1, 1),
  kappa_ell = c(1, 1),
  kappa_ell_periodic = c(1, 1),
  K_phi = c(1, 1),
  phi = c(0.4, 0.3),
  alpha_O_L = 1,
  psi_beta_O_L = 1,
  mu_beta_O_L = 1,
  gamma_O_L = 1,
  iota_O_L = 1,
  kappa_O_L = 1,
  epsilon_O_L = 1,
  verbose = TRUE
) {
  # length checks
  expected_lengths <- list(
    psi_bar = 2L,
    mu_bar = 2L,
    psi_W = 3L,
    mu_W = 3L,
    psi_theta = 2L,
    mu_theta = 2L,
    iota_ell = 2L,
    kappa_ell = 2L,
    kappa_ell_periodic = 2L,
    K_phi = 2L,
    phi = 2L,
    alpha_O_L = 1L,
    psi_beta_O_L = 1L,
    mu_beta_O_L = 1L,
    gamma_O_L = 1L,
    iota_O_L = 1L,
    kappa_O_L = 1L,
    epsilon_O_L = 1L
  )

  # length and positivity
  args <- mget(names(expected_lengths))
  purrr::iwalk(expected_lengths, \(expected_len, nm) {
    if (length(args[[nm]]) != expected_len) {
      cli::cli_abort(
        "{.arg {nm}} must have length {expected_len}, not {length(args[[nm]])}."
      )
    }
    if (any(args[[nm]] < 0)) {
      cli::cli_abort("All values in {.arg {nm}} must not be negative.")
    }
  })

  # set priors
  priors <- structure(
    list(
      psi_bar_beta = psi_bar,
      mu_bar_gamma = mu_bar,
      psi_W_t = psi_W,
      mu_W_t = mu_W,
      psi_theta_gamma = psi_theta,
      mu_theta_gamma = mu_theta,
      iota_ell_inv_gamma = iota_ell,
      kappa_ell_inv_gamma = kappa_ell,
      kappa_ell_periodic_inv_gamma = kappa_ell_periodic,
      K_phi_dirichlet = K_phi,
      phi_inv_gamma = phi,
      alpha_O_L_LKJ = alpha_O_L,
      psi_beta_O_L_LKJ = psi_beta_O_L,
      mu_beta_O_L_LKJ = mu_beta_O_L,
      gamma_O_L_LKJ = gamma_O_L,
      iota_O_L_LKJ = iota_O_L,
      kappa_O_L_LKJ = kappa_O_L,
      epsilon_O_L_LKJ = epsilon_O_L
    ),
    class = "occARU_priors"
  )
  if (verbose) {
    print(priors)
  }
  invisible(priors)
}


#' Print method for occARU_priors objects
#'
#' @param x A `occARU_priors` object.
#' @param ... Ignored.
#' @keywords internal
#' @export
print.occARU_priors <- function(x, ...) {
  cli::cli_h1("occARU priors")
  cli::cli_dl(c(
    "psi_bar" = "Beta({x$psi_bar_beta[1]}, {x$psi_bar_beta[2]})",
    "mu_bar" = "Gamma({x$mu_bar_gamma[1]}, {x$mu_bar_gamma[2]})",
    "psi_W" = "Student-t+({x$psi_W_t[1]}, {x$psi_W_t[2]}, {x$psi_W_t[3]})",
    "mu_W" = "Student-t+({x$mu_W_t[1]}, {x$mu_W_t[2]}, {x$mu_W_t[3]})",
    "psi_theta" = "Gamma({x$psi_theta_gamma[1]}, {x$psi_theta_gamma[2]})",
    "mu_theta" = "Gamma({x$mu_theta_gamma[1]}, {x$mu_theta_gamma[2]})",
    "iota_ell" = "InvGamma({x$iota_ell_inv_gamma[1]}, \\
                  {x$iota_ell_inv_gamma[2]})",
    "kappa_ell (exp. quad.)" = "InvGamma({x$kappa_ell_inv_gamma[1]}, \\
                                         {x$kappa_ell_inv_gamma[2]})",
    "kappa_ell (periodic)" = "InvGamma({x$kappa_ell_periodic_inv_gamma[1]}, \\
                             {x$kappa_ell_periodic_inv_gamma[2]})",
    "K_phi" = "Dirichlet({x$K_phi_dirichlet[1]}, {x$K_phi_dirichlet[2]})",
    "phi" = "InvGamma({x$phi_inv_gamma[1]}, {x$phi_inv_gamma[2]})",
    "alpha_O_L" = "LKJ({x$alpha_O_L_LKJ})",
    "psi_beta_O_L" = "LKJ({x$psi_beta_O_L_LKJ})",
    "mu_beta_O_L" = "LKJ({x$mu_beta_O_L_LKJ})",
    "gamma_O_L" = "LKJ({x$gamma_O_L_LKJ})",
    "iota_O_L" = "LKJ({x$iota_O_L_LKJ})",
    "kappa_O_L" = "LKJ({x$kappa_O_L_LKJ})",
    "epsilon_O_L" = "LKJ({x$epsilon_O_L_LKJ})"
  ))
}
