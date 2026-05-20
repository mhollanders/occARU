#' Set priors for the occARU model
#'
#' Constructs and validates a named list of prior hyperparameters for use with
#' [occARU()]. Any unspecified priors use the defaults listed below.
#'
#' @param psi_bar Numeric vector of length 2. `c(a, b)` for a Beta(a, b) prior
#'   on mean occupancy probability. Default: `c(1, 1)`.
#' @param q_bar Numeric vector of length 2. `c(shape, rate)` for
#'   Gamma(shape, rate) priors on the mean annual colonisation and emigration
#'   rates for multiseason models. Default: `c(1, 3)`.
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
#'   InvGamma(alpha, beta) prior on the detection spatial GP length scale.
#'   Default: `c(1, 1)`.
#' @param kappa_ell Numeric vector of length 2. `c(alpha, beta)` for an
#'   InvGamma(alpha, beta) prior on the temporal GP length scale. Default:
#'   `c(1, 1)`.
#' @param kappa_ell_periodic Numeric vector of length 2. `c(alpha, beta)` for an
#'   InvGamma(alpha, beta) prior on the periodic temporal GP length scale.
#'   Default: `c(1, 1)`.
#' @param K_phi Numeric vector of length 2. `c(alpha[1], alpha[2])` for a
#'   Dirichlet(alpha) prior on the temporal GP variance partitions of the
#'   Matern and periodic kernels. Default: `c(1, 1)`.
#' @param nu_ell Numeric vector of length 2. `c(alpha, beta)` for an
#'   InvGamma(alpha, beta) prior on the detection season GP length scale.
#'   Default: `c(1, 1)`.
#' @param psi_iota_ell Numeric vector of length 2. `c(alpha, beta)` for an
#'   InvGamma(alpha, beta) prior on the occupancy spatial GP length scale.
#'   Default: `c(1, 1)`.
#' @param psi_nu_ell Numeric vector of length 2. `c(alpha, beta)` for an
#'   InvGamma(alpha, beta) prior on the occupancy season GP length scale.
#'   Default: `c(1, 1)`.
#' @param phi Numeric vector of length 2. `c(alpha, beta)` for an
#'   InvGamma(alpha, beta) prior on species-specific negative binomial
#'   overdispersion parameters. Default: `c(0.4, 0.3)`.
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
#'   matrix of species-specific site effects on detection. Default: `1`.
#' @param kappa_O_L Positive scalar. LKJ prior on the \eqn{[S, S]} correlation
#'   matrix of species-specific survey effects. Default: `1`.
#' @param nu_O_L Positive scalar. LKJ prior on the \eqn{[S, S]} correlation
#'   matrix of species-specific season effects on detection. Default: `1`.
#' @param psi_iota_O_L Positive scalar. LKJ prior on the \eqn{[S, S]}
#'   correlation matrix of species-specific season effects on detection.
#'   Default: `1`.
#' @param psi_nu_O_L Positive scalar. LKJ prior on the \eqn{[S, S]} correlation
#'   matrix of species-specific season effects on occupancy. Default: `1`.
#' @param epsilon_O_L Positive scalar. LKJ prior on the \eqn{[S, S]}
#'   correlation matrix of species-specific OLRE residuals. Default: `1`.
#' @param verbose Logical. If `TRUE` (default), prints list of priors.
#'
#' @return An `occARU_priors` object (a named list) for use with [occARU()].
#'
#' @seealso [occARU()]
#' @export
set_priors <- function(
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
  psi_iota_ell = c(1, 1),
  psi_nu_ell = c(1, 1),
  phi = c(0.4, 0.3),
  alpha_O_L = 1,
  psi_beta_O_L = 1,
  mu_beta_O_L = 1,
  gamma_O_L = 1,
  iota_O_L = 1,
  kappa_O_L = 1,
  nu_O_L = 1,
  psi_iota_O_L = 1,
  psi_nu_O_L = 1,
  epsilon_O_L = 1,
  verbose = TRUE
) {
  # length checks
  expected_lengths <- list(
    psi_bar = 2L,
    q_bar = 2L,
    mu_bar = 2L,
    psi_W = 3L,
    mu_W = 3L,
    psi_theta = 2L,
    mu_theta = 2L,
    iota_ell = 2L,
    kappa_ell = 2L,
    kappa_ell_periodic = 2L,
    nu_ell = 2L,
    psi_iota_ell = 2L,
    psi_nu_ell = 2L,
    K_phi = 2L,
    phi = 2L,
    alpha_O_L = 1L,
    psi_beta_O_L = 1L,
    mu_beta_O_L = 1L,
    gamma_O_L = 1L,
    iota_O_L = 1L,
    kappa_O_L = 1L,
    nu_O_L = 1L,
    psi_iota_O_L = 1L,
    psi_nu_O_L = 1L,
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
      q_bar_gamma = q_bar,
      mu_bar_gamma = mu_bar,
      psi_W_t = psi_W,
      mu_W_t = mu_W,
      psi_theta_gamma = psi_theta,
      mu_theta_gamma = mu_theta,
      iota_ell_inv_gamma = iota_ell,
      kappa_ell_inv_gamma = kappa_ell,
      kappa_ell_periodic_inv_gamma = kappa_ell_periodic,
      K_phi_dirichlet = K_phi,
      nu_ell_inv_gamma = nu_ell,
      psi_iota_ell_inv_gamma = psi_iota_ell,
      psi_nu_ell_inv_gamma = psi_nu_ell,
      phi_inv_gamma = phi,
      alpha_O_L_LKJ = alpha_O_L,
      psi_beta_O_L_LKJ = psi_beta_O_L,
      mu_beta_O_L_LKJ = mu_beta_O_L,
      gamma_O_L_LKJ = gamma_O_L,
      iota_O_L_LKJ = iota_O_L,
      kappa_O_L_LKJ = kappa_O_L,
      nu_O_L_LKJ = nu_O_L,
      psi_iota_O_L_LKJ = psi_iota_O_L,
      psi_nu_O_L_LKJ = psi_nu_O_L,
      epsilon_O_L_LKJ = epsilon_O_L
    ),
    class = "occARU_priors"
  )
  if (verbose) {
    print(priors)
  }
  invisible(priors)
}

#' Automatically determine inverse gamma hyperparameters, used for Gaussian
#' process length scales in occARU models.
#'
#' Finds inverse gamma shape (alpha) and rate (beta) parameters such that
#' \code{tail_prob} prior probability mass falls below \code{bounds[1]} and
#' \code{tail_prob} falls above \code{bounds[2]}, placing
#' \code{1 - 2 * tail_prob} of the mass within the bounds.
#'
#' @param bounds Numeric vector of length 2. Lower and upper bounds on the
#'   length scale. For a spatial GP, the bounds could be the minimum and maximum
#'   distance between sites
#'
#'   For a periodic kernel with a fixed annual cycle, the lower
#'   bound should be the minimum temporal spacing between observations (e.g.
#'   \code{1/52} for weekly data) and the upper bound should be the period
#'   (e.g. \code{1} for an annual cycle).
#' @param tail_prob Numeric scalar in (0, 0.5). The probability mass to place
#'   in each tail. For example, \code{0.05} places 5% of the prior below
#'   \code{bounds[1]} and 5% above \code{bounds[2]}.
#' @return Numeric vector of length 2, \code{c(alpha, beta)}, giving the shape
#'   and rate of the fitted inverse gamma prior.
#' @keywords internal
auto_inv_gamma <- function(bounds, tail_prob) {
  # checks
  if (tail_prob < 0 | tail_prob > 1) {
    cli::cli_abort(
      "{.arg tails} must be between 0 and 1."
    )
  }
  if (any(bounds < 0)) {
    cli::cli_abort(
      "{.arg bounds} must be greater than 0."
    )
  }

  # inverse gamma CDF
  pinvgamma <- function(q, alpha, beta) {
    stats::pgamma(1 / q, alpha, beta, lower.tail = FALSE)
  }

  # initial values
  mid <- sqrt(bounds[1] * bounds[2])
  alpha0 <- 3
  beta0 <- mid * (alpha0 + 1)
  init <- log(c(alpha0, beta0))

  fit <- stats::optim(
    init,
    \(params) {
      alpha <- exp(params[1])
      beta <- exp(params[2])
      r1 <- pinvgamma(bounds[1], alpha, beta) - tail_prob
      r2 <- 1 - pinvgamma(bounds[2], alpha, beta) - tail_prob
      r1^2 + r2^2
    },
    method = "Nelder-Mead"
  )
  fit$par
}

#' Print method for occARU_priors objects
#'
#' @param x A `occARU_priors` object.
#' @param ... Ignored.
#' @keywords internal
#' @export
print.occARU_priors <- function(x, ...) {
  cli::cli_h1("occARU priors")
  K <- attr(x, "K")
  cli::cli_dl(
    c(
      "psi_bar" = "Beta({x$psi_bar_beta[1]}, {x$psi_bar_beta[2]})",
      "q_bar" = "Gamma({x$q_bar_gamma[1]}, {x$q_bar_gamma[2]})",
      "mu_bar" = "Gamma({x$mu_bar_gamma[1]}, {x$mu_bar_gamma[2]})",
      "psi_W" = "Student-t+({x$psi_W_t[1]}, {x$psi_W_t[2]}, {x$psi_W_t[3]})",
      "mu_W" = "Student-t+({x$mu_W_t[1]}, {x$mu_W_t[2]}, {x$mu_W_t[3]})",
      "psi_theta" = "Gamma({x$psi_theta_gamma[1]}, {x$psi_theta_gamma[2]})",
      "mu_theta" = "Gamma({x$mu_theta_gamma[1]}, {x$mu_theta_gamma[2]})",
      "iota_ell" = "InvGamma({x$iota_ell_inv_gamma[1]}, \\
                  {x$iota_ell_inv_gamma[2]})",
      "kappa_ell" = "InvGamma({x$kappa_ell_inv_gamma[1]}, \\
                                         {x$kappa_ell_inv_gamma[2]})",
      "kappa_ell (periodic)" = "InvGamma({x$kappa_ell_periodic_inv_gamma[1]}, \\
                             {x$kappa_ell_periodic_inv_gamma[2]})",
      "K_phi" = "Dirichlet({x$K_phi_dirichlet[1]}, {x$K_phi_dirichlet[2]})",
      "nu_ell" = "InvGamma({x$nu_ell_inv_gamma[1]}, {x$nu_ell_inv_gamma[2]})",
      "psi_iota_ell" = "InvGamma({x$psi_iota_ell_inv_gamma[1]}, {x$psi_iota_ell_inv_gamma[2]})",
      "psi_nu_ell" = "InvGamma({x$psi_nu_ell_inv_gamma[1]}, {x$psi_nu_ell_inv_gamma[2]})",
      "phi" = "InvGamma({x$phi_inv_gamma[1]}, {x$phi_inv_gamma[2]})",
      "alpha_O_L" = "LKJ({x$alpha_O_L_LKJ})",
      "psi_beta_O_L" = "LKJ({x$psi_beta_O_L_LKJ})",
      "mu_beta_O_L" = "LKJ({x$mu_beta_O_L_LKJ})",
      "gamma_O_L" = "LKJ({x$gamma_O_L_LKJ})",
      "iota_O_L" = "LKJ({x$iota_O_L_LKJ})",
      "kappa_O_L" = "LKJ({x$kappa_O_L_LKJ})",
      "nu_O_L" = "LKJ({x$nu_O_L_LKJ})",
      "psi_iota_O_L" = "LKJ({x$psi_iota_O_L_LKJ})",
      "psi_nu_O_L" = "LKJ({x$psi_nu_O_L_LKJ})",
      "epsilon_O_L" = "LKJ({x$epsilon_O_L_LKJ})"
    )
  )
}
