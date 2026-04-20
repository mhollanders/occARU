functions {
  #include occARU.stanfunctions
  #include util.stanfunctions
}

data {
  // dimensions, recording effort, site UTMs (km), and det. history
  int<lower=2> I, J;
  matrix<lower=0, upper=1>[I, J] Delta;
  array[I] vector[2] XY;
  array[I, J] int<lower=0> y;

  // continuous preds for site occ. and site and survey det.
  array[3] int<lower=0> P;
  matrix[I, P[1]] X1;
  matrix[I, P[2]] X2;
  array[I] matrix[J, P[3]] X3;

  // categorical preds for site occ. and site and survey det.
  array[3] int<lower=0> P_cat;
  array[I, P_cat[1]] int<lower=1> X_cat1;
  array[I, P_cat[2]] int<lower=1> X_cat2;
  array[I, J, P_cat[3]] int<lower=1> X_cat3;

  // ordinal preds for site occ. and site and survey det.
  array[3] int<lower=0> P_ord;
  array[I, P_ord[1]] int<lower=0> X_ord1;
  array[I, P_ord[2]] int<lower=0> X_ord2;
  array[I, J, P_ord[3]] int<lower=0> X_ord3;

  // scalars and indicators
  real<lower=0> period;  // period for periodic kernel
  int<lower=0, upper=I> grainsize;  // grainsize for threading
  int<lower=0> D;  // Monte Carlo draws for LOO
  int<lower=0, upper=2> spatial,  // site effects (GP, MVN, none)
                        temporal,  // survey effects (GP, MVN, none)
                        OD;  // overdispersion (Poisson, OLRE, negbin)
  int<lower=0, upper=1> dirichlet,  // var. decomp. (logi-norm or Dirichlet)
                        project_kappa,  // orth. proj. survey effects
                        latent,  // recover latent occ. states
                        PPC_y,  // posterior predictions of y
                        PPC_Q;  // posterior predictions of Q
  // priors
  vector<lower=0>[2] psi_bar_beta,  // mean occ. (prob. scale)
                     mu_bar_gamma,  // mean det. (rate scale)
                     psi_theta_gamma,  // occ. var. partition scale
                     mu_theta_gamma,  // det. var. partition scale
                     iota_ell_inv_gamma,  // spat. GP l-scale
                     kappa_ell_inv_gamma,  // temp. GP l-scale
                     kappa_ell_periodic_inv_gamma,  // periodic temp. GP l-scale
                     K_phi_dirichlet,  // temp. GPs var. partitions
                     phi_inv_gamma;  // negbin overdispersion
  vector<lower=0>[3] psi_W_t,  // occ. log odds var.
                     mu_W_t;  // log det. var.
}

transformed data {
  // transformed indicators
  int SP = spatial > 0, TE = temporal > 0, I_GP = spatial == 2,
      J_GP = temporal == 2, periodic = period > 0 && J_GP, OLRE = OD == 1,
      NB = OD == 2, MC = (D > 0) * (SP || OLRE), Im1 = I - 1;
  real log_D = log(D);
  array[3] int P_sum;
  for (d in 1:3) {
    P_sum[d] = P[d] + P_cat[d] + P_ord[d];
  }

  // survey indices, total surveys, and offsets
  tuple(array[I, J] int, array[I] int) indices = survey_indices(Delta);
  array[I, J] int j_idx = indices.1;
  array[I] int J_i = indices.2;
  int J_sum = sum(J_i);
  matrix[J, I] log_Delta = log(Delta');

  // aggregated counts
  array[I] int Q;
  for (i in 1:I) {
    Q[i] = sum(y[i, j_idx[i, :J_i[i]]]);
  }

  // site and survey sequences
  array[I] int sites = linspaced_int_array(I, 1, I);
  array[J] real surveys = linspaced_array(J, 1, J);

  // var. partitions for occ. and det.
  int psi_V = P_sum[1],
      mu_V = sum(P_sum[2:3]) + (SP + TE + OLRE);

  // categorical levels
  array[3, max(P_cat)] int C = rep_array(0, 3, max(P_cat));
  for (p in 1:P_cat[1]) {
    C[1, p] = max(X_cat1[:, p]);
  }
  for (p in 1:P_cat[2]) {
    C[2, p] = max(X_cat2[:, p]);
  }
  for (p in 1:P_cat[3]) {
    for (i in 1:I) {
      C[3, p] = max(C[3, p], max(X_cat3[i, :, p]));
    }
  }
  array[3] int C_max = zeros_int_array(3), C_sum = zeros_int_array(3);
  for (d in 1:3) {
    if (P_cat[d]) {
      C_max[d] = max(C[d, :P_cat[d]]);
      C_sum[d] = sum(C[d, :P_cat[d]]);
    }
  }

  // ordinal levels
  array[3, max(P_ord)] int O = rep_array(0, 3, max(P_ord));
  for (p in 1:P_ord[1]) {
    O[1, p] = max(X_ord1[:, p]);
  }
  for (p in 1:P_ord[2]) {
    O[2, p] = max(X_ord2[:, p]);
  }
  for (p in 1:P_ord[3]) {
    for (i in 1:I) {
      O[3, p] = max(O[3, p], max(X_ord3[i, :, p]));
    }
  }
  array[3] int O_max = zeros_int_array(3), O_sum = zeros_int_array(3);
  for (d in 1:3) {
    if (P_ord[d]) {
      O_max[d] = max(O[d, :P_ord[d]]);
      O_sum[d] = sum(O[d, :P_ord[d]]);
    }
  }

  // site-level design matrix pseudo-inverses for orthogonal projection
  array[2] int X_plus_cols;
  for (d in 2:3) {
    X_plus_cols[d - 1] = P[d] + C_sum[d] - P_cat[d] + P_ord[d];
  }
  matrix[I, X_plus_cols[1]] X2_aug;
  matrix[X_plus_cols[1], I] X2_plus;
  if (P_sum[2]) {
    X2_aug[:, :P[2]] = X2;
    if (P_cat[2]) {
      int idx = P[2];
      for (p in 1:P_cat[2]) {
        int C_p = C[2, p], C_pm1 = C_p - 1;
        for (i in 1:I) {
          X2_aug[i, idx + 1:idx + C_pm1] =
            one_hot_row_vector(C_p, X_cat2[i, p])[:C_pm1];
        }
        idx += C_pm1;
      }
    }
    if (P_ord[2]) {
      int idx = X_plus_cols[1] - P_ord[2] + 1;
      X2_aug[:, idx:] = to_matrix(X_ord2);
    }
    X2_plus = pseudo_inverse(X2_aug);
  }

  // site-by-survey design matrix pseudo-inverses for orthogonal projection
  matrix[J, X_plus_cols[2]] X3_mean = rep_matrix(0, J, X_plus_cols[2]);
  matrix[X_plus_cols[2], J] X3_plus;
  if (P_sum[3] && project_kappa) {
    if (P[3]) {
      for (i in 1:I) {
        X3_mean[:, :P[3]] += X3[i];
      }
    }
    if (P_cat[3]) {
      int idx = P[3];
      for (p in 1:P_cat[3]) {
        int C_p = C[3, p], C_pm1 = C_p - 1;
        for (i in 1:I) {
          for (j in 1:J) {
            X3_mean[j, idx + 1:idx + C_pm1] +=
              one_hot_row_vector(C_p, X_cat3[i, j, p])[:C_pm1];
          }
        }
        idx += C_pm1;
      }
    }
    if (P_ord[3]) {
      int idx = X_plus_cols[2] - P_ord[3] + 1;
      for (i in 1:I) {
        X3_mean[:, idx:] += to_matrix(X_ord3[i]);
      }
    }
    X3_mean /= I;
    X3_plus = pseudo_inverse(X3_mean);
  }
}

parameters {
  // intercepts
  real<lower=0, upper=1> psi_bar;
  real<lower=0> mu_bar;

  // total variances, partition sparsity, and unconstrained variance partitions
  vector<lower=0>[psi_V > 0] psi_W;
  vector<lower=0>[psi_V > 1] psi_theta;
  array[psi_V > 1] sum_to_zero_vector[psi_V > 1 ? psi_V : 1] psi_phi_z;
  vector<lower=0>[mu_V > 0] mu_W;
  vector<lower=0>[mu_V > 1] mu_theta;
  array[mu_V > 1] sum_to_zero_vector[mu_V > 1 ? mu_V : 1] mu_phi_z;

  // coefficients for continuous predictors
  vector[P[1]] psi_beta_z;
  vector[P[2]] mu_beta_z;
  vector[P[3]] gamma_z;

  // unconstrained coefficients for categorical predictors
  vector[C_sum[1] - P_cat[1]] psi_beta_cat_u;
  vector[C_sum[2] - P_cat[2]] mu_beta_cat_u;
  vector[C_sum[3] - P_cat[3]] gamma_cat_u;

  // coefficients for ordinal predictors
  vector[P_ord[1]] psi_beta_ord_z;
  vector[O_sum[1] - P_ord[1]] psi_beta_ord_c_u;
  vector[P_ord[2]] mu_beta_ord_z;
  vector[O_sum[2] - P_ord[2]] mu_beta_ord_c_u;
  vector[P_ord[3]] gamma_ord_z;
  vector[O_sum[3] - P_ord[3]] gamma_ord_c_u;

  // site and survey effects and GP length-scales
  array[SP] sum_to_zero_vector[I] iota_z;
  vector<lower=0>[I_GP] iota_ell;
  array[TE] sum_to_zero_vector[J] kappa_z;
  vector<lower=0>[J_GP + periodic] kappa_ell;
  array[periodic] simplex[2] K_phi;

  // OLRE residuals or negbin overdispersion
  array[OLRE] sum_to_zero_vector[J_sum] epsilon_z;
  vector<lower=0>[NB] phi;
}

transformed parameters {
  // occupancy log odds variance partitions
  vector[psi_V] psi_phi, psi_tau;
  if (psi_V > 1) {
    psi_phi = dirichlet ?
              sum_to_zero_simplex_jacobian(psi_phi_z[1])
              : softmax(psi_theta[1] * psi_phi_z[1]);
    psi_tau = sqrt(psi_W[1] * psi_phi);
  } else if (psi_V == 1) {
    psi_phi[1] = 1;
    psi_tau[1] = sqrt(psi_W[1]);
  }

  // log detection variance partitions
  vector[mu_V] mu_phi, mu_tau;
  if (mu_V > 1) {
    mu_phi = dirichlet ?
             sum_to_zero_simplex_jacobian(mu_phi_z[1])
             : softmax(mu_theta[1] * mu_phi_z[1]);
    mu_tau = sqrt(mu_W[1] * mu_phi);
  } else if (psi_V == 1) {
    mu_phi[1] = 1;
    mu_tau[1] = sqrt(mu_W[1]);
  }

  // intercepts
  vector[2] alpha = [ logit(psi_bar), log(mu_bar) ]';

  // continuous predictor coefficients
  vector[P[1]] psi_beta;
  vector[P[2]] mu_beta;
  vector[P[3]] gamma;

  // categorical predictor coefficients
  matrix[C_max[1], P_cat[1]] psi_beta_cat_z, psi_beta_cat;
  matrix[C_max[2], P_cat[2]] mu_beta_cat_z, mu_beta_cat;
  matrix[C_max[3], P_cat[3]] gamma_cat_z, gamma_cat;

  // ordinal predictor coefficients
  vector[P_ord[1]] psi_beta_ord;
  matrix[O_max[1], P_ord[1]] psi_beta_ord_cs;
  vector[P_ord[2]] mu_beta_ord;
  matrix[O_max[2], P_ord[2]] mu_beta_ord_cs;
  vector[P_ord[3]] gamma_ord;
  matrix[O_max[3], P_ord[3]] gamma_ord_cs;

  // conditional and unconditional site and survey effects
  vector[SP * I] iota;
  vector[SP * P_sum[2] ? I : 0] iota2;
  vector[TE * J] kappa;
  vector[TE * P_sum[3] * project_kappa ? J : 0] kappa2;
  {
    int psi_idx = 1, mu_idx = 1;

    // occupancy continuous coefficients
    if (P[1]) {
      psi_beta = segment(psi_tau, psi_idx, P[1]) .* psi_beta_z;
      psi_idx += P[1];
    }

    // occupancy categorical coefficients
    if (P_cat[1]) {
      int P_c = P_cat[1], C_m = C_max[1];
      tuple(matrix[C_m, P_c], matrix[C_m, P_c]) coef =
        coef_cat_jacobian(segment(psi_tau, psi_idx, P_c), psi_beta_cat_u,
                          C[1, :P_c]);
      psi_beta_cat = coef.1;
      psi_beta_cat_z = coef.2;
      psi_idx += P_c;
    }

    // occupancy ordinal coefficients
    if (P_ord[1]) {
      int P_o = P_ord[1];
      tuple(vector[P_o], matrix[O_max[1], P_o]) coef =
        coef_ord_jacobian(segment(psi_tau, psi_idx, P_o), psi_beta_ord_z,
                          psi_beta_ord_c_u, O[1, :P_o]);
      psi_beta_ord = coef.1;
      psi_beta_ord_cs = coef.2;
      psi_idx += P_o;
    }

    // detection continuous coefficients (site-level)
    if (P[2]) {
      mu_beta = segment(mu_tau, mu_idx, P[2]) .* mu_beta_z;
      mu_idx += P[2];
    }

    // detection categorical coefficients (site-level)
    if (P_cat[2]) {
      int P_c = P_cat[2], C_m = C_max[2];
      tuple(matrix[C_m, P_c], matrix[C_m, P_c]) coef =
        coef_cat_jacobian(segment(mu_tau, mu_idx, P_c), mu_beta_cat_u,
                          C[2, :P_c]);
      mu_beta_cat = coef.1;
      mu_beta_cat_z = coef.2;
      mu_idx += P_c;
    }

    // detection ordinal coefficients (site-level)
    if (P_ord[2]) {
      int P_o = P_ord[2];
      tuple(vector[P_o], matrix[O_max[1], P_o]) coef =
        coef_ord_jacobian(segment(mu_tau, mu_idx, P_o), mu_beta_ord_z,
                          mu_beta_ord_c_u, O[2, :P_o]);
      mu_beta_ord = coef.1;
      mu_beta_ord_cs = coef.2;
      mu_idx += P_o;
    }

    // detection continuous coefficients (site-by-survey level)
    if (P[3]) {
      gamma = segment(mu_tau, mu_idx, P[3]) .* gamma_z;
      mu_idx += P[3];
    }

    // detection categorical coefficients (site-by-survey level)
    if (P_cat[3]) {
      int P_c = P_cat[3], C_m = C_max[3];
      tuple(matrix[C_m, P_c], matrix[C_m, P_c]) coef =
        coef_cat_jacobian(segment(mu_tau, mu_idx, P_c), gamma_cat_u,
                          C[3, :P_c]);
      gamma_cat = coef.1;
      gamma_cat_z = coef.2;
      mu_idx += P_c;
    }

    // detection ordinal coefficients (site-by-survey level)
    if (P_ord[3]) {
      int P_o = P_ord[3];
      tuple(vector[P_o], matrix[O_max[1], P_o]) coef =
        coef_ord_jacobian(segment(mu_tau, mu_idx, P_o), gamma_ord_z,
                          gamma_ord_c_u, O[3, :P_o]);
      gamma_ord = coef.1;
      gamma_ord_cs = coef.2;
      mu_idx += P_o;
    }

    // site effects
    if (SP) {
      iota = mu_tau[mu_idx] * iota_z[1];
      if (I_GP) {
        matrix[I, I] iota_K = gp_exp_quad_cov(XY, 1, iota_ell[1]),
                     iota_L = cholesky_decompose(add_diag(iota_K, 1e-9))';
        iota = iota_L * iota;
      }

      // increment orthogonal projection
      if (P_sum[2]) {
        iota2 = iota;
        iota = orthogonalise(iota, X2_aug, X2_plus);
      }
      mu_idx += 1;
    }

    // survey effects
    if (TE) {
      kappa = mu_tau[mu_idx] * kappa_z[1];
      if (J_GP) {
        matrix[J, J] kappa_K, kappa_L;
        vector[periodic * 2] kappa_t;
        if (periodic) {
          kappa_t = sqrt(K_phi[1]);
          kappa_K = gp_exp_quad_cov(surveys, kappa_t[1], kappa_ell[1])
                    + gp_periodic_cov(surveys, kappa_t[2], kappa_ell[2],
                                      period);
        } else {
          kappa_K = gp_exp_quad_cov(surveys, 1, kappa_ell[1]);
        }
        kappa_L = cholesky_decompose(add_diag(kappa_K, 1e-9));
        kappa = kappa_L * kappa;
      }

      // increment orthogonal projection
      if (P_sum[3] && project_kappa) {
        kappa2 = kappa;
        kappa = orthogonalise(kappa, X3_mean, X3_plus);
      }
    }
  }

  // occupancy
  vector[I] logit_psi = alpha[1] + X1 * psi_beta;
  for (p in 1:P_cat[1]) {
    logit_psi += psi_beta_cat[X_cat1[:, p], p];
  }
  if (P_ord[1]) {
    logit_psi += increment_ordinal(X_ord1, coef_ord_realised(O[1],
                                   psi_beta_ord, psi_beta_ord_cs));
  }

  // priors
  real lprior = beta_lpdf(psi_bar | psi_bar_beta[1], psi_bar_beta[2])
                + gamma_lpdf(mu_bar | mu_bar_gamma[1], mu_bar_gamma[2]);
  if (psi_V) {
    lprior += student_t_lpdf(psi_W[1] | psi_W_t[1], psi_W_t[2], psi_W_t[3]);
    if (psi_V > 1) {
      lprior += gamma_lpdf(psi_theta[1] | psi_theta_gamma[1],
                                          psi_theta_gamma[2]);
    }
  }
  if (mu_V) {
    lprior += student_t_lpdf(mu_W[1] | mu_W_t[1], mu_W_t[2], mu_W_t[3]);
    if (mu_V > 1) {
      lprior += gamma_lpdf(mu_theta[1] | mu_theta_gamma[1], mu_theta_gamma[2]);
    }
  }
  if (I_GP) {
    lprior += inv_gamma_lpdf(iota_ell[1] | iota_ell_inv_gamma[1],
                                           iota_ell_inv_gamma[2]);
  }
  if (J_GP) {
    lprior += inv_gamma_lpdf(kappa_ell[1] | kappa_ell_inv_gamma[1],
                                            kappa_ell_inv_gamma[2]);
    if (periodic) {
      lprior += inv_gamma_lpdf(kappa_ell[2] | kappa_ell_periodic_inv_gamma[1],
                                              kappa_ell_periodic_inv_gamma[2])
                + dirichlet_lpdf(K_phi | K_phi_dirichlet);
    }
  }
  if (NB) {
    lprior += inv_gamma_lpdf(phi[1] | phi_inv_gamma[1], phi_inv_gamma[2]);
  }
}

model {
  target += lprior
            + std_normal_lupdf(psi_beta_z)
            + std_normal_lupdf(psi_beta_ord_z)
            + std_normal_lupdf(mu_beta_z)
            + std_normal_lupdf(mu_beta_ord_z)
            + std_normal_lupdf(gamma_z)
            + std_normal_lupdf(gamma_ord_z);
  if (psi_V > 1) {
    target += dirichlet ?
              dirichlet_lupdf(psi_phi | rep_vector(inv(psi_theta[1]), psi_V))
              : std_normal_lupdf(psi_phi_z[1]);
  }
  if (mu_V > 1) {
    target += dirichlet ?
              dirichlet_lupdf(mu_phi | rep_vector(inv(mu_theta[1]), mu_V))
              : std_normal_lupdf(mu_phi_z[1]);
  }
  for (p in 1:P_cat[1]) {
    target += std_normal_lupdf(psi_beta_cat_z[:C[1, p], p]);
  }
  for (p in 1:P_cat[2]) {
    target += std_normal_lupdf(mu_beta_cat_z[:C[2, p], p]);
  }
  for (p in 1:P_cat[3]) {
    target += std_normal_lupdf(gamma_cat_z[:C[3, p], p]);
  }
  if (SP) {
    target += std_normal_lupdf(iota_z[1]);
  }
  if (TE) {
    target += std_normal_lupdf(kappa_z[1]);
  }

  // Poisson OLREs
  matrix[OLRE * J, I] epsilon;
  if (OLRE) {
    target += std_normal_lupdf(epsilon_z[1]);
    epsilon = fill_epsilon(j_idx, J_i, mu_tau[mu_V] * epsilon_z[1]);
  }

  // likelihood
  target += grainsize ?
            reduce_sum(partial_occARU_lupmf, sites, grainsize, y, Q, j_idx, J_i,
                       log_Delta, X2, X_cat2, X_ord2, X3, X_cat3, X_ord3,
                       logit_psi, alpha[2], mu_beta, mu_beta_cat,
                       coef_ord_realised(O[2], mu_beta_ord, mu_beta_ord_cs),
                       gamma, gamma_cat, coef_ord_realised(O[3], gamma_ord,
                       gamma_ord_cs), iota, kappa, epsilon, phi)
            : occARU_lupmf(y | Q, j_idx, J_i, log_Delta, X2, X_cat2, X_ord2, X3,
                           X_cat3, X_ord3, logit_psi, alpha[2], mu_beta,
                           mu_beta_cat, coef_ord_realised(O[2], mu_beta_ord,
                           mu_beta_ord_cs), gamma, gamma_cat,
                           coef_ord_realised(O[3], gamma_ord, gamma_ord_cs),
                           iota, kappa, epsilon, phi);
}

generated quantities {
  // unconditional site coefficients
  vector[SP * P[2]] mu_beta2;
  vector[SP * P_ord[2]] mu_beta_ord2;
  if (SP) {
    if (P[2]) {
      mu_beta2 = mu_bar - X2_plus[:P[2]] * iota2;
    }
    if (P_ord[2]) {
      mu_beta_ord2 = mu_beta_ord
                     - X2_plus[X_plus_cols[1] - P_ord[2] + 1:] * iota2;
    }
  }

  // unconditional survey coefficients
  vector[TE * P[3] * project_kappa] gamma2;
  vector[TE * P_ord[3] * project_kappa] gamma_ord2;
  if (TE * project_kappa) {
    if (P[3]) {
      gamma2 = gamma - X3_plus[:P[3]] * kappa2;
    }
    if (P_ord[3]) {
      gamma_ord2 = gamma_ord
                   - X3_plus[:, X_plus_cols[2] - P_ord[3] + 1:] * kappa2;
    }
  }

  // log likelihood, latent occupancy, and posterior predictions
  vector[I] log_lik;
  vector[MC * I] log_lik2;
  array[latent * I] int z = ones_int_array(latent * I);
  array[PPC_y * I, J] int yrep = rep_array(0, PPC_y * I, J);
  array[PPC_Q * I] int Qrep = zeros_int_array(PPC_Q * I);

  {
    // reconstruct log likelihood
    matrix[OLRE * J, I] epsilon;
    if (OLRE) {
      epsilon = fill_epsilon(j_idx, J_i, mu_tau[mu_V] * epsilon_z[1]);
    }
    tuple(vector[I], matrix[2, I], matrix[J, I]) lp =
      occARU(y, Q, j_idx, J_i, log_Delta, X2, X_cat2, X_ord2, X3, X_cat3, X_ord3,
             logit_psi, alpha[2], mu_beta, mu_beta_cat,
             coef_ord_realised(O[2], mu_beta_ord, mu_beta_ord_cs), gamma,
             gamma_cat, coef_ord_realised(O[3], gamma_ord, gamma_ord_cs), iota,
             kappa, epsilon, phi);
    log_lik = lp.1;

    // recover latent occupancy
    if (latent) {
      for (i in 1:I) {
        if (!Q[i]) {
          z[i] = bernoulli_logit_rng(lp.2[2, i] - log_lik[i]);
        }
      }
    }

    // posterior predictions
    if (PPC_y || PPC_Q) {
      if (OLRE) {
        epsilon = pp_epsilon_rng(j_idx, J_i, mu_tau[mu_V]);
      }
      tuple(array[I, J] int, array[I] int) pp =
        pp_occARU_rng(j_idx, J_i, X3, logit_psi, lp.3, gamma, epsilon, phi);
      if (PPC_y) {
        yrep = pp.1;
      }
      if (PPC_Q) {
        Qrep = pp.2;
      }
    }

    // Monte Carlo integration of observation-level parameters for loo
    if (MC) {
      matrix[I, D] log_lik_k;

      // prepare site effect covariances
      matrix[SP * I, I] iota_L;
      if (I_GP) {
        iota_L =
          cholesky_decompose(add_diag(gp_exp_quad_cov(XY, 1, iota_ell[1]),
                                      1e-9));
      }

      // posterior predict site effects and OLREs and compute log likelihood
      vector[Im1] zeros = zeros_vector(Im1);
      vector[SP * Im1] iota_u;
      vector[SP * I] iota_rep;
      int mu_idx = sum(P_sum[2:3]) + 1;
      for (d in 1:D) {
        if (SP) {
          iota_u = to_vector(normal_rng(zeros, mu_tau[mu_idx]));
          iota_rep = sum_to_zero_constrain(iota_u);
          if (I_GP) {
            iota_rep = iota_L * iota_rep;
          }
        }
        if (OLRE) {
          epsilon = pp_epsilon_rng(j_idx, J_i, mu_tau[mu_V]);
        }
        lp = occARU(y, Q, j_idx, J_i, log_Delta, X2, X_cat2, X_ord2, X3, X_cat3,
                    X_ord3, logit_psi, alpha[2], mu_beta, mu_beta_cat,
                    coef_ord_realised(O[2], mu_beta_ord, mu_beta_ord_cs),
                    gamma, gamma_cat, coef_ord_realised(O[3], gamma_ord,
                    gamma_ord_cs), iota_rep, kappa, epsilon, phi);
        log_lik_k[:, d] = lp.1;
      }
      log_lik2 = rep_vector(-log_D, I);
      for (i in 1:I) {
        log_lik2[i] += log_sum_exp(log_lik_k[i]);
      }
    }
  }
}
