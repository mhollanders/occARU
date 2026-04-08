functions {
  #include occARU.stanfunctions
  #include util.stanfunctions
}

data {
  // dimensions, recording effort, site UTMs (km), and det. history
  int<lower=2> I, J, S;
  matrix<lower=0>[I, J] Delta;
  array[I] vector[2] XY;
  array[I, J, S] int<lower=0> y;

  // continuous preds for site occ. and site and survey det.
  array[3] int<lower=0> P;
  matrix[P[1], I] X1;
  matrix[P[2], I] X2;
  array[I] matrix[P[3], J] X3;

  // categorical preds for site occ. and site and survey det.
  array[3] int<lower=0> P_cat;
  array[P_cat[1], I] int<lower=1> X_cat1;
  array[P_cat[2], I] int<lower=1> X_cat2;
  array[I, P_cat[3], J] int<lower=1> X_cat3;

  // ordinal preds for site occ. and site and survey det.
  array[3] int<lower=0> P_ord;
  array[P_ord[1], I] int<lower=0> X_ord1;
  array[P_ord[2], I] int<lower=0> X_ord2;
  array[I, P_ord[3], J] int<lower=0> X_ord3;

  // scalars and indicators
  real<lower=0> period;  // period for periodic kernel
  int<lower=0, upper=I> grainsize;  // grainsize for threading
  int<lower=0> D;  // Monte Carlo draws for LOO
  int<lower=0, upper=2> spatial,  // site effects (GP, MVN, none)
                        temporal,  // survey effects (GP, MVN, none)
                        OD;  // overdispersion (Poisson, OLRE, negbin)
  int<lower=0, upper=1> dirichlet,  // var. decomp. (logi-norm or Dirichlet)
                        SS,  // species-specific l-scales
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
                     kappa_v_dirichlet,  // temp. GPs var. partitions
                     phi_inv_gamma;  // negbin overdispersion
  vector<lower=0>[3] psi_W_t,  // occ. log odds var.
                     mu_W_t;  // log det. var.
  real<lower=0> alpha_O_L_LKJ, psi_beta_O_L_LKJ,
                mu_beta_O_L_LKJ, gamma_O_L_LKJ,
                iota_O_L_LKJ, kappa_O_L_LKJ,
                epsilon_O_L_LKJ;
}

transformed data {
  // transformed indicators
  int SP = spatial > 0, TE = temporal > 0, I_GP = spatial == 2,
      J_GP = temporal == 2, periodic = period > 0 && J_GP, OLRE = OD == 1,
      NB = OD == 2, MC = (D > 0) * (SP || OLRE), Im1 = I - 1, Sm1 = S - 1,
      Sp1 = S + 1;
  real log_D = log(D);
  array[3] int P_sum;
  for (d in 1:3) {
    P_sum[d] = P[d] + P_cat[d] + P_ord[d];
  }

  // deployment start and end and survey offsets (0-1)
  array[I, 2] int f_l = first_last_survey(Delta);
  matrix[J, I] log_Delta = log(Delta');

  // aggregated counts and total surveys
  array[I, S] int Q = rep_array(0, I, S);
  int N = 0;
  for (i in 1:I) {
    int f = f_l[i, 1], l = f_l[i, 2];
    for (j in f:l) {
      if (!is_inf(log_Delta[j, i])) {
        N += 1;
        for (s in 1:S) {
          Q[i, s] += y[i, j, s];
        }
      }
    }
  }

  // site and survey sequences
  array[I] int sites = linspaced_int_array(I, 1, I);
  array[J] real surveys = linspaced_array(J, 1, J);

  // var. partitions for occ. and det.
  int psi_V = 1 + 2 * P_sum[1],
      mu_V = 1 + 2 * sum(P_sum[2:3]) + (SP + TE + OLRE) * Sp1;

  // categorical levels
  array[3, max(P_cat)] int C = rep_array(0, 3, max(P_cat));
  for (p in 1:P_cat[1]) {
    C[1, p] = max(X_cat1[p]);
  }
  for (p in 1:P_cat[2]) {
    C[2, p] = max(X_cat2[p]);
  }
  for (p in 1:P_cat[3]) {
    for (i in 1:I) {
      C[3, p] = max(C[3, p], max(X_cat3[i, p]));
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
    O[1, p] = max(X_ord1[p]);
  }
  for (p in 1:P_ord[2]) {
    O[2, p] = max(X_ord2[p]);
  }
  for (p in 1:P_ord[3]) {
    for (i in 1:I) {
      O[3, p] = max(O[3, p], max(X_ord3[i, p]));
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
  matrix[X_plus_cols[1], I] X2_aug;
  matrix[I, X_plus_cols[1]] X2_plus;
  if (P_sum[2]) {
    X2_aug[:P[2]] = X2;
    if (P_cat[2]) {
      int idx = P[2];
      for (p in 1:P_cat[2]) {
        int C_p = C[2, p], C_pm1 = C_p - 1;
        for (i in 1:I) {
          X2_aug[idx + 1:idx + C_pm1, i] =
            one_hot_vector(C_p, X_cat2[p, i])[:C_pm1];
        }
        idx += C_pm1;
      }
    }
    if (P_ord[2]) {
      int idx = X_plus_cols[1] - P_ord[2] + 1;
      X2_aug[idx:] = to_matrix(X_ord2);
    }
    X2_plus = pseudo_inverse(X2_aug);
  }

  // site-by-survey design matrix pseudo-inverses for orthogonal projection
  matrix[X_plus_cols[2], J] X3_mean = rep_matrix(0, X_plus_cols[2], J);
  matrix[J, X_plus_cols[2]] X3_plus;
  if (P_sum[3]) {
    if (P[3]) {
      for (i in 1:I) {
        X3_mean[:P[3]] += X3[i];
      }
    }
    if (P_cat[3]) {
      int idx = P[3];
      for (p in 1:P_cat[3]) {
        int C_p = C[3, p], C_pm1 = C_p - 1;
        for (i in 1:I) {
          for (j in 1:J) {
            X3_mean[idx + 1:idx + C_pm1, j] +=
              one_hot_vector(C_p, X_cat3[i, p, j])[:C_pm1];
          }
        }
        idx += C_pm1;
      }
    }
    if (P_ord[3]) {
      int idx = X_plus_cols[2] - P_ord[3] + 1;
      for (i in 1:I) {
        X3_mean[idx:] += to_matrix(X_ord3[i]);
      }
    }
    X3_mean /= I;
    X3_plus = pseudo_inverse(X3_mean);
  }
}

parameters {
  // intercepts and interspecific correlations
  real<lower=0, upper=1> psi_bar;
  real<lower=0> mu_bar;
  cholesky_factor_corr[2] alpha_O_L;
  array[2] sum_to_zero_vector[S] alpha_z;

  // total variances, partition sparsity, and unconstrained variance partitions
  real<lower=0> psi_W, mu_W;
  vector<lower=0>[psi_V > 1] psi_theta;
  array[psi_V > 1] sum_to_zero_vector[psi_V] psi_phi_z;
  vector<lower=0>[mu_V > 1] mu_theta;
  array[mu_V > 1] sum_to_zero_vector[mu_V] mu_phi_z;

  // correlation matrices for species-level predictors
  cholesky_factor_corr[P_sum[1] ? S : 0] psi_beta_O_L;
  cholesky_factor_corr[P_sum[2] ? S : 0] mu_beta_O_L;
  cholesky_factor_corr[P_sum[3] ? S : 0] gamma_O_L;

  // coefficients for continuous predictors
  row_vector[P[1]] psi_beta_bar_z;
  array[P[1]] sum_to_zero_vector[S] psi_beta_z;
  row_vector[P[2]] mu_beta_bar_z;
  array[P[2]] sum_to_zero_vector[S] mu_beta_z;
  row_vector[P[3]] gamma_bar_z;
  array[P[3]] sum_to_zero_vector[S] gamma_z;

  // unconstrained coefficients for categorical predictors
  vector[C_sum[1] - P_cat[1]] psi_beta_cat_bar_u;
  matrix[Sm1, C_sum[1] - P_cat[1]] psi_beta_cat_u;
  vector[C_sum[2] - P_cat[2]] mu_beta_cat_bar_u;
  matrix[Sm1, C_sum[2] - P_cat[2]] mu_beta_cat_u;
  vector[C_sum[3] - P_cat[3]] gamma_cat_bar_u;
  matrix[Sm1, C_sum[3] - P_cat[3]] gamma_cat_u;

  // coefficients for ordinal predictors
  row_vector[P_ord[1]] psi_beta_ord_bar_z;
  array[P_ord[1]] sum_to_zero_vector[S] psi_beta_ord_z;
  vector[O_sum[1] - P_ord[1]] psi_beta_ord_c_u;
  row_vector[P_ord[2]] mu_beta_ord_bar_z;
  array[P_ord[2]] sum_to_zero_vector[S] mu_beta_ord_z;
  vector[O_sum[2] - P_ord[2]] mu_beta_ord_c_u;
  row_vector[P_ord[3]] gamma_ord_bar_z;
  array[P_ord[3]] sum_to_zero_vector[S] gamma_ord_z;
  vector[O_sum[3] - P_ord[3]] gamma_ord_c_u;

  // site and survey effects and GP length-scales
  array[SP] sum_to_zero_vector[I] iota_bar_z;
  array[SP] sum_to_zero_matrix[S, I] iota_z;
  cholesky_factor_corr[SP * S] iota_O_L;
  vector<lower=0>[I_GP * (1 + SS * S)] iota_ell;
  array[TE] sum_to_zero_vector[J] kappa_bar_z;
  array[TE] sum_to_zero_matrix[S, J] kappa_z;
  cholesky_factor_corr[TE * S] kappa_O_L;
  matrix<lower=0>[1 + (SS * S), J_GP + periodic] kappa_ell;
  array[periodic] simplex[2] kappa_v;

  // OLRE residuals or negbin overdispersion
  array[OLRE] sum_to_zero_vector[N] epsilon_bar_z;
  array[OLRE] sum_to_zero_matrix[S, N] epsilon_z;
  cholesky_factor_corr[OLRE * S] epsilon_O_L;
  vector<lower=0>[NB * S] phi;
}

transformed parameters {
  // occupancy log odds variance partitions
  row_vector[psi_V] psi_phi, psi_tau;
  if (psi_V > 1) {
    psi_phi = dirichlet ?
              sum_to_zero_simplex_jacobian(psi_phi_z[1])'
              : softmax(psi_theta[1] * psi_phi_z[1])';
    psi_tau = sqrt(psi_W * psi_phi);
  } else {
    psi_phi[1] = 1;
    psi_tau[1] = sqrt(psi_W);
  }

  // log detection variance partitions
  row_vector[mu_V] mu_phi, mu_tau;
  if (mu_V > 1) {
    mu_phi = dirichlet ?
             sum_to_zero_simplex_jacobian(mu_phi_z[1])'
             : softmax(mu_theta[1] * mu_phi_z[1])';
    mu_tau = sqrt(mu_W * mu_phi);
  } else {
    mu_phi[1] = 1;
    mu_tau[1] = sqrt(mu_W);
  }

  // bivariate normal intercepts
  matrix[S, 2] alpha = rep_matrix([ logit(psi_bar), log(mu_bar) ], S)
                       + append_col(alpha_z[1], alpha_z[2])
                         * diag_post_multiply(alpha_O_L',
                                              [ psi_tau[1], mu_tau[1] ]);

  // mean and species-level continuous predictor coefficients
  row_vector[P[1]] psi_beta_bar;
  matrix[S, P[1]] psi_beta;
  row_vector[P[2]] mu_beta_bar;
  matrix[S, P[2]] mu_beta;
  row_vector[P[3]] gamma_bar;
  matrix[S, P[3]] gamma;

  // mean and species-level categorical predictor coefficients
  matrix[P_cat[1], C_max[1]] psi_beta_cat_bar_z, psi_beta_cat_bar;
  array[P_cat[1]] matrix[S, C_max[1]] psi_beta_cat_z, psi_beta_cat;
  matrix[P_cat[2], C_max[2]] mu_beta_cat_bar_z, mu_beta_cat_bar;
  array[P_cat[2]] matrix[S, C_max[2]] mu_beta_cat_z, mu_beta_cat;
  matrix[P_cat[3], C_max[3]] gamma_cat_bar_z, gamma_cat_bar;
  array[P_cat[3]] matrix[S, C_max[3]] gamma_cat_z, gamma_cat;

  // mean and species-level ordinal predictor coefficients
  row_vector[P_ord[1]] psi_beta_ord_bar;
  matrix[S, P_ord[1]] psi_beta_ord;
  matrix[P_ord[1], O_max[1]] psi_beta_ord_cs;
  row_vector[P_ord[2]] mu_beta_ord_bar;
  matrix[S, P_ord[2]] mu_beta_ord;
  matrix[P_ord[2], O_max[2]] mu_beta_ord_cs;
  row_vector[P_ord[3]] gamma_ord_bar;
  matrix[S, P_ord[3]] gamma_ord;
  matrix[P_ord[3], O_max[3]] gamma_ord_cs;

  // species-level conditional and unconditional site and survey effects
  row_vector[SP * I] iota_bar;
  row_vector[SP * P_sum[2] ? I : 0] iota_bar_unc;
  matrix[S, SP * I] iota;
  matrix[S, SP * P_sum[2] ? I : 0] iota_unc;
  row_vector[TE * J] kappa_bar;
  row_vector[TE * P_sum[3] ? J : 0] kappa_bar_unc;
  matrix[S, TE * J] kappa;
  matrix[S, TE * P_sum[3] ? J : 0] kappa_unc;
  {
    int psi_idx = 2, mu_idx = 2;

    // occupancy continuous coefficients
    if (P[1]) {
      tuple(row_vector[P[1]], matrix[S, P[1]]) coef =
        coef_cont(segment(psi_tau, psi_idx, 2 * P[1]), psi_beta_O_L,
                  psi_beta_bar_z, psi_beta_z);
      psi_beta_bar = coef.1;
      psi_beta = coef.2;
      psi_idx += 2 * P[1];
    }

    // occupancy categorical coefficients
    if (P_cat[1]) {
      int P_c = P_cat[1], C_m = C_max[1];
      tuple(matrix[P_c, C_m], matrix[P_c, C_m], array[P_c] matrix[S, C_m],
            array[P_c] matrix[S, C_m]) coef =
        coef_cat_jacobian(segment(psi_tau, psi_idx, 2 * P_c), psi_beta_O_L,
                          psi_beta_cat_bar_u, psi_beta_cat_u, C[1, :P_c]);
      psi_beta_cat_bar = coef.1;
      psi_beta_cat_bar_z = coef.2;
      psi_beta_cat = coef.3;
      psi_beta_cat_z = coef.4;
      psi_idx += 2 * P_c;
    }

    // occupancy ordinal coefficients
    if (P_ord[1]) {
      int P_o = P_ord[1];
      tuple(row_vector[P_o], matrix[S, P_o], matrix[P_o, O_max[1]]) coef =
        coef_ord_jacobian(segment(psi_tau, psi_idx, 2 * P_o), psi_beta_O_L,
                          psi_beta_ord_bar_z, psi_beta_ord_z, psi_beta_ord_c_u,
                          O[1, :P_o]);
      psi_beta_ord_bar = coef.1;
      psi_beta_ord = coef.2;
      psi_beta_ord_cs = coef.3;
      psi_idx += 2 * P_o;
    }

    // detection continuous coefficients (site-level)
    if (P[2]) {
      tuple(row_vector[P[2]], matrix[S, P[2]]) coef =
        coef_cont(segment(mu_tau, mu_idx, 2 * P[2]), mu_beta_O_L, mu_beta_bar_z,
                  mu_beta_z);
      mu_beta_bar = coef.1;
      mu_beta = coef.2;
      mu_idx += 2 * P[2];
    }

    // detection categorical coefficients (site-level)
    if (P_cat[2]) {
      int P_c = P_cat[2], C_m = C_max[2];
      tuple(matrix[P_c, C_m], matrix[P_c, C_m],
            array[P_c] matrix[S, C_m], array[P_c] matrix[S, C_m]) coef =
        coef_cat_jacobian(segment(mu_tau, mu_idx, 2 * P_c), mu_beta_O_L,
                          mu_beta_cat_bar_u, mu_beta_cat_u, C[2, :P_c]);
      mu_beta_cat_bar = coef.1;
      mu_beta_cat_bar_z = coef.2;
      mu_beta_cat = coef.3;
      mu_beta_cat_z = coef.4;
      mu_idx += 2 * P_c;
    }

    // detection ordinal coefficients (site-level)
    if (P_ord[2]) {
      int P_o = P_ord[2];
      tuple(row_vector[P_o], matrix[S, P_o], matrix[P_o, O_max[2]]) coef =
        coef_ord_jacobian(segment(mu_tau, mu_idx, 2 * P_o), mu_beta_O_L,
                          mu_beta_ord_bar_z, mu_beta_ord_z, mu_beta_ord_c_u,
                          O[2, :P_o]);
      mu_beta_ord_bar = coef.1;
      mu_beta_ord = coef.2;
      mu_beta_ord_cs = coef.3;
      mu_idx += 2 * P_o;
    }

    // detection continuous coefficients (site-by-survey level)
    if (P[3]) {
      tuple(row_vector[P[3]], matrix[S, P[3]]) coef =
        coef_cont(segment(mu_tau, mu_idx, 2 * P[3]), gamma_O_L, gamma_bar_z,
                  gamma_z);
      gamma_bar = coef.1;
      gamma = coef.2;
      mu_idx += 2 * P[3];
    }

    // detection categorical coefficients (site-by-survey level)
    if (P_cat[3]) {
      int P_c = P_cat[3], C_m = C_max[3];
      tuple(matrix[P_c, C_m], matrix[P_c, C_m],
            array[P_c] matrix[S, C_m], array[P_c] matrix[S, C_m]) coef =
        coef_cat_jacobian(segment(mu_tau, mu_idx, 2 * P_c), gamma_O_L,
                          gamma_cat_bar_u, gamma_cat_u, C[3, :P_c]);
      gamma_cat_bar = coef.1;
      gamma_cat_bar_z = coef.2;
      gamma_cat = coef.3;
      gamma_cat_z = coef.4;
      mu_idx += 2 * P_c;
    }

    // detection ordinal coefficients (site-by-survey level)
    if (P_ord[3]) {
      int P_o = P_ord[3];
      tuple(row_vector[P_o], matrix[S, P_o], matrix[P_o, O_max[3]]) coef =
        coef_ord_jacobian(segment(mu_tau, mu_idx, 2 * P_o), gamma_O_L,
                          gamma_ord_bar_z, gamma_ord_z, gamma_ord_c_u,
                          O[3, :P_o]);
      gamma_ord_bar = coef.1;
      gamma_ord = coef.2;
      gamma_ord_cs = coef.3;
      mu_idx += 2 * P_o;
    }

    // mean and species-specific site effects
    if (SP) {
      iota_bar = mu_tau[mu_idx] * iota_bar_z[1]';
      matrix[S, I] iota_s = diag_pre_multiply(segment(mu_tau, mu_idx + 1, S),
                                              iota_O_L) * iota_z[1];

      // increment spatial GP
      if (I_GP) {
        matrix[I, I] iota_K = gp_exp_quad_cov(XY, 1, iota_ell[1]),
                     iota_U = cholesky_decompose(add_diag(iota_K, 1e-9))';
        iota_bar *= iota_U;
        if (SS) {
          for (s in 1:S) {
            iota_K = gp_exp_quad_cov(XY, 1, iota_ell[1 + s]);
            iota_U = cholesky_decompose(add_diag(iota_K, 1e-9))';
            iota_s[s] *= iota_U;
          }
        } else {
          iota_s *= iota_U;
        }
      }
      iota = rep_matrix(iota_bar, S) + iota_s;

      // increment orthogonal projection
      if (P_sum[2]) {
        iota_bar_unc = iota_bar;
        iota_bar = orthogonalise(iota_bar_unc, X2_aug, X2_plus);
        iota_unc = iota;
        iota = orthogonalise(iota_unc, X2_aug, X2_plus);
      }
      mu_idx += Sp1;
    }

    // mean and species-specific survey effects
    if (TE) {
      kappa_bar = mu_tau[mu_idx] * kappa_bar_z[1]';
      matrix[S, J] kappa_s = diag_pre_multiply(segment(mu_tau, mu_idx + 1, S),
                                               kappa_O_L) * kappa_z[1];

      // increment temporal GP
      if (J_GP) {
        matrix[J, J] kappa_K, kappa_U;
        vector[periodic * 2] kappa_t;
        if (periodic) {
          kappa_t = sqrt(kappa_v[1]);
          kappa_K = gp_exp_quad_cov(surveys, kappa_t[1], kappa_ell[1, 1])
                    + gp_periodic_cov(surveys, kappa_t[2], kappa_ell[1, 2],
                                      period);
        } else {
          kappa_K = gp_exp_quad_cov(surveys, 1, kappa_ell[1, 1]);
        }
        kappa_U = cholesky_decompose(add_diag(kappa_K, 1e-9))';
        kappa_bar *= kappa_U;
        if (SS) {
          for (s in 1:S) {
            int sp1 = s + 1;
            if (periodic) {
              kappa_K = gp_exp_quad_cov(surveys, kappa_t[1], kappa_ell[sp1, 1])
                        + gp_periodic_cov(surveys, kappa_t[2],
                                          kappa_ell[sp1, 2], period);
            } else {
              kappa_K = gp_exp_quad_cov(surveys, 1, kappa_ell[sp1, 1]);
            }
            kappa_U = cholesky_decompose(add_diag(kappa_K, 1e-9))';
            kappa_s[s] *= kappa_U;
          }
        } else {
          kappa_s *= kappa_U;
        }
      }
      kappa = rep_matrix(kappa_bar, S) + kappa_s;

      // increment orthogonal projection
      if (P_sum[3]) {
        kappa_bar_unc = kappa_bar;
        kappa_bar = orthogonalise(kappa_bar_unc, X3_mean, X3_plus);
        kappa_unc = kappa;
        kappa = orthogonalise(kappa_unc, X3_mean, X3_plus);
      }
    }
  }

  // occupancy
  matrix[S, I] logit_psi = rep_matrix(alpha[:, 1], I) + psi_beta * X1;
  for (p in 1:P_cat[1]) {
    logit_psi += psi_beta_cat[p, :, X_cat1[p]];
  }
  if (P_ord[1]) {
    logit_psi += increment_ordinal(X_ord1, coef_ord_realised(O[1],
                                   psi_beta_ord, psi_beta_ord_cs));
  }

  // priors
  real lprior = beta_lpdf(psi_bar | psi_bar_beta[1], psi_bar_beta[2])
                + gamma_lpdf(mu_bar | mu_bar_gamma[1], mu_bar_gamma[2])
                + student_t_lpdf(psi_W | psi_W_t[1], psi_W_t[2], psi_W_t[3])
                + student_t_lpdf(mu_W | mu_W_t[1], mu_W_t[2], mu_W_t[3])
                + lkj_corr_cholesky_lpdf(alpha_O_L | alpha_O_L_LKJ);
  if (psi_V > 1) {
    lprior += gamma_lpdf(psi_theta[1] | psi_theta_gamma[1], psi_theta_gamma[2]);
  }
  if (mu_V > 1) {
    lprior += gamma_lpdf(mu_theta[1] | mu_theta_gamma[1], mu_theta_gamma[2]);
  }
  if (P_sum[1]) {
    lprior += lkj_corr_cholesky_lpdf(psi_beta_O_L | psi_beta_O_L_LKJ);
  }
  if (P_sum[2]) {
    lprior += lkj_corr_cholesky_lpdf(mu_beta_O_L | mu_beta_O_L_LKJ);
  }
  if (P_sum[3]) {
    lprior += lkj_corr_cholesky_lpdf(gamma_O_L | gamma_O_L_LKJ);
  }
  if (SP) {
    lprior += lkj_corr_cholesky_lpdf(iota_O_L | iota_O_L_LKJ);
    if (I_GP) {
      lprior += inv_gamma_lpdf(iota_ell | iota_ell_inv_gamma[1],
                                          iota_ell_inv_gamma[2]);
    }
  }
  if (TE) {
    lprior += lkj_corr_cholesky_lpdf(kappa_O_L | kappa_O_L_LKJ);
    if (J_GP) {
      lprior += inv_gamma_lpdf(kappa_ell[:, 1] | kappa_ell_inv_gamma[1],
                                                 kappa_ell_inv_gamma[2]);
      if (periodic) {
        lprior += inv_gamma_lpdf(kappa_ell[:, 2] |
                                   kappa_ell_periodic_inv_gamma[1],
                                   kappa_ell_periodic_inv_gamma[2])
                  + dirichlet_lpdf(kappa_v[1] | kappa_v_dirichlet);
      }
    }
  }
  if (OLRE) {
    lprior += lkj_corr_cholesky_lpdf(epsilon_O_L | epsilon_O_L_LKJ);
  } else if (NB) {
    lprior += inv_gamma_lpdf(phi | phi_inv_gamma[1], phi_inv_gamma[2]);
  }
}

model {
  target += lprior
            + std_normal_lupdf(append_row(alpha_z[1], alpha_z[2]))
            + std_normal_lupdf(psi_beta_bar_z)
            + std_normal_lupdf(psi_beta_ord_bar_z)
            + std_normal_lupdf(mu_beta_bar_z)
            + std_normal_lupdf(mu_beta_ord_bar_z)
            + std_normal_lupdf(gamma_bar_z)
            + std_normal_lupdf(gamma_ord_bar_z);
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
  for (p in 1:P[1]) {
    target += std_normal_lupdf(psi_beta_z[p]);
  }
  for (p in 1:P_cat[1]) {
    int C_p = C[1, p];
    target += std_normal_lupdf(psi_beta_cat_bar_z[p, :C_p]) +
              std_normal_lupdf(to_vector(psi_beta_cat_z[p, :, :C_p]));
  }
  for (p in 1:P_ord[1]) {
    target += std_normal_lupdf(psi_beta_ord_z[p]);
  }
  for (p in 1:P[2]) {
    target += std_normal_lupdf(mu_beta_z[p]);
  }
  for (p in 1:P_cat[2]) {
    int C_p = C[2, p];
    target += std_normal_lupdf(mu_beta_cat_bar_z[p, :C_p]) +
              std_normal_lupdf(to_vector(mu_beta_cat_z[p, :, :C_p]));
  }
  for (p in 1:P_ord[2]) {
    target += std_normal_lupdf(mu_beta_ord_z[p]);
  }
  for (p in 1:P[3]) {
    target += std_normal_lupdf(gamma_z[p]);
  }
  for (p in 1:P_cat[3]) {
    int C_p = C[3, p];
    target += std_normal_lupdf(gamma_cat_bar_z[p, :C_p]) +
              std_normal_lupdf(to_vector(gamma_cat_z[p, :, :C_p]));
  }
  for (p in 1:P_ord[3]) {
    target += std_normal_lupdf(gamma_ord_z[p]);
  }
  if (SP) {
    target += std_normal_lupdf(iota_bar_z[1])
              + std_normal_lupdf(to_vector(iota_z[1]));
  }
  if (TE) {
    target += std_normal_lupdf(kappa_bar_z[1])
              + std_normal_lupdf(to_vector(kappa_z[1]));
  }

  // Poisson OLREs
  array[OLRE * I] matrix[S, J] epsilon;
  if (OLRE) {
    target += std_normal_lupdf(epsilon_bar_z[1])
              + std_normal_lupdf(to_vector(epsilon_z[1]));
    row_vector[N] epsilon_bar = mu_tau[mu_V - S] * epsilon_bar_z[1]';
    matrix[S, N] epsilon_mat = rep_matrix(epsilon_bar, S)
                               + diag_pre_multiply(tail(mu_tau, S), epsilon_O_L)
                                 * epsilon_z[1];
    epsilon = fill_epsilon(f_l, log_Delta, epsilon_mat);
  }

  // likelihood
  target += grainsize ?
            reduce_sum(partial_occARU_lupmf, sites, grainsize, y, Q, f_l,
                       log_Delta, X2, X_cat2, X_ord2, X3, X_cat3, X_ord3,
                       logit_psi, alpha[:, 2], mu_beta, mu_beta_cat,
                       coef_ord_realised(O[2], mu_beta_ord, mu_beta_ord_cs),
                       gamma, gamma_cat, coef_ord_realised(O[3], gamma_ord,
                       gamma_ord_cs), iota, kappa, epsilon, phi)
            : occARU_lupmf(y | Q, f_l, log_Delta, X2, X_cat2, X_ord2, X3,
                           X_cat3, X_ord3, logit_psi, alpha[:, 2], mu_beta,
                           mu_beta_cat, coef_ord_realised(O[2], mu_beta_ord,
                           mu_beta_ord_cs), gamma, gamma_cat,
                           coef_ord_realised(O[3], gamma_ord, gamma_ord_cs),
                           iota, kappa, epsilon, phi);
}

generated quantities {
  // correlations
  real alpha_rho = multiply_lower_tri_self_transpose(alpha_O_L)[1, 2];
  corr_matrix[P_sum[1] ? S : 0] psi_beta_O =
    multiply_lower_tri_self_transpose(psi_beta_O_L);
  corr_matrix[P_sum[2] ? S : 0] mu_beta_O =
    multiply_lower_tri_self_transpose(mu_beta_O_L);
  corr_matrix[P_sum[3] ? S : 0] gamma_O =
    multiply_lower_tri_self_transpose(gamma_O_L);
  corr_matrix[SP * S] iota_O = multiply_lower_tri_self_transpose(iota_O_L);
  corr_matrix[TE * S] kappa_O = multiply_lower_tri_self_transpose(kappa_O_L);
  corr_matrix[OLRE * S] epsilon_O =
    multiply_lower_tri_self_transpose(epsilon_O_L);

  // unconditional site coefficients
  row_vector[SP * P[2]] mu_beta_bar_unc;
  matrix[S, SP * P[2]] mu_beta_unc;
  matrix[SP * P_cat[2], C_max[2]] mu_beta_cat_bar_unc;
  array[SP * P_cat[2]] matrix[S, C_max[2]] mu_beta_cat_unc;
  row_vector[SP * P_ord[2]] mu_beta_ord_bar_unc;
  matrix[S, SP * P_ord[2]] mu_beta_ord_unc;
  if (SP) {
    if (P[2]) {
      mu_beta_bar_unc = mu_beta_bar - iota_bar_unc * X2_plus[:, :P[2]];
      mu_beta_unc = mu_beta - iota_unc * X2_plus[:, :P[2]];
    }
    if (P_cat[2]) {
      int idx = P[2];
      for (p in 1:P_cat[2]) {
        int C_p = C[2, p];
        array[C_p] int seq = linspaced_int_array(C_p, idx + 1, idx + C_p);
        mu_beta_cat_bar_unc[p, :C_p] = mu_beta_cat_bar[p, :C_p]
                                       - iota_bar_unc * X2_plus[:, seq];
        mu_beta_cat_unc[p, :, :C_p] = mu_beta_cat[p, :, :C_p]
                                      - iota_unc * X2_plus[:, seq];
        idx += C_p;
      }
    }
    if (P_ord[2]) {
      int idx = X_plus_cols[1] - P_ord[2] + 1;
      mu_beta_ord_bar_unc = mu_beta_ord_bar - iota_bar_unc * X2_plus[:, idx:];
      mu_beta_ord_unc = mu_beta_ord - iota_unc * X2_plus[:, idx:];
    }
  }

  // unconditional survey coefficients
  row_vector[TE * P[3]] gamma_bar_unc;
  matrix[S, TE * P[3]] gamma_unc;
  matrix[TE * P_cat[3], C_max[3]] gamma_cat_bar_unc;
  array[TE * P_cat[3]] matrix[S, C_max[3]] gamma_cat_unc;
  row_vector[TE * P_ord[3]] gamma_ord_bar_unc;
  matrix[S, TE * P_ord[3]] gamma_ord_unc;
  if (TE) {
    if (P[3]) {
      gamma_bar_unc = gamma_bar - kappa_bar_unc * X3_plus[:, :P[3]];
      gamma_unc = gamma - kappa_unc * X3_plus[:, :P[3]];
    }
    if (P_cat[3]) {
      int idx = P[3];
      for (p in 1:P_cat[3]) {
        int C_p = C[3, p];
        array[C_p] int seq = linspaced_int_array(C_p, idx + 1, idx + C_p);
        gamma_cat_bar_unc[p, :C_p] = gamma_cat_bar[p, :C_p]
                                     - kappa_bar_unc * X3_plus[:, seq];
        gamma_cat_unc[p, :, :C_p] = gamma_cat[p, :, :C_p]
                                    - kappa_unc * X3_plus[:, seq];
        idx += C_p;
      }
    }
    if (P_ord[3]) {
      int idx = X_plus_cols[2] - P_ord[3] + 1;
      gamma_ord_bar_unc = gamma_ord_bar - kappa_bar_unc * X3_plus[:, idx:];
      gamma_ord_unc = gamma_ord - kappa_unc * X3_plus[:, idx:];
    }
  }

  // log likelihood, latent occupancy, and posterior predictions
  matrix[S, I] log_lik;
  matrix[MC * S, I] log_lik2;
  array[latent * I, S] int z = rep_array(1, latent * I, S);
  array[PPC_y * I, J, S] int yrep = rep_array(0, PPC_y * I, J, S);
  array[PPC_Q * I, S] int Qrep = rep_array(0, PPC_Q * I, S);

  {
    // reconstruct log likelihood
    row_vector[N] epsilon_bar;
    matrix[S, S] epsilon_S_L;
    array[OLRE * I] matrix[S, J] epsilon;
    if (OLRE) {
      epsilon_bar = mu_tau[mu_V - S] * epsilon_bar_z[1]';
      epsilon_S_L = diag_pre_multiply(tail(mu_tau, S), epsilon_O_L);
      matrix[S, N] epsilon_mat = rep_matrix(epsilon_bar, S)
                                 + epsilon_S_L * epsilon_z[1];
      epsilon = fill_epsilon(f_l, log_Delta, epsilon_mat);
    }
    tuple(matrix[S, I], array[I] matrix[S, 2], array[I] matrix[S, J]) lp =
      occARU(y, Q, f_l, log_Delta, X2, X_cat2, X_ord2, X3, X_cat3, X_ord3,
             logit_psi, alpha[:, 2], mu_beta, mu_beta_cat,
             coef_ord_realised(O[2], mu_beta_ord, mu_beta_ord_cs), gamma,
             gamma_cat, coef_ord_realised(O[3], gamma_ord, gamma_ord_cs), iota,
             kappa, epsilon, phi);
    log_lik = lp.1;

    // recover latent occupancy
    if (latent) {
      for (i in 1:I) {
        for (s in 1:S) {
          if (!Q[i, s]) {
            z[i, s] = bernoulli_logit_rng(lp.2[i, s, 2] - log_lik[s, i]);
          }
        }
      }
    }

    // posterior predictions
    if (PPC_y || PPC_Q) {
      if (OLRE) {
        epsilon = pp_epsilon_rng(f_l, log_Delta, epsilon_bar, epsilon_S_L);
      }
      tuple(array[I, J, S] int, array[I, S] int) pp =
        pp_occARU_rng(f_l, log_Delta, logit_psi, lp.3, epsilon, phi);
      if (PPC_y) {
        yrep = pp.1;
      }
      if (PPC_Q) {
        Qrep = pp.2;
      }
    }

    // Monte Carlo integration of observation-level parameters for loo
    if (MC) {
      array[D] matrix[S, I] log_lik_k;

      // prepare site effect covariances
      int mu_idx = 2 + 2 * sum(P_sum[2:3]) + 1;
      array[I_GP * (SS ? S : 1)] matrix[I, I] iota_U;
      matrix[S, S] iota_S_L;
      if (SP) {
        iota_S_L = diag_pre_multiply(segment(mu_tau, mu_idx, S), iota_O_L);
        if (I_GP) {
          if (SS) {
            for (s in 1:S) {
              iota_U[s] =
                cholesky_decompose(add_diag(gp_exp_quad_cov(XY, 1,
                                            iota_ell[1 + s]), 1e-9))';
            }
          } else {
            iota_U[1] =
              cholesky_decompose(add_diag(gp_exp_quad_cov(XY, 1, iota_ell[1]),
                                 1e-9))';
          }
        }
      }

      // posterior predict site effects and OLREs and compute log likelihood
      int SIm1 = Sm1 * Im1;
      array[SIm1] int zeros = zeros_int_array(SIm1);
      matrix[Sm1, SP * Im1] iota_u;
      matrix[S, SP * I] iota_s, iota_rep;
      for (d in 1:D) {
        if (SP) {
          iota_u = to_matrix(normal_rng(zeros, 1), Sm1, Im1);
          iota_s = iota_S_L * sum_to_zero_constrain(iota_u);
          if (I_GP) {
            if (SS) {
              for (s in 1:S) {
                iota_s[s] *= iota_U[s];
              }
            } else {
              iota_s *= iota_U[1];
            }
          }
          if (P_sum[2]) {
            iota_s = orthogonalise(iota_s, X2_aug, X2_plus);
          }
          iota_rep = rep_matrix(iota_bar, S) + iota_s;
        }
        if (OLRE) {
          epsilon = pp_epsilon_rng(f_l, log_Delta, epsilon_bar, epsilon_S_L);
        }
        lp = occARU(y, Q, f_l, log_Delta, X2, X_cat2, X_ord2, X3, X_cat3,
                    X_ord3, logit_psi, alpha[:, 2], mu_beta, mu_beta_cat,
                    coef_ord_realised(O[2], mu_beta_ord, mu_beta_ord_cs),
                    gamma, gamma_cat, coef_ord_realised(O[3], gamma_ord,
                    gamma_ord_cs), iota_rep, kappa, epsilon, phi);
        log_lik_k[d] = lp.1;
      }
      log_lik2 = rep_matrix(-log_D, S, I);
      for (i in 1:I) {
        for (s in 1:S) {
          log_lik2[s, i] += log_sum_exp(log_lik_k[:, s, i]);
        }
      }
    }
  }
}
