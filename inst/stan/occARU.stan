functions {
  #include occARU.stanfunctions
  #include partial.stanfunctions
  #include util.stanfunctions
  #include predictors.stanfunctions
  #include coefficients.stanfunctions
  #include kernels.stanfunctions
  #include predict.stanfunctions
}

data {
  // dimensions, recording effort, and detection history
  int<lower=2> I, J;
  int<lower=1> S, R, K;
  matrix<lower=0, upper=1>[I, J] Delta;
  array[I, J, K, S] int<lower=0> y;

  // site UTMs (km) and regions (assumed ordered)
  array[I] vector[2] XY;
  array[I] int<lower=1, upper=R> region;

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

  // random effects
  array[2] int<lower=0, upper=1> random,  // indicators
                                 SS,  // species length scales
                                 project;  // orthogonal projection
  array[2] int<lower=0, upper=3> kernel;  // GP kernel type
  real<lower=0> period;  // period for periodic kernel

  // indicators
  int<lower=0, upper=2> OD;  // overdispersion (Poisson, OLRE, negbin)
  int<lower=0, upper=1> dirichlet,  // var. decomp. (logi-norm or Dirichlet)
                        latent,  // recover latent occ. states
                        PPC_y,  // posterior predictions of y
                        PPC_Q;  // posterior predictions of Q
  int<lower=0, upper=I> grainsize;  // grainsize for threading
  int<lower=0> D;  // Monte Carlo draws for LOO

  // prior hyperparameters
  vector<lower=0>[2] psi_bar_beta,  // mean occupancy (probability scale)
                     mu_bar_gamma,  // mean detetection (rate scale)
                     psi_theta_gamma,  // occupancy variance partition scale
                     mu_theta_gamma,  // detection variance partition scale
                     iota_ell_inv_gamma,  // site GP length scale
                     kappa_ell_inv_gamma,  // survey GP length scale
                     kappa_ell_periodic_inv_gamma,  // periodic GP length scale
                     K_phi_dirichlet,  // survey GP variance partitions
                     phi_inv_gamma;  // negative binomial overdispersion
  vector<lower=0>[3] psi_W_t,  // occupancy log odds variance
                     mu_W_t;  // log detection rate variance
  real<lower=0> alpha_O_L_LKJ, psi_beta_O_L_LKJ, mu_beta_O_L_LKJ, gamma_O_L_LKJ,
                iota_O_L_LKJ, kappa_O_L_LKJ, epsilon_O_L_LKJ;
}

transformed data {
  // transformed indicators
  int MS = S > 1, MR = R > 1, periodic = kernel[2] > 0 && period > 0,
      OLRE = OD == 1, NB = OD == 2, MC = (D > 0) && (random[1] || OLRE),
      Im1 = I - 1, Sm1 = S - 1;
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

  // number of sites per region
  array[R] int I_r = zeros_int_array(R);
  for (i in 1:I) {
    I_r[region[i]] += 1;
  }
  int I_max = max(I_r);

  // aggregated counts
  array[I, S] int Q;
  for (i in 1:I) {
    int J_ii = J_i[i];
    array[J_ii] int idx = j_idx[i, :J_ii];
    for (s in 1:S) {
      Q[i, s] = sum(y[i, idx, 1, s]);
    }
  }

  // sequences
  array[I] int sites = linspaced_int_array(I, 1, I);
  array[R] int regions = linspaced_int_array(R, 1, R);
  array[J] real surveys = linspaced_array(J, 1, J);

  // variance partitions
  int psi_V = MS + (1 + MS) * (MR + P_sum[1]),
      mu_V = MS + (1 + MS) * (MR + sum(P_sum[2:3]) + sum(random) + OLRE);

  // categorical levels
  tuple(array[3, max(P_cat)] int, array[3] int, array[3] int) cat_levels =
    levels(X_cat1, X_cat2, X_cat3);
  array[3, max(P_cat)] int C = cat_levels.1;
  array[3] int C_max = cat_levels.2, C_sum = cat_levels.3;

  // ordinal levels
  tuple(array[3, max(P_ord)] int, array[3] int, array[3] int) ord_levels =
    levels(X_ord1, X_ord2, X_ord3);
  array[3, max(P_ord)] int O = ord_levels.1;
  array[3] int O_max = ord_levels.2, O_sum = ord_levels.3;

  // add 1 to ordinal arrays for faster incrementation
  array[I, P_ord[1]] int X_ord1_p1 = X_ord1;
  array[I, P_ord[2]] int X_ord2_p1 = X_ord2;
  array[I, J, P_ord[3]] int X_ord3_p1 = X_ord3;
  for (i in 1:I) {
    for (p in 1:P_ord[1]) {
      X_ord1_p1[i, p] += 1;
    }
    for (p in 1:P_ord[2]) {
      X_ord2_p1[i, p] += 1;
    }
    for (p in 1:P_ord[3]) {
      for (j in 1:J) {
        X_ord3_p1[i, j, p] += 1;
      }
    }
  }

  // augmented site-level design matrix and pseudo-inverse for projection
  array[3] int P_aug;
  for (d in 2:3) {
    P_aug[d] = P[d] + C_sum[d] - P_cat[d] + P_ord[d];
  }
  matrix[I, project[1] * P_aug[2]] X2_aug;
  matrix[project[1] * P_aug[2], I] X2_plus;
  if (project[1]) {
    X2_aug = augment_design_matrix(X2, X_cat2, X_ord2, C[2, :P_cat[2]],
                                   O[2, :P_ord[2]]);
    X2_plus = pseudo_inverse(X2_aug);
  }

  // site-by-survey design matrix pseudo-inverses for orthogonal projection
  array[project[2] * I] matrix[J, P_aug[3]] X3_aug;
  matrix[J, project[2] * P_aug[3]] X3_mean;
  matrix[project[2] * P_aug[3], J] X3_plus;
  if (project[2]) {
    X3_aug = augment_design_matrix(X3, X_cat3, X_ord3, C[3, :P_cat[3]],
                                   O[3, :P_ord[3]]);
    X3_mean = average_survey_design_matrix(X3_aug, j_idx);
    X3_plus = pseudo_inverse(X3_mean);
  }
}

parameters {
  // intercepts by species and region
  real<lower=0, upper=1> psi_bar;
  real<lower=0> mu_bar;
  cholesky_factor_corr[MS * 2] alpha_O_L;
  array[MS * 2] sum_to_zero_vector[S] alpha_s_bar_z;
  array[MR * 2] sum_to_zero_vector[R] alpha_r_bar_z;
  array[MS * MR * 2] sum_to_zero_matrix[R, S] alpha_r_s_z;

  // total variances, partition sparsity, and unconstrained variance partitions
  vector<lower=0>[psi_V > 0] psi_W;
  vector<lower=0>[psi_V > 1] psi_theta;
  array[psi_V > 1] sum_to_zero_vector[psi_V > 1 ? psi_V : 1] psi_phi_z;
  vector<lower=0>[mu_V > 0] mu_W;
  vector<lower=0>[mu_V > 1] mu_theta;
  array[mu_V > 1] sum_to_zero_vector[mu_V > 1 ? mu_V : 1] mu_phi_z;

  // correlation matrices
  cholesky_factor_corr[MS && P_sum[1] ? S : 0] psi_beta_O_L;
  cholesky_factor_corr[MS && P_sum[2] ? S : 0] mu_beta_O_L;
  cholesky_factor_corr[MS && P_sum[3] ? S : 0] gamma_O_L;
  cholesky_factor_corr[MS && random[1] ? S : 0] iota_O_L;
  cholesky_factor_corr[MS && random[2] ? S : 0] kappa_O_L;
  cholesky_factor_corr[MS && OLRE ? S : 0] epsilon_O_L;

  // coefficients for continuous predictors
  vector[P[1]] psi_beta_bar_z;
  array[MS * P[1]] sum_to_zero_vector[S] psi_beta_s_z;
  vector[P[2]] mu_beta_bar_z;
  array[MS * P[2]] sum_to_zero_vector[S] mu_beta_s_z;
  vector[P[3]] gamma_bar_z;
  array[MS * P[3]] sum_to_zero_vector[S] gamma_s_z;

  // unconstrained coefficients for categorical predictors
  vector[C_sum[1] - P_cat[1]] psi_beta_cat_bar_u;
  matrix[C_sum[1] - P_cat[1], MS * Sm1] psi_beta_cat_s_u;
  vector[C_sum[2] - P_cat[2]] mu_beta_cat_bar_u;
  matrix[C_sum[2] - P_cat[2], MS * Sm1] mu_beta_cat_s_u;
  vector[C_sum[3] - P_cat[3]] gamma_cat_bar_u;
  matrix[C_sum[3] - P_cat[3], MS * Sm1] gamma_cat_s_u;

  // coefficients for ordinal predictors
  vector[P_ord[1]] psi_beta_ord_bar_z;
  array[P_ord[1]] sum_to_zero_vector[S] psi_beta_ord_s_z;
  vector[O_sum[1] - P_ord[1]] psi_beta_ord_c_u;
  vector[P_ord[2]] mu_beta_ord_bar_z;
  array[P_ord[2]] sum_to_zero_vector[S] mu_beta_ord_s_z;
  vector[O_sum[2] - P_ord[2]] mu_beta_ord_c_u;
  vector[P_ord[3]] gamma_ord_bar_z;
  array[P_ord[3]] sum_to_zero_vector[S] gamma_ord_s_z;
  vector[O_sum[3] - P_ord[3]] gamma_ord_c_u;

  // site effects
  array[random[1]] sum_to_zero_vector[I] iota_bar_z;
  array[MS * random[1]] sum_to_zero_matrix[I, S] iota_s_z;
  array[MS * random[1]] simplex[S] iota_phi;
  vector<lower=0>[(kernel[1] > 0) * (1 + SS[1] * S)] iota_ell;

  // survey effects
  array[random[2]] sum_to_zero_vector[J] kappa_bar_z;
  array[MS * random[2]] sum_to_zero_matrix[J, S] kappa_s_z;
  array[MS * random[2]] simplex[S] kappa_phi;
  matrix<lower=0>[1 + (SS[2] * S), (kernel[2] > 0) + periodic] kappa_ell;
  array[periodic] simplex[2] K_phi;

  // negative binomial overdispersion or OLRE residuals
  vector<lower=0>[NB * S] phi;
  array[OLRE] sum_to_zero_vector[J_sum] epsilon_bar_z;
  array[MS * OLRE] sum_to_zero_matrix[J_sum, S] epsilon_s_z;
  array[MS * OLRE] simplex[S] epsilon_phi;
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
  vector[2] alpha_bar = [ logit(psi_bar), log(mu_bar) ]';
  matrix[2, MS * S] alpha_s_bar;
  matrix[R, MR * 2] alpha_r_bar;
  array[MS * MR * 2] matrix[R, S] alpha_r_s;
  array[2] matrix[R, S] alpha;
  for (d in 1:2) {
    alpha[d] = rep_matrix(alpha_bar[d], R, S);
  }

  // continuous predictor coefficients
  vector[P[1]] psi_beta_bar;
  matrix[P[1], MS * S] psi_beta_s;
  matrix[P[1], S] psi_beta;
  vector[P[2]] mu_beta_bar;
  matrix[P[2], MS * S] mu_beta_s;
  matrix[P[2], S] mu_beta;
  vector[P[3]] gamma_bar;
  matrix[P[3], MS * S] gamma_s;
  matrix[P[3], S] gamma;

  // categorical predictor coefficients
  matrix[C_max[1], P_cat[1]] psi_beta_cat_bar_z, psi_beta_cat_bar;
  array[MS * P_cat[1]] matrix[C_max[1], S] psi_beta_cat_s_z, psi_beta_cat_s;
  array[P_cat[1]] matrix[C_max[1], S] psi_beta_cat;
  matrix[C_max[2], P_cat[2]] mu_beta_cat_bar_z, mu_beta_cat_bar;
  array[MS * P_cat[2]] matrix[C_max[2], S] mu_beta_cat_s_z, mu_beta_cat_s;
  array[P_cat[2]] matrix[C_max[2], S] mu_beta_cat;
  matrix[C_max[3], P_cat[3]] gamma_cat_bar_z, gamma_cat_bar;
  array[MS * P_cat[3]] matrix[C_max[3], S] gamma_cat_s_z, gamma_cat_s;
  array[P_cat[3]] matrix[C_max[3], S] gamma_cat;

  // ordinal predictor coefficients
  vector[P_ord[1]] psi_beta_ord_bar;
  matrix[P_ord[1], MS * S] psi_beta_ord_s;
  matrix[O_max[1], P_ord[1]] psi_beta_ord_cs;
  matrix[P_ord[1], S] psi_beta_ord;
  array[P_ord[1]] matrix[O_max[1] + 1, S] psi_beta_ord_realised;
  vector[P_ord[2]] mu_beta_ord_bar;
  matrix[P_ord[2], MS * S] mu_beta_ord_s;
  matrix[O_max[2], P_ord[2]] mu_beta_ord_cs;
  matrix[P_ord[2], S] mu_beta_ord;
  array[P_ord[2]] matrix[O_max[2] + 1, S] mu_beta_ord_realised;
  vector[P_ord[3]] gamma_ord_bar;
  matrix[P_ord[3], MS * S] gamma_ord_s;
  matrix[O_max[3], P_ord[3]] gamma_ord_cs;
  matrix[P_ord[3], S] gamma_ord;
  array[P_ord[3]] matrix[O_max[3] + 1, S] gamma_ord_realised;

  // conditional and unconditional survey effects
  vector[random[2] * J] kappa_bar;
  vector[project[2] ? J : 0] kappa_bar2;
  matrix[random[2] * J, S] kappa;
  matrix[project[2] ? J : 0, S] kappa2;

  // random scales
  vector[random[1]] iota_bar_t;
  vector[MS * random[1] * S] iota_t;
  vector[random[2]] kappa_bar_t;
  vector[MS * random[2] * S] kappa_t;
  vector[OLRE] epsilon_bar_t;
  vector[MS * OLRE * S] epsilon_t;

  profile("transformed_parameters") {

    // intercepts
    int idx = 1;
    if (MS) {
      alpha_s_bar = alpha_O_L * append_col(psi_tau[idx] * alpha_s_bar_z[1],
                                           mu_tau[idx] * alpha_s_bar_z[2])';
      for (d in 1:2) {
        alpha[d] += rep_matrix(alpha_s_bar[d], R);
      }
      idx += 1;
    }
    if (MR) {
      alpha_r_bar = append_col(psi_tau[idx] * alpha_r_bar_z[1],
                               mu_tau[idx] * alpha_r_bar_z[2]);
      idx += 1;
      if (MS) {
        for (d in 1:2) {
          alpha_r_s[d] = (d == 1 ? psi_tau[idx] : mu_tau[idx]) * alpha_r_s_z[d];
          alpha[d] += rep_matrix(alpha_r_bar[:, d], S) + alpha_r_s[d];
        }
        idx += 1;
      } else {
        for (d in 1:2) {
          alpha[d, :, 1] += alpha_r_bar[:, d];
        }
      }
    }
    int psi_idx = idx, mu_idx = idx;

    // occupancy continuous coefficients
    if (P[1]) {
      int N = P[1] * (1 + MS);
      tuple(vector[P[1]], matrix[P[1], MS * S]) coef =
        coef_cont(segment(psi_tau, psi_idx, N), psi_beta_O_L, psi_beta_bar_z,
                  psi_beta_s_z);
      psi_beta_bar = coef.1;
      if (MS) {
        psi_beta_s = coef.2;
        psi_beta = rep_matrix(psi_beta_bar, S) + psi_beta_s;
      } else {
        psi_beta[:, 1] = psi_beta_bar;
      }
      psi_idx += N;
    }

    // occupancy categorical coefficients
    if (P_cat[1]) {
      int P_c = P_cat[1], C_m = C_max[1], N = P_c * (1 + MS);
      tuple(matrix[C_m, P_c], matrix[C_m, P_c], array[P_c] matrix[C_m, MS * S],
            array[P_c] matrix[C_m, MS * S]) coef =
        coef_cat_jacobian(segment(psi_tau, psi_idx, N), psi_beta_O_L,
                          psi_beta_cat_bar_u, psi_beta_cat_s_u, C[1, :P_c]);
      psi_beta_cat_bar = coef.1;
      psi_beta_cat_bar_z = coef.2;
      if (MS) {
        psi_beta_cat_s = coef.3;
        psi_beta_cat_s_z = coef.4;
      }
      for (p in 1:P_c) {
        int C_p = C[1, p];
        if (MS) {
          psi_beta_cat[p, :C_p] = rep_matrix(psi_beta_cat_bar[:C_p, p], S)
                                  + psi_beta_cat_s[p, :C_p];
        } else {
          psi_beta_cat[p, :C_p, 1] = psi_beta_cat_bar[:C_p, p];
        }
      }
      psi_idx += N;
    }

    // occupancy ordinal coefficients
    if (P_ord[1]) {
      int P_o = P_ord[1], N = P_o * (1 + MS);
      tuple(vector[P_o], matrix[P_o, MS * S], matrix[O_max[1], P_o]) coef =
        coef_ord_jacobian(segment(psi_tau, psi_idx, N), psi_beta_O_L,
                          psi_beta_ord_bar_z, psi_beta_ord_s_z,
                          psi_beta_ord_c_u, O[1, :P_o]);
      psi_beta_ord_bar = coef.1;
      if (MS) {
        psi_beta_ord_s = coef.2;
        psi_beta_ord = rep_matrix(psi_beta_ord_bar, S) + psi_beta_ord_s;
      } else {
        psi_beta_ord[:, 1] = psi_beta_ord_bar;
      }
      psi_beta_ord_cs = coef.3;
      psi_idx += N;
    }
    psi_beta_ord_realised = coef_ord_realised(O[1], psi_beta_ord,
                                              psi_beta_ord_cs);

    // detection continuous coefficients (site-level)
    if (P[2]) {
      int N = P[2] * (1 + MS);
      tuple(vector[P[2]], matrix[P[2], MS * S]) coef =
        coef_cont(segment(mu_tau, mu_idx, N), mu_beta_O_L, mu_beta_bar_z,
                  mu_beta_s_z);
      mu_beta_bar = coef.1;
      if (MS) {
        mu_beta_s = coef.2;
        mu_beta = rep_matrix(mu_beta_bar, S) + mu_beta_s;
      } else {
        mu_beta[:, 1] = mu_beta_bar;
      }
      mu_idx += N;
    }

    // detection categorical coefficients (site-level)
    if (P_cat[2]) {
      int P_c = P_cat[2], C_m = C_max[2], N = P_c * (1 + MS);
      tuple(matrix[C_m, P_c], matrix[C_m, P_c], array[P_c] matrix[C_m, MS * S],
            array[P_c] matrix[C_m, MS * S]) coef =
        coef_cat_jacobian(segment(mu_tau, mu_idx, N), mu_beta_O_L,
                          mu_beta_cat_bar_u, mu_beta_cat_s_u, C[2, :P_c]);
      mu_beta_cat_bar = coef.1;
      mu_beta_cat_bar_z = coef.2;
      if (MS) {
        mu_beta_cat_s = coef.3;
        mu_beta_cat_s_z = coef.4;
      }
      for (p in 1:P_c) {
        int C_p = C[2, p];
        if (MS) {
          mu_beta_cat[p, :C_p] = rep_matrix(mu_beta_cat_bar[:C_p, p], S)
                                 + mu_beta_cat_s[p, :C_p];
        } else {
          mu_beta_cat[p, :C_p, 1] = mu_beta_cat_bar[:C_p, p];
        }
      }
      mu_idx += N;
    }

    // detection ordinal coefficients (site-level)
    if (P_ord[2]) {
      int P_o = P_ord[2], N = P_o * (1 + MS);
      tuple(vector[P_o], matrix[P_o, MS * S], matrix[O_max[2], P_o]) coef =
        coef_ord_jacobian(segment(mu_tau, mu_idx, N), mu_beta_O_L,
                          mu_beta_ord_bar_z, mu_beta_ord_s_z, mu_beta_ord_c_u,
                          O[2, :P_o]);
      mu_beta_ord_bar = coef.1;
      if (MS) {
        mu_beta_ord_s = coef.2;
        mu_beta_ord = rep_matrix(mu_beta_ord_bar, S) + mu_beta_ord_s;
      } else {
        mu_beta_ord[:, 1] = mu_beta_ord_bar;
      }
      mu_beta_ord_cs = coef.3;
      mu_idx += N;
    }
    mu_beta_ord_realised = coef_ord_realised(O[2], mu_beta_ord, mu_beta_ord_cs);

    // detection continuous coefficients (site-by-survey level)
    if (P[3]) {
      int N = P[3] * (1 + MS);
      tuple(vector[P[3]], matrix[P[3], MS * S]) coef =
        coef_cont(segment(mu_tau, mu_idx, N), gamma_O_L, gamma_bar_z,
                  gamma_s_z);
      gamma_bar = coef.1;
      if (MS) {
        gamma_s = coef.2;
        gamma = rep_matrix(gamma_bar, S) + gamma_s;
      } else {
        gamma[:, 1] = gamma_bar;
      }
      mu_idx += N;
    }

    // detection categorical coefficients (site-by-survey level)
    if (P_cat[3]) {
      int P_c = P_cat[3], C_m = C_max[3], N = P_c * (1 + MS);
      tuple(matrix[C_m, P_c], matrix[C_m, P_c], array[P_c] matrix[C_m, MS * S],
            array[P_c] matrix[C_m, MS * S]) coef =
        coef_cat_jacobian(segment(mu_tau, mu_idx, N), gamma_O_L,
                          gamma_cat_bar_u, gamma_cat_s_u, C[3, :P_c]);
      gamma_cat_bar = coef.1;
      gamma_cat_bar_z = coef.2;
      if (MS) {
        gamma_cat_s = coef.3;
        gamma_cat_s_z = coef.4;
      }
      for (p in 1:P_c) {
        int C_p = C[3, p];
        if (MS) {
          gamma_cat[p, :C_p] = rep_matrix(gamma_cat_bar[:C_p, p], S)
                               + gamma_cat_s[p, :C_p];
        } else {
          gamma_cat[p, :C_p, 1] = gamma_cat_bar[:C_p, p];
        }
      }
      mu_idx += N;
    }

    // detection ordinal coefficients (site-by-survey level)
    if (P_ord[3]) {
      int P_o = P_ord[3], N = P_o * (1 + MS);
      tuple(vector[P_o], matrix[P_o, MS * S], matrix[O_max[3], P_o]) coef =
        coef_ord_jacobian(segment(mu_tau, mu_idx, N), gamma_O_L,
                          gamma_ord_bar_z, gamma_ord_s_z, gamma_ord_c_u,
                          O[3, :P_o]);
      gamma_ord_bar = coef.1;
      if (MS) {
        gamma_ord_s = coef.2;
        gamma_ord = rep_matrix(gamma_ord_bar, S) + gamma_ord_s;
      } else {
        gamma_ord[:, 1] = gamma_ord_bar;
      }
      gamma_ord_cs = coef.3;
      mu_idx += N;
    }
    gamma_ord_realised = coef_ord_realised(O[3], gamma_ord, gamma_ord_cs);

    // mean and species-specific site effects
    if (random[1]) {
      iota_bar_t[1] = mu_tau[mu_idx];
      if (MS) {
        iota_t = mu_tau[mu_idx + 1] * sqrt(iota_phi[1]);
      }
      mu_idx += 1 + MS;
    }

    // mean and species-specific survey effects
    if (random[2]) {
      kappa_bar_t[1] = mu_tau[mu_idx];
      kappa_bar = kappa_bar_t[1] * kappa_bar_z[1];
      matrix[MS * J, S] kappa_s;
      if (MS) {
        kappa_t = mu_tau[mu_idx + 1] * sqrt(kappa_phi[1]);
        kappa_s = kappa_s_z[1] * diag_post_multiply(kappa_O_L', kappa_t);
      }

      // increment temporal GP
      if (kernel[2]) {
        vector[periodic * 2] kappa_v;
        if (periodic) {
          kappa_v = sqrt(K_phi[1]);
        }
        array[1 + SS[2] * S] matrix[J, J] kappa_L =
          survey_kernels(surveys, kernel[2], period, kappa_v, kappa_ell);
        kappa_bar = kappa_L[1] * kappa_bar;
        if (MS) {
          if (SS[2]) {
            for (s in 1:S) {
              kappa_s[:, s] = kappa_L[s + 1] * kappa_s[:, s];
            }
          } else {
            kappa_s = kappa_L[1] * kappa_s;
          }
        }
      }
      if (MS) {
        kappa = rep_matrix(kappa_bar, S) + kappa_s;
      } else {
        kappa[:, 1] = kappa_bar;
      }

      // increment orthogonal projection
      if (project[2]) {
        kappa_bar2 = kappa_bar;
        kappa_bar = orthogonalise(X3_mean, X3_plus, kappa_bar);
        if (MS) {
          kappa2 = kappa;
          kappa = orthogonalise(X3_mean, X3_plus, kappa);
        }
      }
      mu_idx += 1 + MS;
    }

    // OLRE scales
    if (OLRE) {
      epsilon_bar_t[1] = mu_tau[mu_idx];
      if (MS) {
        epsilon_t = mu_tau[mu_idx + 1] * sqrt(epsilon_phi[1]);
      }
    }
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
  if (kernel[1]) {
    lprior += inv_gamma_lpdf(iota_ell | iota_ell_inv_gamma[1],
                                        iota_ell_inv_gamma[2]);
  }
  if (kernel[2]) {
    lprior += inv_gamma_lpdf(kappa_ell[:, 1] | kappa_ell_inv_gamma[1],
                                               kappa_ell_inv_gamma[2]);
    if (periodic) {
      lprior += inv_gamma_lpdf(kappa_ell[:, 2] |
                                 kappa_ell_periodic_inv_gamma[1],
                                 kappa_ell_periodic_inv_gamma[2])
                + dirichlet_lpdf(K_phi[1] | K_phi_dirichlet);
    }
  }
  if (MS) {
    lprior += lkj_corr_cholesky_lpdf(alpha_O_L | alpha_O_L_LKJ);
    if (P_sum[1]) {
      lprior += lkj_corr_cholesky_lpdf(psi_beta_O_L | psi_beta_O_L_LKJ);
    }
    if (P_sum[2]) {
      lprior += lkj_corr_cholesky_lpdf(mu_beta_O_L | mu_beta_O_L_LKJ);
    }
    if (P_sum[3]) {
      lprior += lkj_corr_cholesky_lpdf(gamma_O_L | gamma_O_L_LKJ);
    }
    if (random[1]) {
      lprior += lkj_corr_cholesky_lpdf(iota_O_L | iota_O_L_LKJ);
    }
    if (random[2]) {
      lprior += lkj_corr_cholesky_lpdf(kappa_O_L | kappa_O_L_LKJ);
    }
    if (OLRE) {
      lprior += lkj_corr_cholesky_lpdf(epsilon_O_L | epsilon_O_L_LKJ);
    }
  }
  if (NB) {
    lprior += inv_gamma_lpdf(phi | phi_inv_gamma[1], phi_inv_gamma[2]);
  }
}

model {
  profile("priors") {
    target += lprior;
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
    if (MS) {
      for (d in 1:2) {
        target += std_normal_lupdf(alpha_s_bar_z[d]);
      }
    }
    if (MR) {
      for (d in 1:2) {
        target += std_normal_lupdf(alpha_r_bar_z[d]);
        if (MS) {
          target += std_normal_lupdf(to_vector(alpha_r_s_z[d]));
        }
      }
    }
    if (P[1]) {
      target += std_normal_lupdf(psi_beta_bar_z);
      if (MS) {
        for (p in 1:P[1]) {
          target += std_normal_lupdf(psi_beta_s_z[p]);
        }
      }
    }
    for (p in 1:P_cat[1]) {
      int C_p = C[1, p];
      target += std_normal_lupdf(psi_beta_cat_bar_z[:C_p, p]);
      if (MS) {
        target += std_normal_lupdf(to_vector(psi_beta_cat_s_z[p, :C_p]));
      }
    }
    if (P_ord[1]) {
      target += std_normal_lupdf(psi_beta_ord_bar_z);
      if (MS) {
        for (p in 1:P_ord[1]) {
          target += std_normal_lupdf(psi_beta_ord_s_z[p]);
        }
      }
    }
    if (P[2]) {
      target += std_normal_lupdf(mu_beta_bar_z);
      if (MS) {
        for (p in 1:P[2]) {
          target += std_normal_lupdf(mu_beta_s_z[p]);
        }
      }
    }
    for (p in 1:P_cat[2]) {
      int C_p = C[2, p];
      target += std_normal_lupdf(mu_beta_cat_bar_z[:C_p, p]);
      if (MS) {
        target += std_normal_lupdf(to_vector(mu_beta_cat_s_z[p, :C_p]));
      }
    }
    if (P_ord[2]) {
      target += std_normal_lupdf(mu_beta_ord_bar_z);
      if (MS) {
        for (p in 1:P_ord[2]) {
          target += std_normal_lupdf(mu_beta_ord_s_z[p]);
        }
      }
    }
    if (P[3]) {
      target += std_normal_lupdf(gamma_bar_z);
      if (MS) {
        for (p in 1:P[3]) {
          target += std_normal_lupdf(gamma_s_z[p]);
        }
      }
    }
    for (p in 1:P_cat[3]) {
      int C_p = C[3, p];
      target += std_normal_lupdf(gamma_cat_bar_z[:C_p, p]);
      if (MS) {
        target += std_normal_lupdf(to_vector(gamma_cat_s_z[p, :C_p]));
      }
    }
    if (P_ord[3]) {
      target += std_normal_lupdf(gamma_ord_bar_z);
      if (MS) {
        for (p in 1:P_ord[3]) {
          target += std_normal_lupdf(gamma_ord_s_z[p]);
        }
      }
    }
    if (random[2]) {
      target += std_normal_lupdf(kappa_bar_z[1]);
      if (MS) {
        target += std_normal_lupdf(to_vector(kappa_s_z[1]));
      }
    }
  }

  // site and OLRE covariances
  matrix[MS * S, S] iota_S_L, epsilon_S_L;
  if (MS) {
    if (random[1]) {
      iota_S_L = diag_pre_multiply(iota_t, iota_O_L);
    }
    if (OLRE) {
      epsilon_S_L = diag_pre_multiply(epsilon_t, epsilon_O_L);
    }
  }

  // parallelise over regions
  if (grainsize && MR) {
    profile("likelihood_partial_regions") {
      target +=
        reduce_sum(partial_regions_occARU_lupmf, regions, grainsize, I_r, XY,
                   kernel[1], X2_aug, X2_plus, y[,,1,], Q, j_idx, J_i, log_Delta,
                   X1, X_cat1, X_ord1_p1, X2, X_cat2, X_ord2_p1, X3, X_cat3,
                   X_ord3_p1, alpha[:, region], psi_beta, psi_beta_cat,
                   psi_beta_ord_realised, mu_beta, mu_beta_cat,
                   mu_beta_ord_realised, gamma, gamma_cat, gamma_ord_realised,
                   iota_bar_z, iota_s_z, iota_bar_t, iota_S_L, iota_ell, kappa,
                   epsilon_bar_z, epsilon_s_z, epsilon_bar_t, epsilon_S_L, phi);
    }
  } else {

    // random site effects with orthogonal projection
    matrix[random[1] * I, S] iota;
    if (random[1]) {
      profile("iota") {
        vector[I] iota_bar = iota_bar_t[1] * iota_bar_z[1];
        matrix[MS * I, S] iota_s;
        if (MS) {
          iota_s = iota_s_z[1] * iota_S_L';
        }
        if (kernel[1]) {
          array[R, 1 + SS[1] * S] matrix[I_max, I_max] iota_L =
            site_kernels(I_r, XY, kernel[1], iota_ell);
            if (MS) {
              if (SS[1]) {
                iota_bar = increment_spatial(I_r, iota_L[:, 1], iota_bar);
                iota_s = increment_spatial(I_r, iota_L, iota_s);
                iota = rep_matrix(iota_bar, S) + iota_s;
              } else {
                iota = increment_spatial(I_r, iota_L,
                                              rep_matrix(iota_bar, S) + iota_s);
              }
            } else {
              iota[:, 1] = increment_spatial(I_r, iota_L[:, 1], iota_bar);
            }
        } else {
          if (MS) {
            iota = rep_matrix(iota_bar, S) + iota_s;
          } else {
            iota[:, 1] = iota_bar;
          }
        }
        if (project[1]) {
          iota = orthogonalise(X2_aug, X2_plus, iota);
        }
      }
    }

    // parallelise over sites
    if (grainsize) {
      profile("likelihood_partial_sites") {
        target +=
          reduce_sum(partial_sites_occARU_lupmf, sites, grainsize, y[,,1,], Q,
                     j_idx, J_i, log_Delta, X1, X_cat1, X_ord1_p1, X2, X_cat2,
                     X_ord2_p1, X3, X_cat3, X_ord3_p1, alpha[:, region],
                     psi_beta, psi_beta_cat, psi_beta_ord_realised, mu_beta,
                     mu_beta_cat, mu_beta_ord_realised, gamma, gamma_cat,
                     gamma_ord_realised, iota, iota_bar_z, iota_s_z, kappa,
                     epsilon_bar_z, epsilon_s_z, epsilon_bar_t, epsilon_S_L,
                     phi);
      }
    } else {

      // no parallelisation
      if (random[1]) {
        target += std_normal_lupdf(iota_bar_z[1]);
        if (MS) {
          target += std_normal_lupdf(to_vector(iota_s_z[1]));
        }
      }

      // Poisson OLREs
      array[OLRE * I] matrix[J, S] epsilon;
      if (OLRE) {
        profile("OLRE") {
          vector[J_sum] epsilon_bar = epsilon_bar_t[1] * epsilon_bar_z[1];
          target += std_normal_lupdf(epsilon_bar_z[1]);
          matrix[J_sum, S] epsilon_mat;
          if (MS) {
            epsilon_mat = rep_matrix(epsilon_bar, S)
                          + epsilon_s_z[1] * epsilon_S_L';
            target += std_normal_lupdf(to_vector(epsilon_s_z[1]));
          } else {
            epsilon_mat[:, 1] = epsilon_bar;
          }
          epsilon = fill_epsilon(j_idx, J_i, epsilon_mat);
        }
      }
      profile("likelihood") {
        target +=
          occARU_lupmf(y[,,1,] | Q, j_idx, J_i, log_Delta, X1, X_cat1, X_ord1_p1,
                       X2, X_cat2, X_ord2_p1, X3, X_cat3, X_ord3_p1,
                       alpha[:, region], psi_beta, psi_beta_cat,
                       psi_beta_ord_realised, mu_beta, mu_beta_cat,
                       mu_beta_ord_realised, gamma, gamma_cat,
                       gamma_ord_realised, iota, kappa, epsilon, phi);
      }
    }
  }
}

generated quantities {
  // correlations
  corr_matrix[MS * 2] alpha_O;
  corr_matrix[MS && P_sum[1] ? S : 0] psi_beta_O;
  corr_matrix[MS && P_sum[2] ? S : 0] mu_beta_O;
  corr_matrix[MS && P_sum[3] ? S : 0] gamma_O;
  corr_matrix[MS && random[1] ? S : 0] iota_O;
  corr_matrix[MS && random[2] ? S : 0] kappa_O;
  corr_matrix[MS && OLRE ? S : 0] epsilon_O;
  if (MS) {
    alpha_O = multiply_lower_tri_self_transpose(alpha_O_L);
    if (P_sum[1]) {
      psi_beta_O = multiply_lower_tri_self_transpose(psi_beta_O_L);
    }
    if (P_sum[2]) {
      mu_beta_O = multiply_lower_tri_self_transpose(mu_beta_O_L);
    }
    if (P_sum[3]) {
      gamma_O = multiply_lower_tri_self_transpose(gamma_O_L);
    }
    if (random[1]) {
      iota_O = multiply_lower_tri_self_transpose(iota_O_L);
    }
    if (random[2]) {
      kappa_O = multiply_lower_tri_self_transpose(kappa_O_L);
    }
    if (OLRE) {
      epsilon_O = multiply_lower_tri_self_transpose(epsilon_O_L);
    }
  }

  // unconditional coefficients
  vector[project[1] * P[2]] mu_beta_bar2;
  matrix[project[1] * MS * P[2], S] mu_beta2;
  vector[project[1] * P_ord[2]] mu_beta_ord_bar2;
  matrix[project[1] * MS * P_ord[2], S] mu_beta_ord2;
  vector[project[2] * P[3]] gamma_bar2;
  matrix[project[2] * MS * P[3], S] gamma2;
  vector[project[2] * P_ord[3]] gamma_ord_bar2;
  matrix[project[2] * MS * P_ord[3], S] gamma_ord2;
  if (project[2]) {
    if (P[3]) {
      gamma_bar2 = uncondition(X3_plus[:P[3]], gamma_bar, kappa_bar2);
      if (MS) {
        gamma2 = uncondition(X3_plus[:P[3]], gamma, kappa2);
      }
    }
    if (P_ord[3]) {
      int idx = P_aug[3] - P_ord[3] + 1;
      gamma_ord_bar2 = uncondition(X3_plus[idx:], gamma_ord_bar, kappa_bar2);
      if (MS) {
        gamma_ord2 = uncondition(X3_plus[idx:], gamma_ord, kappa2);
      }
    }
  }

  // log likelihood, latent occupancy, and posterior predictions
  vector[random[1] * I] iota_bar;
  vector[project[1] * I] iota_bar2;
  matrix[random[1] * I, S] iota;
  matrix[project[1] * I, S] iota2;
  matrix[I, S] logit_psi;
  matrix[S, I] log_lik;
  matrix[MC * S, I] log_lik2;
  array[latent * I, S] int z;
  array[PPC_y * I, J, S] int yrep;
  array[PPC_Q * I, S] int Qrep;

  profile("generated_quantities") {

    // site effects and unconditional coefficients
    matrix[MS * S, S] iota_S_L;
    array[(kernel[1] > 0) * R, 1 + SS[1] * S] matrix[I_max, I_max] iota_L;
    matrix[MS * random[1] * I, S] iota_s;
    if (random[1]) {
      iota_bar = iota_bar_t[1] * iota_bar_z[1];
      if (MS) {
        iota_S_L = diag_pre_multiply(iota_t, iota_O_L);
        iota_s = iota_s_z[1] * iota_S_L';
      }
      if (kernel[1]) {
        iota_L = site_kernels(I_r, XY, kernel[1], iota_ell);
        iota_bar = increment_spatial(I_r, iota_L[:, 1], iota_bar);
        if (MS) {
          iota_s = increment_spatial(I_r, iota_L, iota_s);
        }
      }
      if (MS) {
        iota = rep_matrix(iota_bar, S) + iota_s;
      } else {
        iota[:, 1] = iota_bar;
      }
      if (project[1]) {
        iota_bar2 = iota_bar;
        iota_bar = orthogonalise(X2_aug, X2_plus, iota_bar);
        if (MS) {
          iota2 = iota;
          iota = orthogonalise(X2_aug, X2_plus, iota);
        }
        if (P[2]) {
          mu_beta_bar2 = uncondition(X2_plus[:P[2]], mu_beta_bar, iota_bar2);
          if (MS) {
            mu_beta2 = uncondition(X2_plus[:P[2]], mu_beta, iota2);
          }
        }
        if (P_ord[2]) {
          int idx = P_aug[2] - P_ord[2] + 1;
          mu_beta_ord_bar2 = uncondition(X2_plus[idx:], mu_beta_ord_bar,
                                         iota_bar2);
          if (MS) {
            mu_beta_ord2 = uncondition(X2_plus[idx:], mu_beta_ord, iota2);
          }
        }
      }
    }

    // reconstruct log likelihood and recover occupancy
    vector[OLRE * J_sum] epsilon_bar;
    matrix[MS * OLRE * S, S] epsilon_S_L;
    matrix[OLRE * J_sum, S] epsilon_mat;
    array[OLRE * I] matrix[J, S] epsilon;
    if (OLRE) {
      epsilon_bar = epsilon_bar_t[1] * epsilon_bar_z[1];
      if (MS) {
        epsilon_S_L = diag_pre_multiply(epsilon_t, epsilon_O_L);
        epsilon_mat = rep_matrix(epsilon_bar, S)
                      + epsilon_s_z[1] * epsilon_S_L';
      } else {
        epsilon_mat[:, 1] = epsilon_bar;
      }
      epsilon = fill_epsilon(j_idx, J_i, epsilon_mat);
    }
    tuple(matrix[S, I], array[I] matrix[2, S], matrix[I, S],
          array[I] matrix[J, S]) lp =
      occARU(y[,,1,], Q, j_idx, J_i, log_Delta, X1, X_cat1, X_ord1_p1, X2, X_cat2,
             X_ord2_p1, X3, X_cat3, X_ord3_p1, alpha[:, region], psi_beta,
             psi_beta_cat, psi_beta_ord_realised, mu_beta, mu_beta_cat,
             mu_beta_ord_realised, gamma, gamma_cat, gamma_ord_realised, iota,
             kappa, epsilon, phi);
    log_lik = lp.1;
    if (latent) {
      z = occARU_rng(Q, lp);
    }
    logit_psi = lp.3;

    // posterior predictions
    if (PPC_y || PPC_Q) {
      tuple(array[I, J, S] int, array[I, S] int) pp =
        pp_occARU_rng(j_idx, J_i, X3, logit_psi, lp.4, gamma, phi);
      if (PPC_y) {
        yrep = pp.1;
      }
      if (PPC_Q) {
        Qrep = pp.2;
      }
    }

    // Monte Carlo integration of observation-level parameters for loo
    if (MC) {
      array[D] matrix[S, I] log_lik_d;

      // posterior predict site effects and OLREs and compute log likelihood
      matrix[Im1, MS ? Sm1 : 1] u;
      vector[!MS * random[1] * Im1] iota_bar_rep;
      matrix[random[1] * I, S] iota_rep;
      for (d in 1:D) {
        if (random[1]) {
          if (MS) {
            u = to_matrix(normal_rng(zeros_vector(Im1 * Sm1), 1), Im1, Sm1);
            iota_s = sum_to_zero_constrain(u) * iota_S_L';
            if (kernel[1]) {
              iota_s = increment_spatial(I_r, iota_L, iota_s);
            }
            if (project[1]) {
              iota_s = orthogonalise(X2_aug, X2_plus, iota_s);
            }
            iota_rep = rep_matrix(iota_bar, S) + iota_s;
          } else {
            u[:, 1] = to_vector(normal_rng(zeros_vector(Im1), 1));
            iota_rep[:, 1] = iota_bar_t[1] * sum_to_zero_constrain(u[:, 1]);
            if (kernel[1]) {
              iota_rep[:, 1] = increment_spatial(I_r, iota_L[:, 1],
                                                 iota_rep[:, 1]);
            }
            if (project[1]) {
              iota_rep[:, 1] = orthogonalise(X2_aug, X2_plus, iota_rep[:, 1]);
            }
          }
        }
        if (OLRE) {
          epsilon = pp_epsilon_rng(j_idx, J_i, epsilon_bar, epsilon_bar_t,
                                   epsilon_S_L);
        }
        log_lik_d[d] =
          occARU(y[,,1,], Q, j_idx, J_i, log_Delta, X1, X_cat1, X_ord1_p1, X2, X_cat2,
                 X_ord2_p1, X3, X_cat3, X_ord3_p1, alpha[:, region], psi_beta,
                 psi_beta_cat, psi_beta_ord_realised, mu_beta, mu_beta_cat,
                 mu_beta_ord_realised, gamma, gamma_cat, gamma_ord_realised,
                 iota_rep, kappa, epsilon, phi).1;
      }

      // average log likelihood over draws
      log_lik2 = rep_matrix(-log_D, S, I);
      for (i in 1:I) {
        for (s in 1:S) {
          log_lik2[s, i] += log_sum_exp(log_lik_d[:, s, i]);
        }
      }
    }
  }
}
