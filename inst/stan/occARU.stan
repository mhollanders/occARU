functions {
  #include occARU.stanfunctions
  #include util.stanfunctions
  #include predictors.stanfunctions
  #include coefficients.stanfunctions
  #include kernels.stanfunctions
  #include predict.stanfunctions
}

data {
  // dimensions, survey effort, and detection history
  int<lower=1> I, J, K, S, R;
  array[K] matrix<lower=0, upper=1>[J, I] Delta;
  array[K, I, J, S] int<lower=0> y;

  // site UTMs (km), regions (ordered), season intervals and averages
  array[I] vector[2] XY;
  array[I] int<lower=1, upper=R> region;
  matrix<lower=0>[K - 1, I] tau;
  array[K] real<lower=0> seasons;

  // continuous preds for site occ. and site and survey det.
  array[3] int<lower=0> P;
  array[K] matrix[I, P[1]] X1;
  array[K] matrix[I, P[2]] X2;
  array[K, I] matrix[J, P[3]] X3;

  // categorical preds for site occ. and site and survey det.
  array[3] int<lower=0> P_cat;
  array[K, I, P_cat[1]] int<lower=1> X_cat1;
  array[K, I, P_cat[2]] int<lower=1> X_cat2;
  array[K, I, J, P_cat[3]] int<lower=1> X_cat3;

  // ordinal preds for site occ. and site and survey det.
  array[3] int<lower=0> P_ord;
  array[K, I, P_ord[1]] int<lower=0> X_ord1;
  array[K, I, P_ord[2]] int<lower=0> X_ord2;
  array[K, I, J, P_ord[3]] int<lower=0> X_ord3;

  // random effects
  array[5] int<lower=0, upper=1> random,  // indicators
                                 SS;  // species length scales
  array[5] int<lower=0, upper=3> kernel;  // GP kernel type
  real<lower=0> period;  // period for periodic kernel

  // indicators
  int<lower=0, upper=2> OD;  // overdispersion (Poisson, OLRE, negbin)
  int<lower=0, upper=1> dirichlet,  // var. decomp. (logi-norm or Dirichlet)
                        latent,  // recover latent occ. states
                        PPC_y,  // posterior predictions of y
                        PPC_Q,  // posterior predictions of Q
                        dyn;  // dynamic occupancy
  int<lower=0, upper=I> grainsize;  // grainsize for threading
  int<lower=0> D;  // Monte Carlo draws for LOO

  // prior hyperparameters
  vector<lower=0>[2] psi_bar_beta,  // mean occupancy (probability scale)
                     q_bar_gamma,  // mean transitions (rate scale)
                     mu_bar_gamma,  // mean detetection (rate scale)
                     psi_theta_gamma,  // occupancy variance partition scale
                     mu_theta_gamma,  // detection variance partition scale
                     iota_ell_inv_gamma,  // site GP length scale
                     kappa_ell_inv_gamma,  // survey GP length scale
                     kappa_ell_periodic_inv_gamma,  // periodic GP length scale
                     K_phi_dirichlet,  // survey GP variance partitions
                     nu_ell_inv_gamma,  // season GP length scale
                     psi_iota_ell_inv_gamma,  // occupancy site GP length scale
                     psi_nu_ell_inv_gamma,  // occupancy season GP length scale
                     phi_inv_gamma;  // negative binomial overdispersion
  vector<lower=0>[3] psi_W_t,  // occupancy log odds variance
                     mu_W_t;  // log detection rate variance
  real<lower=0> alpha_O_L_LKJ, psi_beta_O_L_LKJ, mu_beta_O_L_LKJ, gamma_O_L_LKJ,
                iota_O_L_LKJ, kappa_O_L_LKJ, nu_O_L_LKJ, psi_iota_O_L_LKJ,
                psi_nu_O_L_LKJ, epsilon_O_L_LKJ;
}

transformed data {
  // transformed indicators
  int MS = S > 1, MR = R > 1, periodic = kernel[2] > 0 && period > 0,
      OLRE = OD == 1, NB = OD == 2, MC = (D > 0) && (random[1] || OLRE),
      Km1 = K - 1, Im1 = I - 1, Sm1 = S - 1;
  real log_D = log(D);
  array[3] int P_sum;
  for (d in 1:3) {
    P_sum[d] = P[d] + P_cat[d] + P_ord[d];
  }

  // survey indices, total surveys, and offsets
  tuple(array[K, I, J] int, array[K, I] int) survey_idx = survey_indices(Delta);
  array[K, I, J] int j_idx = survey_idx.1;
  array[K, I] int J_i = survey_idx.2;
  int J_sum = 0;
  for (k in 1:K) {
    J_sum += sum(J_i[k]);
  }
  array[K] matrix[J, I] log_Delta = log(Delta);

  // first and last season per site
  array[I, 2] int f_l = first_last(J_i);

  // season indices
  tuple(array[K, I] int, array[K] int) site_idx = site_indices(J_i);
  array[K, I] int i_idx = site_idx.1;
  array[K] int I_k = site_idx.2;

  // number of surveys per season
  array[K] int J_k = zeros_int_array(K);
  for (k in 1:K) {
    for (i in 1:I) {
      if (J_i[k, i]) {
        J_k[k] = max(J_k[k], max(j_idx[k, i, :J_i[k, i]]));
      }
    }
  }
  int J_ksum = sum(J_k);

  // number of sites per region
  array[R] int I_r = zeros_int_array(R);
  for (i in 1:I) {
    I_r[region[i]] += 1;
  }
  int I_max = max(I_r);

  // aggregated counts
  array[K, I, S] int Q;
  for (k in 1:K) {
    for (i in i_idx[k, :I_k[k]]) {
      int J_ik = J_i[k, i];
      array[J_ik] int idx = j_idx[k, i, :J_ik];
      for (s in 1:S) {
        Q[k, i, s] = sum(y[k, i, idx, s]);
      }
    }
  }

  // site, region, and survey sequences
  array[I] int sites = linspaced_int_array(I, 1, I);
  array[R] int regions = linspaced_int_array(R, 1, R);
  array[K, J] real surveys;
  for (k in 1:K) {
    int J_kk = J_k[k];
    surveys[k, :J_kk] = linspaced_array(J_kk, 1, J_kk);
  }

  // variance partitions
  int psi_V = MS + (1 + MS) * (MR + P_sum[1] + sum(random[4:5])),
      mu_V = MS + (1 + MS) * (MR + sum(P_sum[2:3]) + sum(random[1:3]) + OLRE);

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
  array[K, I, P_ord[1]] int X_ord1_p1 = X_ord1;
  array[K, I, P_ord[2]] int X_ord2_p1 = X_ord2;
  array[K, I, J, P_ord[3]] int X_ord3_p1 = X_ord3;
  for (k in 1:K) {
    for (i in i_idx[k, :I_k[k]]) {
      for (p in 1:P_ord[1]) {
        X_ord1_p1[k, i, p] += 1;
      }
      for (p in 1:P_ord[2]) {
        X_ord2_p1[k, i, p] += 1;
      }
      for (j in j_idx[k, i, :J_i[k, i]]) {
        for (p in 1:P_ord[3]) {
          X_ord3_p1[k, i, j, p] += 1;
        }
      }
    }
  }

  // augmented site-level design matrices and pseudo-inverses
  array[3] int P_aug;
  for (d in 1:3) {
    P_aug[d] = P[d] + C_sum[d] - P_cat[d] + P_ord[d];
  }

  // occupancy site-level predictors
  array[K] matrix[I, P_aug[1]] X1_aug;
  matrix[I, P_aug[1]] X1_mean;
  matrix[P_aug[1], I] X1_plus;
  if (P_sum[1]) {
    X1_aug = augment_design_matrix(X1, X_cat1, X_ord1, C[1, :P_cat[1]],
                                   O[1, :P_ord[1]]);
    X1_mean = average_design_matrix(X1_aug, i_idx, I_k);
    X1_plus = pseudo_inverse(X1_mean);
  }

  // detection site-level predictors
  array[K] matrix[I, P_aug[2]] X2_aug;
  matrix[I, P_aug[2]] X2_mean;
  matrix[P_aug[2], I] X2_plus;
  if (P_sum[2]) {
    X2_aug = augment_design_matrix(X2, X_cat2, X_ord2, C[2, :P_cat[2]],
                                   O[2, :P_ord[2]]);
    X2_mean = average_design_matrix(X2_aug, i_idx, I_k);
    X2_plus = pseudo_inverse(X2_mean);
  }

  // detection survey-level predictors
  array[K, I] matrix[J, P_aug[3]] X3_aug;
  array[K] matrix[J, P_aug[3]] X3_mean;
  matrix[P_aug[3], J_ksum] X3_plus;
  if (P_sum[3]) {
    X3_aug = augment_design_matrix(X3, X_cat3, X_ord3, C[3, :P_cat[3]],
                                   O[3, :P_ord[3]]);
    X3_mean = average_design_matrix(X3_aug, i_idx, I_k, j_idx, J_i);
    matrix[J_ksum, P_aug[3]] X3_long;
    int idx = 0;
    for (k in 1:K) {
      X3_long[sum(J_k[:k - 1]) + 1:sum(J_k[:k])] = X3_mean[k, :J_k[k]];
    }
    X3_plus = pseudo_inverse(X3_long);
  }
}

parameters {
  // intercepts by species and region and factor loadings
  real<lower=0, upper=1> psi_bar;
  vector<lower=0>[dyn * 2] q_bar;
  real<lower=0> mu_bar;
  cholesky_factor_corr[MS * 2] alpha_O_L;
  array[MS * 2] sum_to_zero_vector[S] alpha_s_bar_z;
  array[MR * 2] sum_to_zero_vector[R] alpha_r_bar_z;
  array[MS * MR * 2] sum_to_zero_matrix[R, S] alpha_r_s_z;
  vector<lower=0>[dyn * 3] lambda_bar;
  array[MS * dyn * 3] sum_to_zero_vector[S] lambda_s;

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
  cholesky_factor_corr[MS && random[3] ? S : 0] nu_O_L;
  cholesky_factor_corr[MS && random[4] ? S : 0] psi_iota_O_L;
  cholesky_factor_corr[MS && random[5] ? S : 0] psi_nu_O_L;
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
  vector[random[2] * (J_ksum - K)] kappa_bar_z_u;
  matrix[random[2] * (J_ksum - K), Sm1] kappa_s_z_u;
  array[MS * random[2]] simplex[S] kappa_phi;
  matrix<lower=0>[1 + (SS[2] * S), (kernel[2] > 0) + periodic] kappa_ell;
  array[periodic] simplex[2] K_phi;

  // season effects
  array[random[3]] sum_to_zero_vector[K] nu_bar_z;
  array[MS * random[3]] sum_to_zero_matrix[K, S] nu_s_z;
  array[MS * random[3]] simplex[S] nu_phi;
  vector<lower=0>[(kernel[3] > 0) * (1 + SS[3] * S)] nu_ell;

  // site and season effects on occupancy
  array[random[4]] sum_to_zero_vector[I] psi_iota_bar_z;
  array[MS * random[4]] sum_to_zero_matrix[I, S] psi_iota_s_z;
  array[MS * random[4]] simplex[S] psi_iota_phi;
  vector<lower=0>[(kernel[4] > 0) * (1 + SS[4] * S)] psi_iota_ell;
  array[random[5]] sum_to_zero_vector[K] psi_nu_bar_z;
  array[MS * random[5]] sum_to_zero_matrix[K, S] psi_nu_s_z;
  array[MS * random[5]] simplex[S] psi_nu_phi;
  vector<lower=0>[(kernel[5] > 0) * (1 + SS[5] * S)] psi_nu_ell;

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
  } else if (mu_V == 1) {
    mu_phi[1] = 1;
    mu_tau[1] = sqrt(mu_W[1]);
  }

  // intercepts
  vector[dyn ? 4 : 2] alpha_bar;
  alpha_bar[1] = logit(psi_bar);
  if (dyn) {
    alpha_bar[2:3] = log(q_bar);
  }
  alpha_bar[dyn ? 4 : 2] = log(mu_bar);
  matrix[2, MS * S] alpha_s_bar;
  matrix[R, MR * 2] alpha_r_bar;
  array[MS * MR * 2] matrix[R, S] alpha_r_s;
  array[2] matrix[R, S] alpha = rep_array(rep_matrix(0, R, S), 2);
  if (dyn) {
    alpha[2] += alpha_bar[4];
  } else {
    for (d in 1:2) {
      alpha[d] += alpha_bar[d];
    }
  }

  // factor loadings
  matrix[dyn * 3, S] lambda;
  if (dyn) {
    if (MS) {
      lambda = rep_matrix(lambda_bar, S);
      for (d in 1:3) {
        lambda[d] += lambda_s[d]';
      }
    } else {
      lambda[:, 1] = lambda_bar;
    }

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
  matrix[random[2] * J, K] kappa_bar_z, kappa_bar;
  array[random[2] * K] matrix[J, S] kappa_s_z, kappa;

  // season effects
  vector[random[3] * K] nu_bar;
  matrix[random[3] * K, S] nu;
  vector[random[5] * K] psi_nu_bar;
  matrix[random[5] * K, S] psi_nu;

  // random scales
  vector[random[1]] iota_bar_t;
  vector[MS * random[1] * S] iota_t;
  vector[random[2]] kappa_bar_t;
  vector[MS * random[2] * S] kappa_t;
  vector[random[3]] nu_bar_t;
  vector[MS * random[3] * S] nu_t;
  vector[random[4]] psi_iota_bar_t;
  vector[MS * random[4] * S] psi_iota_t;
  vector[random[5]] psi_nu_bar_t;
  vector[MS * random[5] * S] psi_nu_t;
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
      matrix[MS * S, S] kappa_S_U;
      array[MS * K] matrix[J, S] kappa_s;
      if (MS) {
        kappa_t = mu_tau[mu_idx + 1] * sqrt(kappa_phi[1]);
        kappa_S_U = diag_post_multiply(kappa_O_L', kappa_t);
      }
      vector[periodic * 2] kappa_v;
      if (periodic) {
        kappa_v = sqrt(K_phi[1]);
      }
      int start = 0;
      for (k in 1:K) {
        int J_kk = J_k[k], J_km1 = J_kk - 1;
        array[J_km1] int idx_j =
          linspaced_int_array(J_km1, start + 1, start + J_km1);
        kappa_bar_z[:J_kk, k] = sum_to_zero_jacobian(kappa_bar_z_u[idx_j]);
        kappa_bar[:J_kk, k] = kappa_bar_t[1] * kappa_bar_z[:J_kk, k];
        if (MS) {
          kappa_s_z[k, :J_kk] = sum_to_zero_jacobian(kappa_s_z_u[idx_j]);
          kappa_s[k, :J_kk] = kappa_s_z[k, :J_kk] * kappa_S_U;
        }
        start += J_km1;

        // increment temporal GP
        if (kernel[2]) {
          array[1 + SS[2] * S] matrix[J_kk, J_kk] kappa_L =
            survey_kernels(surveys[k, :J_kk], kernel[2], period, kappa_v,
                           kappa_ell);
          kappa_bar[:J_kk, k] = kappa_L[1] * kappa_bar[:J_kk, k];
          if (MS) {
            if (SS[2]) {
              for (s in 1:S) {
                kappa_s[k, :J_kk, s] = kappa_L[1 + s] * kappa_s[k, :J_kk, s];
              }
            } else {
              kappa_s[k, :J_kk] = kappa_L[1] * kappa_s[k, :J_kk];
            }
          }
        }
        if (MS) {
          kappa[k, :J_kk] = rep_matrix(kappa_bar[:J_kk, k], S)
                            + kappa_s[k, :J_kk];
        } else {
          kappa[k, :J_kk, 1] = kappa_bar[:J_kk, k];
        }
      }
      mu_idx += 1 + MS;
    }

    // random season effects
    if (random[3]) {
      nu_bar = mu_tau[mu_idx] * nu_bar_z[1];
      matrix[MS * K, S] nu_s;
      if (MS) {
        nu_t = mu_tau[mu_idx + 1] * sqrt(nu_phi[1]);
        nu_s = nu_s_z[1] * diag_post_multiply(nu_O_L', nu_t);
      }
      if (kernel[3]) {
        array[1 + SS[3] * S] matrix[K, K] nu_L =
          season_kernels(seasons, kernel[3], nu_ell);
        nu_bar = nu_L[1] * nu_bar;
        if (MS) {
          if (SS[3]) {
            for (s in 1:S) {
              nu_s[:, s] = nu_L[1 + s] * nu_s[:, s];
            }
          } else {
            nu_s = nu_L[1] * nu_s;
          }
        }
      }
      if (MS) {
        nu = rep_matrix(nu_bar, S) + nu_s;
      } else {
        nu[:, 1] = nu_bar;
      }
      mu_idx += 1 + MS;
    }

    // mean and species-specific occupancy site and season effects
    if (random[4]) {
      psi_iota_bar_t[1] = psi_tau[psi_idx];
      if (MS) {
        psi_iota_t = psi_tau[psi_idx + 1] * sqrt(psi_iota_phi[1]);
      }
      psi_idx += 1 + MS;
    }
    if (random[5]) {
      psi_nu_bar = psi_tau[psi_idx] * psi_nu_bar_z[1];
      matrix[MS * K, S] psi_nu_s;
      if (MS) {
        psi_nu_t = psi_tau[psi_idx + 1] * sqrt(psi_nu_phi[1]);
        psi_nu_s = psi_nu_s_z[1] * diag_post_multiply(psi_nu_O_L', psi_nu_t);
      }
      if (kernel[5]) {
        array[1 + SS[5] * S] matrix[K, K] psi_nu_L =
          season_kernels(seasons, kernel[5], psi_nu_ell);
        psi_nu_bar = psi_nu_L[1] * psi_nu_bar;
        if (MS) {
          if (SS[5]) {
            for (s in 1:S) {
              psi_nu_s[:, s] = psi_nu_L[1 + s] * psi_nu_s[:, s];
            }
          } else {
            psi_nu_s = psi_nu_L[1] * psi_nu_s;
          }
        }
      }
      if (MS) {
        psi_nu = rep_matrix(psi_nu_bar, S) + psi_nu_s;
      } else {
        psi_nu[:, 1] = psi_nu_bar;
      }
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
  if (dyn) {
    lprior += gamma_lpdf(q_bar | q_bar_gamma[1], q_bar_gamma[2])
              + std_normal_lpdf(lambda_bar);
    if (MS) {
      for (d in 1:3) {
        lprior += std_normal_lpdf(lambda_s[d]);
      }
    }
  }
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
  if (kernel[3]) {
    lprior += inv_gamma_lpdf(nu_ell | nu_ell_inv_gamma[1], nu_ell_inv_gamma[2]);
  }
  if (kernel[4]) {
    lprior += inv_gamma_lpdf(psi_iota_ell | psi_iota_ell_inv_gamma[1],
                                            psi_iota_ell_inv_gamma[2]);
  }
  if (kernel[5]) {
    lprior += inv_gamma_lpdf(psi_nu_ell | psi_nu_ell_inv_gamma[1],
                                          psi_nu_ell_inv_gamma[2]);
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
    if (random[3]) {
      lprior += lkj_corr_cholesky_lpdf(nu_O_L | nu_O_L_LKJ);
    }
    if (random[4]) {
      lprior += lkj_corr_cholesky_lpdf(psi_iota_O_L | psi_iota_O_L_LKJ);
    }
    if (random[5]) {
      lprior += lkj_corr_cholesky_lpdf(psi_nu_O_L | psi_nu_O_L_LKJ);
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
      for (k in 1:K) {
        int J_kk = J_k[k];
        target += std_normal_lupdf(kappa_bar_z[:J_kk, k]);
        if (MS) {
          target += std_normal_lupdf(to_vector(kappa_s_z[k, :J_kk]));
        }
      }
    }
    if (random[3]) {
      target += std_normal_lupdf(nu_bar_z[1]);
      if (MS) {
        target += std_normal_lupdf(to_vector(nu_s_z[1]));
      }
    }
    if (random[5]) {
      target += std_normal_lupdf(psi_nu_bar_z[1]);
      if (MS) {
        target += std_normal_lupdf(to_vector(psi_nu_s_z[1]));
      }
    }
  }

  // site and OLRE covariances
  matrix[MS * S, S] iota_S_U, psi_iota_S_U, epsilon_S_U;
  if (MS) {
    if (random[1]) {
      iota_S_U = diag_post_multiply(iota_O_L', iota_t);
    }
    if (random[4]) {
      psi_iota_S_U = diag_post_multiply(psi_iota_O_L', psi_iota_t);
    }
    if (OLRE) {
      epsilon_S_U = diag_post_multiply(epsilon_O_L', epsilon_t);
    }
  }

  // Poisson OLREs
  array[OLRE * K, I] matrix[J, S] epsilon;
  if (OLRE) {
    profile("OLRE") {
      vector[J_sum] epsilon_bar = epsilon_bar_t[1] * epsilon_bar_z[1];
      target += std_normal_lupdf(epsilon_bar_z[1]);
      matrix[J_sum, S] epsilon_mat;
      if (MS) {
        epsilon_mat = rep_matrix(epsilon_bar, S)
                      + epsilon_s_z[1] * epsilon_S_U;
        target += std_normal_lupdf(to_vector(epsilon_s_z[1]));
      } else {
        epsilon_mat[:, 1] = epsilon_bar;
      }
      epsilon = fill_epsilon(j_idx, J_i, epsilon_mat);
    }
  }

  // parallelise over regions
  if (grainsize && MR) {
    profile("likelihood_partial_regions") {
      target +=
        reduce_sum(partial_regions_occARU_lupmf, regions, grainsize, I_r, y, Q,
                   f_l, i_idx, I_k, j_idx, J_i, log_Delta, tau, XY, region, P,
                   P_cat, P_ord, random, kernel, SS, X1, X_cat1, X_ord1_p1, X2,
                   X_cat2, X_ord2_p1, X3, X_cat3, X_ord3_p1, alpha_bar, alpha,
                   lambda, psi_beta, psi_beta_cat, psi_beta_ord_realised,
                   mu_beta, mu_beta_cat, mu_beta_ord_realised, gamma, gamma_cat,
                   gamma_ord_realised, iota_bar_t, iota_bar_z, iota_S_U,
                   iota_s_z, iota_ell, kappa, nu, psi_iota_bar_t,
                   psi_iota_bar_z, psi_iota_S_U, psi_iota_s_z, psi_iota_ell,
                   psi_nu, epsilon, phi);
    }
  } else {

    // random site effects
    matrix[random[1] * I, S] iota;
    if (random[1]) {
      profile("iota") {
        vector[I] iota_bar = iota_bar_t[1] * iota_bar_z[1];
        matrix[MS * I, S] iota_s;
        if (MS) {
          iota_s = iota_s_z[1] * iota_S_U;
        }
        if (kernel[1]) {
          array[R, 1 + SS[1] * S] matrix[I_max, I_max] L =
            site_kernels(I_r, XY, kernel[1], iota_ell);
          if (MS) {
            if (SS[1]) {
              iota_bar = increment_spatial(I_r, L[:, 1], iota_bar);
              iota_s = increment_spatial(I_r, L, iota_s);
              iota = rep_matrix(iota_bar, S) + iota_s;
            } else {
              iota = increment_spatial(I_r, L,
                                       rep_matrix(iota_bar, S) + iota_s);
            }
          } else {
            iota[:, 1] = increment_spatial(I_r, L[:, 1], iota_bar);
          }
        } else {
          if (MS) {
            iota = rep_matrix(iota_bar, S) + iota_s;
          } else {
            iota[:, 1] = iota_bar;
          }
        }
      }
    }

    // occupancy random site effects
    matrix[random[4] * I, S] psi_iota;
    if (random[4]) {
      profile("psi_iota") {
        vector[I] psi_iota_bar = psi_iota_bar_t[1] * psi_iota_bar_z[1];
        matrix[MS * I, S] psi_iota_s;
        if (MS) {
          psi_iota_s = psi_iota_s_z[1] * psi_iota_S_U;
        }
        if (kernel[4]) {
          array[R, 1 + SS[4] * S] matrix[I_max, I_max] L =
            site_kernels(I_r, XY, kernel[4], psi_iota_ell);
          if (MS) {
            if (SS[4]) {
              psi_iota_bar = increment_spatial(I_r, L[:, 1], psi_iota_bar);
              psi_iota_s = increment_spatial(I_r, L, psi_iota_s);
              psi_iota = rep_matrix(psi_iota_bar, S) + psi_iota_s;
            } else {
              psi_iota = increment_spatial(I_r, L,
                                           rep_matrix(psi_iota_bar, S)
                                           + psi_iota_s);
            }
          } else {
            psi_iota[:, 1] = increment_spatial(I_r, L[:, 1], psi_iota_bar);
          }
        } else {
          if (MS) {
            psi_iota = rep_matrix(psi_iota_bar, S) + psi_iota_s;
          } else {
            psi_iota[:, 1] = psi_iota_bar;
          }
        }
      }
    }

    // parallelise over sites
    if (grainsize) {
      profile("likelihood_partial_sites") {
        target +=
          reduce_sum(partial_sites_occARU_lupmf, sites, grainsize, y, Q, f_l,
                     i_idx, I_k, j_idx, J_i, log_Delta, tau, region, P, P_cat,
                     P_ord, random, X1, X_cat1, X_ord1_p1, X2, X_cat2,
                     X_ord2_p1, X3, X_cat3, X_ord3_p1, alpha_bar, alpha, lambda,
                     psi_beta, psi_beta_cat, psi_beta_ord_realised, mu_beta,
                     mu_beta_cat, mu_beta_ord_realised, gamma, gamma_cat,
                     gamma_ord_realised, iota, iota_bar_z, iota_s_z, kappa, nu,
                     psi_iota, psi_iota_bar_z, psi_iota_s_z, psi_nu, epsilon,
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
      if (random[4]) {
        target += std_normal_lupdf(psi_iota_bar_z[1]);
        if (MS) {
          target += std_normal_lupdf(to_vector(psi_iota_s_z[1]));
        }
      }
      profile("likelihood") {
        target +=
          occARU_lupmf(y | Q, f_l, i_idx, I_k, j_idx, J_i, log_Delta, tau,
                       region, P, P_cat, P_ord, random, X1, X_cat1, X_ord1_p1,
                       X2, X_cat2, X_ord2_p1, X3, X_cat3, X_ord3_p1, alpha_bar,
                       alpha, lambda, psi_beta, psi_beta_cat,
                       psi_beta_ord_realised, mu_beta, mu_beta_cat,
                       mu_beta_ord_realised, gamma, gamma_cat,
                       gamma_ord_realised, iota, kappa, nu, psi_iota, psi_nu,
                       epsilon, phi);

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
  corr_matrix[MS && random[3] ? S : 0] nu_O;
  corr_matrix[MS && random[4] ? S : 0] psi_iota_O;
  corr_matrix[MS && random[5] ? S : 0] psi_nu_O;
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
    if (random[3]) {
      nu_O = multiply_lower_tri_self_transpose(nu_O_L);
    }
    if (random[4]) {
      psi_iota_O = multiply_lower_tri_self_transpose(psi_iota_O_L);
    }
    if (random[5]) {
      psi_nu_O = multiply_lower_tri_self_transpose(psi_nu_O_L);
    }
    if (OLRE) {
      epsilon_O = multiply_lower_tri_self_transpose(epsilon_O_L);
    }
  }

  // site effects
  vector[random[1] * I] iota_bar;
  matrix[random[1] * I, S] iota;
  vector[random[4] * I] psi_iota_bar;
  matrix[random[4] * I, S] psi_iota;

  // restricted coefficients
  vector[random[1] * P[2]] mu_beta_bar2;
  matrix[random[1] * P[2], S] mu_beta2;
  vector[random[1] * P_ord[2]] mu_beta_ord_bar2;
  matrix[random[1] * P_ord[2], S] mu_beta_ord2;
  vector[random[2] * P[3]] gamma_bar2;
  matrix[random[2] * P[3], S] gamma2;
  vector[random[2] * P_ord[3]] gamma_ord_bar2;
  matrix[random[2] * P_ord[3], S] gamma_ord2;
  if (random[2]) {
    vector[J_ksum] kappa_bar_long;
    matrix[MS * J_ksum, S] kappa_long;
    for (k in 1:K) {
      int J_kk = J_k[k], idx = sum(J_k[:k - 1]);
      kappa_bar_long[idx + 1:idx + J_kk] = kappa_bar[:J_kk, k];
      kappa_long[idx + 1:idx + J_kk] = kappa[k, :J_kk];
    }
    if (P[3]) {
      gamma_bar2 = gamma_bar + X3_plus[:P[3]] * kappa_bar_long;
      if (MS) {
        gamma2 = gamma + X3_plus[:P[3]] * kappa_long;
      } else {
        gamma2[:, 1] = gamma_bar2;
      }
    }
    if (P_ord[3]) {
      int idx = P_aug[3] - P_ord[3] + 1;
      gamma_ord_bar2 = gamma_ord_bar + X3_plus[idx:] * kappa_bar_long;
      if (MS) {
        gamma_ord2 = gamma_ord + X3_plus[idx:] * kappa_long;
      } else {
        gamma_ord2[:, 1] = gamma_ord_bar2;
      }
    }
  }
  vector[random[4] * P[1]] psi_beta_bar2;
  matrix[random[4] * P[1], S] psi_beta2;
  vector[random[4] * P_ord[1]] psi_beta_ord_bar2;
  matrix[random[4] * P_ord[1], S] psi_beta_ord2;

  // log likelihood, latent occupancy, and posterior predictions
  array[dyn * K] matrix[I, S] eta;
  matrix[I, S] logit_psi;
  matrix[S, I] log_lik;
  matrix[MC * S, I] log_lik2;
  array[latent * K, I, S] int z;
  array[PPC_y * K, I, J, S] int yrep;
  array[PPC_Q * K, I, S] int Qrep;

  profile("generated_quantities") {

    // site effects and unconditional coefficients
    matrix[MS * S, S] iota_S_U;
    array[(kernel[1] > 0) * R, 1 + SS[1] * S] matrix[I_max, I_max] iota_L;
    matrix[MS * random[1] * I, S] iota_s;
    if (random[1]) {
      iota_bar = iota_bar_t[1] * iota_bar_z[1];
      if (MS) {
        iota_S_U = diag_post_multiply(iota_O_L', iota_t);
        iota_s = iota_s_z[1] * iota_S_U;
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
      if (P[2]) {
        mu_beta_bar2 = mu_beta_bar + X2_plus[:P[2]] * iota_bar;
        if (MS) {
          mu_beta2 = mu_beta + X2_plus[:P[2]] * iota;
        } else {
          mu_beta2[:, 1] = mu_beta_bar2;
        }
      }
      if (P_ord[2]) {
        int idx = P_aug[2] - P_ord[2] + 1;
        mu_beta_ord_bar2 = mu_beta_ord_bar + X2_plus[idx:] * iota_bar;
        if (MS) {
          mu_beta_ord2 = mu_beta_ord + X2_plus[idx:] * iota;
        } else {
          mu_beta_ord2[:, 1] = mu_beta_ord_bar2;
        }
      }
    }

    // occupancy site effects and coefficients
    matrix[MS * S, S] psi_iota_S_U;
    array[(kernel[4] > 0) * R, 1 + SS[4] * S] matrix[I_max, I_max] psi_iota_L;
    matrix[MS * random[4] * I, S] psi_iota_s;
    if (random[4]) {
      psi_iota_bar = psi_iota_bar_t[1] * psi_iota_bar_z[1];
      if (MS) {
        psi_iota_S_U = diag_post_multiply(psi_iota_O_L', psi_iota_t);
        psi_iota_s = psi_iota_s_z[1] * psi_iota_S_U;
      }
      if (kernel[4]) {
        psi_iota_L = site_kernels(I_r, XY, kernel[4], psi_iota_ell);
        psi_iota_bar = increment_spatial(I_r, psi_iota_L[:, 1], psi_iota_bar);
        if (MS) {
          psi_iota_s = increment_spatial(I_r, psi_iota_L, psi_iota_s);
        }
      }
      if (MS) {
        psi_iota = rep_matrix(psi_iota_bar, S) + psi_iota_s;
      } else {
        psi_iota[:, 1] = psi_iota_bar;
      }
      if (P[1]) {
        psi_beta_bar2 = psi_beta_bar + X1_plus[:P[1]] * psi_iota_bar;
        if (MS) {
          psi_beta2 = psi_beta + X1_plus[:P[1]] * psi_iota;
        } else {
          psi_beta2[:, 1] = psi_beta_bar2;
        }
      }
      if (P_ord[1]) {
        int idx = P_aug[1] - P_ord[1] + 1;
        psi_beta_ord_bar2 = psi_beta_ord_bar + X1_plus[idx:] * psi_iota_bar;
        if (MS) {
          psi_beta_ord2 = psi_beta_ord + X1_plus[idx:] * psi_iota;
        } else {
          psi_beta_ord2[:, 1] = psi_beta_ord_bar2;
        }
      }
    }

    // reconstruct log likelihood and recover occupancy
    vector[OLRE * J_sum] epsilon_bar;
    matrix[MS * OLRE * S, S] epsilon_S_U;
    matrix[OLRE * J_sum, S] epsilon_mat;
    array[OLRE * K, I] matrix[J, S] epsilon;
    if (OLRE) {
      epsilon_bar = epsilon_bar_t[1] * epsilon_bar_z[1];
      if (MS) {
        epsilon_S_U = diag_post_multiply(epsilon_O_L', epsilon_t);
        epsilon_mat = rep_matrix(epsilon_bar, S)
                      + epsilon_s_z[1] * epsilon_S_U;
      } else {
        epsilon_mat[:, 1] = epsilon_bar;
      }
      epsilon = fill_epsilon(j_idx, J_i, epsilon_mat);
    }
    tuple(matrix[S, I], array[I, S] matrix[2, K], array[K] matrix[I, S],
          matrix[I, S], array[dyn * I, S, Km1] matrix[2, 2],
          array[K, I] matrix[J, S]) lp =
      occARU(y, Q, f_l, i_idx, I_k, j_idx, J_i, log_Delta, tau, region, P,
             P_cat, P_ord, random, X1, X_cat1, X_ord1_p1, X2, X_cat2, X_ord2_p1,
             X3, X_cat3, X_ord3_p1, alpha_bar, alpha, lambda, psi_beta,
             psi_beta_cat, psi_beta_ord_realised, mu_beta, mu_beta_cat,
             mu_beta_ord_realised, gamma, gamma_cat, gamma_ord_realised, iota,
             kappa, nu, psi_iota, psi_nu, epsilon, phi);
    log_lik = lp.1;
    if (latent) {
      z = occARU_rng(Q, f_l, lp.2, lp.5);
    }
    if (dyn) {
      eta = lp.3;
    }
    logit_psi = lp.4;

    // posterior predictions
    if (PPC_y || PPC_Q) {
      tuple(array[K, I, J, S] int, array[K, I, S] int) pp =
        pp_occARU_rng(f_l, j_idx, J_i, logit_psi, lp.5, lp.6, phi);
      if (PPC_y) {
        yrep = pp.1;
      }
      if (PPC_Q) {
        Qrep = pp.2;
      }
    }

    // Monte Carlo integration of observation-level parameters for loo
    if (MC) {
      array[K] matrix[I, S] eta_d;
      matrix[I, S] logit_psi_d;
      row_vector[S] log_psi;
      array[Km1] matrix[2, S] log_q;
      array[Km1] matrix[2, 2] log_P;
      matrix[2, K] Omega;
      array[K] matrix[S, I] lp_y;
      array[D] matrix[S, I] log_lik_d;

      // posterior predict site effects and OLREs
      matrix[Im1, MS ? Sm1 : 1] u;
      vector[random[1] * I] iota_bar_rep;
      matrix[random[1] * I, S] iota_rep;
      vector[random[4] * I] psi_iota_bar_rep;
      matrix[random[4] * I, S] psi_iota_rep;
      vector[OLRE * J_sum] epsilon_bar_rep;
      array[OLRE * K, I] matrix[J, S] epsilon_rep;
      for (d in 1:D) {
        if (random[1]) {
          u[:, 1] = to_vector(normal_rng(zeros_vector(Im1), 1));
          iota_bar_rep = iota_bar_t[1] * sum_to_zero_constrain(u[:, 1]);
          if (MS) {
            u = to_matrix(normal_rng(zeros_vector(Im1 * Sm1), 1), Im1, Sm1);
            iota_s = sum_to_zero_constrain(u) * iota_S_U;
            iota_rep = rep_matrix(iota_bar, S) + iota_s;
          } else {
            iota_rep[:, 1] = iota_bar_rep;
          }
          if (kernel[1]) {
            iota_rep = increment_spatial(I_r, iota_L, iota_rep);
          }
        }
        if (random[4]) {
          u[:, 1] = to_vector(normal_rng(zeros_vector(Im1), 1));
          psi_iota_bar_rep = psi_iota_bar_t[1] * sum_to_zero_constrain(u[:, 1]);
          if (MS) {
            u = to_matrix(normal_rng(zeros_vector(Im1 * Sm1), 1), Im1, Sm1);
            psi_iota_s = sum_to_zero_constrain(u) * psi_iota_S_U;
            psi_iota_rep = rep_matrix(psi_iota_bar, S) + psi_iota_s;
          } else {
            psi_iota_rep[:, 1] = psi_iota_bar;
          }
          if (kernel[4]) {
            psi_iota_rep = increment_spatial(I_r, psi_iota_L, psi_iota_rep);
          }
        }
        if (OLRE) {
          int N = J_sum - 1;
          matrix[N, MS ? Sm1 : 1] epsilon_u;
          epsilon_u[:, 1] = to_vector(normal_rng(zeros_vector(N), 1));
          epsilon_bar_rep = epsilon_bar_t[1]
                            * sum_to_zero_constrain(epsilon_u[:, 1]);
          if (MS) {
            epsilon_u = to_matrix(normal_rng(zeros_vector(N), 1), N, Sm1);
            epsilon_mat = rep_matrix(epsilon_bar_rep, S)
                          + sum_to_zero_constrain(epsilon_u) * epsilon_S_U;
          } else {
            epsilon_mat[:, 1] = epsilon_bar_rep;
          }
          epsilon_rep = fill_epsilon(j_idx, J_i, epsilon_mat);
        }

        // compute likelihood
        eta_d = lp.3;
        for (k in 1:K) {
          int I_kk = I_k[k];
          array[I_kk] int idx_i = i_idx[k, :I_kk];
          if (random[4]) {
            eta_d[k, idx_i] += psi_iota_rep[idx_i] - psi_iota[idx_i];
          }
          for (i in idx_i) {
            int J_ki = J_i[k, i];
            array[J_ki] int idx_j = j_idx[k, i, :J_ki];
            matrix[J_ki, S] log_mu = lp.6[k, i, idx_j];
            if (random[1]) {
              log_mu += rep_matrix(iota_rep[i] - iota[i], J_ki);
            }
            if (OLRE) {
              log_mu += epsilon_rep[k, i, idx_j] - epsilon[k, i, idx_j];
            }
            array[J_ki, S] int y_i = y[k, i, idx_j];
            for (s in 1:S) {
              lp_y[k, s, i] = NB ?
                neg_binomial_2_log_lpmf(y_i[:, s] | log_mu[:, s], phi[s])
                : poisson_log_lpmf(y_i[:, s] | log_mu[:, s]);
            }
          }
        }
        if (K == 1) {
          logit_psi_d = eta_d[1];
        }
        for (i in 1:I) {
          int f = f_l[i, 1], l = f_l[i, 2];
          if (dyn) {
            logit_psi_d[i] = alpha_bar[1] + lambda[1] .* eta_d[f, i];
            for (k in f + 1:l) {
              int km1 = k - 1;
              log_q[km1, 1] = alpha_bar[2] + lambda[2] .* eta_d[k, i];
              log_q[km1, 2] = alpha_bar[3] - lambda[3] .* eta_d[k, i];
            }
          } else if (K > 1) {
            logit_psi_d[i] = eta_d[f, i];
          }
          log_psi = log_inv_logit(logit_psi_d[i]);
          for (s in 1:S) {
            Omega[:, f] = [ Q[f, i, s] ?
                            negative_infinity()
                            : log_psi[s] - logit_psi_d[i, s],
                            log_psi[s] + lp_y[f, s, i] ]';
            for (k in f + 1:l) {
              int km1 = k - 1;
              log_P[km1] = log_matrix_exp(log_q[km1, :, s], tau[km1, i]);
              Omega[:, k] = Q[k, i, s] ?
                            [ negative_infinity(),
                              Q[km1, i, s] ?
                              Omega[2, km1] + log_P[km1, 2, 2]
                              : log_sum_exp(Omega[:, km1] + log_P[km1, :, 2]) ]'
                            : Q[km1, i, s] ?
                              Omega[2, km1] + log_P[km1, 2]'
                              : log_prod_exp(log_P[km1]', Omega[:, km1]);
              if (J_i[k, i]) {
                Omega[2, k] += lp_y[k, s, i];
              }
            }
            log_lik_d[d, s, i] = Q[l, i, s] ?
                                 Omega[2, l] : log_sum_exp(Omega[:, l]);
          }
        }
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
