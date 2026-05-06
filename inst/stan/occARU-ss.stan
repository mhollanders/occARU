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
  matrix<lower=0, upper=1>[I, J] Delta;
  array[I, J] int<lower=0> y;

  // site UTMs (km) and regions
  array[I] vector[2] XY;
  int<lower=1> R;
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

  // indicators and scalars
  array[2] int<lower=0, upper=1> random,  // random effects
                                 SS;  // species length scales
  array[2] int<lower=0, upper=3> kernel;  // GP kernels
  real<lower=0> period;  // period for periodic kernel
  int<lower=0, upper=2> OD;  // overdispersion (Poisson, OLRE, negbin)
  int<lower=0, upper=1> dirichlet,  // var. decomp. (logi-norm or Dirichlet)
                        project_kappa,  // orth. proj. survey effects
                        latent,  // recover latent occ. states
                        PPC_y,  // posterior predictions of y
                        PPC_Q;  // posterior predictions of Q
  int<lower=0, upper=I> grainsize;  // grainsize for threading
  int<lower=0> D;  // Monte Carlo draws for LOO

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
  int MR = R > 1, periodic = period > 0 && kernel[2] > 0, OLRE = OD == 1,
      NB = OD == 2, MC = (D > 0) * (kernel[1] || OLRE), Im1 = I - 1;
  real log_D = log(D);
  array[3] int P_sum;
  for (d in 1:3) {
    P_sum[d] = P[d] + P_cat[d] + P_ord[d];
  }

  // survey indices, total surveys, and offsets
  tuple(array[I, J] int, array[I] int) surv_indices = survey_indices(Delta);
  array[I, J] int j_idx = surv_indices.1;
  array[I] int J_i = surv_indices.2;
  int J_sum = sum(J_i);
  matrix[J, I] log_Delta = log(Delta');

  // region indices
  tuple(array[R, I] int, array[R] int) r_indices = region_indices(region);
  array[R] int I_r = r_indices.2;
  int I_max = max(I_r);
  array[R, I_max] int r_idx = r_indices.1[:, :I_max];

  // aggregated counts
  array[I] int Q;
  for (i in 1:I) {
    Q[i] = sum(y[i, j_idx[i, :J_i[i]]]);
  }

  // site, region, and survey sequences
  array[I] int sites = linspaced_int_array(I, 1, I);
  array[R] int regions = linspaced_int_array(R, 1, R);
  array[J] real surveys = linspaced_array(J, 1, J);

  // var. partitions for occ. and det.
  int psi_V = MR + P_sum[1],
      mu_V = MR + sum(P_sum[2:3]) + sum(random) + OLRE;

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

  // augmented site-level design matrix and pseudo-inverse for projection
  array[3] int P_aug;
  for (d in 2:3) {
    P_aug[d] = P[d] + C_sum[d] - P_cat[d] + P_ord[d];
  }
  matrix[I, P_aug[2]] X2_aug =
    augment_design_matrix(X2, X_cat2, X_ord2, C[2, :P_cat[2]], O[2, :P_ord[2]]);
  matrix[P_aug[2], I] X2_plus = pseudo_inverse(X2_aug);

  // site-by-survey design matrix pseudo-inverses for orthogonal projection
  array[I] matrix[J, P_aug[3]] X3_aug =
    augment_design_matrix(X3, X_cat3, X_ord3, C[3, :P_cat[3]], O[3, :P_ord[3]]);
  matrix[J, P_aug[3]] X3_mean = average_survey_design_matrix(X3_aug, j_idx);
  matrix[P_aug[3], J] X3_plus = pseudo_inverse(X3_mean);
}

parameters {
  // intercepts
  real<lower=0, upper=1> psi_bar;
  real<lower=0> mu_bar;
  cholesky_factor_corr[MR * 2] alpha_r_O_L;
  array[MR * 2] sum_to_zero_vector[R] alpha_r_z;

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
  array[random[1]] sum_to_zero_vector[I] iota_z;
  vector<lower=0>[kernel[1] > 0] iota_ell;
  array[random[2]] sum_to_zero_vector[J] kappa_z;
  vector<lower=0>[(kernel[2] > 0) + periodic] kappa_ell;
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
  row_vector[2] alpha = [ logit(psi_bar), log(mu_bar) ];
  matrix[R, 2] alpha_r;
  if (MR) {
    alpha_r = append_col(alpha_r_z[1], alpha_r_z[2])
              * diag_post_multiply(alpha_r_O_L', [ psi_tau[1], mu_tau[1] ]);
  }

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

  // conditional and unconditional survey effects
  vector[random[2] * J] kappa;
  vector[random[2] * P_sum[3] * project_kappa ? J : 0] kappa2;

  // random scales
  vector[random[1]] iota_t;
  vector[random[2] + periodic] kappa_t;
  vector[OLRE] epsilon_t;
  if (OLRE) {
    epsilon_t[1] = mu_tau[mu_V];
  }
  {
    int psi_idx = MR + 1, mu_idx = MR + 1;

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
      tuple(vector[P_o], matrix[O_max[2], P_o]) coef =
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
      tuple(vector[P_o], matrix[O_max[3], P_o]) coef =
        coef_ord_jacobian(segment(mu_tau, mu_idx, P_o), gamma_ord_z,
                          gamma_ord_c_u, O[3, :P_o]);
      gamma_ord = coef.1;
      gamma_ord_cs = coef.2;
      mu_idx += P_o;
    }

    // site effects
    if (random[1]) {
      iota_t[1] = mu_tau[mu_idx];
      mu_idx += 1;
    }

    // survey effects
    if (random[2]) {
      kappa_t[1] = mu_tau[mu_idx];
      kappa = kappa_t[1] * kappa_z[1];
      if (kernel[2]) {
        vector[periodic * 2] kappa_v;
        if (periodic) {
          kappa_v = sqrt(K_phi[1]);
        }
        kappa = survey_kernels(surveys, kernel[2], period, kappa_v, kappa_ell)
                * kappa;
      }

      // increment orthogonal projection
      if (P_sum[3] && project_kappa) {
        kappa2 = kappa;
        kappa = orthogonalise(X3_mean, X3_plus, kappa);
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
    lprior += inv_gamma_lpdf(iota_ell[1] | iota_ell_inv_gamma[1],
                                           iota_ell_inv_gamma[2]);
  }
  if (kernel[2]) {
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
  if (MR) {
    target += std_normal_lupdf(append_row(alpha_r_z[1], alpha_r_z[2]));
  }
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
  if (random[2]) {
    target += std_normal_lupdf(kappa_z[1]);
  }

  // intercepts
  matrix[I, 2] alpha_mat = rep_matrix(alpha, I);
  if (MR) {
    alpha_mat += alpha_r[region];
  }

  // parallelise over regions
  if (grainsize && MR) {
    target +=
      reduce_sum(partial_regions_occARU_lupmf, regions, grainsize, r_idx, I_r,
                 XY, kernel[1], X2_aug, X2_plus, y, Q, j_idx, J_i, log_Delta,
                 X1, X_cat1, X_ord1, X2, X_cat2, X_ord2, X3, X_cat3, X_ord3,
                 alpha_mat, psi_beta, psi_beta_cat,
                 coef_ord_realised(O[1], psi_beta_ord, psi_beta_ord_cs),
                 mu_beta, mu_beta_cat, coef_ord_realised(O[2], mu_beta_ord,
                 mu_beta_ord_cs), gamma, gamma_cat, coef_ord_realised(O[3],
                 gamma_ord, gamma_ord_cs), iota_z, iota_t, iota_ell, kappa,
                 epsilon_z, epsilon_t, phi);
  } else {

    // random site effects with orthogonal projection
    vector[random[1] * I] iota;
    if (random[1]) {
      iota = iota_t[1] * iota_z[1];
      if (kernel[1]) {
        array[R] matrix[I_max, I_max] iota_L =
          site_kernels(r_idx, I_r, XY, kernel[1], iota_ell[1]);
        iota = increment_spatial(r_idx, I_r, iota_L, iota);
      }
      if (P_sum[2]) {
        iota = orthogonalise(X2_aug, X2_plus, iota);
      }
    }

    // parallelise over sites
    if (grainsize) {
      target += reduce_sum(partial_sites_occARU_lupmf, sites, grainsize, y, Q,
                           j_idx, J_i, log_Delta, X1, X_cat1, X_ord1, X2,
                           X_cat2, X_ord2, X3, X_cat3, X_ord3, alpha_mat,
                           psi_beta, psi_beta_cat, coef_ord_realised(O[1],
                           psi_beta_ord, psi_beta_ord_cs), mu_beta, mu_beta_cat,
                           coef_ord_realised(O[2], mu_beta_ord, mu_beta_ord_cs),
                           gamma, gamma_cat, coef_ord_realised(O[3], gamma_ord,
                           gamma_ord_cs), iota, iota_z, kappa, epsilon_z,
                           epsilon_t, phi);

      // no parallelisation
    } else {
      if (random[1]) {
        target += std_normal_lupdf(iota_z[1]);
      }

      // Poisson OLREs
      matrix[OLRE * J, I] epsilon;
      if (OLRE) {
        epsilon = fill_epsilon(j_idx, J_i, epsilon_t[1] * epsilon_z[1]);
        target += std_normal_lupdf(epsilon_z[1]);
      }
      target += occARU_lupmf(y | Q, j_idx, J_i, log_Delta, X1, X_cat1, X_ord1,
                             X2, X_cat2, X_ord2, X3, X_cat3, X_ord3, alpha_mat,
                             psi_beta, psi_beta_cat, coef_ord_realised(O[1],
                             psi_beta_ord, psi_beta_ord_cs), mu_beta,
                             mu_beta_cat, coef_ord_realised(O[2], mu_beta_ord,
                             mu_beta_ord_cs), gamma, gamma_cat,
                             coef_ord_realised(O[3], gamma_ord, gamma_ord_cs),
                             iota, kappa, epsilon, phi);
    }
  }
}

generated quantities {
  // unconditional coefficients
  vector[random[1] * P[2]] mu_beta2;
  vector[random[2] * P_ord[2]] mu_beta_ord2;
  vector[random[2] * P[3] * project_kappa] gamma2;
  vector[random[2] * P_ord[3] * project_kappa] gamma_ord2;
  if (random[2] * project_kappa) {
    if (P[3]) {
      gamma2 = uncondition(X3_plus[:P[3]], gamma, kappa2);
    }
    if (P_ord[3]) {
      gamma_ord2 = uncondition(X3_plus[P_aug[3] - P_ord[3] + 1:], gamma_ord,
                               kappa2);
    }
  }

  // log likelihood, latent occupancy, and posterior predictions
  vector[random[1] * I] iota, iota2;
  vector[I] log_lik, logit_psi;
  vector[MC * I] log_lik2;
  array[latent * I] int z;
  array[PPC_y * I, J] int yrep;
  array[PPC_Q * I] int Qrep;

  {
    // site effects and unconditional coefficients
    array[(kernel[1] > 0) * R] matrix[I_max, I_max] iota_L;
    if (random[1]) {
      iota = iota_t[1] * iota_z[1];
      if (kernel[1]) {
        iota_L = site_kernels(r_idx, I_r, XY, kernel[1], iota_ell[1]);
        iota = increment_spatial(r_idx, I_r, iota_L, iota);
      }
      if (P_sum[2]) {
        iota2 = iota;
        iota = orthogonalise(X2_aug, X2_plus, iota);
        if (P[2]) {
          mu_beta2 = uncondition(X2_plus[:P[2]], mu_beta, iota2);
        }
        if (P_ord[2]) {
          mu_beta_ord2 = uncondition(X2_plus[P_aug[2] - P_ord[2] + 1:],
                                     mu_beta_ord, iota2);
        }
      }
    }

    // reconstruct log likelihood and recover occupancy
    matrix[I, 2] alpha_mat = rep_matrix(alpha, I);
    if (MR) {
      alpha_mat += alpha_r[region];
    }
    matrix[OLRE * J, I] epsilon;
    if (OLRE) {
      epsilon = fill_epsilon(j_idx, J_i, epsilon_t[1] * epsilon_z[1]);
    }
    tuple(vector[I], matrix[2, I], vector[I], matrix[J, I]) lp =
      occARU(y, Q, j_idx, J_i, log_Delta, X1, X_cat1, X_ord1, X2, X_cat2,
             X_ord2, X3, X_cat3, X_ord3, alpha_mat, psi_beta, psi_beta_cat,
             coef_ord_realised(O[1], psi_beta_ord, psi_beta_ord_cs), mu_beta,
             mu_beta_cat, coef_ord_realised(O[2], mu_beta_ord, mu_beta_ord_cs),
             gamma, gamma_cat, coef_ord_realised(O[3], gamma_ord, gamma_ord_cs),
             iota, kappa, epsilon, phi);
    log_lik = lp.1;
    if (latent) {
      z = occARU_rng(Q, lp);
    }
    logit_psi = lp.3;

    // posterior predictions
    if (PPC_y || PPC_Q) {
      if (OLRE) {

  matrix[N, S] epsilon_mat;
  if (S > 1) {
    int Sm1 = S - 1, NSm1 = Nm1 * Sm1;
    matrix[Nm1, Sm1] u = to_matrix(normal_rng(zeros_vector(NSm1), 1), Nm1, Sm1);
    epsilon_mat = rep_matrix(bar, S) + sum_to_zero_constrain(u) * S_L';
  } else {
    vector[Nm1] u = to_vector(normal_rng(zeros_vector(Nm1), 1));
    epsilon_mat = rep_matrix(sum_to_zero_constrain(u) * tau, S);
  }
  return fill_epsilon(j_idx, J_i, epsilon_mat);
        epsilon = pp_epsilon_rng(j_idx, J_i, mu_tau[mu_V]);
      }
      tuple(array[I, J] int, array[I] int) pp =
        pp_occARU_rng(j_idx, J_i, X3, lp.3, lp.4, gamma, epsilon, phi);
      if (PPC_y) {
        yrep = pp.1;
      }
      if (PPC_Q) {
        Qrep = pp.2;
      }
    }

    // Monte Carlo integration of observation-level parameters for loo
    if (MC) {
      matrix[I, D] log_lik_d;

      // posterior predict site effects and OLREs and compute log likelihood
      vector[Im1] zeros = zeros_vector(Im1);
      vector[random[1] * I] iota_rep;
      for (d in 1:D) {
        if (random[1]) {
          iota_rep = iota_t[1]
                     * sum_to_zero_constrain(to_vector(normal_rng(zeros, 1)));
          if (kernel[1]) {
            iota_rep = increment_spatial(r_idx, I_r, iota_L, iota_rep);
          }
        }
        if (OLRE) {
          epsilon = pp_epsilon_rng(j_idx, J_i, epsilon_t[1]);
        }
        lp = occARU(y, Q, j_idx, J_i, log_Delta, X1, X_cat1, X_ord1, X2, X_cat2,
                    X_ord2, X3, X_cat3, X_ord3, alpha_mat, psi_beta,
                    psi_beta_cat, coef_ord_realised(O[1], psi_beta_ord,
                    psi_beta_ord_cs), mu_beta, mu_beta_cat,
                    coef_ord_realised(O[2], mu_beta_ord, mu_beta_ord_cs), gamma,
                    gamma_cat, coef_ord_realised(O[3], gamma_ord, gamma_ord_cs),
                    iota_rep, kappa, epsilon, phi);
        log_lik_d[:, d] = lp.1;
      }
      log_lik2 = rep_vector(-log_D, I);
      for (i in 1:I) {
        log_lik2[i] += log_sum_exp(log_lik_d[i]);
      }
    }
  }
}
