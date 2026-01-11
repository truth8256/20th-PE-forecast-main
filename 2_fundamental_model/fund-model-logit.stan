data {
  int<lower=2> P;                 // number of provinces
  int<lower=1> TT;                // number of observed elections (fit period)
  int<lower=1> N;                 // N = P * TT
  int<lower=1> K;                 // # national-level coeffs
  int<lower=1> L;                 // # province-level coeffs

  matrix[N, K] X;                 // stacked by time blocks of size P
  matrix[N, L] Z;

  matrix[P, K] X_pred;            // predictors for TT+1
  matrix[P, L] Z_pred;

  vector<lower=0, upper=1>[N] Y;  // vote share (two-party), fit period only

  vector[P] Alpha_init;           // initial alpha mean for t=1 (free P-1 used)
  matrix[P, TT + 1] Pop_weight;   // weights for constraint (fit + pred)
}

transformed data {
  vector[N] Y_logit;
  real eps;

  eps = 1e-6;

  for (n in 1:N) {
    real y_adj = fmin(1 - eps, fmax(eps, Y[n]));
    Y_logit[n] = logit(y_adj);
  }
}

parameters {
  // -------- state: province effects (RW, non-centered) --------
  // [CHANGE] only fit period t=1..TT (no TT+1 state in parameters)
  matrix[P - 1, TT] z_alpha;

  // -------- regression coefficients --------
  vector[K] beta;
  vector[L] gamma;

  // -------- national shock (RW, anchored) --------
  // delta[1]=0, then increments for t=2..TT
  vector[TT - 1] z_delta;

  // -------- scales --------
  real<lower=0> sigma_alpha;       // RW innovation sd for alpha
  real<lower=0> sigma_delta;       // RW innovation sd for delta
  real<lower=0> sigma_epsilon;     // observation noise on logit scale
}

transformed parameters {
  matrix[P - 1, TT] alpha0;
  matrix[P, TT] alpha;

  vector[TT] delta;

  // ---- alpha RW (non-centered), fit period only ----
  alpha0[, 1] = Alpha_init[1:(P - 1)] + sigma_alpha * z_alpha[, 1];
  for (t in 2:TT) {
    alpha0[, t] = alpha0[, t - 1] + sigma_alpha * z_alpha[, t];
  }

  // ---- apply weighted-mean constraint: sum_p w_{p,t} * alpha_{p,t} = 0 ----
  alpha[1:(P - 1), ] = alpha0;
  for (t in 1:TT) {
    alpha[P, t] =
      -dot_product(Pop_weight[1:(P - 1), t], alpha0[, t]) / Pop_weight[P, t];
  }

  // ---- delta RW (anchored), fit period only ----
  delta[1] = 0;
  for (t in 2:TT) {
    delta[t] = delta[t - 1] + sigma_delta * z_delta[t - 1];
  }
}

model {
  // ---------------- priors ----------------
  to_vector(z_alpha) ~ normal(0, 1);
  z_delta ~ normal(0, 1);

  beta  ~ normal(0, 0.5);
  gamma ~ normal(0, 0.5);

  sigma_alpha   ~ normal(0, 0.2);
  sigma_delta   ~ normal(0, 0.2);
  sigma_epsilon ~ normal(0, 0.3);

  // --------------- likelihood (fit period only) ---------------
  for (t in 1:TT) {
    int start = 1 + (t - 1) * P;
    int end_  = t * P;

    matrix[P, K] Xt = X[start:end_, ];
    matrix[P, L] Zt = Z[start:end_, ];

    Y_logit[start:end_] ~ normal(
      alpha[, t] + Xt * beta + Zt * gamma + rep_vector(delta[t], P),
      sigma_epsilon
    );
  }
}

generated quantities {
  // ---- pointwise log-lik for LOO etc. (fit period only) ----
  vector[N] log_lik;

  // ---- 1-step ahead forecast states ----
  vector[P] alpha_TT1;
  real delta_TT1;

  // ---- prediction ----
  vector[P] y_pred_logit;
  vector[P] y_pred_share;

  // log-lik
  {
    vector[N] xb = X * beta;
    vector[N] zg = Z * gamma;

    for (t in 1:TT) {
      int start = 1 + (t - 1) * P;
      int end_  = t * P;
      for (n in start:end_) {
        log_lik[n] = normal_lpdf(
          Y_logit[n] |
          alpha[n - start + 1, t] + xb[n] + zg[n] + delta[t],
          sigma_epsilon
        );
      }
    }
  }

  // ---- alpha forecast at TT+1 ----
  // free P-1 provinces
  {
    vector[P - 1] alpha0_TT1_free;
    vector[P - 1] z_alpha_pred;

    // draw new RW shocks for prediction time
    for (p in 1:(P - 1)) z_alpha_pred[p] = normal_rng(0, 1);

    // 1-step RW: alpha0_free(TT+1) = alpha0_free(TT) + sigma_alpha * shock
    for (p in 1:(P - 1)) {
      alpha0_TT1_free[p] = alpha0[p, TT] + sigma_alpha * z_alpha_pred[p];
      alpha_TT1[p] = alpha0_TT1_free[p];
    }

    // constrained province P to satisfy weighted-mean 0 at time TT+1
    alpha_TT1[P] =
      -dot_product(Pop_weight[1:(P - 1), TT + 1], alpha0_TT1_free)
      / Pop_weight[P, TT + 1];
  }

  // ---- delta forecast at TT+1 ----
  // delta(TT+1) = delta(TT) + sigma_delta * shock
  delta_TT1 = delta[TT] + sigma_delta * normal_rng(0, 1);

  // ---- prediction mean and back-transform ----
  y_pred_logit =
    alpha_TT1
    + X_pred * beta
    + Z_pred * gamma
    + rep_vector(delta_TT1, P);

  for (p in 1:P)
    y_pred_share[p] = inv_logit(y_pred_logit[p]);
}
