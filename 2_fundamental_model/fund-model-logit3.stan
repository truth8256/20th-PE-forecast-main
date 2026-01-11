data {
  int<lower=2> P;
  int<lower=1> TT;
  int<lower=1> N;                 // N = P * TT
  int<lower=1> K;
  int<lower=1> L;

  matrix[N, K] X;                 // stacked by time blocks of size P
  matrix[N, L] Z;

  matrix[P, K] X_pred;
  matrix[P, L] Z_pred;

  vector<lower=0, upper=1>[N] Y;

  vector[P] Alpha_init;
  matrix[P, TT + 1] Pop_weight;   // provided, but we will build w_bar for fit constraint
}

transformed data {
  vector[N] Y_logit;
  real eps;
  vector[P] w_bar;                // ✅ fixed weights for constraint in fit period

  eps = 1e-4;

  // continuity correction: avoid logit(0)/logit(1)
  for (n in 1:N) {
    real y_adj = fmin(1 - eps, fmax(eps, Y[n]));
    Y_logit[n] = logit(y_adj);
  }

  // ✅ compute mean weights over fit period (t=1..TT)
  // (pred period TT+1 excluded from w_bar by design)
  w_bar = rep_vector(0.0, P);
  for (t in 1:TT) {
    w_bar += Pop_weight[, t];
  }
  w_bar /= TT;
}

parameters {
  // non-centered RW increments
  matrix[P - 1, TT + 1] z_alpha;

  vector[K] beta;
  vector[L] gamma;

  vector[TT] delta;

  real<lower=0> sigma_alpha;
  real<lower=0> sigma_epsilon;    // (오차 단순화)
}

transformed parameters {
  matrix[P - 1, TT + 1] alpha0;
  matrix[P, TT + 1] alpha;

  // non-centered random walk
  alpha0[, 1] = Alpha_init[1:(P - 1)] + sigma_alpha * z_alpha[, 1];
  for (t in 2:(TT + 1)) {
    alpha0[, t] = alpha0[, t - 1] + sigma_alpha * z_alpha[, t];
  }

  // apply weighted-mean constraint
  alpha[1:(P - 1), ] = alpha0;

  // ✅ Fit period constraint uses fixed weights w_bar (no rotation across t)
  for (t in 1:TT) {
    alpha[P, t] =
      -dot_product(w_bar[1:(P - 1)], alpha0[, t]) / w_bar[P];
  }

  // ✅ Prediction period constraint: choose ONE option below

  // ---- Option A (recommended for stability): use same fixed weights w_bar
  alpha[P, TT + 1] =
    -dot_product(w_bar[1:(P - 1)], alpha0[, TT + 1]) / w_bar[P];

  // ---- Option B (original intention): use Pop_weight at TT+1
  // alpha[P, TT + 1] =
  //   -dot_product(Pop_weight[1:(P - 1), TT + 1], alpha0[, TT + 1]) / Pop_weight[P, TT + 1];
}

model {
  // priors (you can tune scales; these are reasonable “regularizing” defaults)
  to_vector(z_alpha) ~ normal(0, 1);

  beta  ~ normal(0, 0.3);
  gamma ~ normal(0, 0.3);
  delta ~ normal(0, 0.15);

  sigma_alpha   ~ normal(0, 0.6);
  sigma_epsilon ~ normal(0, 0.15);

  // likelihood (fit period only)
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
  vector[P] y_pred_logit;
  vector[P] y_pred_share;
  vector[N] log_lik;

  real delta_TTp1 = 0;

  y_pred_logit =
    alpha[, TT + 1] + X_pred * beta + Z_pred * gamma + rep_vector(delta_TTp1, P);

  for (p in 1:P)
    y_pred_share[p] = inv_logit(y_pred_logit[p]);

  // log-lik (fit period)
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
}
