data {
  int<lower=2> P;                 // number of provinces
  int<lower=1> TT;                // number of times (not including the prediction)
  int<lower=1> N;                 // N = P * TT
  int<lower=1> K;                 // # national-level coeffs
  int<lower=1> L;                 // # province-level coeffs

  matrix[N, K] X;                 // stacked by time blocks of size P
  matrix[N, L] Z;

  matrix[P, K] X_pred;            // predictors for TT+1
  matrix[P, L] Z_pred;

  vector<lower=0, upper=1>[N] Y;  // vote share (two-party), fit period only

  vector[P] Alpha_init;           // initial alpha mean for t=1
  matrix[P, TT + 1] Pop_weight;   // weights for constraint (fit + pred)
}

transformed data {
  vector[N] Y_logit;
  real eps;                       // continuity correction amount

  // continuity correction: avoid logit(0) / logit(1)
  // (you may tune eps; 1e-6 is usually safe if Y never equals 0/1 exactly)
  eps = 1e-6;

  for (n in 1:N) {
    real y_adj = fmin(1 - eps, fmax(eps, Y[n]));
    Y_logit[n] = logit(y_adj);
  }
}

parameters {
  // non-centered RW increments
  matrix[P - 1, TT + 1] z_alpha;

  vector[K] beta;
  vector[L] gamma;
  vector[TT] delta;                // ❗ TT까지만

  real<lower=0> sigma_alpha;
  // real<lower=0> sigma_beta;
  // real<lower=0> sigma_gamma;
  // real<lower=0> sigma_delta;

  vector<lower=0>[TT] sigma_epsilon;
}


transformed parameters {
  matrix[P - 1, TT + 1] alpha0;
  matrix[P, TT + 1] alpha;

  // non-centered random walk construction
  alpha0[, 1] = Alpha_init[1:(P - 1)] + sigma_alpha * z_alpha[, 1];
  for (t in 2:(TT + 1)) {
    alpha0[, t] = alpha0[, t - 1] + sigma_alpha * z_alpha[, t];
  }

  // apply weighted-mean constraint
  alpha[1:(P - 1), ] = alpha0;
  for (t in 1:(TT + 1)) {
    alpha[P, t] =
      -dot_product(Pop_weight[1:(P - 1), t], alpha0[, t]) / Pop_weight[P, t];
  }
}

model {
  // priors
  to_vector(z_alpha) ~ normal(0,1);
  // delta ~ normal(0, sigma_delta);
  // beta  ~ normal(0, sigma_beta);
  // gamma ~ normal(0, sigma_gamma);
  beta  ~ normal(0, 0.5);
  gamma ~ normal(0, 0.5);
  delta ~ normal(0, 0.5);   // delta가 벡터면 동일


  sigma_alpha ~ normal(0, 0.2);
  // sigma_delta ~ normal(0, 0.5);
  // sigma_beta  ~ normal(0, 1.0);
  // sigma_gamma ~ normal(0, 1.0);
  sigma_epsilon ~ normal(0, 0.3);

  // likelihood (fit period only)
  for (t in 1:TT) {
    int start = 1 + (t - 1) * P;
    int end_  = t * P;

    matrix[P, K] Xt = X[start:end_, ];
    matrix[P, L] Zt = Z[start:end_, ];

    Y_logit[start:end_] ~ normal(
      alpha[, t] + Xt * beta + Zt * gamma + rep_vector(delta[t], P),
      sigma_epsilon[t]
    );
  }
}



generated quantities {
  vector[P] y_pred_logit;
  vector[P] y_pred_share;
  vector[N] log_lik;

  // 예측용 delta (prior draw)
  // real delta_TTp1 = normal_rng(0, sigma_delta);
  real delta_TTp1 = 0;


  // 예측 평균
  y_pred_logit =
    alpha[, TT + 1]
    + X_pred * beta
    + Z_pred * gamma
    + rep_vector(delta_TTp1, P);

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
          sigma_epsilon[t]
        );
      }
    }
  }
}

