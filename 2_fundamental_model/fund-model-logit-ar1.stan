data {
  int<lower=2> P;                 // provinces
  int<lower=1> TT;                // time points (fit only)
  int<lower=1> N;                 // N = P * TT
  int<lower=1> K;
  int<lower=1> L;

  matrix[N, K] X;
  matrix[N, L] Z;

  matrix[P, K] X_pred;
  matrix[P, L] Z_pred;

  vector<lower=0, upper=1>[N] Y;

  vector[P] Alpha_init;           // logit-scale initial mean
  matrix[P, TT + 1] Pop_weight;   // weights
}

transformed data {
  vector[N] Y_logit;
  real eps = 1e-6;

  for (n in 1:N) {
    real y_adj = fmin(1 - eps, fmax(eps, Y[n]));
    Y_logit[n] = logit(y_adj);
  }
}

parameters {
  // non-centered AR(1) innovations
  matrix[P - 1, TT + 1] z_alpha;

  real<lower=0, upper=1> rho;     // AR(1)
  real<lower=0> sigma_alpha;

  vector[K] beta;
  vector[L] gamma;
  vector[TT] delta;

  vector<lower=0>[TT] sigma_epsilon;
}

transformed parameters {
  matrix[P - 1, TT + 1] alpha0;
  matrix[P, TT + 1] alpha;

  // AR(1), non-centered
  alpha0[,1] =
    Alpha_init[1:(P - 1)]
    + sigma_alpha / sqrt(1 - square(rho)) * z_alpha[,1];

  for (t in 2:(TT + 1)) {
    alpha0[,t] =
      rho * alpha0[,t - 1]
      + sigma_alpha * z_alpha[,t];
  }

  // weighted-mean constraint (logit scale)
  alpha[1:(P - 1), ] = alpha0;
  for (t in 1:(TT + 1)) {
    alpha[P, t] =
      -dot_product(Pop_weight[1:(P - 1), t], alpha0[, t])
      / Pop_weight[P, t];
  }
}

model {
  // priors
  to_vector(z_alpha) ~ normal(0, 1);

  rho ~ beta(8, 2);                // weakly informative, persistent
  sigma_alpha ~ normal(0, 0.2);

  beta  ~ normal(0, 0.5);
  gamma ~ normal(0, 0.5);
  delta ~ normal(0, 0.5);

  sigma_epsilon ~ normal(0, 0.3);

  // likelihood
  for (t in 1:TT) {
    int start = 1 + (t - 1) * P;
    int end_  = t * P;

    matrix[P, K] Xt = X[start:end_, ];
    matrix[P, L] Zt = Z[start:end_, ];

    Y_logit[start:end_] ~ normal(
      alpha[, t]
      + Xt * beta
      + Zt * gamma
      + rep_vector(delta[t], P),
      sigma_epsilon[t]
    );
  }
}

generated quantities {
  vector[P] y_pred_logit;
  vector[P] y_pred_share;

  y_pred_logit =
    alpha[, TT + 1]
    + X_pred * beta
    + Z_pred * gamma;

  for (p in 1:P)
    y_pred_share[p] = inv_logit(y_pred_logit[p]);
}
