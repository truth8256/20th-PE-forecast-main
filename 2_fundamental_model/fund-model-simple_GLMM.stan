data {
  int<lower=2> P;
  int<lower=1> TT;
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> L;

  matrix[N, K] X;
  matrix[N, L] Z;

  array[N] int<lower=0> Y_count;
  array[N] int<lower=0> N_count;

  matrix[P, K] X_pred;
  matrix[P, L] Z_pred;

  matrix[P, TT + 1] Pop_weight;
}

parameters {
  // province effects INCLUDING TT+1
  matrix[P - 1, TT + 1] alpha0;

  // election shocks INCLUDING TT+1
  vector[TT + 1] delta;

  // fixed effects
  vector[K] beta;
  vector[L] gamma;

  // scales
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_delta;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;

  // pseudo outcome for TT+1
  vector[P] y_pred;
}

transformed parameters {
  matrix[P, TT + 1] alpha;
  vector[N] eta;

  // weighted-mean-zero constraint
  alpha[1:(P - 1), ] = alpha0;
  for (t in 1:(TT + 1)) {
    alpha[P, t] =
      -sum(Pop_weight[1:(P - 1), t] .* alpha0[, t]) / Pop_weight[P, t];
  }

  // linear predictor (observed elections only)
  for (t in 1:TT) {
    int start = 1 + (t - 1) * P;
    int end_  = t * P;

    eta[start:end_] =
      alpha[, t]
      + (X * beta)[start:end_]
      + (Z * gamma)[start:end_]
      + rep_vector(delta[t], P);
  }
}

model {
  // --- state evolution INCLUDING TT+1 ---
  for (t in 1:(TT + 1)) {
    if (t == 1)
      alpha[, 1] ~ normal(0, sigma_alpha);
    else
      alpha[, t] ~ normal(alpha[, t - 1], sigma_alpha);
  }

  // --- likelihood (observed only) ---
  Y_count ~ binomial_logit(N_count, eta);

  // --- pseudo-likelihood for TT+1 ---
  y_pred ~ binomial_logit(
    rep_array(1, P),  // dummy denominator
    alpha[, TT + 1] + X_pred * beta + Z_pred * gamma + delta[TT + 1]
  );

  // --- priors ---
  delta ~ normal(0, sigma_delta);
  beta  ~ normal(0, sigma_beta);
  gamma ~ normal(0, sigma_gamma);

  sigma_alpha ~ normal(0, 1);
  sigma_delta ~ normal(0, 1);
  sigma_beta  ~ normal(0, 1);
  sigma_gamma ~ normal(0, 1);
}
