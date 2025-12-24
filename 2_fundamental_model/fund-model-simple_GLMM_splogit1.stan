functions {
  real inv_splogit(real x, real r) {
    if (r <= 1.0) {
      return pow(inv_logit(x / r), r);
    } else {
      return 1.0 - pow(inv_logit(-r * x), 1.0 / r);
    }
  }

  real binomial_splogit_lpmf(int y, int n, real eta, real r) {
    real p = inv_splogit(eta, r);
    p = fmax(1e-10, fmin(1.0 - 1e-10, p));
    return binomial_lpmf(y | n, p);
  }
}

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
  // initial alpha (t = 1)
  vector[P - 1] alpha1_raw;

  // innovations for t = 2..TT
  matrix[P - 1, TT - 1] alpha_innov_raw;

  // national swing (observed elections only)
  vector[TT] delta;

  // fixed effects
  vector[K] beta;
  vector[L] gamma;

  // splogit shape parameters
  vector<lower=0.1>[P] r_p;
  real<lower=0> mu_r;
  real<lower=0> sigma_r;

  // scales
  real<lower=0> sigma_alpha0;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_delta;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;
}

transformed parameters {
  matrix[P, TT] alpha;
  vector[N] eta;

  // t = 1
  alpha[1:(P - 1), 1] = alpha1_raw;
  alpha[P, 1] =
    -sum(Pop_weight[1:(P - 1), 1] .* alpha1_raw) / Pop_weight[P, 1];

  // t = 2..TT
  for (t in 2:TT) {
    alpha[1:(P - 1), t] =
      alpha[1:(P - 1), t - 1] + sigma_alpha * alpha_innov_raw[, t - 1];

    alpha[P, t] =
      -sum(Pop_weight[1:(P - 1), t] .* alpha[1:(P - 1), t])
      / Pop_weight[P, t];
  }

  // linear predictor
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
  // priors
  alpha1_raw ~ normal(0, sigma_alpha0);
  to_vector(alpha_innov_raw) ~ normal(0, 1);

  delta ~ normal(0, sigma_delta);

  beta  ~ normal(0, sigma_beta);
  gamma ~ normal(0, sigma_gamma);

  mu_r ~ gamma(2, 2);
  sigma_r ~ normal(0, 0.5);
  r_p ~ lognormal(log(mu_r), sigma_r);

  sigma_alpha0 ~ normal(0, 1);
  sigma_alpha  ~ normal(0, 1);
  sigma_delta  ~ normal(0, 1);
  sigma_beta   ~ normal(0, 1);
  sigma_gamma  ~ normal(0, 1);

  // likelihood (observed elections only)
  for (n in 1:N)
    Y_count[n] ~ binomial_splogit(N_count[n], eta[n], r_p[(n - 1) % P + 1]);
}

generated quantities {
  vector[P] alpha_TTp1;
  real delta_TTp1;

  vector[P] eta_pred;
  vector[P] pi_pred;

  vector[N] log_lik;

  // 1) alpha one-step-ahead prediction
  {
    vector[P - 1] alpha_free;
    for (p in 1:(P - 1)) {
      alpha_free[p] = normal_rng(alpha[p, TT], sigma_alpha);
      alpha_TTp1[p] = alpha_free[p];
    }
    alpha_TTp1[P] =
      -dot_product(Pop_weight[1:(P - 1), TT + 1], alpha_free)
      / Pop_weight[P, TT + 1];
  }

  // 2) national swing prediction
  delta_TTp1 = normal_rng(0, sigma_delta);

  // 3) linear predictor and probability
  eta_pred =
      alpha_TTp1
    + X_pred * beta
    + Z_pred * gamma
    + rep_vector(delta_TTp1, P);

  for (p in 1:P)
    pi_pred[p] = inv_splogit(eta_pred[p], r_p[p]);

  // 4) log-likelihood (observed elections only)
  for (n in 1:N)
    log_lik[n] =
      binomial_splogit_lpmf(
        Y_count[n] | N_count[n],
        eta[n],
        r_p[(n - 1) % P + 1]
      );
}
