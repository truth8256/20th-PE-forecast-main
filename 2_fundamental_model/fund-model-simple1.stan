// Baseline fundamentals model (Kang-style) rewritten so that the (TT+1) election
// is generated from the posterior (posterior predictive), not treated as parameters.
//
// Key changes vs. your original:
// 1) Remove TT+1 from parameters: alpha0 becomes [P-1, TT], delta becomes [TT].
// 2) Remove y_pred from parameters (no pseudo-likelihood for unobserved outcome).
// 3) Generate alpha_{.,TT+1}, delta_{TT+1}, and y_tilde_{.,TT+1} in generated quantities
//    using RNG, conditional on posterior draws of parameters and states.
// 4) Keep the weighted-mean-to-0.5 constraint for each t, including TT+1 (using Pop_weight[,TT+1]).

data {
  int<lower=2> P;                 // number of provinces
  int<lower=1> TT;                // number of observed elections (prediction is TT+1)
  int<lower=1> N;                 // number of observations = P*TT
  int<lower=1> K;                 // number of national-level coeffs
  int<lower=1> L;                 // number of province-level coeffs

  matrix[N, K] X;                 // stacked by time: (t-1)*P + p
  matrix[N, L] Z;                 // stacked by time: (t-1)*P + p

  // covariates for the (TT+1) election, one row per province
  matrix[P, K] X_pred;
  matrix[P, L] Z_pred;

  vector<lower=0, upper=1>[N] Y;  // observed vote share (14..TT)

  vector[P] Alpha_init;           // prior mean for alpha at t=1 (length P)

  matrix[P, TT + 1] Pop_weight;   // turnout weights, includes prediction column TT+1
}

parameters {
  // Province effects for observed elections only (t = 1..TT)
  matrix[P - 1, TT] alpha0;       // last province computed via constraint

  vector[K] beta;
  vector[L] gamma;

  // election-level intercepts for observed elections only
  vector[TT] delta;

  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;
  real<lower=0> sigma_delta;

  array[TT] real<lower=0> sigma_epsilon;

  real<lower=0> sigma;            // global scale for half-normal priors
}

transformed parameters {
  // weighted mean to 0.5 constraint, for observed elections only (1..TT)
  matrix[P, TT] alpha;

  alpha[1:(P - 1), ] = alpha0;
  for (t in 1:TT) {
    alpha[P, t] =
      (0.5 - sum(Pop_weight[1:(P - 1), t] .* alpha0[, t])) / Pop_weight[P, t];
  }
}

model {
  // --- State evolution for observed elections only ---
  for (t in 1:TT) {
    if (t == 1) {
      alpha[, 1] ~ normal(Alpha_init, sigma_alpha);
    } else {
      alpha[, t] ~ normal(alpha[, t - 1], sigma_alpha);
    }
  }

  // --- Likelihood: only observed elections (1..TT) ---
  for (t in 1:TT) {
    int start = 1 + (t - 1) * P;
    int end_  = t * P;

    Y[start:end_] ~ normal(
      alpha[, t]
      + (X * beta)[start:end_]
      + (Z * gamma)[start:end_]
      + delta[t],
      sigma_epsilon[t]
    );
  }

  // --- Priors ---
  delta ~ normal(0, sigma_delta);

  // beta ~ cauchy(0, sigma_beta);
  // gamma ~ cauchy(0, sigma_gamma);
  beta  ~ normal(0, sigma_beta);
  gamma ~ normal(0, sigma_gamma);


  // sigma_alpha ~ normal(0, sigma);
  // sigma_delta ~ normal(0, sigma);
  // sigma_beta  ~ normal(0, sigma);
  // sigma_gamma ~ normal(0, sigma);
  // sigma_epsilon ~ normal(0, sigma);
  
  sigma_alpha ~ normal(0, 0.05);   // 예시: α는 선거 간 변동이 작다는 가정
  sigma_delta ~ normal(0, 0.10);   // 선거-level intercept
  sigma_beta  ~ normal(0, 0.20);   // X 계수 스케일
  sigma_gamma ~ normal(0, 0.20);   // Z 계수 스케일
  for (t in 1:TT) sigma_epsilon[t] ~ normal(0, 0.03); // 관측오차

}

generated quantities {
  // One-step-ahead prediction for the (TT+1) election

  vector[P] alpha_TTp1;          // alpha for prediction election (TT+1)
  real delta_TTp1;               // delta for prediction election (TT+1)

  vector[P] mu_pred;             // conditional mean for y at (TT+1)
  vector[P] y_tilde;             // posterior predictive draw (vote share)

  // 1) Draw next election-level intercept from its prior (consistent with delta ~ N(0, sigma_delta))
  delta_TTp1 = normal_rng(0, sigma_delta);

  // 2) Draw next alpha using the RW evolution for provinces 1..(P-1)
  {
    vector[P - 1] alpha0_TTp1_free;
    for (p in 1:(P - 1)) {
      alpha0_TTp1_free[p] = normal_rng(alpha[p, TT], sigma_alpha);
      alpha_TTp1[p] = alpha0_TTp1_free[p];
    }

    // enforce weighted-mean-to-0.5 constraint at TT+1
    alpha_TTp1[P] =
      (0.5 - dot_product(Pop_weight[1:(P - 1), TT + 1], alpha0_TTp1_free))
      / Pop_weight[P, TT + 1];
  }

  // 3) Compute prediction mean and draw posterior predictive replicate
  mu_pred = alpha_TTp1 + X_pred * beta + Z_pred * gamma + rep_vector(delta_TTp1, P);

  // use sigma_epsilon[TT] as in your original code (carry-forward assumption)
  for (p in 1:P) {
    y_tilde[p] = normal_rng(mu_pred[p], sigma_epsilon[TT]);

    // optional safety clamp to [0,1] because Y is a vote share.
    // If you prefer not to truncate, comment this out.
    if (y_tilde[p] < 0) y_tilde[p] = 0;
    if (y_tilde[p] > 1) y_tilde[p] = 1;
  }
}
