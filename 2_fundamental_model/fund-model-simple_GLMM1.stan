data {
  int<lower=2> P;                 // # regions
  int<lower=1> TT;                // # elections used for fitting
  int<lower=1> N;                 // N = P * TT
  int<lower=1> K;                 // # national covariates
  int<lower=1> L;                 // # regional covariates

  matrix[N, K] X;
  matrix[N, L] Z;

  array[N] int<lower=0> Y_count;  // Dem votes
  array[N] int<lower=0> N_count;  // Dem+Con votes

  matrix[P, K] X_pred;
  matrix[P, L] Z_pred;

  matrix[P, TT + 1] Pop_weight;   // 각 선거별 지역 인구 비중
}

parameters {
  // non-centered RW for alpha (P-1 regions)
  vector[P-1] alpha0;                 // initial state at t=1
  matrix[P-1, TT-1] z_alpha;          // innovations for t=2..TT

  // non-centered RW for delta
  real delta0;                        // initial delta at t=1
  vector[TT-1] z_delta;               // innovations for t=2..TT

  // fixed effects
  vector[K] beta;
  vector[L] gamma;

  // scales
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_delta;
}


transformed parameters {
  matrix[P-1, TT] alpha_raw;
  vector[TT] delta;

  matrix[P, TT + 1] alpha;
  vector[N] eta;

  // reconstruct alpha_raw (RW)
  alpha_raw[, 1] = alpha0;
  for (t in 2:TT)
    alpha_raw[, t] = alpha_raw[, t-1] + sigma_alpha * z_alpha[, t-1];

  // reconstruct delta (RW)
  delta[1] = delta0;
  for (t in 2:TT)
    delta[t] = delta[t-1] + sigma_delta * z_delta[t-1];

  // alpha init
  for (t in 1:(TT + 1))
    alpha[, t] = rep_vector(0.0, P);

  // weighted-sum constraint (t=1..TT)
  for (t in 1:TT) {
    alpha[1:(P - 1), t] = alpha_raw[, t];
    alpha[P, t] =
      -sum(Pop_weight[1:(P - 1), t] .* alpha_raw[, t]) / Pop_weight[P, t];
  }

  // eta (fit period)
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
  // innovations
  to_vector(z_alpha) ~ normal(0, 1);
  z_delta ~ normal(0, 1);

  // initial states
  alpha0 ~ normal(0, 1);
  delta0 ~ normal(0, 1);

  // scales (조금 더 강하게 주는 게 대체로 안정적)
  sigma_alpha ~ normal(0, 0.3);
  sigma_delta ~ normal(0, 0.3);

  // fixed effects (이미 sigma 제거 버전 기준)
  beta  ~ normal(0, 0.5);
  gamma ~ normal(0, 0.5);

  // likelihood
  Y_count ~ binomial_logit(N_count, eta);
}

generated quantities {
  vector[P] eta_pred;
  vector[P] pi_pred;
  vector[N] log_lik;

  vector[P] alpha_TTp1;
  vector[P-1] alpha_TTp1_free;
  real delta_TTp1;

  array[P] int y_pred_count; // ✅ 여기서 선언 (루프 밖)

  // alpha_{TT+1}
  for (p in 1:(P - 1))
    alpha_TTp1_free[p] = normal_rng(alpha_raw[p, TT], sigma_alpha);

  alpha_TTp1[1:(P - 1)] = alpha_TTp1_free;
  alpha_TTp1[P] =
    -sum(Pop_weight[1:(P - 1), TT + 1] .* alpha_TTp1[1:(P - 1)])
    / Pop_weight[P, TT + 1];

  // delta_{TT+1}
  delta_TTp1 = normal_rng(delta[TT], sigma_delta);

  // prediction
  eta_pred = alpha_TTp1 + X_pred * beta + Z_pred * gamma + rep_vector(delta_TTp1, P);
  pi_pred  = inv_logit(eta_pred);

  // log-lik
  for (n in 1:N)
    log_lik[n] = binomial_logit_lpmf(Y_count[n] | N_count[n], eta[n]);

  // ✅ 예측 득표수 생성 (현재는 마지막 선거의 분모를 임시 사용)
  for (p in 1:P)
    y_pred_count[p] = binomial_rng(N_count[(TT - 1) * P + p], pi_pred[p]);
}

