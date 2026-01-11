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
  matrix[P, TT + 1] Pop_weight;
}

transformed data {
  vector[N] Y_logit;
  real eps;
  eps = 1e-4;  // (안정화 목적) 1e-6보다 크게 시작 권장

  for (n in 1:N) {
    real y_adj = fmin(1 - eps, fmax(eps, Y[n]));
    Y_logit[n] = logit(y_adj);
  }
}

parameters {
  // non-centered RW increments (RW 유지)
  matrix[P - 1, TT + 1] z_alpha;

  vector[K] beta;
  vector[L] gamma;

  // 전국 공통 충격(절편 역할). 유지하되 수축 prior로 안정화
  vector[TT] delta;

  // scales
  real<lower=0> sigma_alpha;

  // ✅ 오차 단순화: 시점별 sigma_epsilon[t] 제거 → 단일 sigma_epsilon
  real<lower=0> sigma_epsilon;
}

transformed parameters {
  matrix[P - 1, TT + 1] alpha0;
  matrix[P, TT + 1] alpha;

  // non-centered random walk
  alpha0[, 1] = Alpha_init[1:(P - 1)] + sigma_alpha * z_alpha[, 1];
  for (t in 2:(TT + 1)) {
    alpha0[, t] = alpha0[, t - 1] + sigma_alpha * z_alpha[, t];
  }

  // weighted-mean constraint: sum_p w_{p,t} * alpha_{p,t} = 0  (logit 중심)
  alpha[1:(P - 1), ] = alpha0;
  for (t in 1:(TT + 1)) {
    alpha[P, t] =
      -dot_product(Pop_weight[1:(P - 1), t], alpha0[, t]) / Pop_weight[P, t];
  }
}

model {
  // -------------------------
  // 1) Stronger / regularizing priors
  // -------------------------
  to_vector(z_alpha) ~ normal(0, 1);

  // (회귀계수) 공변량 스케일이 크면 더 줄이세요. 표준화되어 있으면 0.3 권장.
  beta  ~ normal(0, 0.3);
  gamma ~ normal(0, 0.3);

  // (전국충격) TT가 짧을수록 delta는 강수축이 안정화에 도움
  // "절편" 역할이므로 너무 자유로우면 alpha/epsilon과 경쟁합니다.
  delta ~ normal(0, 0.15);

  // (상태 변화 크기) half-normal: 데이터 변화가 크더라도 우선은 중간 스케일부터 시작
  sigma_alpha ~ normal(0, 0.6);

  // ✅ (관측오차) 단일 sigma_epsilon + 강수축
  // 로짓-정규에서 sigma_epsilon이 커지면 변동을 오차가 먹어 ridge가 생깁니다.
  sigma_epsilon ~ normal(0, 0.15);

  // -------------------------
  // 2) Likelihood (fit period only)
  // -------------------------
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

  // 예측용 delta (선택 1: 0으로 고정 / 선택 2: prior draw)
  real delta_TTp1 = 0;
  // real delta_TTp1 = normal_rng(0, 0.15);  // delta prior 스케일과 일치시키고 싶으면 사용

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
