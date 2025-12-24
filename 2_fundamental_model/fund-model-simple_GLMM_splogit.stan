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
}

parameters {
  vector[P] alpha1;
  matrix[P, TT-1] alpha_innov_raw;
  vector[P] alpha_pred_innov_raw;
  vector[TT] delta;
  real delta_pred;
  vector[K] beta;
  vector[L] gamma;
  
  // --- 지역별 r_p 추정을 위한 계층적 파라미터 ---
  vector<lower=0.1>[P] r_p;      // 각 지역의 비대칭 매개변수
  real<lower=0> mu_r;            // 전국 평균 r의 중심
  real<lower=0> sigma_r;         // 지역 간 r의 변동성
  
  real<lower=0> sigma_alpha0;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_delta;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;
}

transformed parameters {
  matrix[P, TT] alpha;
  vector[N] eta;
  vector[P] eta_pred;
  vector[P] pi_pred;
  
  alpha[,1] = alpha1;
  for (t in 2:TT)
    alpha[,t] = alpha[,t-1] + sigma_alpha * alpha_innov_raw[, t-1];
  
  for (t in 1:TT) {
    int start = 1 + (t - 1) * P;
    int end_  = t * P;
    eta[start:end_] = alpha[,t] + (X * beta)[start:end_] + (Z * gamma)[start:end_] + rep_vector(delta[t], P);
  }
  
  vector[P] alpha_TTp1 = alpha[,TT] + sigma_alpha * alpha_pred_innov_raw;
  eta_pred = alpha_TTp1 + X_pred * beta + Z_pred * gamma + rep_vector(delta_pred, P);
  
  for (p in 1:P) {
    pi_pred[p] = inv_splogit(eta_pred[p], r_p[p]); // 지역별 r_p 적용
  }
}

model {
  alpha1 ~ normal(0, sigma_alpha0);
  delta  ~ normal(0, sigma_delta);
  delta_pred ~ normal(0, sigma_delta);
  to_vector(alpha_innov_raw) ~ normal(0, 1);
  alpha_pred_innov_raw ~ normal(0, 1);
  beta  ~ normal(0, sigma_beta);
  gamma ~ normal(0, sigma_gamma);
  
  // --- r_p에 대한 계층적 사전분포 (Lognormal 구조) ---
  mu_r ~ gamma(2, 2);           // r의 중심값은 1 근처로 유도
  sigma_r ~ normal(0, 0.5);      // 지역 간 차이가 과도하게 벌어지지 않도록 제어
  r_p ~ lognormal(log(mu_r), sigma_r); 
  
  sigma_alpha0 ~ normal(0, 1);
  sigma_alpha  ~ normal(0, 1);
  sigma_delta  ~ normal(0, 1);
  sigma_beta   ~ normal(0, 1);
  sigma_gamma  ~ normal(0, 1);
  
  // --- Likelihood: n번째 데이터의 지역 p를 매핑하여 r_p[p] 사용 ---
  for (t in 1:TT) {
    for (p in 1:P) {
      int n = p + (t - 1) * P;
      Y_count[n] ~ binomial_splogit(N_count[n], eta[n], r_p[p]);
    }
  }
}

generated quantities {
  vector[N] log_lik;
  for (t in 1:TT) {
    for (p in 1:P) {
      int n = p + (t - 1) * P;
      log_lik[n] = binomial_splogit_lpmf(Y_count[n] | N_count[n], eta[n], r_p[p]);
    }
  }
}