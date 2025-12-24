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
  
  matrix[P, TT + 1] Pop_weight;   // 각 선거별 지역 인구 비중 (추가됨)
}

parameters {
  // P-1 지역에 대해서만 자유도(degrees of freedom) 부여
  matrix[P - 1, TT] alpha_raw;    
  
  // election shock (national swing)
  vector[TT] delta;
  real delta_pred;
  
  // fixed effects
  vector[K] beta;
  vector[L] gamma;
  
  // scales
  real<lower=0> sigma_alpha;      // alpha의 시계열 변화량 (필요 시)
  real<lower=0> sigma_delta;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;
}

transformed parameters {
  matrix[P, TT + 1] alpha;        // 예측(TT+1) 포함을 위해 행렬 확장
  vector[N] eta;
  
  // 1. 가중 평균 0 제약을 활용한 alpha 구성
  for (t in 1:TT) {
    alpha[1:(P - 1), t] = alpha_raw[, t];
    // 마지막 P번째 지역은 가중 합이 0이 되도록 결정 (로짓 스케일)
    alpha[P, t] = -sum(Pop_weight[1:(P - 1), t] .* alpha_raw[, t]) / Pop_weight[P, t];
  }

  // 2. Linear Predictor (eta) 계산
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
  // alpha_raw에 대한 시계열/계층적 사전분포 설정
  to_vector(alpha_raw) ~ normal(0, 5); // 필요시 sigma_alpha를 활용한 Random Walk로 변경 가능
  
  delta  ~ normal(0, sigma_delta);
  delta_pred ~ normal(0, sigma_delta);
  
  beta  ~ normal(0, sigma_beta);
  gamma ~ normal(0, sigma_gamma);
  
  sigma_delta ~ normal(0, 1);
  sigma_beta  ~ normal(0, 1);
  sigma_gamma ~ normal(0, 1);
  
  // likelihood
  Y_count ~ binomial_logit(N_count, eta);
}

generated quantities {
  vector[P] eta_pred;
  vector[P] pi_pred;
  vector[N] log_lik;
  vector[P] alpha_TTp1;

  // (1) 예측 시점(TT+1)의 alpha를 "새 변수"로 구성
  //     여기서는 단순히 마지막 시점(TT)의 alpha_raw를 그대로 가져오는 규칙
  alpha_TTp1[1:(P - 1)] = alpha_raw[, TT];
  alpha_TTp1[P] = -sum(Pop_weight[1:(P - 1), TT + 1] .* alpha_TTp1[1:(P - 1)])
                  / Pop_weight[P, TT + 1];

  // (2) 예측 선형예측자 및 확률
  eta_pred = alpha_TTp1 + X_pred * beta + Z_pred * gamma + rep_vector(delta_pred, P);
  pi_pred  = inv_logit(eta_pred);

  // (3) pointwise log-likelihood (fit 구간만)
  for (n in 1:N)
    log_lik[n] = binomial_logit_lpmf(Y_count[n] | N_count[n], eta[n]);
}

