// ------------------------------------------------------------------------------
// [Modified Baseline Model]
// 1. Delta (National Trend) Removed:
//    The national swing is now solely explained by predictors (X * beta).
//    Alpha carries the baseline intercept (anchored at 0.5).
//
// 2. Hard Constraint Maintained:
//    Alpha[P] is deterministically calculated so that sum(Weight * Alpha) = 0.5.
// ------------------------------------------------------------------------------
data {
  int<lower=1> P;                 // number of provinces
  int<lower=1> TT;                // number of observed elections
  int<lower=1> N;                 // number of observations = P*TT
  int<lower=1> K;                 // number of national-level coeffs
  int<lower=1> L;                 // number of province-level coeffs

  matrix[N, K] X;                 // data matrix of national effects
  matrix[N, L] Z;                 // data matrix of province effects

  // covariates for the (TT+1) election
  matrix[P, K] X_pred;
  matrix[P, L] Z_pred;

  vector<lower=0, upper=1>[N] Y;
  vector[P] Alpha_init;
  matrix[P, TT + 1] Pop_weight;   
}

parameters {
  // 비중심화(Non-Centered): alpha0_raw (Standard Normal)
  matrix[P - 1, TT] alpha0_raw;    
   
  vector[K] beta;
  vector[L] gamma;
   
  // [삭제됨] Delta 관련 파라미터 제거
  // vector[TT] delta;                
   
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;
  
  // [삭제됨] sigma_delta 제거
  // real<lower=0> sigma_delta;
   
  array[TT] real<lower=0> sigma_epsilon;
   
  real<lower=0> sigma;            // global scale
}

transformed parameters {
  matrix[P - 1, TT] alpha0;
  matrix[P, TT] alpha;

  // 1. Alpha0 복원 (Random Walk)
  // t=1
  alpha0[, 1] = Alpha_init[1:(P - 1)] + sigma_alpha * alpha0_raw[, 1];
  // t=2~TT
  for (t in 2:TT) {
    alpha0[, t] = alpha0[, t - 1] + sigma_alpha * alpha0_raw[, t];
  }

  // 2. 전체 Alpha 행렬 구성 (1 ~ P-1)
  alpha[1:(P - 1), ] = alpha0;

  // 3. Hard Constraint (P번째 지역 계산)
  // "모든 지역 Alpha의 가중평균은 0.5여야 한다"
  // (Delta가 없으므로 Alpha가 50% 기준선을 잡아줌)
  for (t in 1:TT) {
    alpha[P, t] = (0.5 - sum(Pop_weight[1:(P - 1), t] .* alpha0[, t])) / Pop_weight[P, t];
  }
}

model {
  // --- Priors ---
  to_vector(alpha0_raw) ~ std_normal();

  // [삭제됨] Delta Prior 제거
  // delta ~ normal(0, sigma_delta);
   
  beta  ~ normal(0, sigma_beta);        
  gamma ~ normal(0, sigma_gamma);       

  // Scale parameters
  sigma_alpha   ~ normal(0, sigma);
  sigma_beta    ~ normal(0, sigma);
  sigma_gamma   ~ normal(0, sigma);
  sigma_epsilon ~ normal(0, sigma);
  
  // [삭제됨] sigma_delta 제거
  // sigma_delta   ~ normal(0, sigma);
  
  sigma ~ normal(0, 1); // Global scale prior

  // --- Likelihood ---
  for (t in 1:TT) {
    int start_idx = 1 + (t - 1) * P;
    int end_idx = t * P;
    
    Y[start_idx:end_idx] ~ normal(
      // [수정] delta[t] 제거됨
      alpha[, t] + 
      (X * beta)[start_idx:end_idx] + 
      (Z * gamma)[start_idx:end_idx], 
      sigma_epsilon[t]
    );
  }
}

generated quantities {
  // --- Prediction for TT+1 ---
  vector[P] alpha_TTp1;
  // [삭제됨] delta_TTp1 제거
  
  vector[P] mu_pred;
  vector[P] y_pred;

  // [삭제됨] Delta 예측 제거
  // delta_TTp1 = normal_rng(0, sigma_delta);

  {
    vector[P - 1] alpha0_TTp1_free;
    for (p in 1:(P - 1)) {
      // Random Walk Forward
      alpha0_TTp1_free[p] = normal_rng(alpha0[p, TT], sigma_alpha);
      alpha_TTp1[p] = alpha0_TTp1_free[p];
    }

    // Future Hard Constraint
    alpha_TTp1[P] = (0.5 - dot_product(Pop_weight[1:(P - 1), TT + 1], alpha0_TTp1_free)) 
                    / Pop_weight[P, TT + 1];
  }

  // [수정] 예측식에서 delta 항 제거
  mu_pred = alpha_TTp1 + X_pred * beta + Z_pred * gamma;

  for (p in 1:P) {
    y_pred[p] = normal_rng(mu_pred[p], sigma_epsilon[TT]);
  }
}
