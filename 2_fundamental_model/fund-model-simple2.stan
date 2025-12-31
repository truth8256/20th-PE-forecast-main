// ------------------------------------------------------------------------------
// [Key Modifications vs. Original Implementation]
// 1. Alpha Evolution Constraint: 
//    The Random Walk prior is applied ONLY to the free parameters (alpha0: provinces 1 to P-1).
//    The P-th province is deterministically calculated to satisfy the weighted-mean-to-0.5 constraint,
//    resolving the logical conflict in the original full-vector specification.
//
// 2. Prediction Strategy (TT+1):
//    Future outcomes are generated in the 'generated quantities' block (Posterior Predictive Checks),
//    instead of treating unobserved Y and parameters as missing data to be imputed.
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
  // [수정 1] 비중심화(Non-Centered): alpha0 대신 표준정규분포를 따르는 raw 변수 선언
  matrix[P - 1, TT] alpha0_raw;   
  
  vector[K] beta;
  vector[L] gamma;
  
  vector[TT] delta;               
  
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;
  real<lower=0> sigma_delta;
  
  array[TT] real<lower=0> sigma_epsilon;
  
  real<lower=0> sigma;            // global scale
}

transformed parameters {
  // [수정 1의 연장] 실제 alpha0와 alpha를 여기서 재조립
  matrix[P - 1, TT] alpha0;
  matrix[P, TT] alpha;

  // 1. Alpha0 복원 (Random Walk 구현: 누적합 방식이 효율적)
  // t=1
  alpha0[, 1] = Alpha_init[1:(P - 1)] + sigma_alpha * alpha0_raw[, 1];
  // t=2~TT
  for (t in 2:TT) {
    alpha0[, t] = alpha0[, t - 1] + sigma_alpha * alpha0_raw[, t];
  }

  // 2. 전체 Alpha 행렬 구성
  alpha[1:(P - 1), ] = alpha0;

  // 3. P번째 지역 제약조건 적용 (논리 유지)
  for (t in 1:TT) {
    alpha[P, t] = (0.5 - sum(Pop_weight[1:(P - 1), t] .* alpha0[, t])) / Pop_weight[P, t];
  }
}

model {
  // --- State evolution Priors ---
  // [수정 1] Raw 변수는 무조건 표준정규분포를 따름 (샘플링 효율 극대화)
  to_vector(alpha0_raw) ~ std_normal();

  // --- Likelihood (Vectorized for speed) ---
  for (t in 1:TT) {
    int start_idx = 1 + (t - 1) * P;
    int end_idx = t * P;
    
    Y[start_idx:end_idx] ~ normal(
      alpha[, t] + 
      (X * beta)[start_idx:end_idx] + 
      (Z * gamma)[start_idx:end_idx] + 
      delta[t], 
      sigma_epsilon[t]
    );
  }

  // --- Priors ---
  delta ~ normal(0, sigma_delta);
  
  // [수정 2] Cauchy -> Normal 변경 (Treedepth 문제 해결의 핵심)
  beta  ~ normal(0, sigma_beta);       
  gamma ~ normal(0, sigma_gamma);      

  // Scale parameters (Weakly informative)
  sigma_alpha   ~ normal(0, sigma);
  sigma_delta   ~ normal(0, sigma);
  sigma_beta    ~ normal(0, sigma);
  sigma_gamma   ~ normal(0, sigma);
  sigma_epsilon ~ normal(0, sigma);
}

generated quantities {
  // --- Prediction for TT+1 ---
  vector[P] alpha_TTp1;
  real delta_TTp1;
  vector[P] mu_pred;
  vector[P] y_pred;

  delta_TTp1 = normal_rng(0, sigma_delta);

  {
    vector[P - 1] alpha0_TTp1_free;
    for (p in 1:(P - 1)) {
      // 마지막 시점(TT)에서 한 발자국 더 Random Walk
      alpha0_TTp1_free[p] = normal_rng(alpha0[p, TT], sigma_alpha);
      alpha_TTp1[p] = alpha0_TTp1_free[p];
    }

    alpha_TTp1[P] = (0.5 - dot_product(Pop_weight[1:(P - 1), TT + 1], alpha0_TTp1_free)) 
                    / Pop_weight[P, TT + 1];
  }

  mu_pred = alpha_TTp1 + X_pred * beta + Z_pred * gamma + delta_TTp1;

  for (p in 1:P) {
    y_pred[p] = normal_rng(mu_pred[p], sigma_epsilon[TT]);
  }
}
