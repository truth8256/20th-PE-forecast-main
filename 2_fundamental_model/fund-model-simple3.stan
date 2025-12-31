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
  // [수정] 모든 주요 파라미터를 "raw" (표준정규분포) 형태로 선언
  matrix[P - 1, TT] alpha0_raw;   
  vector[K] beta_raw;             // beta 비중심화
  vector[L] gamma_raw;            // gamma 비중심화
  vector[TT] delta_raw;           // delta 비중심화
  
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;
  real<lower=0> sigma_delta;
  
  array[TT] real<lower=0> sigma_epsilon;
  
  real<lower=0> sigma;            // global scale
}

transformed parameters {
  // [변환] raw 값에 sigma를 곱해서 실제 파라미터로 복원
  matrix[P - 1, TT] alpha0;
  matrix[P, TT] alpha;
  
  vector[K] beta;
  vector[L] gamma;
  vector[TT] delta;

  // 1. 파라미터 복원 (Matt Trick)
  // beta ~ normal(0, sigma_beta)와 수학적으로 동일하지만 샘플링은 훨씬 빠름
  beta = beta_raw * sigma_beta;
  gamma = gamma_raw * sigma_gamma;
  delta = delta_raw * sigma_delta;

  // 2. Alpha0 복원 (Random Walk)
  alpha0[, 1] = Alpha_init[1:(P - 1)] + sigma_alpha * alpha0_raw[, 1];
  for (t in 2:TT) {
    alpha0[, t] = alpha0[, t - 1] + sigma_alpha * alpha0_raw[, t];
  }

  // 3. 전체 Alpha 행렬 구성 및 제약조건 적용
  alpha[1:(P - 1), ] = alpha0;
  for (t in 1:TT) {
    alpha[P, t] = (0.5 - sum(Pop_weight[1:(P - 1), t] .* alpha0[, t])) / Pop_weight[P, t];
  }
}

model {
  // --- Priors (All Standard Normal) ---
  // Stan이 가장 좋아하는 분포(N(0,1))에서 샘플링하게 함
  to_vector(alpha0_raw) ~ std_normal();
  beta_raw  ~ std_normal();
  gamma_raw ~ std_normal();
  delta_raw ~ std_normal();

  // --- Likelihood ---
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

  // --- Hyper-priors ---
  // Scale parameters
  sigma_alpha   ~ normal(0, sigma);
  sigma_delta   ~ normal(0, sigma);
  sigma_beta    ~ normal(0, sigma);
  sigma_gamma   ~ normal(0, sigma);
  sigma_epsilon ~ normal(0, sigma);
  
  // Global scale prior (조금 더 구체적으로 잡아주는 것이 좋음)
  sigma ~ normal(0, 1); 
}

generated quantities {
  vector[P] alpha_TTp1;
  real delta_TTp1;
  vector[P] mu_pred;
  vector[P] y_pred;

  // 예측 시에는 sigma_delta를 사용하여 직접 생성
  delta_TTp1 = normal_rng(0, sigma_delta);

  {
    vector[P - 1] alpha0_TTp1_free;
    for (p in 1:(P - 1)) {
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
