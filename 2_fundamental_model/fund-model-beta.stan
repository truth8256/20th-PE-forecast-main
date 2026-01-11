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

  // [주의] Y는 0과 1을 포함하지 않는 (0, 1) 구간이어야 함 (베타분포 특성)
  vector<lower=0, upper=1>[N] Y;
  vector[P] Alpha_init;
  matrix[P, TT + 1] Pop_weight;   
}

parameters {
  // [수정 1] 비중심화(Non-Centered)
  matrix[P - 1, TT] alpha0_raw;   
   
  vector[K] beta;
  vector[L] gamma;
   
  vector[TT] delta;                
   
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;
  real<lower=0> sigma_delta;
   
  // [변경] Normal의 sigma_epsilon 대신 Beta 분포의 정밀도(precision) 파라미터 phi 사용
  real<lower=0> phi;            
   
  real<lower=0> sigma;            // global scale
}

transformed parameters {
  matrix[P - 1, TT] alpha0;
  matrix[P, TT] alpha;

  // 1. Alpha0 복원
  alpha0[, 1] = Alpha_init[1:(P - 1)] + sigma_alpha * alpha0_raw[, 1];
  for (t in 2:TT) {
    alpha0[, t] = alpha0[, t - 1] + sigma_alpha * alpha0_raw[, t];
  }

  // 2. 전체 Alpha 행렬 구성
  alpha[1:(P - 1), ] = alpha0;

  // 3. [중요 변경] P번째 지역 제약조건 수정
  // 로짓 스케일에서 0은 확률 0.5를 의미합니다.
  // 따라서 가중평균이 '0.5'가 아니라 로짓 스케일의 중심인 '0'이 되도록 맞춰줍니다.
  for (t in 1:TT) {
    alpha[P, t] = (0.0 - sum(Pop_weight[1:(P - 1), t] .* alpha0[, t])) / Pop_weight[P, t];
  }
}

model {
  // --- State evolution Priors ---
  to_vector(alpha0_raw) ~ std_normal();

  // --- Likelihood (Beta Regression with Logit Link) ---
  // mu는 [0, 1] 사이의 확률값으로 변환됨 (inv_logit 적용)
  vector[N] mu;
  vector[N] linear_predictor;
  
  for (t in 1:TT) {
    int start_idx = 1 + (t - 1) * P;
    int end_idx = t * P;
    
    // 선형 예측식 (Linear Predictor) - 범위: (-inf, +inf)
    linear_predictor[start_idx:end_idx] = 
      alpha[, t] + 
      (X * beta)[start_idx:end_idx] + 
      (Z * gamma)[start_idx:end_idx] + 
      delta[t];
  }
  
  // Link Function: Logit^-1
  mu = inv_logit(linear_predictor);

  // Beta Likelihood
  // Stan의 beta_proportion(mu, phi) 함수 사용: mu=평균, phi=정밀도
  Y ~ beta_proportion(mu, phi);

  // --- Priors ---
  delta ~ normal(0, sigma_delta);
   
  beta  ~ normal(0, sigma_beta);        
  gamma ~ normal(0, sigma_gamma);       

  // Scale parameters
  sigma_alpha  ~ normal(0, sigma);
  sigma_delta  ~ normal(0, sigma);
  sigma_beta   ~ normal(0, sigma);
  sigma_gamma  ~ normal(0, sigma);
  sigma ~ normal(0, 1);
  
  // [신규] Beta precision prior (값이 클수록 분산이 작음)
  phi ~ gamma(0.1, 0.1); 
}

generated quantities {
  // --- Prediction for TT+1 ---
  vector[P] alpha_TTp1;
  real delta_TTp1;
  vector[P] mu_pred_prob; // 확률 스케일
  vector[P] y_pred_share;

  delta_TTp1 = normal_rng(0, sigma_delta);

  {
    vector[P - 1] alpha0_TTp1_free;
    for (p in 1:(P - 1)) {
      alpha0_TTp1_free[p] = normal_rng(alpha0[p, TT], sigma_alpha);
      alpha_TTp1[p] = alpha0_TTp1_free[p];
    }
    // 예측 시점에서도 로짓 스케일 합이 0이 되도록 제약
    alpha_TTp1[P] = (0.0 - dot_product(Pop_weight[1:(P - 1), TT + 1], alpha0_TTp1_free)) 
                    / Pop_weight[P, TT + 1];
  }

  // 예측값 생성 (Linear -> Inverse Logit -> Beta RNG)
  {
      vector[P] eta_pred;
      eta_pred = alpha_TTp1 + X_pred * beta + Z_pred * gamma + delta_TTp1;
      mu_pred_prob = inv_logit(eta_pred);
  }

  for (p in 1:P) {
    // beta_proportion_rng를 사용하여 샘플링
    y_pred_share[p] = beta_proportion_rng(mu_pred_prob[p], phi);
  }
}
