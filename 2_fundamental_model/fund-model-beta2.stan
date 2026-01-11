data {
  int<lower=2> P;                 
  int<lower=1> TT;                
  int<lower=1> N;                 
  int<lower=1> K; 
  int<lower=1> L; 

  matrix[N, K] X; 
  matrix[N, L] Z; 

  matrix[P, K] X_pred;
  matrix[P, L] Z_pred;

  // [주의] Y는 0과 1을 포함하지 않는 구간 (0 < Y < 1)
  vector<lower=0, upper=1>[N] Y;
  
  vector[P] Alpha_init;
  matrix[P, TT + 1] Pop_weight;   
}

parameters {
  // 1. 지역별 추세 (Alpha): Random Walk
  matrix[P, TT] z_alpha; 

  // [삭제됨] Delta 관련 파라미터 전부 제거
  // real delta_start;
  // vector[TT - 1] z_delta;
  // real<lower=0> sigma_delta;

  // 2. 회귀 계수
  vector[K] beta;
  vector[L] gamma;
  
  // 3. Scales
  real<lower=0> sigma_alpha; 
  
  // 4. Beta Precision (높을수록 분산 작음)
  real<lower=0> phi; 
}

transformed parameters {
  matrix[P, TT] alpha;
  // [삭제됨] vector[TT] delta;

  // [삭제됨] Delta RW 로직 제거

  // [Alpha RW] 지역별 고유 변동 (Random Walk)
  // NCP: alpha는 0 주변에서 움직임
  alpha[, 1] = Alpha_init + sigma_alpha * z_alpha[, 1];
  for (t in 2:TT) {
    alpha[, t] = alpha[, t - 1] + sigma_alpha * z_alpha[, t];
  }
}

model {
  // --- Priors ---
  
  // [삭제됨] Delta Priors 제거
  // z_delta ~ std_normal(); ...

  // 1. 지역 추세 (Alpha)
  to_vector(z_alpha) ~ std_normal();
  sigma_alpha ~ normal(0, 0.1); // 지역 구도 안정화 (타이트하게)

  // 2. Coefficients
  beta  ~ normal(0, 0.5);
  gamma ~ normal(0, 0.5);

  // 3. Beta Precision Prior
  // 유연한 분산 허용
  phi ~ gamma(2, 0.1); 

  // --- Soft Constraint (매우 중요) ---
  // Delta가 없으므로 "지역 편차의 합 = 0" 조건이 
  // "X*beta가 곧 전국 평균이다"라는 정의를 완성합니다.
  for (t in 1:TT) {
     target += normal_lpdf(dot_product(Pop_weight[, t], alpha[, t]) | 0, 0.001);
  }

  // --- Likelihood (Beta Regression) ---
  {
    vector[N] mu;
    vector[N] eta; 
    
    for (t in 1:TT) {
      int start = 1 + (t - 1) * P;
      int end_  = t * P;
      
      // 선형 결합: Delta 항 제거됨
      // 오직 X*beta(전국요인) + Z*gamma(지역요인) + alpha(지역잔차)
      eta[start:end_] = 
        alpha[, t] 
        + X[start:end_] * beta 
        + Z[start:end_] * gamma;
    }
    
    // Link Function: Logit^-1 -> Probability Scale
    mu = inv_logit(eta);
    
    // Beta Likelihood
    Y ~ beta_proportion(mu, phi);
  }
}

generated quantities {
  vector[P] alpha_TT1;
  // [삭제됨] real delta_TT1;
  
  vector[P] mu_pred; 
  vector[P] y_pred_share;

  // [삭제됨] Forecast Delta

  // Forecast Alpha (RW)
  {
    vector[P] z_alpha_pred;
    for (p in 1:P) z_alpha_pred[p] = normal_rng(0, 1);
    
    // 어제 값에서 출발 (Random Walk)
    alpha_TT1 = alpha[, TT] + sigma_alpha * z_alpha_pred;
    
    // Hard Centering for Forecast
    real weighted_sum = dot_product(Pop_weight[, TT + 1], alpha_TT1);
    real total_weight = sum(Pop_weight[, TT + 1]);
    alpha_TT1 = alpha_TT1 - (weighted_sum / total_weight);
  }

  // Calculate Prediction
  {
    vector[P] eta_pred;
    
    // 예측식에서도 Delta 제거
    eta_pred = alpha_TT1 
               + X_pred * beta 
               + Z_pred * gamma;
               
    mu_pred = inv_logit(eta_pred);
  }

  // Beta RNG
  for (p in 1:P) {
    y_pred_share[p] = beta_proportion_rng(mu_pred[p], phi);
  }
}
