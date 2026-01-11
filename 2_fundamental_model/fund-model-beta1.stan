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

  // [주의] Y는 0과 1을 포함하면 안 됨 (0 < Y < 1)
  // 만약 0이나 1이 있다면 R에서 0.001 / 0.999 등으로 보정해서 넣어야 함
  vector<lower=0, upper=1>[N] Y;
  
  vector[P] Alpha_init;
  matrix[P, TT + 1] Pop_weight;   
}

parameters {
  // 1. 지역별 추세 (Alpha): P개 모두 자유롭게 선언 (Soft Constraint용)
  matrix[P, TT] z_alpha; 

  // 2. 전국 추세 (Delta): 시작점 자유화 & RW
  real delta_start;
  vector[TT - 1] z_delta;

  // 3. 회귀 계수
  vector[K] beta;
  vector[L] gamma;
  
  // 4. Scales (Reality-Based Strategy 적용)
  real<lower=0> sigma_alpha; 
  real<lower=0> sigma_delta; 
  
  // 5. Beta Precision (높을수록 분산 작음)
  real<lower=0> phi; 
}

transformed parameters {
  matrix[P, TT] alpha;
  vector[TT] delta;

  // [Delta RW] 전국 파도
  delta[1] = delta_start;
  for (t in 2:TT) {
    delta[t] = delta[t - 1] + sigma_delta * z_delta[t - 1];
  }

  // [Alpha RW] 지역별 고유 변동 (NCP)
  // t=1
  alpha[, 1] = Alpha_init + sigma_alpha * z_alpha[, 1];
  // t=2~TT
  for (t in 2:TT) {
    alpha[, t] = alpha[, t - 1] + sigma_alpha * z_alpha[, t];
  }
}

model {
  // --- Priors (Reality-Based) ---
  
  // 1. 전국 추세: "크게 흔들릴 수 있다"
  z_delta ~ std_normal();
  delta_start ~ normal(0, 5); 
  sigma_delta ~ normal(0, 0.5); // 0.9 변동 허용

  // 2. 지역 추세: "지역 구도는 안정적이다"
  to_vector(z_alpha) ~ std_normal();
  sigma_alpha ~ normal(0, 0.1); // 지역별 이탈 억제

  // 3. Coefficients
  beta  ~ normal(0, 0.5);
  gamma ~ normal(0, 0.5);

  // 4. Beta Precision Prior
  // phi가 너무 작으면(분산 큼) 예측 범위가 넓어짐.
  // 적당히 큰 값(분산 작음)을 선호하도록 감마 분포 설정
  phi ~ gamma(2, 0.1); // Mean=20, Variance=200 (유연함)

  // --- Soft Constraint ---
  // "Alpha(지역편차)의 인구 가중 합은 0이어야 한다" (식별성 확보)
  for (t in 1:TT) {
     target += normal_lpdf(dot_product(Pop_weight[, t], alpha[, t]) | 0, 0.0001);
  }

  // --- Likelihood (Beta Regression) ---
  {
    vector[N] mu;
    vector[N] eta; // Linear Predictor
    
    for (t in 1:TT) {
      int start = 1 + (t - 1) * P;
      int end_  = t * P;
      
      // 선형 결합 (Linear Predictor) -> Logit Scale
      eta[start:end_] = 
        rep_vector(delta[t], P) 
        + alpha[, t] 
        + X[start:end_] * beta 
        + Z[start:end_] * gamma;
    }
    
    // Link Function: Logit^-1 -> Probability Scale (0~1)
    mu = inv_logit(eta);
    
    // Stan Beta Likelihood (Mean-Precision parametrization)
    Y ~ beta_proportion(mu, phi);
  }
}

generated quantities {
  vector[P] alpha_TT1;
  real delta_TT1;
  vector[P] mu_pred; // 확률 평균
  vector[P] y_pred_share;

  // Forecast Delta
  delta_TT1 = delta[TT] + sigma_delta * normal_rng(0, 1);

  // Forecast Alpha
  {
    vector[P] z_alpha_pred;
    for (p in 1:P) z_alpha_pred[p] = normal_rng(0, 1);
    
    alpha_TT1 = alpha[, TT] + sigma_alpha * z_alpha_pred;
    
    // Hard Centering for Forecast
    real weighted_sum = dot_product(Pop_weight[, TT + 1], alpha_TT1);
    real total_weight = sum(Pop_weight[, TT + 1]);
    alpha_TT1 = alpha_TT1 - (weighted_sum / total_weight);
  }

  // Calculate Prediction
  {
    vector[P] eta_pred;
    eta_pred = rep_vector(delta_TT1, P) 
               + alpha_TT1 
               + X_pred * beta 
               + Z_pred * gamma;
               
    mu_pred = inv_logit(eta_pred);
  }

  // Beta RNG
  for (p in 1:P) {
    y_pred_share[p] = beta_proportion_rng(mu_pred[p], phi);
  }
}
