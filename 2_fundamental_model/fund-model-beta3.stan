// Likelihood: beta_proportion 유지 (0~1 구간의 엄밀함).
// Alpha: AR(1) 구조 도입 (기존 베타 모델은 RW였음 $\to$ m4처럼 rho 도입).
// Delta: Random Walk를 쓰되, sigma를 0.2~0.3 정도로 넉넉하게 줍니다. (0.5는 너무 넓고, 0.1은 너무 좁음. 21대 총선을 잡으려면 0.25 정도가 적절)
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

  // [주의] Y는 (0, 1) 구간 (0과 1 포함 X)
  vector<lower=0, upper=1>[N] Y;
  
  vector[P] Alpha_init;
  matrix[P, TT + 1] Pop_weight;   
}

parameters {
  // 1. Alpha (지역별 추세): AR(1) 적용을 위한 raw parameters
  matrix[P, TT] z_alpha; 
  real<lower=0, upper=0.99> rho;   // AR(1) 지속성 계수

  // 2. Delta (전국 추세): Random Walk
  real delta_start;
  vector[TT - 1] z_delta;

  // 3. Coefficients
  vector[K] beta;
  vector[L] gamma;
  
  // 4. Scales
  real<lower=0> sigma_alpha; 
  real<lower=0> sigma_delta; // 21대 총선을 잡기 위한 핵심 키
  
  // 5. Beta Precision
  real<lower=0> phi; 
}

transformed parameters {
  matrix[P, TT] alpha;
  vector[TT] delta;

  // [Delta RW]
  delta[1] = delta_start;
  for (t in 2:TT) {
    delta[t] = delta[t - 1] + sigma_delta * z_delta[t - 1];
  }

  // [Alpha AR(1)] m4 모델의 성공 요인 이식
  // t=1
  alpha[, 1] = Alpha_init + sigma_alpha * z_alpha[, 1];
  // t=2~TT
  for (t in 2:TT) {
    alpha[, t] = rho * alpha[, t - 1] + sigma_alpha * z_alpha[, t];
  }
}

model {
  // --- Priors ---
  
  // 1. 전국 추세 (Delta): "파도를 허용하라"
  z_delta ~ std_normal();
  delta_start ~ normal(0, 2);
  
  // [핵심 튜닝] 
  // 0.1(Tight) -> 실패, 0.5(Loose) -> Width 폭발
  // 0.25 (약 6~7%p 변동 허용): 21대 총선 편향을 잡으면서도 Width를 제어하는 타협점
  sigma_delta ~ normal(0, 0.25); 

  // 2. 지역 추세 (Alpha): "지역 구도는 꽉 잡아라"
  to_vector(z_alpha) ~ std_normal();
  rho ~ beta(10, 1);             // 강한 관성 (0.9 이상)
  sigma_alpha ~ normal(0, 0.1);  // 지역 이탈 최소화

  // 3. Coefficients
  beta  ~ normal(0, 0.5);
  gamma ~ normal(0, 0.5);

  // 4. Beta Precision
  phi ~ gamma(50, 1.0); 

  // --- Soft Constraint ---
  // "Alpha 가중합 = 0" -> X*beta가 전국 기준선이 됨
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
      
      // 선형 예측식
      eta[start:end_] = 
        rep_vector(delta[t], P) 
        + alpha[, t] 
        + X[start:end_] * beta 
        + Z[start:end_] * gamma;
    }
    
    // Link Function: Logit^-1
    mu = inv_logit(eta);
    
    // Beta Likelihood
    Y ~ beta_proportion(mu, phi);
  }
}

generated quantities {
  vector[P] alpha_TT1;
  real delta_TT1;
  vector[P] mu_pred; 
  vector[P] y_pred_share;

  // Forecast Delta (RW)
  delta_TT1 = delta[TT] + sigma_delta * normal_rng(0, 1);

  // Forecast Alpha (AR1)
  {
    vector[P] z_alpha_pred;
    for (p in 1:P) z_alpha_pred[p] = normal_rng(0, 1);
    
    // AR(1) Logic applied to Forecast
    alpha_TT1 = rho * alpha[, TT] + sigma_alpha * z_alpha_pred;
    
    // Hard Centering
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
