
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

  vector<lower=0, upper=1>[N] Y;
  vector[P] Alpha_init;            
  matrix[P, TT + 1] Pop_weight;    
}

transformed data {
  vector[N] Y_logit;
  Y_logit = logit(Y);
}

parameters {
  // 1. 지역별 추세 (AR1)
  matrix[P, TT] z_alpha; 
  real<lower=0, upper=0.99> rho;

  // 2. 전국 추세 (Robust Random Walk)
  real delta_start;
  vector[TT - 1] z_delta; // ★핵심: 정규분포가 아님

  // 3. Coefficients
  vector[K] beta;
  vector[L] gamma;
  
  // 4. Scales
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_delta; 
  real<lower=0> sigma_epsilon; 
  
  // ★[NEW] 자유도 (Degrees of Freedom)
  // nu가 작을수록(3~4) "대격변"을 잘 허용함
  real<lower=2> nu; 
}

transformed parameters {
  matrix[P, TT] alpha;
  vector[TT] delta;

  // [Delta RW] 
  delta[1] = delta_start;
  for (t in 2:TT) {
    // z_delta가 t-분포에서 나왔으므로, 가끔 큰 충격이 와도 됨
    delta[t] = delta[t - 1] + sigma_delta * z_delta[t - 1];
  }

  // [Alpha AR1]
  alpha[, 1] = Alpha_init + sigma_alpha * z_alpha[, 1];
  for (t in 2:TT) {
    alpha[, t] = rho * alpha[, t - 1] + sigma_alpha * z_alpha[, t];
  }
}

model {
  // --- Priors (Robust Strategy) ---
  
  // 1. 전국 추세 (Delta): Student-t 분포 적용
  // nu(자유도)가 4 근처면 "가끔 큰 파도가 친다"는 뜻
  nu ~ gamma(2, 0.1); 
  
  // ★핵심: z_delta가 t-분포를 따름
  for (t in 1:(TT-1)) {
      z_delta[t] ~ student_t(nu, 0, 1); 
  }
  
  delta_start ~ normal(0, 2);
  
  // Sigma는 타이트하게(0.1) 유지! 
  // (평소엔 좁게, 비상시엔 nu 덕분에 넓게)
  sigma_delta ~ normal(0, 0.1); 

  // 2. 지역 추세 (Alpha): 정규분포 (지역 구도는 안정적)
  to_vector(z_alpha) ~ std_normal();
  rho ~ beta(10, 1);            
  sigma_alpha ~ normal(0, 0.1); 

  // 3. Coefficients & Noise
  beta  ~ normal(0, 0.5); 
  gamma ~ normal(0, 0.5);
  sigma_epsilon ~ normal(0, 0.05); 

  // --- Soft Constraint ---
  for (t in 1:TT) {
     target += normal_lpdf(dot_product(Pop_weight[, t], alpha[, t]) | 0, 0.001);
  }

  // --- Likelihood ---
  for (t in 1:TT) {
    int start = 1 + (t - 1) * P;
    int end_  = t * P;
    
    Y_logit[start:end_] ~ normal(
      rep_vector(delta[t], P)      
      + alpha[, t]                 
      + X[start:end_] * beta 
      + Z[start:end_] * gamma,
      sigma_epsilon
    );
  }
}

generated quantities {
  vector[P] alpha_TT1;
  real delta_TT1;
  
  vector[P] y_pred_logit;
  vector[P] y_pred_share; 

  // Forecast Delta (Robust RW)
  // 예측 시점에서도 "혹시 모를 대격변" 가능성을 열어둠
  delta_TT1 = delta[TT] + sigma_delta * student_t_rng(nu, 0, 1);

  // Forecast Alpha (AR1)
  {
     vector[P] z_alpha_pred;
     for (p in 1:P) z_alpha_pred[p] = normal_rng(0, 1);
     
     alpha_TT1 = rho * alpha[, TT] + sigma_alpha * z_alpha_pred;
     
     real weighted_sum = dot_product(Pop_weight[, TT + 1], alpha_TT1);
     real total_weight = sum(Pop_weight[, TT + 1]);
     alpha_TT1 = alpha_TT1 - (weighted_sum / total_weight);
  }

  y_pred_logit = rep_vector(delta_TT1, P) 
                 + alpha_TT1 
                 + X_pred * beta 
                 + Z_pred * gamma;
                 
  for (p in 1:P) {
    y_pred_share[p] = inv_logit(y_pred_logit[p]);
  }
}
