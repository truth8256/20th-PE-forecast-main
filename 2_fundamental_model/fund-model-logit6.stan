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
  // 1. Alpha (지역별 편차): Random Walk
  matrix[P, TT] z_alpha; 

  // [삭제됨] Delta 관련 파라미터 전부 제거
  // real delta_start;         
  // vector[TT - 1] z_delta;   
  // real<lower=0> sigma_delta;   

  // 2. Coefficients
  vector[K] beta;
  vector[L] gamma;

  // 3. Scales
  real<lower=0> sigma_alpha;      
  real<lower=0> sigma_epsilon;    
}

transformed parameters {
  matrix[P, TT] alpha;
  // [삭제됨] vector[TT] delta;

  // [삭제됨] Delta RW 계산 로직 제거

  // [Alpha RW] 
  // NCP: alpha는 0 주변에서 움직임 (Random Walk)
  alpha[, 1] = sigma_alpha * z_alpha[, 1];
  for (t in 2:TT) {
    alpha[, t] = alpha[, t - 1] + sigma_alpha * z_alpha[, t];
  }
}

model {
  // --- Priors ---
  to_vector(z_alpha) ~ std_normal();
  
  // [삭제됨] Delta Priors 제거
  // z_delta ~ std_normal();
  // delta_start ~ normal(0, 5); 
  // sigma_delta   ~ normal(0, 0.5);

  beta  ~ normal(0, 1);
  gamma ~ normal(0, 1);

  sigma_alpha   ~ normal(0, 0.1); // [추천] 0.5 -> 0.1 (지역 구도 안정화)
  sigma_epsilon ~ normal(0, 0.5); 

  // --- Soft Constraint (필수) ---
  // Delta가 없으므로 이 제약이 "X*beta가 전국 평균임"을 정의해줍니다.
  for (t in 1:TT) {
     target += normal_lpdf(dot_product(Pop_weight[, t], alpha[, t]) | 0, 0.001);
  }

  // --- Likelihood ---
  for (t in 1:TT) {
    int start = 1 + (t - 1) * P;
    int end_  = t * P;
    
    Y_logit[start:end_] ~ normal(
      // [수정] delta[t] 제거 -> X*beta가 기준선이 됨
      alpha[, t] + X[start:end_]*beta + Z[start:end_]*gamma,
      sigma_epsilon
    );
  }
}

generated quantities {
  vector[P] alpha_TT1;
  // [삭제됨] real delta_TT1;
  
  vector[P] y_pred_logit;
  vector[P] y_pred_share;

  // [삭제됨] Delta 예측 제거
  // delta_TT1 = delta[TT] + sigma_delta * normal_rng(0, 1);

  // Alpha 예측 (Random Walk)
  {
     vector[P] z_alpha_pred;
     for (p in 1:P) z_alpha_pred[p] = normal_rng(0, 1);
     
     // Random Walk: 어제 값에서 출발
     alpha_TT1 = alpha[, TT] + sigma_alpha * z_alpha_pred;

     // Hard Centering
     real weighted_sum = dot_product(Pop_weight[, TT + 1], alpha_TT1);
     real total_weight = sum(Pop_weight[, TT + 1]);
     alpha_TT1 = alpha_TT1 - (weighted_sum / total_weight);
  }

  // [수정] 예측식에서 delta_TT1 제거
  y_pred_logit = alpha_TT1 + X_pred * beta + Z_pred * gamma;
  
  for (p in 1:P) y_pred_share[p] = inv_logit(y_pred_logit[p]);
}
