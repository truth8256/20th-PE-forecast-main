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
  
  // Alpha_init은 여기서 직접 쓰지 않고 추정하게 둡니다
  vector[P] Alpha_init;           
  matrix[P, TT + 1] Pop_weight;   
}

transformed data {
  vector[N] Y_logit;
  Y_logit = logit(Y); 
}

parameters {
  // 1. Alpha (지역별 편차): P개 모두 자유롭게 선언
  matrix[P, TT] z_alpha; 

  // 2. Delta (전국 추세): 시작점 자유화
  real delta_start;         
  vector[TT - 1] z_delta;   

  // 3. Coefficients
  vector[K] beta;
  vector[L] gamma;

  // 4. Scales
  real<lower=0> sigma_alpha;      
  real<lower=0> sigma_delta;      
  real<lower=0> sigma_epsilon;    
}

transformed parameters {
  matrix[P, TT] alpha;
  vector[TT] delta;

  // [Delta RW]
  delta[1] = delta_start;
  for (t in 2:TT) {
    delta[t] = delta[t - 1] + sigma_delta * z_delta[t - 1];
  }

  // [Alpha RW] 
  // NCP: alpha는 0 주변에서 움직임
  alpha[, 1] = sigma_alpha * z_alpha[, 1];
  for (t in 2:TT) {
    alpha[, t] = alpha[, t - 1] + sigma_alpha * z_alpha[, t];
  }
}

model {
  // --- Priors ---
  to_vector(z_alpha) ~ std_normal();
  z_delta ~ std_normal();
  delta_start ~ normal(0, 5); 

  beta  ~ normal(0, 1);
  gamma ~ normal(0, 1);

  sigma_alpha   ~ normal(0, 0.5); 
  sigma_delta   ~ normal(0, 0.5);
  sigma_epsilon ~ normal(0, 0.5); 

  // ----------------------------------------------------------------
  // [핵심 처방: Soft Constraint 부활]
  // "전국 추세(Delta)와 지역 편차(Alpha)를 구분해라!"
  // 매 시점마다 인구 가중 평균이 0이 되도록 '강력하게(sd=0.001)' 유도합니다.
  // 나눗셈을 하지 않으므로 '제주도 지렛대 효과'는 발생하지 않습니다.
  // ----------------------------------------------------------------
  for (t in 1:TT) {
     target += normal_lpdf(dot_product(Pop_weight[, t], alpha[, t]) | 0, 0.001);
  }

  // --- Likelihood (Normal) ---
  for (t in 1:TT) {
    int start = 1 + (t - 1) * P;
    int end_  = t * P;
    
    Y_logit[start:end_] ~ normal(
      rep_vector(delta[t], P) + alpha[, t] + X[start:end_]*beta + Z[start:end_]*gamma,
      sigma_epsilon
    );
  }
}

generated quantities {
  vector[P] alpha_TT1;
  real delta_TT1;
  vector[P] y_pred_logit;
  vector[P] y_pred_share;

  // 예측 단계
  delta_TT1 = delta[TT] + sigma_delta * normal_rng(0, 1);

  {
     vector[P] z_alpha_pred;
     for (p in 1:P) z_alpha_pred[p] = normal_rng(0, 1);
     
     // 1. Alpha 예측
     alpha_TT1 = alpha[, TT] + sigma_alpha * z_alpha_pred;

     // 2. 예측값 보정 (Hard Centering)
     // Soft Constraint로 이미 거의 0이지만, 깔끔하게 0으로 딱 맞춤
     real weighted_sum = dot_product(Pop_weight[, TT + 1], alpha_TT1);
     real total_weight = sum(Pop_weight[, TT + 1]);
     
     alpha_TT1 = alpha_TT1 - (weighted_sum / total_weight);
  }

  y_pred_logit = rep_vector(delta_TT1, P) + alpha_TT1 + X_pred * beta + Z_pred * gamma;
  for (p in 1:P) y_pred_share[p] = inv_logit(y_pred_logit[p]);
}
