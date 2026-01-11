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
  // 데이터가 안전하므로 단순 변환
  Y_logit = logit(Y);
}

parameters {
  // 1. 지역별 추세 (AR1)
  matrix[P, TT] z_alpha; 
  real<lower=0, upper=0.99> rho;   // 지속성 (Persistence)

  // 2. [부활] 전국 추세 (Random Walk) - 거대한 파도를 설명
  real delta_start;
  vector[TT - 1] z_delta;

  // 3. 회귀 계수
  vector[K] beta;
  vector[L] gamma;
  
  // 4. Scales (핵심 전략 적용 대상)
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_delta;
  real<lower=0> sigma_epsilon; 
}

transformed parameters {
  matrix[P, TT] alpha;
  vector[TT] delta;

  // [Delta RW] 전국 판세는 자유롭게 흐름
  delta[1] = delta_start;
  for (t in 2:TT) {
    delta[t] = delta[t - 1] + sigma_delta * z_delta[t - 1];
  }

  // [Alpha AR1] 지역 구도는 끈적하게 유지됨 (Time Dependency)
  // t=1: 초기값에서 시작
  alpha[, 1] = Alpha_init + sigma_alpha * z_alpha[, 1];

  // t=2~TT: AR(1) Process
  // alpha(t) = rho * alpha(t-1) + shock
  for (t in 2:TT) {
    alpha[, t] = rho * alpha[, t - 1] + sigma_alpha * z_alpha[, t];
  }
}

model {
  // --- Priors (Reality-Based Strategy) ---
  
  // 1. 전국 추세: "무슨 일이든 일어날 수 있다" (Loose)
  z_delta ~ std_normal();
  delta_start ~ normal(0, 2);
  sigma_delta ~ normal(0, 0.1); 

  // 2. 지역 추세: "지역 민심은 쉽게 안 바뀐다" (Tight)
  to_vector(z_alpha) ~ std_normal();
  rho ~ beta(10, 1);             // 높은 지속성 선호 (0.8 근처)
  sigma_alpha ~ normal(0, 0.1); // 지역별 돌발행동 억제 (핵심!)

  // 3. 기타
  beta  ~ normal(0, 0.5); 
  gamma ~ normal(0, 0.5);
  sigma_epsilon ~ normal(0, 0.05); // 관측 오차는 작게

  // --- Soft Constraint (Identifiability) ---
  // "Alpha는 전국 평균(Delta)을 뺀 나머지여야 한다"
  for (t in 1:TT) {
     target += normal_lpdf(dot_product(Pop_weight[, t], alpha[, t]) | 0, 0.001);
  }

  // --- Likelihood ---
  for (t in 1:TT) {
    int start = 1 + (t - 1) * P;
    int end_  = t * P;
    
    Y_logit[start:end_] ~ normal(
      rep_vector(delta[t], P)      // 전국 파도
      + alpha[, t]                 // 지역 배
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

  // 예측 (Forecast)
  
  // 1. Delta: RW로 한 발자국 전진 (큰 불확실성 반영)
  delta_TT1 = delta[TT] + sigma_delta * normal_rng(0, 1);

  // 2. Alpha: AR1 논리로 한 발자국 전진 (작은 불확실성 반영)
  {
     vector[P] z_alpha_pred;
     for (p in 1:P) z_alpha_pred[p] = normal_rng(0, 1);
     
     // alpha(TT+1) = rho * alpha(TT) + shock
     alpha_TT1 = rho * alpha[, TT] + sigma_alpha * z_alpha_pred;
     
     // 예측 시점에도 중심 잡기 (Hard Centering)
     real weighted_sum = dot_product(Pop_weight[, TT + 1], alpha_TT1);
     real total_weight = sum(Pop_weight[, TT + 1]);
     alpha_TT1 = alpha_TT1 - (weighted_sum / total_weight);
  }

  // 3. 최종 결합
  y_pred_logit = rep_vector(delta_TT1, P) 
                 + alpha_TT1 
                 + X_pred * beta 
                 + Z_pred * gamma;
                 
  for (p in 1:P) {
    y_pred_share[p] = inv_logit(y_pred_logit[p]);
  }
}
