data {
  int<lower=2> P;                 // provinces
  int<lower=1> TT;                // time points (observed)
  int<lower=1> N;                 // N = P * TT
  int<lower=1> K;
  int<lower=1> L;

  matrix[N, K] X;
  matrix[N, L] Z;
  
  // [추가] 예측(TT+1 시점)을 위한 설명변수 데이터
  matrix[P, K] X_pred; 
  matrix[P, L] Z_pred;
  
  vector<lower=0, upper=1>[N] Y;
  vector[P] Alpha_init;           
  matrix[P, TT + 1] Pop_weight;   
}

transformed data {
  vector[N] Y_logit;
  real eps = 1e-6;
  for (n in 1:N) {
    real y_adj = fmin(1 - eps, fmax(eps, Y[n]));
    Y_logit[n] = logit(y_adj);
  }
}

parameters {
  matrix[P, TT + 1] z_alpha; 

  real<lower=0, upper=0.99> rho;  
  real<lower=0> sigma_alpha;

  vector[K] beta;
  vector[L] gamma;
  
  // [수정] 예측 시점(TT+1)의 시간 효과도 포함하기 위해 크기를 TT+1로 확장
  vector[TT + 1] delta; 

  real<lower=0> sigma_epsilon; 
}

transformed parameters {
  matrix[P, TT + 1] alpha;

  // t=1
  alpha[, 1] = Alpha_init + sigma_alpha * z_alpha[, 1];

  // t=2 ~ TT+1 (여기에 이미 TT+1 시점의 지역 추세가 계산됨)
  for (t in 2:(TT + 1)) {
    alpha[, t] = rho * alpha[, t - 1] + sigma_alpha * z_alpha[, t];
  }
}

model {
  // Priors
  to_vector(z_alpha) ~ std_normal();

  rho ~ beta(8, 2); 
  sigma_alpha ~ exponential(5); 

  beta  ~ normal(0, 0.5); 
  gamma ~ normal(0, 0.5);
  delta ~ normal(0, 0.5); // TT+1번째 delta는 Prior에서 샘플링됨

  sigma_epsilon ~ exponential(3); 

  // Soft Constraint
  for (t in 1:(TT + 1)) {
    sum(Pop_weight[, t] .* alpha[, t]) ~ normal(0, 0.001 * P); 
  }

  // Likelihood (관측된 TT 시점까지만 학습)
  for (t in 1:TT) {
    int start = 1 + (t - 1) * P;
    int end_  = t * P;
    
    Y_logit[start:end_] ~ normal(
      alpha[, t] 
      + X[start:end_] * beta 
      + Z[start:end_] * gamma 
      + delta[t], 
      sigma_epsilon
    );
  }
}

// [추가] 예측 블록
generated quantities {
  vector[P] y_pred_logit; // 로짓 스케일 예측값
  vector[P] y_pred_share; // 0~1 득표율 예측값

  for (p in 1:P) {
    // 1. 선형 예측식 계산 (TT+1 시점의 alpha, delta 사용)
    real mu = alpha[p, TT + 1] 
              + dot_product(X_pred[p, ], beta) 
              + dot_product(Z_pred[p, ], gamma) 
              + delta[TT + 1];
    
    // 2. 오차항(sigma_epsilon)을 포함하여 샘플링 (불확실성 반영)
    y_pred_logit[p] = normal_rng(mu, sigma_epsilon);
    
    // 3. 로짓 -> 확률(득표율) 변환
    y_pred_share[p] = inv_logit(y_pred_logit[p]);
  }
}
