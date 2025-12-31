data {
  int<lower=2> P;                 // provinces (16 or 17)
  int<lower=1> TT;                // time points (observed)
  int<lower=1> N;                 // N = P * TT
  int<lower=1> K;                 // number of national predictors (X)
  int<lower=1> L;                 // number of regional predictors (Z)

  matrix[N, K] X;
  matrix[N, L] Z;
  
  // 예측용 데이터 (TT+1 시점)
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
  // 지역별 추세 (Soft Constraint 적용)
  matrix[P, TT + 1] z_alpha; 

  real<lower=0, upper=0.99> rho;  
  real<lower=0> sigma_alpha;

  // 회귀 계수
  vector[K] beta;
  vector[L] gamma;
  
  // [삭제됨] delta (시간 고정 효과) 제거
  // vector[TT] delta; 

  // 오차 분산
  real<lower=0> sigma_epsilon; 
}

transformed parameters {
  matrix[P, TT + 1] alpha;

  // AR(1) 과정 생성
  // t=1
  alpha[, 1] = Alpha_init + sigma_alpha * z_alpha[, 1];

  // t=2 ~ TT+1
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
  
  // [삭제됨] delta에 대한 Prior 제거
  // delta ~ normal(0, 0.5); 

  sigma_epsilon ~ exponential(3); 

  // Soft Constraint: "알파의 가중 평균은 0이어야 한다"
  // delta가 사라졌으므로 이 제약조건은 더욱 중요해집니다 (전국 레벨은 X*beta가 담당)
  for (t in 1:(TT + 1)) {
    sum(Pop_weight[, t] .* alpha[, t]) ~ normal(0, 0.001 * P); 
  }

  // Likelihood
  for (t in 1:TT) {
    int start = 1 + (t - 1) * P;
    int end_  = t * P;
    
    Y_logit[start:end_] ~ normal(
      alpha[, t] 
      + X[start:end_] * beta 
      + Z[start:end_] * gamma,
      // [삭제됨] + delta[t] 제거 -> 이제 X*beta가 전국 스윙을 설명해야 함
      sigma_epsilon
    );
  }
}

generated quantities {
  vector[P] y_pred_logit; 
  vector[P] y_pred_share; 

  for (p in 1:P) {
    // 예측 식에서도 delta 제거
    // 미래의 전국 판세는 오직 X_pred(경제/정치 변수)로만 예측됨
    real mu = alpha[p, TT + 1] 
              + dot_product(X_pred[p, ], beta) 
              + dot_product(Z_pred[p, ], gamma);
    
    // 샘플링
    y_pred_logit[p] = normal_rng(mu, sigma_epsilon);
    y_pred_share[p] = inv_logit(y_pred_logit[p]);
  }
}