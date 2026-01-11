data {
  int<lower=2> P;                 // number of provinces
  int<lower=1> TT;                // fit period
  int<lower=1> N;                 // N = P * TT
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
  real eps = 1e-6;

  for (n in 1:N) {
    real y_adj = fmin(1 - eps, fmax(eps, Y[n]));
    Y_logit[n] = logit(y_adj);
  }
}

parameters {
  // [Alpha] NCP Random Walk
  matrix[P - 1, TT] z_alpha;

  // [Beta, Gamma] Coefficients
  vector[K] beta;
  vector[L] gamma;

  // [처방 1] 시작점(Intercept)을 자유롭게 해방!
  real delta_start;         // t=1 시점의 초기값
  vector[TT - 1] z_delta;   // t=2~TT 변화량

  // [Scales]
  real<lower=0> sigma_alpha;      
  real<lower=0> sigma_delta;      
  real<lower=0> sigma_epsilon;    
  
  // [처방 2] 이상치 대응을 위한 Student-t 자유도 (선택 사항)
  // real<lower=2> nu; 
}

transformed parameters {
  matrix[P - 1, TT] alpha0;
  matrix[P, TT] alpha;
  vector[TT] delta;

  // ---- Alpha RW (동일) ----
  alpha0[, 1] = Alpha_init[1:(P - 1)] + sigma_alpha * z_alpha[, 1];
  for (t in 2:TT) {
    alpha0[, t] = alpha0[, t - 1] + sigma_alpha * z_alpha[, t];
  }

  alpha[1:(P - 1), ] = alpha0;
  for (t in 1:TT) {
    alpha[P, t] = -dot_product(Pop_weight[1:(P - 1), t], alpha0[, t]) / Pop_weight[P, t];
  }

  // ---- [처방 1 적용] Delta RW (시작점 자유화) ----
  delta[1] = delta_start; // 0이 아니라 파라미터로 받음
  for (t in 2:TT) {
    delta[t] = delta[t - 1] + sigma_delta * z_delta[t - 1];
  }
}

model {
  // ---------------- Priors (조금 더 넓게) ----------------
  to_vector(z_alpha) ~ std_normal(); // normal(0,1)과 동일
  z_delta ~ std_normal();
  
  // [처방 1] 초기값 Prior (넓게)
  delta_start ~ normal(0, 5); 

  beta  ~ normal(0, 2); // 0.5 -> 2.0 (로짓 스케일 고려)
  gamma ~ normal(0, 2);

  // [처방 3] Sigma Prior 완화 (0.2 -> 1.0)
  sigma_alpha   ~ normal(0, 1.0); // 로짓 스케일에서 0.2는 너무 빡빡함
  sigma_delta   ~ normal(0, 1.0);
  sigma_epsilon ~ normal(0, 1.0);

  // --------------- Likelihood ---------------
  for (t in 1:TT) {
    int start = 1 + (t - 1) * P;
    int end_  = t * P;

    matrix[P, K] Xt = X[start:end_, ];
    matrix[P, L] Zt = Z[start:end_, ];

    // [처방 2] Normal 대신 Student-t 사용 권장 (df=7 정도면 적당)
    // 데이터가 이상치(0이나 1 근처)에 덜 민감해짐
    Y_logit[start:end_] ~ student_t(20, 
      alpha[, t] + Xt * beta + Zt * gamma + rep_vector(delta[t], P),
      sigma_epsilon
    );
  }
}

generated quantities {
  // (예측 부분 코드는 기존 로직과 동일하게 유지하되, 
  // delta_TT1 계산 시 delta[TT]가 이제 올바른 레벨을 가지고 있으므로 그대로 사용 가능)
  
  vector[P] alpha_TT1;
  real delta_TT1;
  vector[P] y_pred_logit;
  vector[P] y_pred_share;

  // Alpha forecast (동일)
  {
     vector[P - 1] alpha0_TT1_free;
     vector[P - 1] z_alpha_pred;
     for (p in 1:(P - 1)) z_alpha_pred[p] = normal_rng(0, 1);
     
     for (p in 1:(P - 1)) {
       alpha0_TT1_free[p] = alpha0[p, TT] + sigma_alpha * z_alpha_pred[p];
       alpha_TT1[p] = alpha0_TT1_free[p];
     }
     alpha_TT1[P] = -dot_product(Pop_weight[1:(P - 1), TT + 1], alpha0_TT1_free) 
                    / Pop_weight[P, TT + 1];
  }

  // Delta forecast
  delta_TT1 = delta[TT] + sigma_delta * normal_rng(0, 1);

  // Prediction
  y_pred_logit = alpha_TT1 + X_pred * beta + Z_pred * gamma + rep_vector(delta_TT1, P);
  
  for (p in 1:P) y_pred_share[p] = inv_logit(y_pred_logit[p]);
}
