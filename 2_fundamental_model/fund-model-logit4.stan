data {
  int<lower=2> P;                 // 지역 수
  int<lower=1> TT;                // 관측 기간 (선거 횟수)
  int<lower=1> N;                 // 총 관측치 = P * TT
  int<lower=1> K;                 // 전국 레벨 변수 개수
  int<lower=1> L;                 // 지역 레벨 변수 개수

  matrix[N, K] X;                 // 전국 공통 변수 (Time-blocks)
  matrix[N, L] Z;                 // 지역 특화 변수

  matrix[P, K] X_pred;            // 예측 시점(TT+1) 공변량
  matrix[P, L] Z_pred;

  vector<lower=0, upper=1>[N] Y;  // 득표율 데이터 (0~1)
  
  // Alpha_init은 이 모델에서 직접 쓰지 않고, 데이터가 스스로 찾게 둡니다.
  // 다만 코드 호환성을 위해 받아두기는 합니다.
  vector[P] Alpha_init;           
  matrix[P, TT + 1] Pop_weight;   // 전국 지지율 계산(사후 층화)용 가중치
}

transformed data {
  vector[N] Y_logit;
  // 데이터가 0과 1 사이 안전한 구간에 있으므로 단순 변환
  Y_logit = logit(Y); 
}

parameters {
  // 1. 지역별 편차 (Alpha): P개 지역 모두 자유롭게 추정
  //    (sum=0 제약 없음. 대신 Prior가 0 근처로 잡아줌)
  matrix[P, TT] z_alpha; 

  // 2. 전국 공통 추세 (Delta)
  real delta_start;         // 시작점 자유화
  vector[TT - 1] z_delta;   // 변화량

  // 3. 회귀 계수
  vector[K] beta;
  vector[L] gamma;

  // 4. 스케일 파라미터 (Scales)
  real<lower=0> sigma_alpha;      // 지역별 편차의 크기
  real<lower=0> sigma_delta;      // 전국 추세의 변동성
  real<lower=0> sigma_epsilon;    // 관측 노이즈
}

transformed parameters {
  matrix[P, TT] alpha;
  vector[TT] delta;

  // [Delta RW] 전국 추세 (Random Walk)
  delta[1] = delta_start;
  for (t in 2:TT) {
    delta[t] = delta[t - 1] + sigma_delta * z_delta[t - 1];
  }

  // [Alpha RW] 지역별 편차 (Random Walk around 0)
  // NCP 적용: alpha는 '전국 추세 대비 편차'이므로 0에서 시작
  alpha[, 1] = sigma_alpha * z_alpha[, 1];
  for (t in 2:TT) {
    alpha[, t] = alpha[, t - 1] + sigma_alpha * z_alpha[, t];
  }
}

model {
  // --- Priors ---
  // 계층적 구조의 핵심: Alpha는 0을 중심으로 분포한다 (Soft Constraint 효과)
  to_vector(z_alpha) ~ std_normal(); 
  
  z_delta ~ std_normal();
  delta_start ~ normal(0, 5); // 시작점은 데이터에 맡김

  beta  ~ normal(0, 1);
  gamma ~ normal(0, 1);

  // Scales
  sigma_alpha   ~ normal(0, 0.5); 
  sigma_delta   ~ normal(0, 0.5);
  sigma_epsilon ~ normal(0, 0.5); 

  // --- Likelihood (Hierarchical Logit-Normal) ---
  // 개별 지역 예측값 = (전국 추세 Delta) + (지역 편차 Alpha) + (공변량 효과)
  for (t in 1:TT) {
    int start = 1 + (t - 1) * P;
    int end_  = t * P;
    
    matrix[P, K] Xt = X[start:end_, ];
    matrix[P, L] Zt = Z[start:end_, ];

    Y_logit[start:end_] ~ normal(
      rep_vector(delta[t], P) + alpha[, t] + Xt * beta + Zt * gamma,
      sigma_epsilon
    );
  }
}

generated quantities {
  // --- Prediction for TT+1 ---
  vector[P] alpha_TT1;
  real delta_TT1;
  
  vector[P] y_pred_logit;
  vector[P] y_pred_share; // 지역별 예측 득표율

  // --- Post-Stratification (사후 층화) ---
  // 모델이 추정한 지역별 결과를 합쳐서 '진짜 전국 지지율'을 계산
  vector[TT] national_share_est;      // 과거 추정치
  real national_share_pred;           // 미래 예측치

  // 1. Forecast Delta (전국 추세)
  delta_TT1 = delta[TT] + sigma_delta * normal_rng(0, 1);

  // 2. Forecast Alpha (지역별 편차)
  {
     vector[P] z_alpha_pred;
     for (p in 1:P) z_alpha_pred[p] = normal_rng(0, 1);
     // 각 지역은 자신의 과거 편차에서 출발해 랜덤하게 움직임
     alpha_TT1 = alpha[, TT] + sigma_alpha * z_alpha_pred;
  }

  // 3. Calculate Province Predictions (Logit & Share)
  y_pred_logit = rep_vector(delta_TT1, P) + alpha_TT1 + X_pred * beta + Z_pred * gamma;
  for (p in 1:P) {
    y_pred_share[p] = inv_logit(y_pred_logit[p]);
  }

  // 4. Calculate National Share (History)
  // 과거 시점(t=1~TT)에 대해서도 모델이 생각하는 전국 지지율을 복원
  for (t in 1:TT) {
    vector[P] current_share;
    int start = 1 + (t - 1) * P;
    int end_  = t * P;
    
    // (Delta + Alpha + Xb + Zg) -> inv_logit
    vector[P] linear_pred = rep_vector(delta[t], P) + alpha[, t] 
                            + X[start:end_, ] * beta + Z[start:end_, ] * gamma;
    
    for(p in 1:P) current_share[p] = inv_logit(linear_pred[p]);
    
    // 인구 가중 평균 (Post-stratification)
    national_share_est[t] = dot_product(Pop_weight[, t], current_share) / sum(Pop_weight[, t]);
  }

  // 5. Calculate National Share (Prediction)
  national_share_pred = dot_product(Pop_weight[, TT + 1], y_pred_share) / sum(Pop_weight[, TT + 1]);
}
