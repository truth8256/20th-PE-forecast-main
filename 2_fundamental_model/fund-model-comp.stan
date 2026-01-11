data {
  int<lower=2> P;                 // provinces
  int<lower=1> TT;                // time points
  int<lower=1> N;                 // N = P * TT
  
  // [핵심 수정] 압축된 차원 수 2개를 받습니다 (K, L 대신)
  int<lower=1> M_K;               // X의 압축 차원 (예: 4)
  int<lower=1> M_L;               // Z의 압축 차원 (예: 3)
  
  // 행렬의 크기도 압축 차원(M_K, M_L)을 따릅니다
  matrix[N, M_K] X;               // 압축된 X
  matrix[N, M_L] Z;               // 압축된 Z
  
  matrix[P, M_K] X_pred;          // 압축된 예측용 X
  matrix[P, M_L] Z_pred;          // 압축된 예측용 Z
  
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

  // 압축된 변수의 계수 (Theta)
  vector[M_K] theta_beta;  
  vector[M_L] theta_gamma; 
  
  real<lower=0> sigma_epsilon; 
}

transformed parameters {
  matrix[P, TT + 1] alpha;
  
  alpha[, 1] = Alpha_init + sigma_alpha * z_alpha[, 1];
  for (t in 2:(TT + 1)) {
    alpha[, t] = rho * alpha[, t - 1] + sigma_alpha * z_alpha[, t];
  }
}

model {
  // Priors
  to_vector(z_alpha) ~ std_normal();
  rho ~ beta(8, 2); 
  sigma_alpha ~ exponential(5); 
  sigma_epsilon ~ exponential(3); 

  // 계수 Prior (압축변수이므로 normal(0,1))
  theta_beta ~ normal(0, 1);
  theta_gamma ~ normal(0, 1);
  
  // Soft Constraint
  for (t in 1:(TT + 1)) {
    sum(Pop_weight[, t] .* alpha[, t]) ~ normal(0, 0.001 * P); 
  }

  // Likelihood
  for (t in 1:TT) {
    int start = 1 + (t - 1) * P;
    int end_  = t * P;
    
    Y_logit[start:end_] ~ normal(
      alpha[, t] 
      + X[start:end_] * theta_beta    // X(압축) * beta(압축)
      + Z[start:end_] * theta_gamma,  // Z(압축) * gamma(압축)
      sigma_epsilon
    );
  }
}

generated quantities {
  vector[P] y_pred_logit; 
  vector[P] y_pred_share; 

  for (p in 1:P) {
    real mu = alpha[p, TT + 1] 
              + dot_product(X_pred[p, ], theta_beta) 
              + dot_product(Z_pred[p, ], theta_gamma);
    
    y_pred_logit[p] = normal_rng(mu, sigma_epsilon);
    y_pred_share[p] = inv_logit(y_pred_logit[p]);
  }
}

