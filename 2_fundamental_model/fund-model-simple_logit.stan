data {
  int P; // number of provinces
  int TT; // number of times (not including the prediction)
  int N; // number of data = P*T
  int K; // number of national-level coeffs
  int L; // number of province-level coeffs
  matrix[N,K] X; // data matrix of national effects
  matrix[N,L] Z; // data matrix of province effects
  matrix[P,K] X_pred;
  matrix[P,L] Z_pred;
  vector<lower=0, upper=1>[N] Y;
  vector[P] Alpha_init;
  matrix[P,TT+1] Pop_weight; // weight of each province (turnout rate in each election)
}

transformed data {
  vector[N] y_logit;
  for (n in 1:N) {
    // 0이나 1이 있으면 logit가 -inf/inf가 되므로 아주 살짝 조정
    real y_adj = fmin(fmax(Y[n], 1e-6), 1 - 1e-6);
    y_logit[n] = logit(y_adj);
  }
}

parameters {
  matrix[P-1,TT+1] alpha0;
  vector[K] beta;
  vector[L] gamma;
  vector[TT+1] delta;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;
  real<lower=0> sigma_delta;
  array[TT] real<lower=0> sigma_epsilon;
  real<lower=0> sigma;
  vector[P] eta_pred;  // '로그오즈' 스케일의 예측값을 위한 파라미터
}

transformed parameters {
  matrix[P,TT+1] alpha;
  vector[P] y_pred;   // 최종적으로 0~1 스케일의 예측 지지율

  // weighted mean to .5 constraint
  alpha[1:(P-1),] = alpha0;
  for (t in 1:(TT+1)) {
    alpha[P,t] = (0.5 - sum(Pop_weight[1:(P-1),t] .* alpha0[,t])) / Pop_weight[P,t];
  }

  // 로그오즈(eta_pred)를 inv_logit으로 되돌려서 0~1 예측 지지율로 변환
  y_pred = inv_logit(eta_pred);
}

model {
  for(t in 1:(TT+1)){
    if(t==1){
      alpha[,1] ~ normal(Alpha_init,sigma_alpha);
    }else{
      alpha[,t] ~ normal(alpha[,t-1],sigma_alpha);
    } 
  }
  // observation equation on logit scale
  for (t in 1:TT) {
    int start = 1 + (t-1)*P;
    int end_  = t*P;

    y_logit[start:end_] ~ normal(
      alpha[,t]
      + (X * beta)[start:end_]
      + (Z * gamma)[start:end_]
      + delta[t],
      sigma_epsilon[t]
    );
  }

  // prediction on logit scale
  eta_pred ~ normal(
    alpha[,TT+1]
    + X_pred * beta
    + Z_pred * gamma
    + delta[TT+1],
    sigma_epsilon[TT]
  );

  delta ~ normal(0, sigma_delta);
  // beta ~ normal(0, sigma_beta);
  // gamma ~ normal(0, sigma_gamma);
  beta ~ cauchy(0, sigma_beta);
  gamma ~ cauchy(0, sigma_gamma);

  sigma_alpha ~ normal(0, sigma);
  sigma_delta ~ normal(0, sigma);
  sigma_beta ~ normal(0, sigma);
  sigma_gamma ~ normal(0, sigma);
  sigma_epsilon ~ normal(0, sigma);
}

