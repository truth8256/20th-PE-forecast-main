data {
  int<lower=2> P;
  int<lower=1> TT;
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> L;

  matrix[N,K] X;
  matrix[N,L] Z;
  matrix[P,K] X_pred;
  matrix[P,L] Z_pred;

  vector<lower=0, upper=1>[N] Y;
  vector[P] Alpha_init;
  matrix[P,TT+1] Pop_weight;
}

parameters {
  matrix[P-1,TT+1] alpha0;

  vector[K] beta;
  vector[L] gamma;
  vector[TT+1] delta;

  real<lower=-1, upper=1> rho;      // ★ AR(1) 핵심 파라미터
  real<lower=0> sigma_alpha;

  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;
  real<lower=0> sigma_delta;
  vector<lower=0>[TT] sigma_epsilon;

  real<lower=0> sigma;

  vector[P] y_pred;
}

transformed parameters {
  matrix[P,TT+1] alpha;

  alpha[1:(P-1),] = alpha0;
  for (t in 1:(TT+1)) {
    alpha[P,t] =
      (0.5 - dot_product(Pop_weight[1:(P-1),t], alpha0[,t]))
      / Pop_weight[P,t];
  }
}

model {
  // --- AR(1) state evolution ---
  alpha0[,1] ~ normal(Alpha_init[1:(P-1)], sigma_alpha);

  for (t in 2:(TT+1)) {
    alpha0[,t] ~ normal(rho * alpha0[,t-1], sigma_alpha);
  }

  // --- likelihood (fit only) ---
  for (t in 1:TT) {
    int start = 1 + (t-1)*P;
    int end_  = t*P;

    Y[start:end_] ~ normal(
      alpha[,t]
      + (X * beta)[start:end_]
      + (Z * gamma)[start:end_]
      + delta[t],
      sigma_epsilon[t]
    );
  }

  // --- prediction ---
  y_pred ~ normal(
    alpha[,TT+1]
    + X_pred * beta
    + Z_pred * gamma
    + delta[TT+1],
    sigma_epsilon[TT]
  );

  // --- priors ---
  rho ~ normal(0, 0.5);              // ★ AR 안정성

  beta  ~ cauchy(0, sigma_beta);
  gamma ~ cauchy(0, sigma_gamma);
  delta ~ normal(0, sigma_delta);

  sigma_alpha  ~ normal(0, sigma);
  sigma_beta   ~ normal(0, sigma);
  sigma_gamma  ~ normal(0, sigma);
  sigma_delta  ~ normal(0, sigma);
  sigma_epsilon ~ normal(0, sigma);
}
