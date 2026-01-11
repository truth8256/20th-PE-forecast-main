# 설정
n_ensemble <- 20
M_K <- 4  # X(12개) -> 4개로 압축
M_L <- 2  # Z(6개라고 가정) -> 2개로 압축

K <- ncol(X_scaled)
L <- ncol(Z_scaled)

pred_list <- list()
beta_list <- list()  # 복원된 X 계수 저장
gamma_list <- list() # 복원된 Z 계수 저장

mod_comp <- cmdstan_model("2_fundamental_model/fund-model-comp.stan")

for (i in 1:n_ensemble) {
  cat(sprintf("\n[Ensemble %d / %d] Compressing both X and Z...\n", i, n_ensemble))
  
  # (1) 무작위 투영 행렬 생성 (Phi_X, Phi_Z)
  Phi_X <- matrix(rnorm(K * M_K, mean=0, sd=1/sqrt(M_K)), nrow=K, ncol=M_K)
  Phi_Z <- matrix(rnorm(L * M_L, mean=0, sd=1/sqrt(M_L)), nrow=L, ncol=M_L)
  
  # (2) 데이터 압축
  X_comp <- X_scaled %*% Phi_X
  X_pred_comp <- X_scaled[(N+1):(N+P),,drop=FALSE] %*% Phi_X
  
  Z_comp <- Z_scaled %*% Phi_Z      # Z 압축
  Z_pred_comp <- Z_scaled[(N+1):(N+P),,drop=FALSE] %*% Phi_Z
  
  # (3) Stan 데이터 준비
  data_list_comp <- list(
    P = P, TT = TT, N = N,
    Y = Y, Pop_weight = Pop_weight,
    Alpha_init = Alpha_init <- make_alpha_init(y0, Pop_weight[,1], "logit"),
    # 수정된 부분: 차원수와 데이터 모두 압축된 것으로 교체
    M_K = M_K, 
    M_L = M_L,
    X = X_comp[1:N,], 
    Z = Z_comp[1:N,],
    X_pred = X_pred_comp, 
    Z_pred = Z_pred_comp
  )

  # 초기값 생성 함수 정의
  init_fun <- function() {
    list(
      # 1. 스케일 파라미터 (0이 아닌 안전한 값)
      sigma_epsilon = 0.05,
      sigma_alpha = 0.05,
      
      # 2. 자기상관 계수
      rho = 0.5,
      
      # 3. 회귀 계수 (압축 차원 크기만큼 0으로 채움)
      theta_beta = rep(0, data_list_comp$M_K),
      theta_gamma = rep(0, data_list_comp$M_L),
      
      # 4. [추가됨] z_alpha (표준정규분포의 잠재변수)
      # 크기: P행 x (TT+1)열
      # 0으로 초기화하면 "튀는 값" 없이 평균에서 시작하므로 가장 안전합니다.
      z_alpha = matrix(0, nrow = data_list_comp$P, ncol = data_list_comp$TT + 1)
    )
  }
  
  # (4) Stan 실행 (mod_comp 모델 필요 - 아래 3번 참조)
  fit_comp <- mod_comp$sample(
    data = data_list_comp,init = init_fun,
    chains = 4, parallel_chains = 4,
    iter_warmup = 500, iter_sampling = 1000,
    refresh = 0
  )
  
  # (5) 예측값 저장
  pred_list[[i]] <- fit_comp$summary(variables = "y_pred_share")$mean
  
  # (6) 계수 복원 (Back-projection)
  # Beta 복원
  theta_beta <- fit_comp$draws("theta_beta", format = "matrix") 
  beta_restored <- theta_beta %*% t(Phi_X) 
  beta_list[[i]] <- colMeans(beta_restored)
  
  # Gamma 복원 (Z의 계수)
  theta_gamma <- fit_comp$draws("theta_gamma", format = "matrix") 
  gamma_restored <- theta_gamma %*% t(Phi_Z)
  gamma_list[[i]] <- colMeans(gamma_restored)
}

# 결과 종합
final_beta <- colMeans(do.call(rbind, beta_list))
final_gamma <- colMeans(do.call(rbind, gamma_list))

print("=== 복원된 Beta (X) ===")
print(round(final_beta, 4))
print("=== 복원된 Gamma (Z) ===")
print(round(final_gamma, 4))



# 리스트를 행렬로 변환 (행: 20번의 시도, 열: 16개 지역)
pred_matrix <- do.call(rbind, pred_list) 

# (선택) 지역 이름 붙이기 (보기 편하게)
colnames(pred_matrix) <- names(cluster)

# 실제 21대 득표율 (16개 지역 벡터)
# 만약 백분율(0~100)이 아니라 비율(0~1)이라면 스케일 맞춰주세요.
actual_y_21 <- pe21_sum[, 1] / rowSums(pe21_sum[, 1:2])

# 20번 시도의 평균값 계산
final_prediction <- colMeans(pred_matrix)

# 1. RMSE (평균 제곱근 오차) - 낮을수록 좋음
rmse <- sqrt(mean((final_prediction - actual_y_21)^2))

# 2. MAE (평균 절대 오차) - 낮을수록 좋음
mae <- mean(abs(final_prediction - actual_y_21))

# 3. 상관계수 (경향성) - 1에 가까울수록 좋음
corr <- cor(final_prediction, actual_y_21)

cat(sprintf("=== 최종 예측 성적표 ===\nRMSE: %.4f\nMAE : %.4f\nCORR: %.4f\n", rmse, mae, corr))

# 오차의 방향성 확인 (Bias)
# 양수(+)면 과대평가, 음수(-)면 과소평가
bias <- mean(final_prediction - actual_y_21)

cat(sprintf("평균 편향(Bias): %.4f (%%p로 환산하면 약 %.2f%%p)\n", bias, bias * 100))

# 3.3%p 오차의 원인이 '한쪽으로 쏠림' 때문인지 확인
if (bias < 0) {
  cat("해석: 모델이 실제보다 득표율을 '낮게' 예측했습니다. (야당/민주당 바람을 덜 반영함)")
} else {
  cat("해석: 모델이 실제보다 득표율을 '높게' 예측했습니다.")
}
