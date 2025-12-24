# setwd("D:\\OneDrive - 에스티아이\\=대학원\\=논문연구\\20th-PE-forecast-main")

# ─────────────────────────────────────
# 📦 준비
# ─────────────────────────────────────
rm(list=ls())
load("pe.RData")

library(cmdstanr)
library(posterior)
library(ggplot2)

standardize <- function(mat){
  mat <- t(mat) - colMeans(mat)
  t(mat / sqrt(rowSums(mat^2)))
}

# ─────────────────────────────────────
# 📊 데이터 구성
# ─────────────────────────────────────
P = length(cluster)
TT = length(15:19)
N = P * TT

Pop_weight <- matrix(unlist(lapply(pe_sum, \(l) rowSums(l[, 1:(ncol(l-2))]))), nrow = P)
Pop_weight <- t(t(Pop_weight) / colSums(Pop_weight))
dimnames(Pop_weight) <- list(names(cluster), paste0("pe", 14:19))

Alpha_init <- pe_twoparty_vote_share[1, -(1+P)] -
  sum(Pop_weight[,1] * pe_twoparty_vote_share[1, -(1+P)]) + 0.5

data_list <- list(
  P = P,
  TT = TT,
  N = N,
  K = ncol(X),
  L = ncol(Z),
  X = standardize(X)[1:N,,drop=FALSE],
  Z = standardize(Z)[1:N,,drop=FALSE],
  X_pred = standardize(X)[(N+1):(N+P),,drop=FALSE],
  Z_pred = standardize(Z)[(N+1):(N+P),,drop=FALSE],
  Y = c(t(pe_twoparty_vote_share[2:(1+TT), -(1+P)])),
  Alpha_init = Alpha_init,
  Pop_weight = Pop_weight[, c(2:6, 6)]
)

# ─────────────────────────────────────
# 🧠 모델 컴파일 및 샘플링
# ─────────────────────────────────────
mod <- cmdstan_model("2_fundamental_model/fund-model-simple_logit.stan")

fit <- mod$sample(
  data = data_list,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 8000,
  thin = 1
)

# ─────────────────────────────────────
# 📈 결과 요약
# ─────────────────────────────────────
fit$summary() |> View()

# posterior package로 변환
draws <- as_draws_df(fit)
head(draws)

# 저장
fit$save_output_files(dir = "2_fundamental_model")

# 특정 파라미터 추출 (예: beta, gamma, y_pred)
beta_draws <- fit$draws("beta")
gamma_draws <- fit$draws("gamma")
y_pred_draws <- fit$draws("y_pred")

# 저장 (필요 시)
save(beta_draws, gamma_draws, y_pred_draws, file = "samples-fundamental.RData")



# ─────────────────────────────────────
# 🧠 (2) 로짓 스케일 모형: 컴파일 및 샘플링
# ─────────────────────────────────────
mod_logit <- cmdstan_model("2_fundamental_model/fund-model-simple_logit.stan")
# ↑ 여기에 바로 전 메시지에서 만든 logit 버전 Stan 코드가 들어있는 파일

fit_logit <- mod_logit$sample(
  data = data_list,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 8000,
  thin = 1
)

# 결과 요약
fit_logit$summary() |> View()

# posterior package로 변환
draws_logit <- as_draws_df(fit_logit)
head(draws_logit)

# output 파일 저장
fit_logit$save_output_files(dir = "2_fundamental_model")

# 특정 파라미터 추출
beta_draws_logit   <- fit_logit$draws("beta")
gamma_draws_logit  <- fit_logit$draws("gamma")
y_pred_draws_logit <- fit_logit$draws("y_pred")  # 이미 0~1 스케일

# 저장
save(
  beta_draws_logit, gamma_draws_logit, y_pred_draws_logit,
  file = "samples-fundamental-logit.RData"
)

