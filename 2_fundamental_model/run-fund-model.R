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
mod <- cmdstan_model("2_fundamental_model/fund-model-simple.stan")

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
