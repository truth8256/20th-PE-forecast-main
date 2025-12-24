# setwd("D:\\OneDrive - 에스티아이\\=대학원\\=논문연구\\20th-PE-forecast-main")

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
# 📊 기본 설정
# ─────────────────────────────────────
P  <- length(cluster)
TT <- length(15:19)
N  <- P * TT

region_order <- names(cluster)

# 민주당 후보가 후보열에서 몇 번째인지 (14~20)
democ_index <- c(2, 2, 2, 1, 2, 1, 1)  # pe14..pe20
names(democ_index) <- paste0("pe", 14:20)

# fitting에 쓰는 선거(15~19)
pe_names_fit <- paste0("pe", 15:19)

# ─────────────────────────────────────
# ✅ Binomial counts 생성 (민주/보수 열을 democ_index로 결정)
#    가정: 양당 후보가 항상 1~2열에 위치
# ─────────────────────────────────────
get_two_party_counts <- function(pe_df, pe_name, region_order, democ_index_vec){
  df <- pe_df[region_order, , drop=FALSE]
  
  dem_col <- democ_index_vec[pe_name]
  if (is.na(dem_col)) stop("democ_index에 ", pe_name, " 값이 없습니다.")
  
  # 보수 열은 (1,2) 중 민주가 아닌 열
  con_col <- setdiff(1:2, dem_col)
  if (length(con_col) != 1) stop(pe_name, ": 보수 열을 결정할 수 없습니다. (1~2열 가정 점검 필요)")
  
  dem <- as.integer(df[, dem_col])
  con <- as.integer(df[, con_col])
  
  list(Y = dem, N = dem + con)
}

Y_mat <- matrix(NA_integer_, nrow=P, ncol=TT)
N_mat <- matrix(NA_integer_, nrow=P, ncol=TT)

for(t in 1:TT){
  pe_name <- pe_names_fit[t]
  obj <- get_two_party_counts(pe_sum[[pe_name]], pe_name, region_order, democ_index)
  Y_mat[, t] <- obj$Y
  N_mat[, t] <- obj$N
}

Y_count <- as.integer(c(t(Y_mat)))
N_count <- as.integer(c(t(N_mat)))

stopifnot(length(Y_count) == N, length(N_count) == N)
stopifnot(all(N_count > 0), all(Y_count >= 0), all(Y_count <= N_count))

# ─────────────────────────────────────
# 📦 공변량 표준화 및 data_list
# ─────────────────────────────────────
X_std <- standardize(X)
Z_std <- standardize(Z)

data_list <- list(
  P = P,
  TT = TT,
  N = N,
  K = ncol(X),
  L = ncol(Z),
  X = X_std[1:N, , drop=FALSE],
  Z = Z_std[1:N, , drop=FALSE],
  X_pred = X_std[(N+1):(N+P), , drop=FALSE],
  Z_pred = Z_std[(N+1):(N+P), , drop=FALSE],
  Y_count = Y_count,
  N_count = N_count
)

# ─────────────────────────────────────
# 🧠 모델 컴파일 및 샘플링
# ─────────────────────────────────────
mod <- cmdstan_model("2_fundamental_model/fund-model-simple_GLMM.stan")

fit <- mod$sample(
  data = data_list,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 8000,
  thin = 1,
  adapt_delta = 0.95,
  max_treedepth = 12
)

# ─────────────────────────────────────
# 📈 결과 요약/저장
# ─────────────────────────────────────
fit$summary() |> View()

fit$save_output_files(dir = "2_fundamental_model")

beta_draws    <- fit$draws("beta")
gamma_draws   <- fit$draws("gamma")
delta_draws   <- fit$draws("delta")
pi_pred_draws <- fit$draws("pi_pred")

save(beta_draws, gamma_draws, delta_draws, pi_pred_draws,
     file = "samples-fundamental-dglmm.RData")
