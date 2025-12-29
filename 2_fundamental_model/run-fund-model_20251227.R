# setwd("D:\\OneDrive - 에스티아이\\=대학원\\=논문연구\\20th-PE-forecast-main")

# ─────────────────────────────────────
# 📦 준비
# ─────────────────────────────────────
# rm(list=ls())
# load("pe.RData")

library(cmdstanr)
library(posterior)
library(dplyr)

# (선택) 결과 저장 폴더
dir.create("stan_outputs", showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────
# ✅ DATA SWITCH (pe20 / pe21) : 아래 둘 중 하나만 주석 해제
# ─────────────────────────────────────────────────────────────

# [A] 20대 예측: fit pe15~pe19 (TT=5), pred pe20
{
  target_label <- "pe20"

  X <- X_20; Z <- Z_20
  pe_sum <- pe_sum_20
  pe_twoparty_vote_share <- pe_twoparty_vote_share_20
  democ_index <- democ_index_20
  Pop_weight <- Pop_weight_20

  pe_names_fit <- paste0("pe", 15:19)
  pe_pred_name <- "pe20"

  cat("\n▶ DATA SWITCH: 20th PE prediction (fit pe15–pe19 → pred pe20)\n")
}

# 
# # [B] 21대 예측: fit pe15~pe20 (TT=6), pred pe21
# {
#   target_label <- "pe21"
# 
#   X <- X_21; Z <- Z_21
#   pe_sum <- pe_sum_21
#   pe_twoparty_vote_share <- pe_twoparty_vote_share_21
#   democ_index <- democ_index_21
#   Pop_weight <- Pop_weight_21
# 
#   pe_names_fit <- paste0("pe", 15:20)
#   pe_pred_name <- "pe21"
# 
#   cat("\n▶ DATA SWITCH: 21st PE prediction (fit pe15–pe20 → pred pe21)\n")
# }

# ─────────────────────────────────────
# ✅ 공통 파생값 + 안전장치
# ─────────────────────────────────────
P  <- length(cluster)
TT <- length(pe_names_fit)
N  <- P * TT

stopifnot(nrow(X) == P * (TT + 1))
stopifnot(nrow(Z) == P * (TT + 1))
stopifnot(all(dim(Pop_weight) == c(P, TT + 1)))

out_dir <- file.path("stan_outputs", target_label)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────
# ✅ Alpha_init 자동 생성 (fit 첫 선거 기준)
#    - alpha[,1] ~ Normal(Alpha_init, sigma_alpha) 형태에서
#      Alpha_init이 "가중평균 0.5 제약"과 정합되도록 구성
# ─────────────────────────────────────
# fit 첫 선거의 지역 득표율 벡터 (길이 P)
y0 <- as.numeric(pe_twoparty_vote_share[1, names(cluster), drop=TRUE])

make_alpha_init <- function(y0, w, scale = c("prob","logit"), eps = 1e-6) {
  scale <- match.arg(scale)
  if (scale == "prob") {
    a <- y0 - sum(w * y0) + 0.5
    stopifnot(abs(sum(w * a) - 0.5) < 1e-10)
    return(a)
  } else {
    y0_adj <- pmin(1 - eps, pmax(eps, y0))
    y0_logit <- qlogis(y0_adj)
    a <- y0_logit - sum(w * y0_logit)
    stopifnot(abs(sum(w * a) - 0.0) < 1e-10)
    return(a)
  }
}

# sanity check
stopifnot(abs(sum(Pop_weight[,1] * Alpha_init) - 0.5) < 1e-10)


out_dir <- file.path("stan_outputs", target_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ─────────────────────────────────────
# ✅ scale: 더미/서열(1,2) 변수 제외하고 연속형만 스케일
# ─────────────────────────────────────
dummy_X <- c("is_current","twice1","twice2","impeach")
dummy_Z <- c("home_dem","homeground_dem","home_con","homeground_con")

X_scaled <- X
Z_scaled <- Z

colX_scale <- setdiff(colnames(X_scaled), dummy_X)
colZ_scale <- setdiff(colnames(Z_scaled), dummy_Z)

X_scaled[, colX_scale] <- scale(X_scaled[, colX_scale, drop=FALSE])
Z_scaled[, colZ_scale] <- scale(Z_scaled[, colZ_scale, drop=FALSE])


# ─────────────────────────────────────
# (1) Baseline용 Y (fit 기간만)
# ─────────────────────────────────────
Y <- as.numeric(c(t(pe_twoparty_vote_share[, names(cluster), drop=FALSE])))

data_list_1 <- list(
  P = P, TT = TT, N = N,
  K = ncol(X_scaled), L = ncol(Z_scaled),
  X = X_scaled[1:N,,drop=FALSE],
  Z = Z_scaled[1:N,,drop=FALSE],
  X_pred = X_scaled[(N+1):(N+P),,drop=FALSE],
  Z_pred = Z_scaled[(N+1):(N+P),,drop=FALSE],
  Y = Y,
  Alpha_init = Alpha_init <- make_alpha_init(y0, Pop_weight[,1], "prob"),
  Pop_weight = Pop_weight
)

data_list_2 <- list(
  P = P, TT = TT, N = N,
  K = ncol(X_scaled), L = ncol(Z_scaled),
  X = X_scaled[1:N,,drop=FALSE],
  Z = Z_scaled[1:N,,drop=FALSE],
  X_pred = X_scaled[(N+1):(N+P),,drop=FALSE],
  Z_pred = Z_scaled[(N+1):(N+P),,drop=FALSE],
  Y = Y,
  Alpha_init = Alpha_init <- make_alpha_init(y0, Pop_weight[,1], "logit"),
  Pop_weight = Pop_weight
)


# # ─────────────────────────────────────
# # (2)(3) GLMM용 counts 생성 (fit 기간만)
# # ─────────────────────────────────────
# Y_mat <- matrix(NA_integer_, nrow=P, ncol=TT)
# N_mat <- matrix(NA_integer_, nrow=P, ncol=TT)
# 
# for(t in 1:TT){
#   pe_name <- pe_names_fit[t]
#   df <- pe_sum[[pe_name]][names(cluster), , drop=FALSE]
#   
#   dem_col <- as.integer(democ_index[pe_name])
#   con_col <- setdiff(1:2, dem_col)
#   
#   dem <- as.integer(df[, dem_col])
#   con <- as.integer(df[, con_col])
#   
#   Y_mat[, t] <- dem
#   N_mat[, t] <- dem + con
# }
# 
# Y_count <- as.integer(c(t(Y_mat)))
# N_count <- as.integer(c(t(N_mat)))
# 
# data_list_23 <- list(
#   P = P, TT = TT, N = N,
#   K = ncol(X_scaled), L = ncol(Z_scaled),
#   X = X_scaled[1:N, , drop = FALSE],
#   Z = Z_scaled[1:N, , drop = FALSE],
#   X_pred = X_scaled[(N+1):(N+P), , drop = FALSE],
#   Z_pred = Z_scaled[(N+1):(N+P), , drop = FALSE],
#   Y_count = Y_count,
#   N_count = N_count,
#   Pop_weight = Pop_weight
# )


# ─────────────────────────────────────
# ✅ 실행저장함수
# ─────────────────────────────────────

run_and_save <- function(stan_file, data_list, model_tag, out_dir,
                         seed=1234, chains=4, parallel_chains=4,
                         iter_warmup=2000, iter_sampling=8000, adapt_delta = 0.98,   
                         max_treedepth = 15, thin=1){
  
  mod <- cmdstan_model(
    stan_file,
    cpp_options = list(STAN_THREADS = FALSE)  # ✅ TBB 링크 회피
  )
  
  fit <- mod$sample(
    data = data_list,
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    thin = thin,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  fit$save_object(file.path(out_dir, paste0("fit_", model_tag, ".RDS")))
  summ <- fit$summary()
  readr::write_csv(summ, file.path(out_dir, paste0("summary_", model_tag, ".csv")))
  
  fit
}



# ─────────────────────────────────────
# 🧠 1,2,3 모델 순차 실행 + 저장
# ─────────────────────────────────────

# # (1) baseline (Kang2024)
# fit1 <- run_and_save(
#   stan_file = "2_fundamental_model/fund-model-simple1.stan",
#   data_list = data_list_1,
#   model_tag = "m1_baseline",
#   out_dir   = out_dir
# )

# (2) logit 
fit2 <- run_and_save(
  stan_file = "2_fundamental_model/fund-model-logit.stan",
  iter_warmup=1000, iter_sampling=1000, 
  data_list = data_list_2,
  model_tag = "m2_logit",
  out_dir   = out_dir
)


cat("\n✅ All models completed. Outputs saved to: ", out_dir, "\n")


