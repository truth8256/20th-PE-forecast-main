library(posterior)
library(purrr)
library(dplyr)
library(tidyr)
library(readr)

# 📌 실제 득표율(지역별, 전국) 데이터 필요
# 지역별 실제 득표율: actual_y_vec (벡터 길이 P x T)
# 전국 실제 득표율: actual_y_nat (벡터 길이 T)
# Pop_weight: P x T 행렬 (지역 가중치)

evaluate_model_final <- function(
    fit,
    actual_y_region,   # length P (최종 시점 지역 실제값)
    actual_y_nat,      # scalar (최종 시점 전국 실제값)
    Pop_weight,        # length P vector OR P×TT matrix (최종 시점 가중치)
    TT = NULL,         # (선택) y_pred가 시간축을 포함할 때 최종시점 인덱스용
    P  = NULL          # (선택) 지역 수(미지정 시 actual_y_region로 추정)
) {
  library(posterior)
  library(tibble)
  
  ## -------------------------------
  ## 0. 기본 체크/정리
  ## -------------------------------
  actual_y_region <- as.numeric(actual_y_region)
  if (is.null(P)) P <- length(actual_y_region)
  if (length(actual_y_region) != P) stop("actual_y_region 길이가 P와 일치하지 않습니다.")
  if (length(actual_y_nat) != 1) stop("actual_y_nat는 스칼라여야 합니다.")
  
  # Pop_weight: vector(P) 또는 matrix(P×TT)
  if (is.matrix(Pop_weight)) {
    if (nrow(Pop_weight) != P) stop("Pop_weight 행 수가 P와 다릅니다.")
    if (is.null(TT)) TT <- ncol(Pop_weight)
    if (ncol(Pop_weight) < TT) stop("Pop_weight의 열 수가 TT보다 작습니다.")
    w_last <- as.numeric(Pop_weight[, TT])
  } else {
    w_last <- as.numeric(Pop_weight)
    if (length(w_last) != P) stop("Pop_weight 벡터 길이가 P와 다릅니다.")
  }
  
  # (권장) 가중치 합이 1인지 점검 (필수는 아님)
  if (abs(sum(w_last) - 1) > 1e-6) {
    warning("Pop_weight의 최종시점 가중치 합이 1이 아닙니다. (sum = ", signif(sum(w_last), 6), ")")
  }
  
  ## -------------------------------
  ## 1. posterior draws 가져오기 + 예측변수 자동 탐지
  ## -------------------------------
  draws <- as_draws_df(fit$draws())
  draw_vars <- names(draws)
  
  if (any(grepl("^y_pred_share\\[", draw_vars))) {
    var_type <- "y_pred_share"
  } else if (any(grepl("^y_pred\\[", draw_vars))) {
    var_type <- "y_pred"
  } else {
    stop("예측 변수(y_pred_share 또는 y_pred)가 존재하지 않습니다.")
  }
  
  ## -------------------------------
  ## 2. 최종 시점 지역 예측 draw 행렬 pred_mat (S×P) 구성
  ##    - 지원 형태:
  ##      (a) var[p]             : vector[P] (이미 최종시점만 저장된 경우)
  ##      (b) var[p, t]          : matrix[P,TT]
  ##      (c) var[idx] (idx=1..P*TT) : vector[P*TT]로 저장된 경우 (TT 필요)
  ## -------------------------------
  S <- nrow(draws)
  pred_mat <- matrix(NA_real_, nrow = S, ncol = P)
  
  has_2d <- any(grepl(paste0("^", var_type, "\\[\\d+,\\d+\\]$"), draw_vars))
  has_1d <- any(grepl(paste0("^", var_type, "\\[\\d+\\]$"), draw_vars))
  
  if (has_2d) {
    if (is.null(TT)) stop("예측변수가 [p,t] 형태인데 TT가 NULL입니다. TT(최종 시점 인덱스)를 지정해주세요.")
    for (p in 1:P) {
      vn <- paste0(var_type, "[", p, ",", TT, "]")
      if (!vn %in% draw_vars) stop("예측 변수 없음: ", vn)
      pred_mat[, p] <- draws[[vn]]
    }
  } else if (has_1d) {
    # 1차원일 때: P개만 있으면 지역만, P*TT개면 마지막 블록을 사용
    one_d_names <- grep(paste0("^", var_type, "\\[\\d+\\]$"), draw_vars, value = TRUE)
    n_1d <- length(one_d_names)
    
    if (n_1d == P) {
      # var[p] : 최종시점 지역값만 들어있다고 해석
      for (p in 1:P) {
        vn <- paste0(var_type, "[", p, "]")
        if (!vn %in% draw_vars) stop("예측 변수 없음: ", vn)
        pred_mat[, p] <- draws[[vn]]
      }
    } else {
      # var[1..P*TT]로 들어있는 경우로 해석
      if (is.null(TT)) stop("예측변수가 1차원이고 P개가 아닙니다. P*TT 형태로 보이며, TT가 필요합니다.")
      if (n_1d < P * TT) stop("예측변수 개수가 P*TT보다 작습니다. (", n_1d, " < ", P * TT, ")")
      
      idx_last <- ((TT - 1) * P + 1):(TT * P)
      for (j in 1:P) {
        vn <- paste0(var_type, "[", idx_last[j], "]")
        if (!vn %in% draw_vars) stop("예측 변수 없음: ", vn)
        pred_mat[, j] <- draws[[vn]]
      }
    }
  } else {
    stop("예측 변수 차원을 해석할 수 없습니다. y_pred_share[...] 또는 y_pred[...]의 이름 규칙을 확인해주세요.")
  }
  
  ## -------------------------------
  ## 3. 지역 메트릭 (최종 시점 단면 P개)
  ##    - 점추정 기반: RMSE/MAE/Bias (posterior mean vs actual)
  ##    - 분포 기반: Coverage95/MeanWidth (95% CrI)
  ## -------------------------------
  pred_mean_region <- colMeans(pred_mat)
  err_region <- pred_mean_region - actual_y_region
  
  Region_RMSE <- sqrt(mean(err_region^2))
  Region_MAE  <- mean(abs(err_region))
  Region_Bias <- mean(err_region)
  
  ci_low_region  <- apply(pred_mat, 2, quantile, 0.025)
  ci_high_region <- apply(pred_mat, 2, quantile, 0.975)
  
  Region_Coverage95 <- mean(actual_y_region >= ci_low_region & actual_y_region <= ci_high_region)
  Region_MeanWidth  <- mean(ci_high_region - ci_low_region)
  
  ## -------------------------------
  ## 4. 전국 메트릭 (최종 시점 1개)
  ##    - 점추정 기반: National_AbsError (= |posterior mean - actual|)
  ##    - 분포(의사결정) 기반: National_PERMSE
  ##      = sqrt( mean_s ( (y_nat_draw_s - actual)^2 ) )
  ##    - 불확실성: National_CIWidth, National_Covered (0/1)
  ## -------------------------------
  y_hat_nat_draws <- as.numeric(pred_mat %*% w_last)  # length S
  
  y_hat_nat_mean  <- mean(y_hat_nat_draws)
  National_AbsError <- abs(y_hat_nat_mean - actual_y_nat)
  National_Bias     <- (y_hat_nat_mean - actual_y_nat)
  
  # Posterior Expected RMSE (PERMSE)
  National_PERMSE <- sqrt(mean((y_hat_nat_draws - actual_y_nat)^2))
  
  ci_nat_low  <- as.numeric(quantile(y_hat_nat_draws, 0.025))
  ci_nat_high <- as.numeric(quantile(y_hat_nat_draws, 0.975))
  National_CIWidth <- ci_nat_high - ci_nat_low
  National_Covered <- as.numeric(actual_y_nat >= ci_nat_low && actual_y_nat <= ci_nat_high)
  
  ## -------------------------------
  ## 5. 진단: Max R-hat
  ## -------------------------------
  summ <- fit$summary()
  max_rhat <- max(summ$rhat, na.rm = TRUE)
  
  ## -------------------------------
  ## 6. 반환
  ## -------------------------------
  tibble(
    Region_RMSE       = Region_RMSE,
    Region_MAE        = Region_MAE,
    Region_Bias       = Region_Bias,
    Region_Coverage95 = Region_Coverage95,
    Region_MeanWidth  = Region_MeanWidth,
    Region_MaxRhat    = max_rhat,
    
    National_AbsError = National_AbsError,
    National_Bias     = National_Bias,
    National_PERMSE   = National_PERMSE,
    National_CIWidth  = National_CIWidth,
    National_Covered  = National_Covered,
    National_MaxRhat  = max_rhat
  )
}


# ============================================================
# 0) (전제) evaluate_model_final() 함수가 이미 로드되어 있어야 합니다.
#    - 위에서 드린 evaluate_model_final 그대로 사용
# ============================================================


# ============================================================
# 1) 여러 모델을 한꺼번에 평가/비교하는 함수
# ============================================================
compare_models_final <- function(
    model_tags,
    out_dir,
    actual_y_region,
    actual_y_nat,
    Pop_weight,
    TT = NULL,
    P  = NULL,
    file_prefix = "compare_table"
) {
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  
  # P 자동 추정
  if (is.null(P)) P <- length(actual_y_region)
  
  # fit 파일 존재 여부 체크
  fit_paths <- file.path(out_dir, paste0("fit_", model_tags, ".RDS"))
  missing_idx <- which(!file.exists(fit_paths))
  if (length(missing_idx) > 0) {
    stop(
      "아래 모델의 fit RDS 파일이 없습니다:\n",
      paste0(" - ", model_tags[missing_idx], " (", fit_paths[missing_idx], ")", collapse = "\n")
    )
  }
  
  # 모델별 평가
  res_list <- map2(
    model_tags, fit_paths,
    ~{
      fit <- readRDS(.y)
      
      metrics <- evaluate_model_final(
        fit = fit,
        actual_y_region = actual_y_region,
        actual_y_nat    = actual_y_nat,
        Pop_weight      = Pop_weight,
        TT              = TT,
        P               = P
      )
      
      metrics %>%
        mutate(Model = .x, .before = 1)
    }
  )
  
  compare_table <- bind_rows(res_list)
  
  # 보기 좋게 정렬 (원하시면 기준 바꾸세요)
  compare_table_sorted <- compare_table %>%
    arrange(Region_RMSE, National_PERMSE, National_AbsError)
  
  # 저장(선택): CSV + RDS
  write_csv(compare_table_sorted, file.path(out_dir, paste0(file_prefix, ".csv")))
  saveRDS(compare_table_sorted, file.path(out_dir, paste0(file_prefix, ".RDS")))
  
  compare_table_sorted
}


# ============================================================
# 2) 사용 예시
# ============================================================
model_tags <- c(
  "m1_baseline3_0109",  "m2_logit5_0109",  "m3_beta1_0109",  "m4_logit_ar1_3_0109",
  "m4_logit_ar1_4_0109",  "m2_logit6_0109",  "m3_beta2_0109",  "m3_beta3_0109",
  "m4_logit_ar1_5_0109" 
)

out_dir <- file.path("stan_outputs", "pe20")
out_dir <- file.path("stan_outputs", "pe21")
# out_dir: run_and_save에서 저장한 디렉토리와 동일해야 함
# 예: out_dir <- "stan_outputs/model_compare"
# actual_y_region: 예) actual_y_21 (길이 16)
# actual_y_nat:    예) national_share_21 (스칼라)
# Pop_weight:      vector(16) 또는 matrix(16×TT)
# TT: y_pred가 시간축을 포함할 때만 필요 (예: TT=마지막시점 인덱스)
# P:  지역 수 (예: 16)

compare_table <- compare_models_final(
  model_tags       = model_tags,
  out_dir          = out_dir,
  actual_y_region  = actual_y_21,
  actual_y_nat     = national_share_21,
  Pop_weight       = Pop_weight,
  TT               = TT,     # 필요 없으면 NULL로 두세요
  P                = 16,
  file_prefix      = "compare_table_final_0109"
)

print(compare_table, width = Inf)


# ============================================================
# 3) (선택) “요약_최상위 N개”만 뽑아 보기
# ============================================================
top5 <- compare_table %>%
  dplyr::arrange(Region_RMSE) %>%
  dplyr::slice(1:5)

print(top5, width = Inf)

