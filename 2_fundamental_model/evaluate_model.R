library(cmdstanr)
library(dplyr)
library(readr)
library(posterior)
library(loo)

# -------------------------------------------------------------------
# 0. 설정 및 실제 데이터 준비
# -------------------------------------------------------------------
# 평가할 모델들의 태그 리스트 (파일 명에 사용된 태그)
# 예: c("m1_baseline", "m2_logit_normal", "m3_hierarchical")

model_tags <- c("m2_logit5_0109","m2_logit6_0109","m3_beta1_delta","m3_beta2_0103","m3_beta3",
                "m4_logit_ar1_3","m4_logit_ar1_4","m4_logit_ar1_5"
                )
out_dir <- "stan_outputs/pe21" # 파일이 저장된 경로

# [중요] 실제 투표 결과 (actual_y_20)
# P개 지역의 실제 득표율 벡터 (순서는 모델의 지역 순서와 일치해야 함)
# 예시용 더미 데이터입니다. 실제 데이터를 여기에 할당하세요.
# actual_y_20 <- c(0.35, 0.42, ...) 
if (!exists("actual_y_20")) {
  stop("actual_y_20 변수가 정의되지 않았습니다. 실제 20대 대선 득표율 벡터를 생성해주세요.")
}

# 2. 전국 득표율 계산 함수
calc_national_share <- function(actual_vec, weight_col_name) {
  weight_vec <- Pop_weight_21[, weight_col_name] 
  
  return(sum(actual_vec * weight_vec))
}

# [사용 예시] 
# actual_y_20, actual_y_21 벡터가 있다고 가정할 때:
national_share_20 <- calc_national_share(actual_y_20, "pe20")
national_share_21 <- calc_national_share(actual_y_21, "pe21")

# print(paste("2020 전국 득표율:", national_share_20))
# print(paste("2021 전국 득표율:", national_share_21))

# -------------------------------------------------------------------
# 1. 평가 함수 정의
# -------------------------------------------------------------------
evaluate_model <- function(tag, dir, actuals, pop_weights) {
  
  # 1. 파일 확인
  rds_path <- file.path(dir, paste0("fit_", tag, ".RDS"))
  if (!file.exists(rds_path)) {
    # 파일이 정말 없는 경우
    warning(paste("File not found:", rds_path))
    return(NULL)
  }
  
  # 2. 파일 로드
  fit <- readRDS(rds_path)
  
  cat(paste0("\n>>> Evaluating Model: ", tag, " <<<\n"))
  
  # -----------------------------------------------------------
  # [수정] 변수 탐지 로직 (유연하게 변경)
  # 메타데이터를 확인하지 않고, 직접 summary를 시도합니다.
  # -----------------------------------------------------------
  target_var <- "y_pred_share"
  
  # 일단 y_pred_share로 시도해봅니다. (silent=TRUE로 에러 메시지 숨김)
  check_try <- try(fit$summary(target_var, "mean"), silent = TRUE)
  
  # 만약 에러가 났다면(변수가 없다면), y_pred로 다시 시도합니다.
  if (inherits(check_try, "try-error")) {
    target_var <- "y_pred"
    check_try <- try(fit$summary(target_var, "mean"), silent = TRUE)
    
    # y_pred도 없다면? 진짜 없는 겁니다.
    if (inherits(check_try, "try-error")) {
      warning(paste0("[Skip] '", tag, "' 모델에 예측 변수가 없습니다."))
      return(NULL)
    }
  }
  
  cat(paste0("   Target variable found: ", target_var, "\n"))
  
  # -----------------------------------------------------------
  # [기존 로직 동일] 요약 통계량 추출
  # -----------------------------------------------------------
  pred_summary <- fit$summary(variables = target_var, 
                              "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975)))
  
  # 가중치 길이 검증
  if (nrow(pred_summary) != length(pop_weights)) {
    stop(paste0("Error: 모델 지역 수(", nrow(pred_summary), 
                ") != 가중치 수(", length(pop_weights), ")"))
  }
  
  # -----------------------------------------------------------
  # [기존 로직 동일] 전국 지표 계산을 위한 8000개 샘플 추출
  # -----------------------------------------------------------
  # Draws 추출 (여기도 에러 날 수 있으니 tryCatch)
  post_draws_try <- try(fit$draws(variables = target_var, format = "draws_matrix"), silent = TRUE)
  
  if (inherits(post_draws_try, "try-error")) {
    # 만약 draws를 못 뽑으면(용량 문제 등), 전국 RMSE는 NA로 처리하고 넘어갑니다.
    nat_rmse_post <- NA
  } else {
    posterior_draws <- post_draws_try
    weight_mat <- as.matrix(pop_weights) 
    nat_preds_dist <- posterior_draws %*% weight_mat 
    nat_actual_val <- sum(actuals * pop_weights)
    nat_rmse_post <- sqrt(mean((nat_preds_dist - nat_actual_val)^2))
  }
  
  # -----------------------------------------------------------
  # 데이터 프레임 정리
  # -----------------------------------------------------------
  results_df <- pred_summary %>%
    rename(pred_mean = mean, pred_lower = `2.5%`, pred_upper = `97.5%`) %>%
    mutate(
      region_idx = row_number(),
      actual = actuals,
      weight = pop_weights,
      error = pred_mean - actual,
      abs_error = abs(pred_mean - actual),
      sq_error = (pred_mean - actual)^2,
      covered = (actual >= pred_lower & actual <= pred_upper),
      width = pred_upper - pred_lower
    )
  
  # 전국 단위(National) 메트릭 계산 (점추정치 기준)
  nat_actual_val <- sum(results_df$actual * results_df$weight)
  nat_pred_val   <- sum(results_df$pred_mean * results_df$weight)
  nat_bias       <- nat_pred_val - nat_actual_val
  
  metrics <- list(
    Model = tag,
    Var_Name = target_var,
    RMSE_Region = sqrt(mean(results_df$sq_error)),
    Bias_Region = mean(results_df$error),
    Width_Region = mean(results_df$width),
    
    Nat_Actual = nat_actual_val,
    Nat_Pred = nat_pred_val,
    Nat_Bias = nat_bias,
    Nat_RMSE_Post = nat_rmse_post, # 여기가 선생님의 아이디어
    
    Max_Rhat = max(fit$summary()$rhat, na.rm=TRUE)
  )
  
  return(list(metrics = metrics, details = results_df))
}
# -------------------------------------------------------------------
# 2. 전체 모델 일괄 평가 실행
# -------------------------------------------------------------------
library(dplyr)

# 모델 목록
model_tags <- c("m2_logit5_delta","m3_beta1_delta","m4_logit_ar1_3","m4_logit_ar1_4","m2_logit6_0103","m3_beta2_0103",
                "m3_beta3","m4_logit_ar1_4")

# 결과 담을 리스트 초기화
all_metrics <- list()
all_details <- list()
weights21 <- Pop_weight_21[, "pe21"]

# 루프 실행
for (tag in model_tags) {
  
  cat(paste0("Processing: ", tag, "... "))
  
  # tryCatch로 감싸서 에러가 나도 멈추지 않게 함
  tryCatch({
    
    # 함수 실행
    res <- evaluate_model(tag, out_dir, actual_y_21, weights21)
    
    if (!is.null(res)) {
      all_metrics[[tag]] <- res$metrics
      all_details[[tag]] <- res$details
      cat("OK\n")
    } else {
      cat("Skipped (NULL returned)\n")
    }
    
  }, error = function(e) {
    # ★ 에러 발생 시 멈추지 않고 원인을 출력 ★
    cat(paste0("FAILED! Error: ", e$message, "\n"))
  })
}

# 결과 합치기 (성공한 것들만)
if (length(all_metrics) > 0) {
  final_df <- bind_rows(all_metrics) %>% arrange(RMSE_Region)
  print(final_df)
}
# -------------------------------------------------------------------
# 3. 결과 요약 및 출력
# -------------------------------------------------------------------

# (1) 모델별 비교표 (Summary Table)
comparison_table <- bind_rows(all_metrics) %>%
  arrange(Nat_RMSE_Post) # Nat_RMSE_Post가 낮은 순서로 정렬

print("=== Model Comparison Table ===")
print(as.data.frame(comparison_table), digits = 4)

# (2) 상세 결과 (Region-wise Comparison)
# 가장 성능이 좋은(RMSE 최소) 모델의 상세 결과를 봅니다.
best_model_tag <- comparison_table$Model[1]
print(paste0("\n=== Detailed Results for Best Model: ", best_model_tag, " ==="))
print(as.data.frame(all_details[[best_model_tag]]), digits = 4)

# (3) 파일로 저장
write_csv(comparison_table, file.path(out_dir, "model_comparison_metrics_0109.csv"))





# 
# evaluate_model <- function(tag, dir, actuals) {
#   
#   rds_path <- file.path(dir, paste0("fit_", tag, ".RDS"))
#   
#   if (!file.exists(rds_path)) {
#     warning(paste("File not found:", rds_path))
#     return(NULL)
#   }
#   fit <- readRDS(rds_path)
#   
#   cat(paste0("\n>>> Evaluating Model: ", tag, " <<<\n"))
#   
#   # -----------------------------------------------------------
#   # [수정] 변수 이름 자동 탐지 로직 추가
#   # -----------------------------------------------------------
#   # 모델이 가지고 있는 모든 변수 목록 확인
#   model_vars <- fit$metadata()$stan_variables
#   
#   if ("y_pred_share" %in% model_vars) {
#     target_var <- "y_pred_share"
#   } else if ("y_pred" %in% model_vars) {
#     target_var <- "y_pred"
#   } else {
#     warning(paste0("[Skip] 모델 '", tag, "'에서 예측 변수(y_pred_share 또는 y_pred)를 찾을 수 없습니다."))
#     return(NULL)
#   }
#   
#   cat(paste0("   Detected prediction variable: ", target_var, "\n"))
#   
#   # -----------------------------------------------------------
#   # 요약 통계량 추출 (찾아낸 target_var 사용)
#   # -----------------------------------------------------------
#   pred_summary <- fit$summary(variables = target_var, 
#                               "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975)))
#   
#   # 데이터 프레임 정리
#   results_df <- pred_summary %>%
#     rename(
#       pred_mean = mean,
#       pred_lower = `2.5%`,
#       pred_upper = `97.5%`
#     ) %>%
#     mutate(
#       region_idx = row_number(),
#       actual = actuals,
#       error = pred_mean - actual,
#       abs_error = abs(pred_mean - actual),
#       sq_error = (pred_mean - actual)^2,
#       covered = (actual >= pred_lower & actual <= pred_upper),
#       width = pred_upper - pred_lower
#     )
#   
#   # 성능 지표 계산
#   metrics <- list(
#     Model = tag,
#     Var_Name = target_var, # 어떤 변수를 썼는지 기록
#     RMSE = sqrt(mean(results_df$sq_error)),
#     MAE = mean(results_df$abs_error),
#     Bias = mean(results_df$error),
#     Coverage_95 = mean(results_df$covered),
#     Mean_Width = mean(results_df$width),
#     Max_Rhat = max(fit$summary()$rhat, na.rm=TRUE)
#   )
#   
#   return(list(metrics = metrics, details = results_df))
# }