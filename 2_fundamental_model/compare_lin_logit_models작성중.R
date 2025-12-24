library(dplyr)
library(tidyr)
library(readr)
library(posterior)
library(stringr)

# ───────────────────────────────────
# 1) truth, 가중치 등 데이터 불러오기
# ───────────────────────────────────
# load("0_data/fundamental_model_data_ulsan_split.RData")
# 포함 내용:
#   cluster, pe_twoparty_vote_share, pe21_truth,
#   Pop_weight_20, Pop_weight_21, P 등

cluster_old <- list("Seoul" = "Seoul",
                "Gyeonggi-do" = "Gyeonggi-do",
                "Incheon" = "Incheon",
                "Daejeon"="Daejeon",
                "Chungcheongbuk-do"="Chungcheongbuk-do",
                "Chungcheongnam-do"=c("Sejong-si","Chungcheongnam-do"),
                "Gwangju"="Gwangju",
                "Jeollabuk-do"="Jeollabuk-do",
                "Jeollanam-do"="Jeollanam-do",
                "Daegu"="Daegu",
                "Gyeongsangbuk-do"="Gyeongsangbuk-do",
                "Busan"="Busan",
                "Gyeongsangnam-do"=c("Ulsan","Gyeongsangnam-do"),
                "Gangwon-do"="Gangwon-do",
                "Jeju-do"="Jeju-do")

# 20대(지역별) 진짜 민주당 양당 득표율 (길이 P)
pe20_truth_region <- pe_twoparty_vote_share["pe20", 1:P]
names(pe20_truth_region) <- colnames(pe_twoparty_vote_share)[1:P]

# 20대 전국 득표율 (pe_twoparty_vote_share의 'national' 열)
pe20_truth_nat <- pe_twoparty_vote_share["pe20", "national"]

# 21대(지역별) 진짜 민주당 양당 득표율 (data.frame → named vector)
pe21_truth_region <- setNames(pe21_truth$dem_share, pe21_truth$cluster)

# 21대 전국 득표율 (없으면 turnout 가중평균으로 정의)
w21_truth <- Pop_weight_21[, "pe21"]
pe21_truth_nat <- sum(w21_truth * pe21_truth_region)

# 20/21대에서 사용할 지역별 turnout 가중치
w20 <- Pop_weight_20[, "pe20"]   # 20대
w21 <- Pop_weight_21[, "pe21"]   # 21대

cluster_names <- rownames(Pop_weight_20)   # 또는 names(cluster)


# ───────────────────────────────────
# 2) posterior draws 불러오기
#    (※ 파일명은 실제와 맞게 수정 필요)
# ───────────────────────────────────

# ① mod      : 선형, 개선 전 20대
load("samples-fundamental.RData")
#   객체: y_pred_draws_20_linear

# ② mod_logit: 로짓, 개선 전 20대
load("samples-fundamental-logit.RData")
#   객체: y_pred_draws_20_logit

# ③ mod20_2      : 선형, 개선 후 20대
load("2_fundamental_model/samples-fundamental-linear-pe20-2.RData")
#   객체: y_pred_draws_20_2_linear

# ④ mod_logit20_2: 로짓, 개선 후 20대
load("2_fundamental_model/samples-fundamental-logit-pe20-2.RData")
#   객체: y_pred_draws_20_2_logit

# ⑤ mod21        : 선형, 개선 후 21대
load("2_fundamental_model/samples-fundamental-linear-pe21.RData")
#   객체: y_pred_draws_21_linear

# ⑥ mod_logit21  : 로짓, 개선 후 21대
load("2_fundamental_model/samples-fundamental-logit-pe21.RData")
#   객체: y_pred_draws_21_logit



# ───────────────────────────────────
# 3) posterior draws 요약 함수
# ───────────────────────────────────
summarise_model <- function(draws, model_name, election,
                            cluster_names, truth_region, w_region, truth_nat) {
  # draws: fit$draws("y_pred") 결과 (draws_array 등)
  mat <- as_draws_matrix(draws)  # iter × P
  P_here <- length(cluster_names)
  stopifnot(ncol(mat) == P_here)
  
  # 지역별 posterior mean
  pred_region <- colMeans(mat)
  
  # 전국 득표율 posterior draws (turnout 가중 평균)
  nat_draws <- as.numeric(mat %*% w_region)
  nat_mean  <- mean(nat_draws)
  nat_ci    <- quantile(nat_draws, c(0.025, 0.975))
  
  tibble(
    model    = model_name,
    election = election,
    cluster  = cluster_names,
    truth    = as.numeric(truth_region[cluster_names]),
    pred     = pred_region,
    error    = pred - truth,
    abs_error = abs(error),
    # 전국 단위 요약은 행마다 붙지만, distinct()로 나중에 1행만 남길 것
    nat_pred       = nat_mean,
    nat_truth      = truth_nat,
    nat_error      = nat_mean - truth_nat,
    nat_abs_error  = abs(nat_mean - truth_nat),
    nat_lower_95   = nat_ci[1],
    nat_upper_95   = nat_ci[2]
  )
}



# ───────────────────────────────────
# 4) 6개 모형 결과 한 번에 결합
# ───────────────────────────────────
results_region <- bind_rows(
  # # 20대, 개선 전
  # summarise_model(
  #   draws         = y_pred_draws,
  #   model_name    = "mod",
  #   election      = "20대",
  #   cluster_names = cluster_names,
  #   truth_region  = pe20_truth_region,
  #   w_region      = w20,
  #   truth_nat     = pe20_truth_nat
  # ),
  # summarise_model(
  #   draws         = y_pred_draws_logit,
  #   model_name    = "mod_logit",
  #   election      = "20대",
  #   cluster_names = cluster_names,
  #   truth_region  = pe20_truth_region,
  #   w_region      = w20,
  #   truth_nat     = pe20_truth_nat
  # ),
  # 20대, 개선 후
  summarise_model(
    draws         = y_pred_draws_20_2_linear,
    model_name    = "mod20_2",
    election      = "20대",
    cluster_names = cluster,
    truth_region  = pe20_truth_region,
    w_region      = w20,
    truth_nat     = pe20_truth_nat
  ),
  summarise_model(
    draws         = y_pred_draws_20_2_logit,
    model_name    = "mod_logit20_2",
    election      = "20대",
    cluster_names = cluster,
    truth_region  = pe20_truth_region,
    w_region      = w20,
    truth_nat     = pe20_truth_nat
  ),
  # 21대, 개선 후
  summarise_model(
    draws         = y_pred_draws_21_linear,
    model_name    = "mod21",
    election      = "21대",
    cluster_names = cluster,
    truth_region  = pe21_truth_region,
    w_region      = w21,
    truth_nat     = pe21_truth_nat
  ),
  summarise_model(
    draws         = y_pred_draws_21_logit,
    model_name    = "mod_logit21",
    election      = "21대",
    cluster_names = cluster,
    truth_region  = pe21_truth_region,
    w_region      = w21,
    truth_nat     = pe21_truth_nat
  )
)


#########5. 요약 테이블 생성
#########5-1. 지역별 예측 오차(평균 MAE, RMSE)


region_summary <- results_region %>%
  group_by(election, model) %>%
  summarise(
    MAE_region  = mean(abs_error, na.rm = TRUE),
    RMSE_region = sqrt(mean(error^2, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(election, MAE_region)

region_summary

#######5-2. 전국 득표율 예측 성능
national_summary <- results_region %>%
  select(election, model,
         nat_pred, nat_truth,
         nat_error, nat_abs_error,
         nat_lower_95, nat_upper_95) %>%
  distinct() %>%
  arrange(election, nat_abs_error)

national_summary

