library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)

# ─────────────────────────────────────
# 1) 요약 통계 (평균, 95% CI)
# ─────────────────────────────────────

pe_twoparty_vote_share <- readRDS("pe_twoparty_vote_share_14to20.rds")

summ_lin <- summarise_draws(
  y_pred_draws_20_2_linear,
  mean,
  sd,
  q5  = ~quantile(.x, 0.05),
  q95 = ~quantile(.x, 0.95)
) |>
  mutate(model = "linear")

summ_logit <- summarise_draws(
  y_pred_draws_20_2_logit,
  mean,
  sd,
  q5  = ~quantile(.x, 0.05),
  q95 = ~quantile(.x, 0.95)
) |>
  mutate(model = "logit")


# y_pred의 인덱스(1..P)를 province 이름으로 붙이기
prov_names <- names(cluster)       # 이미 쓰고 있는 이름이라고 가정
P <- length(prov_names)

summ_lin <- summ_lin |>
  mutate(province = rep(prov_names, times = 1))  # y_pred[1], ..., y_pred[P]

summ_logit <- summ_logit |>
  mutate(province = rep(prov_names, times = 1))

# 하나로 합치기 (wide 형태로 비교)
cmp <- full_join(
  summ_lin  |> select(province, mean_lin  = mean, q5_lin  = "5%", q95_lin  = "95%"),
  summ_logit|> select(province, mean_log = mean, q5_log = "5%", q95_log  = "95%"),
  by = "province"
)

cmp

cmp2 <- cmp |>
  mutate(
    diff_mean   = mean_log - mean_lin,              # 평균 차이 (로짓 - 선형)
    width_lin   = q95_lin - q5_lin,                 # 선형 90% CI 폭
    width_log   = q95_log - q5_log,                 # 로짓 90% CI 폭
    ratio_width = width_log / width_lin             # 로짓/선형 구간폭 비율
  ) |>
  arrange(desc(abs(diff_mean)))                     # 차이 큰 순서 보기

cmp2


### (1) 산점도로 비교
ggplot(cmp, aes(x = mean_lin, y = mean_log, label = province)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point() +
  geom_text(vjust = -0.5, size = 3) +
  coord_equal() +
  labs(
    x = "선형 모델 y_pred (mean, 0~1 스케일)",
    y = "로짓 모델 y_pred (mean, 0~1 스케일)",
    title = "선형 vs 로짓 펀더멘털 예측 비교 (주별)"
  )

### (2) 관측값과의 RMSE/MAE 비교

# 예: pe_twoparty_vote_share에서 마지막 행이 실제 결과라고 가정
# (실제 구조에 맞게 row 인덱스만 조정해 주세요)
y_true <- as.numeric(pe_twoparty_vote_share[TT, -(1+P)])  # 길이 P
names(y_true) <- prov_names

cmp2 <- cmp |>
  mutate(
    y_true = y_true[province],
    err_lin  = mean_lin  - y_true,
    err_log  = mean_log  - y_true,
    abs_lin  = abs(err_lin),
    abs_log  = abs(err_log),
    sq_lin   = err_lin^2,
    sq_log   = err_log^2
  )

data.frame(
  metric = c("MAE", "RMSE"),
  linear = c(mean(cmp2$abs_lin), sqrt(mean(cmp2$sq_lin))),
  logit  = c(mean(cmp2$abs_log), sqrt(mean(cmp2$sq_log)))
)

### 지역별 성능 비교 테이블 만들기

region_perf <- cmp |>
  mutate(
    # 1) 실제값 매칭
    y_true = y_true[province],
    
    # 2) 선형/로짓 예측 오차
    err_lin = mean_lin - y_true,
    err_log = mean_log - y_true,
    
    # 3) 절대오차, 제곱오차
    abs_lin = abs(err_lin),
    abs_log = abs(err_log),
    sq_lin  = err_lin^2,
    sq_log  = err_log^2,
    
    # 4) RMSE, MAE (지역별, 한 점이지만 이름 통일 차원에서)
    MAE_lin  = abs_lin,
    MAE_log  = abs_log,
    RMSE_lin = sqrt(sq_lin),
    RMSE_log = sqrt(sq_log),
    
    # 5) 예측 평균 차이, 구간폭 등
    diff_mean   = mean_log - mean_lin,          # 로짓 - 선형
    width_lin   = q95_lin - q5_lin,
    width_log   = q95_log - q5_log,
    ratio_width = width_log / width_lin
  ) |>
  select(
    province,
    y_true,
    MAE_lin, MAE_log,
    RMSE_lin, RMSE_log,
    diff_mean, width_lin, width_log, ratio_width
  )

region_perf


region_improve <- region_perf |>
  mutate(
    MAE_improvement  = MAE_lin  - MAE_log,   # +면 로짓이 더 나음
    RMSE_improvement = RMSE_lin - RMSE_log
  ) |>
  arrange(desc(MAE_improvement))

region_improve
