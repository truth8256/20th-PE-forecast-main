# setwd("D:\\OneDrive - 에스티아이\\=대학원\\=논문연구\\20th-PE-forecast-main")


# ============================================================
# 📘 데이터 전처리: 14~20대 (학습용) + 21대 (예측용)
#   - 울산 분리 버전
#   - fundamentals model용 Stan data_list_20 / data_list_21 생성
# ============================================================

rm(list = ls())
load("pe.RData")
library(readr)
library(dplyr)
library(purrr)
library(stringr)

# ------------------------------------------------------------
# 0. 선거 데이터 불러오기 (14~21)
# ------------------------------------------------------------

url <- "https://raw.githubusercontent.com/seungwoo-stat/south-korea-election/main/csv/"

# pe14~pe20: 학습용
pe14 <- read_csv(paste0(url, "14pe.csv"), show_col_types = FALSE)
pe15 <- read_csv(paste0(url, "15pe.csv"), show_col_types = FALSE)
pe16 <- read_csv(paste0(url, "16pe.csv"), show_col_types = FALSE)
pe17 <- read_csv(paste0(url, "17pe.csv"), show_col_types = FALSE)
pe18 <- read_csv(paste0(url, "18pe.csv"), show_col_types = FALSE)
pe19 <- read_csv(paste0(url, "19pe.csv"), show_col_types = FALSE)
pe20 <- read_csv("0_data/20pe.csv", show_col_types = FALSE)

# pe21: 예측용
pe21 <- read_csv("0_data/21pe.csv", show_col_types = FALSE)

# ------------------------------------------------------------
# 1. 14대 데이터에서 울산 분리 (경남 → 울산으로 재분류)
# ------------------------------------------------------------

pe14$primary_cluster[
  pe14$primary_cluster == "Gyeongsangnam-do" &
    grepl("Ulsan", pe14$secondary_cluster)
] <- "Ulsan"

# ------------------------------------------------------------
# 2. 광역자치단체 클러스터 정의 (울산 분리 버전)
# ------------------------------------------------------------

cluster <- list(
  "Seoul"             = "Seoul",
  "Gyeonggi-do"       = "Gyeonggi-do",
  "Incheon"           = "Incheon",
  "Daejeon"           = "Daejeon",
  "Chungcheongbuk-do" = "Chungcheongbuk-do",
  "Chungcheongnam-do" = c("Sejong-si", "Chungcheongnam-do"),
  "Gwangju"           = "Gwangju",
  "Jeollabuk-do"      = "Jeollabuk-do",
  "Jeollanam-do"      = "Jeollanam-do",
  "Daegu"             = "Daegu",
  "Gyeongsangbuk-do"  = "Gyeongsangbuk-do",
  "Busan"             = "Busan",
  "Ulsan"             = "Ulsan",
  "Gyeongsangnam-do"  = "Gyeongsangnam-do",
  "Gangwon-do"        = "Gangwon-do",
  "Jeju-do"           = "Jeju-do"
)

P <- length(cluster)   # 16개 지역

# ------------------------------------------------------------
# 3. 14~20대 데이터 목록화 및 클러스터별 집계 (pe_sum)
# ------------------------------------------------------------

pe <- list(pe14, pe15, pe16, pe17, pe18, pe19, pe20)
names(pe) <- paste0("pe", 14:20)

pe_sum <- vector("list", length = length(pe))
names(pe_sum) <- names(pe)

for (i in seq_along(pe)) {
  df <- pe[[i]]
  pe_sum[[i]] <- matrix(nrow = P, ncol = ncol(df) - 2)
  for (j in seq_along(cluster)) {
    rows <- df$primary_cluster %in% cluster[[j]]
    pe_sum[[i]][j, ] <- colSums(df[rows, -(1:2)], na.rm = TRUE)
  }
  rownames(pe_sum[[i]]) <- names(cluster)
  colnames(pe_sum[[i]]) <- colnames(df)[-(1:2)]
}

# ------------------------------------------------------------
# 4. 민주당 양당 득표율 계산 (14~20대)
# ------------------------------------------------------------

# 민주당 후보가 후보열에서 몇 번째인지 (14~20)
democ_index <- c(2, 2, 2, 1, 2, 1, 1)

twoparty_vote_share <- \(mat, index) {
  if (is.vector(mat)) return(mat[index] / sum(mat[1:2]))
  mat[, index] / rowSums(mat[, 1:2])
}

pe_twoparty_vote_share <- matrix(
  ncol = P + 1,           # 16개 지역 + national
  nrow = length(pe)       # pe14~pe20
)

for (i in seq_along(pe)) {
  pe_twoparty_vote_share[i, 1:P] <-
    twoparty_vote_share(pe_sum[[i]], democ_index[i])
}

rownames(pe_twoparty_vote_share) <- names(pe)                    # "pe14"~"pe20"
colnames(pe_twoparty_vote_share) <- c(names(cluster), "national")

for (i in seq_along(pe)) {
  pe_twoparty_vote_share[i, P + 1] <-
    twoparty_vote_share(colSums(pe[[i]][, -(1:2)]), democ_index[i])
}

cat("\n✅ 민주당 양당 득표율 계산 완료 (14~20대)\n")

# ------------------------------------------------------------
# 5. 21대 예측용: 클러스터 합계 + 민주당 양당 득표율 (진짜 결과)
# ------------------------------------------------------------

pe21_sum <- matrix(nrow = P, ncol = ncol(pe21) - 2)

for (j in seq_along(cluster)) {
  rows <- pe21$primary_cluster %in% cluster[[j]]
  pe21_sum[j, ] <- colSums(pe21[rows, -(1:2)], na.rm = TRUE)
}
rownames(pe21_sum) <- names(cluster)
colnames(pe21_sum) <- colnames(pe21)[-(1:2)]

# 민주당 열 인덱스 (필요시 조정)
democ_index_21 <- 1  # 예: 1열이 민주당

pe21_share <- pe21_sum[, democ_index_21] / rowSums(pe21_sum[, 1:2])

# 나중에 21대 예측 성능 비교할 때 사용할 truth
pe21_truth <- data.frame(
  cluster   = names(cluster),
  dem_share = as.numeric(pe21_share)
)

write_csv(pe21_truth, "0_data/21pe_twoparty_share.csv")
cat("\n💾 21대 실제 양당 득표율 저장 완료 → 0_data/21pe_twoparty_share.csv\n")

# ------------------------------------------------------------
# 6. 14~20대 투표가중치 Pop_weight_14_20 (지역별 turnout 비율)
# ------------------------------------------------------------

# 각 선거(pe14~pe20)에서 지역별 총 득표수 합산
Pop_weight_14_20_raw <- sapply(
  pe_sum,
  \(m) rowSums(m, na.rm = TRUE)
)
Pop_weight_14_20 <- t(t(Pop_weight_14_20_raw) / colSums(Pop_weight_14_20_raw))

rownames(Pop_weight_14_20) <- names(cluster)     # 지역명
colnames(Pop_weight_14_20) <- names(pe_sum)      # "pe14" ~ "pe20"

# ------------------------------------------------------------
# 7. 20대/21대용 Pop_weight, Alpha_init 계산
# ------------------------------------------------------------

# 20대 예측용: t=1..5 → pe15~pe19, t=6 → pe20
Pop_weight_20 <- Pop_weight_14_20[, c("pe15","pe16","pe17","pe18","pe19","pe20")]

# 21대 예측용: t=1..6 → pe15~pe20, t=7 → pe21
turnout21_raw <- rowSums(pe21_sum, na.rm = TRUE)
w21 <- turnout21_raw / sum(turnout21_raw)

Pop_weight_21 <- cbind(
  Pop_weight_14_20[, c("pe15","pe16","pe17","pe18","pe19","pe20")],
  "pe21" = w21
)

# Alpha_init (14대 baseline, 두 모델에서 공통 사용)
Alpha_init <- {
  y14 <- pe_twoparty_vote_share["pe14", 1:P]
  w14 <- Pop_weight_14_20[, "pe14"]
  y14 - sum(w14 * y14) + 0.5
}

# ============================================================
# 8. 펀더멘털 디자인 매트릭스 X, Z 불러오기
#    (⚠️ 여기는 본인 프로젝트에서 사용 중인 파일명으로 수정)
# ============================================================

# 예시: X, Z가 15~21대 × 16개 지역 순서로 112행에 들어있다고 가정
# load("pe.RData")  # 여기에 X, Z, cluster 등이 저장되어 있다면 이렇게 사용
# 또는:
# load("0_data/fundamental_design_pe15_21.RData")  # X, Z 포함 RData

# 여기서는 X, Z가 이미 존재한다고 가정
# dim(X); dim(Z)  # 112 x K, 112 x L 이어야 함

standardize <- function(mat){
  mat <- t(mat) - colMeans(mat)
  t(mat / sqrt(rowSums(mat^2)))
}

X_std <- standardize(X)
Z_std <- standardize(Z)

# ------------------------------------------------------------
# 9. 선거별 인덱스 정의 (pe15~pe21 → X/Z row index 매핑)
# ------------------------------------------------------------

pe_vec <- 15:21
n_pe   <- length(pe_vec)  # 7
P      <- length(cluster)

idx_by_pe <- setNames(
  lapply(seq_along(pe_vec), \(k) {
    start <- (k - 1) * P + 1
    end   <- k * P
    start:end
  }),
  paste0("pe", pe_vec)
)
# 예: idx_by_pe[["pe15"]] == 1:16, idx_by_pe[["pe20"]] == 81:96, idx_by_pe[["pe21"]] == 97:112

# ------------------------------------------------------------
# 10. 20대 예측용 Stan data_list_20  (15~19 학습, 20 예측)
# ------------------------------------------------------------

train_pes_20 <- 15:19
TT_20 <- length(train_pes_20)      # 5
N_20  <- P * TT_20                 # 16 * 5 = 80

idx_train_20 <- unlist(idx_by_pe[paste0("pe", train_pes_20)])  # pe15~pe19
idx_pred_20  <- idx_by_pe[["pe20"]]                            # pe20

# Y_20: pe15~pe19의 민주당 양당 득표율 (지역 16개만)
Y_20 <- c(t(
  pe_twoparty_vote_share[
    paste0("pe", train_pes_20),   # 행: pe15~pe19
    1:P                           # 열: national 제외
  ]
))

data_list_20 <- list(
  P = P,
  TT = TT_20,
  N = N_20,
  K = ncol(X_std),           # ← 여기
  L = ncol(Z_std),           # ← 여기
  X = X_std[idx_train_20, , drop = FALSE],
  Z = Z_std[idx_train_20, , drop = FALSE],
  X_pred = X_std[idx_pred_20, , drop = FALSE],
  Z_pred = Z_std[idx_pred_20, , drop = FALSE],
  Y = Y_20,
  Alpha_init = Alpha_init,
  Pop_weight = Pop_weight_20
)


# ------------------------------------------------------------
# 11. 21대 예측용 Stan data_list_21 (15~20 학습, 21 예측)
# ------------------------------------------------------------

train_pes_21 <- 15:20
TT_21 <- length(train_pes_21)      # 6
N_21  <- P * TT_21                 # 16 * 6 = 96

idx_train_21 <- unlist(idx_by_pe[paste0("pe", train_pes_21)])  # pe15~pe20
idx_pred_21  <- idx_by_pe[["pe21"]]                            # pe21

# Y_21: pe15~pe20의 민주당 양당 득표율
Y_21 <- c(t(
  pe_twoparty_vote_share[
    paste0("pe", train_pes_21),   # 행: pe15~pe20
    1:P
  ]
))

data_list_21 <- list(
  P = P,
  TT = TT_21,
  N = N_21,
  K = ncol(X_std),           # ← 여기
  L = ncol(Z_std),           # ← 여기
  X = X_std[idx_train_21, , drop = FALSE],
  Z = Z_std[idx_train_21, , drop = FALSE],
  X_pred = X_std[idx_pred_21, , drop = FALSE],
  Z_pred = Z_std[idx_pred_21, , drop = FALSE],
  Y = Y_21,
  Alpha_init = Alpha_init,
  Pop_weight = Pop_weight_21
)


cat("\n✅ Stan 모델용 data_list_20, data_list_21 생성 완료\n")

# 원하면 한 번에 저장
save(
  cluster, pe_sum, pe_twoparty_vote_share, pe21_sum, pe21_truth,
  Pop_weight_14_20, Pop_weight_20, Pop_weight_21,
  Alpha_init,
  X, Z, X_std, Z_std,
  data_list_20, data_list_21,
  file = "0_data/fundamental_model_data_ulsan_split.RData"
)
cat("💾 fundamental_model_data_ulsan_split.RData 로 저장 완료\n")


# ─────────────────────────────────────
# 🧠 모델 컴파일 및 샘플링
# ─────────────────────────────────────
library(cmdstanr)
library(posterior)

# ─────────────────────────────────────
# 🔧 Stan 모델 컴파일 (선형 / 로짓)
# ─────────────────────────────────────
mod_linear <- cmdstan_model("2_fundamental_model/fund-model-simple.stan")
mod_logit  <- cmdstan_model("2_fundamental_model/fund-model-simple_logit0.stan")

# 공통 샘플링 옵션
sampling_opts <- list(
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 8000,
  thin = 1
)

# ─────────────────────────────────────
# (1) mod20_2 : 15~19 학습, 20대 예측 (선형)
# ─────────────────────────────────────
fit20_2 <- do.call(
  mod_linear$sample,
  c(list(data = data_list_20), sampling_opts)
)

fit20_2$summary() |> View()

# posterior → data.frame
draws20_2 <- as_draws_df(fit20_2)

# 주요 파라미터 추출
beta_draws_20_2_linear   <- fit20_2$draws("beta")
gamma_draws_20_2_linear  <- fit20_2$draws("gamma")
y_pred_draws_20_2_linear <- fit20_2$draws("y_pred")  # pe20 예측값 (0~1 스케일)

# Stan raw output 저장 (원하면)
fit20_2$save_output_files(dir = "2_fundamental_model")

# RData 저장
save(
  beta_draws_20_2_linear,
  gamma_draws_20_2_linear,
  y_pred_draws_20_2_linear,
  file = "2_fundamental_model/samples-fundamental-linear-pe20-2.RData"
)

# ─────────────────────────────────────
# (2) mod_logit20_2 : 15~19 학습, 20대 예측 (로짓)
# ─────────────────────────────────────
fit_logit20_2 <- do.call(
  mod_logit$sample,
  c(list(data = data_list_20), sampling_opts)
)

fit_logit20_2$summary() |> View()

draws_logit20_2 <- as_draws_df(fit_logit20_2)

beta_draws_20_2_logit   <- fit_logit20_2$draws("beta")
gamma_draws_20_2_logit  <- fit_logit20_2$draws("gamma")
y_pred_draws_20_2_logit <- fit_logit20_2$draws("y_pred")  # 이미 inv_logit까지 적용된 0~1 스케일

fit_logit20_2$save_output_files(dir = "2_fundamental_model")

save(
  beta_draws_20_2_logit,
  gamma_draws_20_2_logit,
  y_pred_draws_20_2_logit,
  file = "2_fundamental_model/samples-fundamental-logit-pe20-2.RData"
)

# ─────────────────────────────────────
# (3) mod21 : 15~20 학습, 21대 예측 (선형)
# ─────────────────────────────────────
fit21 <- do.call(
  mod_linear$sample,
  c(list(data = data_list_21), sampling_opts)
)

fit21$summary() |> View()

draws21 <- as_draws_df(fit21)

beta_draws_21_linear   <- fit21$draws("beta")
gamma_draws_21_linear  <- fit21$draws("gamma")
y_pred_draws_21_linear <- fit21$draws("y_pred")  # pe21 예측값 (0~1 스케일)

fit21$save_output_files(dir = "2_fundamental_model")

save(
  beta_draws_21_linear,
  gamma_draws_21_linear,
  y_pred_draws_21_linear,
  file = "2_fundamental_model/samples-fundamental-linear-pe21.RData"
)

# ─────────────────────────────────────
# (4) mod_logit21 : 15~20 학습, 21대 예측 (로짓)
# ─────────────────────────────────────
fit_logit21 <- do.call(
  mod_logit$sample,
  c(list(data = data_list_21), sampling_opts)
)

fit_logit21$summary() |> View()

draws_logit21 <- as_draws_df(fit_logit21)

beta_draws_21_logit   <- fit_logit21$draws("beta")
gamma_draws_21_logit  <- fit_logit21$draws("gamma")
y_pred_draws_21_logit <- fit_logit21$draws("y_pred")  # 0~1 스케일

fit_logit21$save_output_files(dir = "2_fundamental_model")

save(
  beta_draws_21_logit,
  gamma_draws_21_logit,
  y_pred_draws_21_logit,
  file = "2_fundamental_model/samples-fundamental-logit-pe21.RData"
)

