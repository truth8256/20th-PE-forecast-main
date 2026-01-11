# 클러스터 울산 분리

# ────────────────────────────────────────────────
# 📘 데이터 전처리: 14~20대 (학습용) + 21대 (예측용)
# Kang(2024) fundamentals model 확장 버전
# ────────────────────────────────────────────────

rm(list = ls())
library(readr)
library(dplyr)

# ────────────────────────────────────────────────
# 1️⃣ 파일 불러오기 (14~21)
# ────────────────────────────────────────────────

### 각 선거별 결과는 시군구 단위로 이루어져있으며, 전체 후보자의 득표 및 무효, 기권이 표시되어 있음
### 1열에는 primary_cluster(시도), 2열에는 secondary_cluster(시군구) 표기되어 있음.

# pe14~pe20: 학습용
url <- "https://raw.githubusercontent.com/seungwoo-stat/south-korea-election/main/csv/"
pe14 <- read_csv(paste0(url, "14pe.csv"), show_col_types = FALSE)
pe15 <- read_csv(paste0(url, "15pe.csv"), show_col_types = FALSE)
pe16 <- read_csv(paste0(url, "16pe.csv"), show_col_types = FALSE)
pe17 <- read_csv(paste0(url, "17pe.csv"), show_col_types = FALSE)
pe18 <- read_csv(paste0(url, "18pe.csv"), show_col_types = FALSE)
pe19 <- read_csv(paste0(url, "19pe.csv"), show_col_types = FALSE)
pe20 <- read_csv("0_data/20pe.csv", show_col_types = FALSE)

# pe21: 예측용
pe21 <- read_csv("0_data/21pe.csv", show_col_types = FALSE)



# ────────────────────────────────────────────────
# 1-1️⃣ 14대 데이터에서 울산 분리 (경남 → 울산으로 재분류)
# ────────────────────────────────────────────────
# primary_cluster == "Gyeongsangnam-do" 이면서
# secondary_cluster에 "Ulsan"이 포함된 4개 행을 "Ulsan"으로 변경

pe14$primary_cluster[
  pe14$primary_cluster == "Gyeongsangnam-do" &
    grepl("Ulsan", pe14$secondary_cluster)
] <- "Ulsan"

# (체크용)
# pe14 |> dplyr::count(primary_cluster)
# pe14 |> dplyr::filter(primary_cluster == "Ulsan") |>
#   dplyr::select(primary_cluster, secondary_cluster) |>
#   dplyr::distinct()

# ────────────────────────────────────────────────
# 2️⃣ 광역자치단체 클러스터 정의  (울산 분리 버전)
# ────────────────────────────────────────────────
cluster <- list(
  "Seoul"            = "Seoul",
  "Incheon"          = "Incheon",
  "Daejeon"          = "Daejeon",
  "Chungcheongbuk-do"= "Chungcheongbuk-do",
  "Chungcheongnam-do"= c("Sejong-si", "Chungcheongnam-do"),
  "Gwangju"          = "Gwangju",
  "Jeollabuk-do"     = "Jeollabuk-do",
  "Jeollanam-do"     = "Jeollanam-do",
  "Daegu"            = "Daegu",
  "Gyeongsangbuk-do" = "Gyeongsangbuk-do",
  "Busan"            = "Busan",
  "Ulsan"            = "Ulsan",              # 🔹 새로 분리
  "Gyeongsangnam-do" = "Gyeongsangnam-do",   # 🔹 울산 제외
  "Gangwon-do"       = "Gangwon-do",
  "Jeju-do"          = "Jeju-do",
  "Gyeonggi-do"      = "Gyeonggi-do"
)

# ────────────────────────────────────────────────
# 3️⃣ 선거 데이터 목록화: 20대/21대 예측용 세트 동시 생성
# ────────────────────────────────────────────────
pe_fit_20 <- list(pe15, pe16, pe17, pe18, pe19)           # fit for pe20
names(pe_fit_20) <- paste0("pe", 15:19)

pe_fit_21 <- list(pe15, pe16, pe17, pe18, pe19, pe20)     # fit for pe21
names(pe_fit_21) <- paste0("pe", 15:20)

# prediction targets
pe_pred20 <- pe20
pe_pred21 <- pe21


# ────────────────────────────────────────────────
# 4️⃣ 집계 함수 + fit/pred 집계
# ────────────────────────────────────────────────
aggregate_by_cluster <- function(df, cluster){
  out <- matrix(nrow = length(cluster), ncol = ncol(df) - 2)
  for (j in seq_along(cluster)) {
    rows <- df$primary_cluster %in% cluster[[j]]
    out[j, ] <- colSums(df[rows, -(1:2)], na.rm = TRUE)
  }
  rownames(out) <- names(cluster)
  colnames(out) <- colnames(df)[-(1:2)]
  out
}

pe_sum_fit_20 <- lapply(pe_fit_20, aggregate_by_cluster, cluster = cluster)  # pe15~pe19
pe_sum_fit_21 <- lapply(pe_fit_21, aggregate_by_cluster, cluster = cluster)  # pe15~pe20

pe20_sum <- aggregate_by_cluster(pe_pred20, cluster)                         # pe20
pe21_sum <- aggregate_by_cluster(pe_pred21, cluster)                         # pe21


# ────────────────────────────────────────────────
# 5️⃣ 민주당 양당 득표율 (fit용): 20대/21대 버전 각각 생성
# ────────────────────────────────────────────────
democ_index_20 <- c(pe15=2, pe16=2, pe17=1, pe18=2, pe19=1, pe20=1)          # 필요시 수정
democ_index_21 <- c(pe15=2, pe16=2, pe17=1, pe18=2, pe19=1, pe20=1, pe21=1)  # 필요시 수정

twoparty_vote_share <- function(mat, index){
  mat[, index] / rowSums(mat[, 1:2])
}

stopifnot(all(names(pe_sum_fit_20) %in% names(democ_index_20)))
stopifnot(all(names(pe_sum_fit_21) %in% names(democ_index_21)))

pe_twoparty_vote_share_20 <- t(sapply(names(pe_sum_fit_20), function(nm){
  twoparty_vote_share(pe_sum_fit_20[[nm]], democ_index_20[nm])
}))
rownames(pe_twoparty_vote_share_20) <- names(pe_sum_fit_20)
colnames(pe_twoparty_vote_share_20) <- names(cluster)

pe_twoparty_vote_share_21 <- t(sapply(names(pe_sum_fit_21), function(nm){
  twoparty_vote_share(pe_sum_fit_21[[nm]], democ_index_21[nm])
}))
rownames(pe_twoparty_vote_share_21) <- names(pe_sum_fit_21)
colnames(pe_twoparty_vote_share_21) <- names(cluster)

cat("\n✅ twoparty vote share created: 20-case(pe15~19), 21-case(pe15~20)\n")


# ────────────────────────────────────────────────
# 6️⃣ Pop_weight 생성 (fit+pred 포함): 20대용/21대용
#    - 각 선거별 지역 규모(총합)를 weight로 사용
#    - pe_sum 행렬의 마지막 2열(무효/기권 등) 제외 후 rowSums
# ────────────────────────────────────────────────
make_pop_weight <- function(pe_sum_list){
  # pe_sum_list: list of matrices, each P x C
  raw <- sapply(pe_sum_list, function(m){
    stopifnot(ncol(m) >= 3)                    # 최소한 후보2 + 기타2 가정
    rowSums(m[, 1:(ncol(m) - 2), drop=FALSE])  # 마지막 2열 제외
  })
  # sapply 결과: P x T 매트릭스(보장)
  w <- t(t(raw) / colSums(raw))                # 각 선거열 합이 1이 되게 정규화
  rownames(w) <- rownames(pe_sum_list[[1]])
  colnames(w) <- names(pe_sum_list)
  w
}

# 20대 예측: pe15~pe19(fit) + pe20(pred) → 6열
pe_sum_for_weight_20 <- c(pe_sum_fit_20, list(pe20 = pe20_sum))
names(pe_sum_for_weight_20)[length(pe_sum_for_weight_20)] <- "pe20"
Pop_weight_20 <- make_pop_weight(pe_sum_for_weight_20)

# 21대 예측: pe15~pe20(fit) + pe21(pred) → 7열
pe_sum_for_weight_21 <- c(pe_sum_fit_21, list(pe21 = pe21_sum))
names(pe_sum_for_weight_21)[length(pe_sum_for_weight_21)] <- "pe21"
Pop_weight_21 <- make_pop_weight(pe_sum_for_weight_21)

cat("\n✅ Pop_weight created: Pop_weight_20(pe15~pe20), Pop_weight_21(pe15~pe21)\n")


# ────────────────────────────────────────────────
# 7️⃣ run-fund-model에서 쓰기 좋게 이름 정리(권장)
# ────────────────────────────────────────────────
pe_sum_20 <- pe_sum_fit_20
pe_sum_21 <- pe_sum_fit_21


# ────────────────────────────────────────────────
# actual_y_20, actual_y_21 생성
# ────────────────────────────────────────────────

# actual_y_20: pe20 기준 실제 양당 득표율
actual_y_20 <- as.numeric(twoparty_vote_share(pe20_sum, 1))

# actual_y_21: pe21 기준 실제 양당 득표율
actual_y_21 <- as.numeric(twoparty_vote_share(pe21_sum, 1))

# 이름 붙이기 (선택)
names(actual_y_20) <- rownames(pe20_sum)
names(actual_y_21) <- rownames(pe21_sum)
