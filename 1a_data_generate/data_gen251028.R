# setwd("D:\\OneDrive - 에스티아이\\=대학원\\=논문연구\\20th-PE-forecast-main")

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
# 2️⃣ 광역자치단체 클러스터 정의
# ────────────────────────────────────────────────
cluster <- list(
  "Seoul" = "Seoul",
  "Gyeonggi-do" = "Gyeonggi-do",
  "Incheon" = "Incheon",
  "Daejeon" = "Daejeon",
  "Chungcheongbuk-do" = "Chungcheongbuk-do",
  "Chungcheongnam-do" = c("Sejong-si", "Chungcheongnam-do"),
  "Gwangju" = "Gwangju",
  "Jeollabuk-do" = "Jeollabuk-do",
  "Jeollanam-do" = "Jeollanam-do",
  "Daegu" = "Daegu",
  "Gyeongsangbuk-do" = "Gyeongsangbuk-do",
  "Busan" = "Busan",
  "Gyeongsangnam-do" = c("Ulsan", "Gyeongsangnam-do"),
  "Gangwon-do" = "Gangwon-do",
  "Jeju-do" = "Jeju-do"
)

# ────────────────────────────────────────────────
# 3️⃣ 14~20대 데이터 목록화
# ────────────────────────────────────────────────
pe <- list(pe14, pe15, pe16, pe17, pe18, pe19, pe20)
names(pe) <- paste0("pe", 14:20)

# ────────────────────────────────────────────────
# 4️⃣ 클러스터별 득표수 합계 (광역단위 집계)
# ────────────────────────────────────────────────
pe_sum <- vector("list", length = length(pe))
names(pe_sum) <- names(pe)

for (i in seq_along(pe)) {
  df <- pe[[i]]
  pe_sum[[i]] <- matrix(nrow = length(cluster), ncol = ncol(df) - 2)
  for (j in seq_along(cluster)) {
    rows <- df$primary_cluster %in% cluster[[j]]
    pe_sum[[i]][j, ] <- colSums(df[rows, -(1:2)], na.rm = TRUE)
  }
  rownames(pe_sum[[i]]) <- names(cluster)
  colnames(pe_sum[[i]]) <- colnames(df)[-(1:2)]
}

# ────────────────────────────────────────────────
# 5️⃣ 민주당 양당 득표율 계산
# ────────────────────────────────────────────────
# 민주당 후보가 후보열에서 몇 번째인지
democ_index <- c(2, 2, 2, 1, 2, 1, 1)  # 14~20

twoparty_vote_share <- \(mat, index) {
  if (is.vector(mat)) return(mat[index] / sum(mat[1:2]))
  mat[, index] / rowSums(mat[, 1:2])
}

pe_twoparty_vote_share <- matrix(
  ncol = length(cluster) + 1,
  nrow = length(pe)
)

for (i in seq_along(pe)) {
  pe_twoparty_vote_share[i, 1:length(cluster)] <-
    twoparty_vote_share(pe_sum[[i]], democ_index[i])
}

rownames(pe_twoparty_vote_share) <- names(pe)
colnames(pe_twoparty_vote_share) <- c(names(cluster), "national")

for (i in seq_along(pe)) {
  pe_twoparty_vote_share[i, length(cluster) + 1] <-
    twoparty_vote_share(colSums(pe[[i]][, -(1:2)]), democ_index[i])
}

cat("\n✅ 민주당 양당 득표율 계산 완료 (14~20대)\n")


# ────────────────────────────────────────────────
# 6️⃣ 21대 예측용 데이터: 클러스터링 + 민주당 득표율 계산
# ────────────────────────────────────────────────
pe21_sum <- matrix(nrow = length(cluster), ncol = ncol(pe21) - 2)

for (j in seq_along(cluster)) {
  rows <- pe21$primary_cluster %in% cluster[[j]]
  pe21_sum[j, ] <- colSums(pe21[rows, -(1:2)], na.rm = TRUE)
}
rownames(pe21_sum) <- names(cluster)
colnames(pe21_sum) <- colnames(pe21)[-(1:2)]

# 민주당 열 지정 (필요시 조정)
democ_index_21 <- 1  # 예시: 첫 번째 열이 민주당일 경우
pe21_share <- pe21_sum[, democ_index_21] / rowSums(pe21_sum[, 1:2])

pe21_output <- data.frame(
  cluster = names(cluster),
  dem_share = round(pe21_share, 4)
)

write_csv(pe21_output, "0_data/21pe_twoparty_share.csv")
cat("\n💾 21대 예측용 득표율 저장 완료 → 0_data/21pe_twoparty_share.csv\n")
