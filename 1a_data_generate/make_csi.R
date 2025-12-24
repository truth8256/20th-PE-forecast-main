setwd("D:\\OneDrive - 에스티아이\\=대학원\\=논문연구\\20th-PE-forecast-main")


library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# ─────────────────────────────────────────────
# 0. 지역명 영문 변환 매핑
# ─────────────────────────────────────────────
cluster_map <- c(
  "서울특별시" = "Seoul",
  "경기도" = "Gyeonggi-do",
  "인천광역시" = "Incheon",
  "대전광역시" = "Daejeon",
  "충청북도" = "Chungcheongbuk-do",
  "충청남도" = "Chungcheongnam-do",
  "광주광역시" = "Gwangju",
  "전라북도" = "Jeollabuk-do",
  "전북특별자치도" = "Jeollabuk-do",
  "전라남도" = "Jeollanam-do",
  "대구광역시" = "Daegu",
  "경상북도" = "Gyeongsangbuk-do",
  "부산광역시" = "Busan",
  "경상남도" = "Gyeongsangnam-do",
  "강원도" = "Gangwon-do",
  "강원특별자치도" = "Gangwon-do",
  "제주도" = "Jeju-do",
  "제주특별자치도" = "Jeju-do",
  "울산광역시" = "Ulsan",
  "세종특별자치시" = "Sejong-si"
)

# 충남/세종 명시
nm_cn <- "충청남도"
nm_sj <- "세종특별자치시"

# ─────────────────────────────────────────────
# 1. 파일 읽기
# ─────────────────────────────────────────────

# 시도별 가구수 (2019~2024)
hh_raw <- read_excel("0_data\\시도별가구수from2019to2024.xlsx")

# 소비자물가지수 (1990~2024)
cpi_raw <- read_excel("0_data\\소비자물가지수_시도.xlsx")


# ─────────────────────────────────────────────
# 2. CPI 롱폼 만들기 + 연도열 numeric 정리
# ─────────────────────────────────────────────

cpi_raw_fixed <- cpi_raw %>%
  rename(region = `시도별(2)`)

# 연도처럼 보이는 열만 선택 (1990~2099 범위)
year_cols <- grep("^19[0-9]{2}$|^20[0-9]{2}$", names(cpi_raw_fixed), value = TRUE)

cpi_raw_clean <- cpi_raw_fixed %>%
  mutate(across(
    all_of(year_cols),
    ~ suppressWarnings(as.numeric(.x))
  ))

cpi_long <- cpi_raw_clean %>%
  pivot_longer(
    cols = all_of(year_cols),
    names_to = "year",
    values_to = "cpi"
  ) %>%
  mutate(year = as.integer(year))

# ─────────────────────────────────────────────
# 3. 가구수 롱폼 만들기 + numeric 정리
# ─────────────────────────────────────────────

hh_long <- hh_raw %>%
  rename(region = `행정구역별(읍면동)`) %>%
  pivot_longer(
    cols = matches("^20(19|2[0-4])$"),   # 2019~2024
    names_to = "year",
    values_to = "hh"
  ) %>%
  mutate(
    year = as.integer(year),
    hh   = readr::parse_number(as.character(hh))  # "12,345" → 12345, "-" → NA
  )

# ─────────────────────────────────────────────
# 4. CPI + hh 롱 조인
# ─────────────────────────────────────────────

cpi_hh_long <- cpi_long %>%
  left_join(hh_long, by = c("region", "year"))
# (세종의 옛 시기, 2019 이전 등은 hh가 NA로 남음)

# ─────────────────────────────────────────────
# 5. 연도별 "충남+세종" 합성 CPI (롱폼 기반)
# ─────────────────────────────────────────────

cn_sj_by_year <- cpi_hh_long %>%
  filter(region %in% c(nm_cn, nm_sj)) %>%
  group_by(year) %>%
  summarise(
    # 각 연도에서 충남/세종 CPI, hh 추출
    cpi_cn = cpi[region == nm_cn][1],
    cpi_sj = cpi[region == nm_sj][1],
    hh_cn  = hh[region == nm_cn][1],
    hh_sj  = hh[region == nm_sj][1],
    
    has_cn = !is.na(cpi_cn),
    has_sj = !is.na(cpi_sj),
    
    cpi_cn_sj = case_when(
      # (1) 세종 CPI 없는 시기 → 충남 값 그대로
      has_cn & !has_sj ~ cpi_cn,
      
      # (2) CPI·hh 모두 존재 → 가구수 가중 평균
      has_cn & has_sj &
        !is.na(hh_cn) & !is.na(hh_sj) ~
        (hh_cn * cpi_cn + hh_sj * cpi_sj) / (hh_cn + hh_sj),
      
      # (3) CPI는 둘 다 있는데 hh가 NA → 단순 평균(백업용)
      has_cn & has_sj ~ (cpi_cn + cpi_sj) / 2,
      
      TRUE ~ NA_real_
    ),
    .groups = "drop"
  )

# ─────────────────────────────────────────────
# 6. 세종 제거 + 충남에 합성 CPI 덮어쓰기 (cpi_long_adj)
# ─────────────────────────────────────────────

cpi_long_adj <- cpi_long %>%
  filter(region != nm_sj) %>%  # 세종은 모델에서 사용 안 함 (충남+세종으로 흡수)
  left_join(
    cn_sj_by_year %>% select(year, cpi_cn_sj),
    by = "year"
  ) %>%
  mutate(
    cpi = if_else(region == nm_cn & !is.na(cpi_cn_sj), cpi_cn_sj, cpi)
  ) %>%
  select(region, year, cpi)

# ─────────────────────────────────────────────
# 7. 연도별 시도별 물가상승률 계산 (csi_raw)
# ─────────────────────────────────────────────

csi_raw <- cpi_long_adj %>%
  group_by(region) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    csi = (cpi - dplyr::lag(cpi)) / dplyr::lag(cpi)
  ) %>%
  ungroup() %>%
  mutate(
    cluster = cluster_map[region]
  )

# 원하시면 NA cluster는 제거 가능
# csi_raw <- csi_raw %>% filter(!is.na(cluster))

# 저장
write_csv(csi_raw, "csi_raw.csv")

# ─────────────────────────────────────────────
# 8. 선거별 csi4 / csi5 만들기 (pe_num, cluster, csi4, csi5)
# ─────────────────────────────────────────────

mapping_pe <- tibble(
  pe_num = 15:21,
  year4  = c(1996, 2001, 2006, 2011, 2015, 2020, 2023),
  year5  = c(1997, 2002, 2007, 2012, 2016, 2021, 2024)
)

clusters <- sort(unique(csi_raw$cluster))

csi <- tidyr::expand_grid(
  cluster = clusters,
  mapping_pe
) %>%
  left_join(
    csi_raw %>%
      select(cluster, year, csi) %>%
      rename(year4 = year, csi4 = csi),
    by = c("cluster", "year4")
  ) %>%
  left_join(
    csi_raw %>%
      select(cluster, year, csi) %>%
      rename(year5 = year, csi5 = csi),
    by = c("cluster", "year5")
  ) %>%
  select(pe_num, cluster, csi4, csi5) %>%
  arrange(cluster, pe_num)

write_csv(csi, "csi21.csv")


# #############################################
# # csi.csv(기존 버전) 과 csi21.csv(21대 선거용 최신 버전) 비교
# 
# # 1) 파일 읽기 ---------------------------------------------------------
# csi_old <- read_csv("0_data/csi.csv")      # cluster, pe_num, csi_increase_4, csi_increase_5
# csi_new <- read_csv("0_data/csi21.csv")    # cluster, pe_num, csi4, csi5
# 
# # 2) 컬럼 이름 맞추기 --------------------------------------------------
# 
# csi_old <- csi_old %>%
#   rename(
#     csi4_old = csi_increase_4,
#     csi5_old = csi_increase_5
#   )
# 
# csi_new <- csi_new %>%
#   rename(
#     csi4_new = csi4,
#     csi5_new = csi5
#   )
# 
# # 3) cluster 목록 비교 -------------------------------------------------
# 
# clusters_old <- sort(unique(csi_old$cluster))
# clusters_new <- sort(unique(csi_new$cluster))
# 
# cat("=== 기존 파일(csi) 지역 목록 ===\n")
# print(clusters_old)
# 
# cat("\n=== 신규 파일(csi21) 지역 목록 ===\n")
# print(clusters_new)
# 
# cat("\n=== 신규에만 있는 지역 ===\n")
# print(setdiff(clusters_new, clusters_old))
# 
# cat("\n=== 기존에만 있는 지역 ===\n")
# print(setdiff(clusters_old, clusters_new))
# 
# # 4) 같은 cluster + pe_num 기준으로 값 비교 ---------------------------
# 
# csi_compare <- csi_old %>%
#   full_join(csi_new, by = c("cluster", "pe_num")) %>%
#   mutate(
#     diff_csi4 = csi4_new - csi4_old,
#     diff_csi5 = csi5_new - csi5_old
#   ) %>%
#   arrange(cluster, pe_num)
# 
# write_csv(csi_compare, "0_data/csi_compare.csv")
# 
# # 5) 차이가 큰 것만 따로 보기 (예: 0.1%p 이상) -----------------------
# 
# csi_diff_big <- csi_compare %>%
#   filter(
#     (!is.na(diff_csi4) & abs(diff_csi4) > 0.1) |
#       (!is.na(diff_csi5) & abs(diff_csi5) > 0.1)
#   )
# 
# write_csv(csi_diff_big, "0_data/csi_diff_big.csv")
# 
# cat("\n=== 차이가 큰 부분 (0.1p 이상) ===\n")
# print(csi_diff_big)



