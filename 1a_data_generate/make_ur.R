library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)

setwd("D:\\OneDrive - 에스티아이\\=대학원\\=논문연구\\20th-PE-forecast-main\\0_data")

# ─────────────────────────────────────────────
# 1️⃣ 데이터 읽기 및 long 변환
# ─────────────────────────────────────────────
read_unemp <- function(path, weeks_label) {
  # 1️⃣ 머리글(2행) 읽기
  hdr <- read_excel(path, n_max = 2, col_names = FALSE)
  header1 <- as.character(unlist(hdr[1, ]))
  header2 <- as.character(unlist(hdr[2, ]))
  names_combined <- ifelse(is.na(header2), header1, paste(header1, header2))
  names_combined[1] <- "region"
  
  # 2️⃣ 본문 읽기 (모든 열을 문자로 강제)
  df <- read_excel(path, skip = 2, col_names = FALSE, col_types = "text")
  colnames(df) <- names_combined
  
  # 3️⃣ 길게 펼치기
  df_long <- df %>%
    pivot_longer(-region, names_to = "col", values_to = "value") %>%
    mutate(
      period = str_extract(col, "\\d{4}\\.[1-4]/4"),
      var = case_when(
        str_detect(col, "경제활동") ~ "labor_force",
        str_detect(col, "실업자")   ~ "unemployed",
        str_detect(col, "실업률")   ~ "ur",
        TRUE ~ NA_character_
      ),
      value = gsub(",", "", value),
      value = ifelse(value %in% c("", "-", "NA", "N/A"), NA, value),
      value = as.numeric(value),
      year  = as.integer(str_extract(period, "^\\d{4}")),
      qtr   = as.integer(str_extract(period, "(?<=\\.)[1-4]"))
    ) %>%
    filter(!is.na(var), !is.na(period))
  
  # 4️⃣ 각 region-year-qtr별로 세 변수 결합
  df_wide <- df_long %>%
    group_by(region, year, qtr) %>%
    summarise(
      labor_force = value[var == "labor_force"][1],
      unemployed  = value[var == "unemployed"][1],
      ur          = value[var == "ur"][1],
      .groups = "drop"
    ) %>%
    mutate(
      week_basis = weeks_label
    )
  
  return(df_wide)
}



unemp_1w <- read_unemp("unemp_1w.xlsx", "1w")
unemp_4w <- read_unemp("unemp_4w.xlsx", "4w")

# ─────────────────────────────────────────────
# 2️⃣ 병합 처리 (울산↔경남, 세종↔충남)
# ─────────────────────────────────────────────
combine_regions <- function(df) {
  # ────────────────────────────────
  # 1️⃣ 기본 지역 매핑
  # ────────────────────────────────
  region_map <- c(
    "울산광역시"      = "경상남도_통합",
    "경상남도"        = "경상남도_통합",
    "세종특별자치시"  = "충청남도_통합",
    "충청남도"        = "충청남도_통합"
  )
  
  # ────────────────────────────────
  # 2️⃣ 경상남도_통합: 항상 합산
  # ────────────────────────────────
  df_gyeongnam <- df %>%
    filter(region %in% c("울산광역시", "경상남도")) %>%
    mutate(region = "경상남도_통합") %>%
    group_by(region, year, qtr, week_basis) %>%
    summarise(
      labor_force = sum(labor_force, na.rm = TRUE),
      unemployed  = sum(unemployed, na.rm = TRUE),
      ur = 100 * unemployed / labor_force,
      .groups = "drop"
    )
  
  # ────────────────────────────────
  # 3️⃣ 충청남도_통합: 조건부 합산
  # ────────────────────────────────
  df_chungnam <- df %>%
    filter(region %in% c("충청남도", "세종특별자치시")) %>%
    group_by(year, qtr, week_basis) %>%
    summarise(
      n_exist = sum(!is.na(ur)),  # ur이 존재하는 지역 수
      labor_force_sum = sum(labor_force, na.rm = TRUE),
      unemployed_sum  = sum(unemployed, na.rm = TRUE),
      ur_chungnam = ur[region == "충청남도"][1],
      ur_sejong   = ur[region == "세종특별자치시"][1],
      .groups = "drop"
    ) %>%
    mutate(
      region = "충청남도_통합",
      ur = case_when(
        n_exist == 2 ~ 100 * unemployed_sum / labor_force_sum,  # 둘 다 존재 → 합산
        n_exist == 1 & !is.na(ur_chungnam) ~ ur_chungnam,       # 충남만 존재 → 그대로
        n_exist == 1 & is.na(ur_chungnam) ~ ur_sejong,          # 세종만 존재 → 그대로
        TRUE ~ NA_real_
      ),
      labor_force = ifelse(n_exist == 2, labor_force_sum, NA),
      unemployed  = ifelse(n_exist == 2, unemployed_sum, NA)
    ) %>%
    select(region, year, qtr, week_basis, labor_force, unemployed, ur)
  
  # ────────────────────────────────
  # 4️⃣ 나머지 지역 그대로 유지
  # ────────────────────────────────
  df_keep <- df %>%
    filter(!region %in% names(region_map))
  
  # ────────────────────────────────
  # 5️⃣ 병합
  # ────────────────────────────────
  df_final <- bind_rows(df_keep, df_gyeongnam, df_chungnam) %>%
    arrange(region, year, qtr)
  
  return(df_final)
}



unemp_1w <- combine_regions(unemp_1w)
unemp_4w <- combine_regions(unemp_4w)

# ─────────────────────────────────────────────
# 3️⃣ 1주→4주 회귀 브리징 (1999.3/4~2014.4/4 구간에서 추정)
# ─────────────────────────────────────────────

# 공통 구간 추출
bridge_train <- unemp_1w %>%
  inner_join(unemp_4w, by = c("region", "year", "qtr"), suffix = c("_1w", "_4w")) %>%
  filter(!(year < 1999 | (year == 1999 & qtr < 3))) %>%
  filter(!(year > 2014 | (year == 2014 & qtr > 4))) %>%
  filter(!is.na(ur_1w), !is.na(ur_4w))

# 지역별 회귀 추정
bridge_params <- bridge_train %>%
  group_by(region) %>%
  do({
    fit <- lm(ur_4w ~ ur_1w, data = .)
    tibble(a = coef(fit)[1], b = coef(fit)[2])
  }) %>%
  ungroup()

# ─────────────────────────────────────────────
# 4️⃣ 4주 기준 완전 시계열 생성
# ─────────────────────────────────────────────
unemp_full <- unemp_1w %>%
  full_join(unemp_4w, by = c("region", "year", "qtr"), suffix = c("_1w", "_4w")) %>%
  left_join(bridge_params, by = "region") %>%
  mutate(
    # 4w 데이터 생성 규칙
    ur_final = case_when(
      # ① 1999.2/4 이전: 4w 결측 → 회귀 예측
      year < 1999 | (year == 1999 & qtr < 3) ~ a + b * ur_1w,
      # ② 1999.3/4 이후: 4w 실측 있으면 그대로 사용
      !is.na(ur_4w) ~ ur_4w,
      # ③ 4w 결측이지만 1w는 존재하는 경우 (예외적)
      is.na(ur_4w) & !is.na(ur_1w) ~ a + b * ur_1w,
      TRUE ~ NA_real_
    ),
    source = case_when(
      year < 1999 | (year == 1999 & qtr < 3) ~ "predicted_from_1w",
      !is.na(ur_4w) ~ "actual_4w",
      is.na(ur_4w) & !is.na(ur_1w) ~ "bridged_after_2014",
      TRUE ~ "missing"
    )
  ) %>%
  arrange(region, year, qtr) %>%
  ungroup()

# ─────────────────────────────────────────────
# 5️⃣ 연평균 계산 (단순평균)
# ─────────────────────────────────────────────
unemp_annual <- unemp_full %>%
  group_by(region, year) %>%
  summarise(unemp_year = mean(ur_final, na.rm = TRUE), .groups = "drop")


# ─────────────────────────────────────────────
# 6️⃣ 연평균 계산
# ─────────────────────────────────────────────
unemp_annual <- unemp_full %>%
  group_by(region, year) %>%
  summarise(unemp_year = mean(ur_final, na.rm = TRUE), .groups = "drop")

# ─────────────────────────────────────────────
# 지역명 영문 변환 매핑
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
  "세종특별자치시" = "Sejong-si",
  # 통합 지역 추가
  "경상남도_통합" = "Gyeongsangnam-do",
  "충청남도_통합" = "Chungcheongnam-do"
)

# ─────────────────────────────────────────────
# 영문 변환 및 저장
# ─────────────────────────────────────────────
unemp_annual_en <- unemp_annual %>%
  mutate(region_en = cluster_map[region]) %>%
  relocate(region_en, .after = region)

# 결과 미리보기
head(unemp_annual_en)

# 파일 저장
write_csv(unemp_annual_en, "unemp_annual_combined_2014on_en.csv", append = FALSE)



# ─────────────────────────────────────────────
# 2️⃣ pe_num별 연도 매핑
# ─────────────────────────────────────────────
period_map <- tribble(
  ~pe_num, ~year4, ~year5,
  15, 1996, 1997,
  16, 2001, 2002,
  17, 2006, 2007,
  18, 2011, 2012,
  19, 2015, 2016,
  20, 2020, 2021,
  21, 2023, 2024
)

# ─────────────────────────────────────────────
# 3️⃣ 매핑 기반으로 unemp4, unemp5 생성
# ─────────────────────────────────────────────
unemp_wide <- period_map %>%
  left_join(
    unemp_annual_en %>% rename(cluster = region_en, unemp4 = unemp_year),
    by = c("year4" = "year")
  ) %>%
  left_join(
    unemp_annual_en %>% rename(cluster = region_en, unemp5 = unemp_year),
    by = c("year5" = "year", "cluster" = "cluster")
  ) %>%
  select(pe_num, cluster, unemp4, unemp5) %>%
  filter(!is.na(cluster))   # ← ★ 여기서 NA 지역 완전 제거!

# ─────────────────────────────────────────────
# 4️⃣ 지역 순서 맞추기
# ─────────────────────────────────────────────
region_order <- c(
  "Seoul", "Gyeonggi-do", "Incheon", "Daejeon",
  "Chungcheongbuk-do", "Chungcheongnam-do",
  "Gwangju", "Jeollabuk-do", "Jeollanam-do",
  "Daegu", "Gyeongsangbuk-do", "Busan",
  "Gyeongsangnam-do", "Gangwon-do", "Jeju-do"
)

unemp_wide <- unemp_wide %>%
  mutate(cluster = factor(cluster, levels = region_order)) %>%
  arrange(pe_num, cluster)

# ─────────────────────────────────────────────
# 5️⃣ 결과 저장
# ─────────────────────────────────────────────
write.csv(unemp_wide, "unemp_panel_by_pe.csv", row.names = FALSE)


