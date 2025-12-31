## design matrix generation - for fundamentals model
## X: national-level predictors
## Z: province-level predictors
## this file should be conducted by `save_RData.R` file
## 업데이트 : 울산, 경남 분리
## 업데이트 : 울산 실업률 94,95년 없음 -> 경남 실업률로 넣었음
## 업데이트 : Z 변수에서 home, homeground를 각 후보별로 분리함.

## design matrix generation - for fundamentals model
## outputs: X_20, Z_20 (pe15~pe20) and X_21, Z_21 (pe15~pe21)

library(readr)
library(dplyr)
library(tibble)

P <- length(cluster)

# ────────────────────────────────────────────────
# 0) 공통 설정: 선거별(national) 공변량을 "pe 번호"로 이름 붙여 관리
# ────────────────────────────────────────────────
pe_all <- 15:21

is_current_by_pe   <- setNames(c(-1,  1,  1, -1, -1,  1, -1), pe_all)
twice1_base_by_pe  <- setNames(c( 1,  0,  1,  0,  1,  0,  0), pe_all)
twice2_base_by_pe  <- setNames(c( 2,  1,  2,  1,  2,  1,  1), pe_all)

approval_by_pe     <- setNames(c( 7, 26, 24, 25, 12, 39, 18), pe_all)
netappr_by_pe      <- setNames(c(-67,-27,-42,-33,-68,-14,-55), pe_all)

gdp4_by_pe         <- setNames(c( 8.0, 4.7, 5.2, 3.7, 2.9,-0.7, 1.6), pe_all)
gdp5_by_pe         <- setNames(c( 6.3, 7.7, 5.8, 2.5, 3.2, 4.6, 2.0), pe_all)
gni4_by_pe         <- setNames(c( 6.8, 3.9, 4.0, 1.7, 6.5, 0.1, 2.0), pe_all)
gni5_by_pe         <- setNames(c( 4.3, 8.8, 5.6, 3.0, 4.6, 3.8, 3.9), pe_all)

D_primary_by_pe    <- setNames(c(0.78,0.72,0.44,0.57,0.57,0.50,0.90), pe_all)
P_primary_by_pe    <- setNames(c(0.60,0.68,0.50,0.84,0.47,0.48,0.57), pe_all)

# impeach_base_by_pe <- setNames(c(0,0,1,0,1,0,1), pe_all)
impeach_base_by_pe <- setNames(c(0,0,-1,0,1,0,1), pe_all)

# ────────────────────────────────────────────────
# 1) X 생성 함수: pe 범위에 따라 pe 단위 행렬 -> P 배 확장
# ────────────────────────────────────────────────
make_X <- function(pe_nums, P){
  pe_chr <- as.character(pe_nums)
  
  is_current   <- as.numeric(is_current_by_pe[pe_chr])
  
  twice1 <- as.numeric(twice1_base_by_pe[pe_chr]) * is_current
  twice2 <- as.numeric(twice2_base_by_pe[pe_chr]) * is_current
  
  approval <- as.numeric(approval_by_pe[pe_chr]) * is_current
  net_approval <- as.numeric(netappr_by_pe[pe_chr]) * is_current
  
  gdp4 <- as.numeric(gdp4_by_pe[pe_chr]) * is_current
  gdp5 <- as.numeric(gdp5_by_pe[pe_chr]) * is_current
  gni4 <- as.numeric(gni4_by_pe[pe_chr]) * is_current
  gni5 <- as.numeric(gni5_by_pe[pe_chr]) * is_current
  
  D_primary <- as.numeric(D_primary_by_pe[pe_chr])
  P_primary <- as.numeric(P_primary_by_pe[pe_chr])
  
  impeach <- as.numeric(impeach_base_by_pe[pe_chr]) * is_current
  
  X_temp <- cbind(
    is_current, twice1, twice2, approval, net_approval,
    gdp4, gdp5, gni4, gni5, D_primary, P_primary, impeach
  )
  colnames(X_temp) <- c("is_current","twice1","twice2","approval","net_approval",
                        "gdp4","gdp5","gni4","gni5","D_primary","P_primary","impeach")
  
  # election(행)마다 P개 지역을 반복해 N=P*TT 행 생성
  X <- matrix(rep(X_temp, each=P), nrow = P * length(pe_nums))
  colnames(X) <- colnames(X_temp)
  X
}

# ────────────────────────────────────────────────
# 2) Z 생성 함수: (cluster, pe_num) 인덱스 만들고 home/homeground + CPI + 실업률
# ────────────────────────────────────────────────
make_Z <- function(pe_nums, P, cluster){
  
  TT <- length(pe_nums)
  Z_penum    <- rep(pe_nums, each = P)
  Z_cluster  <- rep(names(cluster), times = TT)
  
  is_current_long <- rep(
    as.numeric(is_current_by_pe[as.character(pe_nums)]),
    each = P
  )
  
  ## ------------------------------------------------
  ## 1) home / homeground (모두 0/1)
  ## ------------------------------------------------
  home_dem        <- integer(P * TT)
  home_con        <- integer(P * TT)
  homeground_dem  <- integer(P * TT)
  homeground_con  <- integer(P * TT)
  
  home_dem_list <- list(
    "15"="Jeollanam-do",
    "16"="Gyeongsangnam-do",
    "17"="Jeollabuk-do",
    "18"="Gyeongsangnam-do",
    "19"="Gyeongsangnam-do",
    "20"="Gyeongsangbuk-do",
    "21"="Gyeongsangbuk-do"
  )
  
  home_con_list <- list(
    "15"="Chungcheongnam-do",
    "16"="Chungcheongnam-do",
    "17"="Gyeongsangbuk-do",
    "18"="Daegu",
    "19"="Gyeongsangnam-do",
    "20"="Seoul",
    "21"="Gyeongsangbuk-do"
  )
  
  homeground_dem_list <- list(
    "20"="Gyeonggi-do",
    "21"="Gyeonggi-do"
  )
  
  homeground_con_list <- list(
    "17"="Seoul",
    "19"="Gyeongsangnam-do",
    "21"="Gyeonggi-do"
  )
  
  for(tt in pe_nums){
    key <- as.character(tt)
    
    if(!is.null(home_dem_list[[key]])){
      home_dem[Z_penum == tt & Z_cluster %in% home_dem_list[[key]]] <- 1
    }
    if(!is.null(home_con_list[[key]])){
      home_con[Z_penum == tt & Z_cluster %in% home_con_list[[key]]] <- 1
    }
    if(!is.null(homeground_dem_list[[key]])){
      homeground_dem[Z_penum == tt & Z_cluster %in% homeground_dem_list[[key]]] <- 1
    }
    if(!is.null(homeground_con_list[[key]])){
      homeground_con[Z_penum == tt & Z_cluster %in% homeground_con_list[[key]]] <- 1
    }
  }
  
  ## ------------------------------------------------
  ## 2) 옵션 B: OR 통합 (연고 존재 여부)
  ## ------------------------------------------------
  dem_link <- as.integer(home_dem == 1 | homeground_dem == 1)
  con_link <- as.integer(home_con == 1 | homeground_con == 1)
  
  ## ------------------------------------------------
  ## 3) CPI / 실업률
  ## ------------------------------------------------
  idx <- tibble(cluster = Z_cluster, pe_num = Z_penum)
  
  cpi <- read_csv("0_data/csi21.csv", show_col_types = FALSE)
  idx_cpi <- idx %>% left_join(cpi, by = c("cluster","pe_num"))
  cpi4 <- idx_cpi$csi4 * is_current_long
  cpi5 <- idx_cpi$csi5 * is_current_long
  
  unemp45 <- read_csv("0_data/unemp_panel_by_pe.csv", show_col_types = FALSE)
  idx_unemp <- idx %>% left_join(unemp45, by = c("cluster","pe_num"))
  unemp4 <- idx_unemp$unemp4 * is_current_long
  unemp5 <- idx_unemp$unemp5 * is_current_long
  
  ## ------------------------------------------------
  ## 4) Z 구성 (home 통합 버전)
  ## ------------------------------------------------
  Z <- cbind(dem_link, con_link, cpi4, cpi5, unemp4, unemp5)
  colnames(Z) <- c("dem_link","con_link","cpi4","cpi5","unemp4","unemp5")
  
  Z
}


# ────────────────────────────────────────────────
# 3) 20대/21대용 생성 (run-fund-model에서 스위치로 선택 사용)
#    - 20대 예측: pe15~pe20 필요
#    - 21대 예측: pe15~pe21 필요
# ────────────────────────────────────────────────
X_20 <- make_X(15:20, P)
Z_20 <- make_Z(15:20, P, cluster)

X_21 <- make_X(15:21, P)
Z_21 <- make_Z(15:21, P, cluster)

# sanity check
stopifnot(nrow(X_20) == P*length(15:20), nrow(Z_20) == P*length(15:20))
stopifnot(nrow(X_21) == P*length(15:21), nrow(Z_21) == P*length(15:21))

cat("\n✅ design matrices created: X_20,Z_20 and X_21,Z_21\n")

