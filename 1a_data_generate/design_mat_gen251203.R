## design matrix generation - for fundamentals model
## X: national-level predictors
## Z: province-level predictors
## this file should be conducted by `save_RData.R` file
## 업데이트 : 울산, 경남 분리
## 업데이트 : 울산 실업률 94,95년 없음 -> 경남 실업률로 넣었음
## 업데이트 : Z 변수에서 home, homeground를 각 후보별로 분리함.

# rm(list=ls())
# load("pe21.RData")

library(readxl)

P <- length(cluster)  # num. of provinces
TT <- length(15:21)   # num. of elections

##============================================================================##
## X: national-level predictors (ncol=length(15:20)*P=90)
##============================================================================##

# is current president democ? from 14th to 20th
is_current <- c(-1,1,1,-1,-1,1,-1) 
# is the incumbent been in the blue house for twice?
twice1 <- c(1,0,1,0,1,0,0)*is_current # yes/no
twice2 <- c(2,1,2,1,2,1,1)*is_current # how many times
# approval rate of the current president at the second quarter of 5th yr
# from Gallup Korea
# (4th quarter of the 4th year for Park Geun-hye, who got impeached)
# (3rd quarter of the 3rd year for Yoon Suk-yeol, who was impeached)
approval <- c(7,26,24,25,12,39,18)*is_current            # positive
net_approval <- c(-67,-27,-42,-33,-68,-14,-55)*is_current # positive-negative
# national poll for the Democrats 100 days before the election 
# (60 days for the 19th PE)
# WE DON'T USE THIS VARIABLE
# poll <- c(29.9/(29.9+18.3),20.4/(20.4+30.2),7.2/(7.2+60.7),
#           15/(15+40),32/(32+1),36/(36+36))
# real GDP and GNI growth of the 4th and 5th year from KOSIS
# https://kosis.kr/statHtml/statHtml.do?orgId=301&tblId=DT_200Y006&conn_path=I2
# gdp4 <- c(7.9,4.9,5.3,3.7,2.8,-0.9)*is_current
# gdp5 <- c(6.2,7.7,5.8,2.4,2.9, 4.0)*is_current
# gni4 <- c(6.7,4.0,4.0,1.6,6.3,-0.2)*is_current
# gni5 <- c(4.0,8.6,5.7,2.9,4.4, 3.5)*is_current
gdp4 <- c(8.0, 4.7, 5.2, 3.7, 2.9, -0.7, 1.6)*is_current
gdp5 <- c(6.3, 7.7, 5.8, 2.5, 3.2, 4.6, 2.0)*is_current
gni4 <- c(6.8, 3.9, 4.0, 1.7, 6.5, 0.1, 2.0)*is_current
gni5 <- c(4.3, 8.8, 5.6, 3.0, 4.6, 3.8, 3.9)*is_current
# primary result
# D_primary <- c(0.78,0.72,0.56,0.72,0.73,0.56)
# P_primary <- c(0.60,0.79,0.51,0.91,0.74,0.54)
D_primary <- c(0.78, 0.72, 0.44, 0.57, 0.57, 0.50, 0.90)
P_primary <- c(0.60, 0.68, 0.50, 0.84, 0.47, 0.48, 0.57)
# impeachment by the national congress?
impeach <- c(0,0,1,0,1,0,1)*is_current

# make matrix X
X_temp <- cbind(is_current,twice1,twice2,approval,net_approval,#poll,
                gdp4,gdp5,gni4,gni5,D_primary,P_primary,impeach)
X <- matrix(rep(X_temp,each=P),nrow=P*length(15:21))
colnames(X) <- colnames(X_temp)
X

##============================================================================##
## Z: province-level predictors (ncol=length(15:21)*P=90)
##============================================================================##

Z_penum <- rep(15:21,each=P)
Z_cluster <- rep(names(cluster),TT)
length(Z_penum)
length(Z_cluster)

home_dem      <- rep(0, P * TT)
home_con      <- rep(0, P * TT)
homeground_dem <- rep(0, P * TT)
homeground_con <- rep(0, P * TT)

## 각 선거별 후보 홈/근거지 정보 리스트화 (예시: 지금 코드 그대로 옮김)

home_dem_list <- list(
  "15" = "Jeollanam-do",       # 김대중
  "16" = "Gyeongsangnam-do",   # 노무현
  "17" = "Jeollabuk-do",       # 정동영
  "18" = "Gyeongsangnam-do",   # 문재인
  "19" = "Gyeongsangnam-do",   # 문재인
  "20" = "Gyeongsangbuk-do",   # (현재 코드 기준)
  "21" = "Gyeongsangbuk-do"    # (현재 코드 기준)
)

home_con_list <- list(
  "15" = "Chungcheongnam-do",  # 이회창
  "16" = "Chungcheongnam-do",  # 이회창
  "17" = "Gyeongsangbuk-do",   # 이명박
  "18" = "Daegu",              # 박근혜
  "19" = "Gyeongsangnam-do",   # 홍준표
  "20" = "Seoul",              # 윤석열
  "21" = "Gyeongsangbuk-do"    # 김문수
)

homeground_dem_list <- list(
  "20" = "Gyeonggi-do",  # 이재명 경기지사
  "21" = "Gyeonggi-do"   # 이재명 경기지사
)

homeground_con_list <- list(
  "17" = "Seoul",            # 이명박 서울시장 
  "19" = "Gyeongsangnam-do", # 홍준표 경남도지사
  "21" = "Gyeonggi-do"       # 김문수 경기지사
)

for (tt in 15:21) {
  key <- as.character(tt)
  
  # 민주당 홈
  if (!is.null(home_dem_list[[key]])) {
    home_dem[Z_penum == tt & Z_cluster %in% home_dem_list[[key]]] <- 1
  }
  
  # 보수 홈
  if (!is.null(home_con_list[[key]])) {
    home_con[Z_penum == tt & Z_cluster %in% home_con_list[[key]]] <- -1
  }
  
  # 민주당 근거지
  if (!is.null(homeground_dem_list[[key]])) {
    homeground_dem[Z_penum == tt & Z_cluster %in% homeground_dem_list[[key]]] <- 1
  }
  
  # 보수 근거지
  if (!is.null(homeground_con_list[[key]])) {
    homeground_con[Z_penum == tt & Z_cluster %in% homeground_con_list[[key]]] <- -1
  }
}


# 1) 인덱스 기반 구조
idx <- tibble(
  cluster = Z_cluster,
  pe_num  = Z_penum
)


## --- CPI 붙이기 ------------------------------------------------
cpi <- readr::read_csv("0_data/csi21.csv", show_col_types = FALSE)
# cpi에는 최소한 cluster, pe_num, csi4, csi5, is_current 열이 있어야 함

idx_cpi <- idx %>%
  dplyr::left_join(cpi, by = c("cluster", "pe_num"))

# 혹시 CPI 매칭 안 되는 조합이 있는지 확인
if (any(is.na(idx_cpi$csi4) | is.na(idx_cpi$csi5))) {
  warning("⚠ CPI 데이터에서 NA가 있습니다. (cluster, pe_num) 매칭을 확인하세요.")
}

cpi4 <- idx_cpi$csi4 * is_current
cpi5 <- idx_cpi$csi5 * is_current

## --- 실업률 (unemp_panel_by_pe) 붙이기 ------------------------
unemp45 <- readr::read_csv("0_data/unemp_panel_by_pe.csv", show_col_types = FALSE)
# 여기에는 pe_num, cluster, unemp4, unemp5 가 있어야 함

idx_unemp <- idx %>%
  dplyr::left_join(unemp45, by = c("cluster", "pe_num"))

# 울산 1996/97 보정이 잘 먹었는지 확인 (NA가 없어야 함)
if (any(is.na(idx_unemp$unemp4) | is.na(idx_unemp$unemp5))) {
  warning("⚠ 실업률 unemp4/unemp5에 NA가 남아 있습니다. 울산 96/97 보정 여부를 확인하세요.")
}

unemp4 <- idx_unemp$unemp4 * is_current
unemp5 <- idx_unemp$unemp5 * is_current

## --- 최종 Z 매트릭스 ------------------------------------------
Z <- cbind(home_dem, homeground_dem, home_con, homeground_con, cpi4, cpi5, unemp4, unemp5)

dim(Z)
summary(Z)
# 
# # CPI growth
# cpi <- read.csv("0_data/csi21.csv")
# all(cpi$pe_num == Z_penum)
# all(cpi$cluster == Z_cluster)
# cpi4 <- cpi$csi4 * is_current
# cpi5 <- cpi$csi5 * is_current
# # unemployment rate
# ## 임기 4,5년차 실업률
# ## 구직기간 1주 기준으로 작성된 자료는 2014년까지 제공되고,
# ## 구직기간 4주 기준으로 작성된 자료는 99~부터 제공되어 
# ## 두 자료의 분기별 실업률을 단순선형회귀로 보정하였음.
# ## 이후 분기별 값을 단순 평균내어 시도별 4년차 z, 5년차 z 생성.
# ## 참고로 17년 선거부터는 15~16년으로 두 해 전 자료 사용하고, 
# ## 그 전 선거까지는 (예. 12년) 11~12년 선거 년도 포함한 두 해 자료 사용.
# unemp45 <- data.frame(read_excel("0_data/unemp_panel_by_pe.csv"))
# all(unemp45$pe_num == Z_penum)
# all(unemp45$cluster == Z_cluster)
# unemp4 <- unemp45$unemp4 * is_current
# unemp5 <- unemp45$unemp5 * is_current
# 
# # make matrix Z
# Z <- cbind(home,homeground,cpi4,cpi5,unemp4,unemp5)
# Z


