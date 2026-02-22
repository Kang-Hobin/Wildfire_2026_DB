## =========================================================
## Wildfire_9_Annual_Lag_Panel.R
## 목적: 보험료 산정을 위한 연간 시차(Lag-1) 패널 데이터 구축
## 구조: t-1년 기상 요약 지표 -> t년 산불 손해액 예측
## =========================================================

library(tidyverse)
library(lubridate)

# 1. 데이터 로드
# ---------------------------------------------------------
# [R7]에서 생성된 '최종 무결성' 데이터를 불러옵니다.
df_cleaned <- readRDS("Wildfire_Data/processed_data/wildfire_model_ready.rds")

# 2. 연간 데이터 요약 (Annual Aggregation)
# ---------------------------------------------------------
message(">>> 일일 데이터를 보험 수리적 연간 지표로 요약 중...")

df_yearly <- df_cleaned %>%
  mutate(year = year(date)) %>%
  group_by(sigungu_cd, sigungu_nm, year) %>%
  summarise(
    # [Outcomes] t년의 결과 (보험금 지급의 근거)
    fire_cnt_year      = sum(fire_cnt, na.rm = TRUE),
    damage_total_year  = sum(damage_area, na.rm = TRUE),
    fire_any_year      = as.integer(any(fire_yn == 1)),
    
    # [Features] t-1년의 위험 지표 (보험료 산정의 근거)
    # 지상(Surface) 요약
    ta_sfc_mean    = mean(sfc_ta, na.rm = TRUE),
    hm_sfc_min     = min(sfc_hm, na.rm = TRUE),   # 최저 습도가 낮을수록 위험
    ws_sfc_max     = max(sfc_ws, na.rm = TRUE),   # 최대 풍속이 높을수록 위험
    prcp_sfc_sum   = sum(sfc_prcp, na.rm = TRUE),
    he_sfc_mean    = mean(sfc_he, na.rm = TRUE),  # 실효습도 연간 평균
    
    # 산악(Mountain) 요약
    ta_mtw_mean    = mean(mtw_ta, na.rm = TRUE),
    hm_mtw_min     = min(mtw_hm, na.rm = TRUE),
    ws_mtw_max     = max(mtw_ws, na.rm = TRUE),
    prcp_mtw_sum   = sum(mtw_prcp, na.rm = TRUE),
    he_mtw_mean    = mean(mtw_he, na.rm = TRUE),
    
    .groups = "drop"
  )

# 3. 1년 시차(Lag-1) 변수 생성
# ---------------------------------------------------------
message(">>> t-1년 기상 변수를 t년 리스크와 매칭 중 (Lag-1)...")

df_annual_lagged <- df_yearly %>%
  arrange(sigungu_cd, year) %>%
  group_by(sigungu_cd) %>%
  mutate(across(
    # 요약된 기상 변수들에 대해서만 시차 적용
    .cols = c(starts_with("ta_"), starts_with("hm_"), 
              starts_with("ws_"), starts_with("prcp_"), starts_with("he_")),
    .fns = ~ lag(.x, 1),
    .names = "{.col}_lag1"
  )) %>%
  ungroup() %>%
  # 시차 데이터가 없는 첫해(예: 2017년) 데이터는 분석에서 제외
  filter(!is.na(ta_sfc_mean_lag1))

# 4. 최종 저장
# ---------------------------------------------------------
saveRDS(df_annual_lagged, "Wildfire_Data/processed_data/df_annual_insurance_panel.rds")

message(">>> [완료] 보험 분석용 연간 패널 생성 완료!")
message(">>> 최종 차원: ", nrow(df_annual_lagged), " 행 x ", ncol(df_annual_lagged), " 열")