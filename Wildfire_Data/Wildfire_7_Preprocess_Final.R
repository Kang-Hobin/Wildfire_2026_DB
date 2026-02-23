## =========================================================
## Wildfire_7_Preprocess_Final.R
## 목적: 결측치 처리 및 파생변수 생성
## =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(zoo)
})

# 1. 데이터 로드
# ---------------------------------------------------------
df_raw <- read_rds("Wildfire_Data/processed_data/final_wildfire_panel.rds")

# 2. 초기 스크리닝 (이상치 제거)
# ---------------------------------------------------------
message(">>> 1단계: 이상치 스크리닝 및 논리적 결측치 처리...")

df_cleaned <- df_raw %>%
  arrange(sigungu_cd, date) %>%
  mutate(
    across(c(sfc_prcp, mtw_prcp), ~replace_na(.x, 0)),
    damage_area = if_else(is.na(damage_area) & fire_cnt == 0, 0, damage_area),
    across(contains("hm"), ~if_else(.x < 0 | .x > 100, NA_real_, .x)),
    across(contains("ta"), ~if_else(.x < -45 | .x > 50, NA_real_, .x)),
    across(contains("ws"), ~if_else(.x < 0 | .x > 60, NA_real_, .x))
  )

# 3. MTW 가용 구간 트리밍
# ---------------------------------------------------------
message(">>> 2단계: MTW 관측 가용 구간 트리밍...")

df_trimmed <- df_cleaned %>%
  group_by(sigungu_cd) %>%
  mutate(
    has_mtw = !is.na(mtw_ta),
    mtw_start = if(any(has_mtw)) min(date[has_mtw]) else as.Date(NA),
    mtw_end   = if(any(has_mtw)) max(date[has_mtw]) else as.Date(NA)
  ) %>%
  filter(!is.na(mtw_start) & date >= mtw_start & date <= mtw_end) %>%
  ungroup() %>%
  dplyr::select(-has_mtw, -mtw_start, -mtw_end)

# 4. [핵심] 3중 계층적 결측치 복구 (Triple-Guard Imputation)
# ---------------------------------------------------------
message(">>> 3단계: 계층적 보간 엔진 가동 (SFC-MTW 상호 보정)...")

eps <- 0.001 

df_imputed <- df_trimmed %>%
  group_by(sigungu_cd) %>%
  arrange(date) %>%
  mutate(
    # [Layer 1] 선형 보간: 짧은 결측(1~3일)을 주변 값으로 메움
    across(c(sfc_ta, sfc_hm, sfc_ws, mtw_ta, mtw_hm, mtw_ws), 
           ~ zoo::na.approx(.x, na.rm = FALSE, maxgap = 3)),
    
    # [Layer 2] 지상-산악 상관성 보정 (Local Ratio)
    ws_ratio = (mean(mtw_ws, na.rm = TRUE) + eps) / (mean(sfc_ws, na.rm = TRUE) + eps),
    hm_ratio = (mean(mtw_hm, na.rm = TRUE) + eps) / (mean(sfc_hm, na.rm = TRUE) + eps),
    ta_diff  = mean(mtw_ta - sfc_ta, na.rm = TRUE),
    
    mtw_ta = if_else(is.na(mtw_ta), sfc_ta + ta_diff, mtw_ta),
    mtw_hm = if_else(is.na(mtw_hm), sfc_hm * hm_ratio, mtw_hm),
    mtw_ws = if_else(is.na(mtw_ws), sfc_ws * ws_ratio, mtw_ws)
  ) %>%
  ungroup() %>%
  
  # [Layer 3] 시군구-월별 평균 대치 (Sigungu-Specific Seasonal Mean)
  # 특정 시군구 데이터가 한 달 내내 비었을 때를 대비
  group_by(sigungu_cd, month = month(date)) %>%
  mutate(across(c(starts_with("sfc_"), starts_with("mtw_")), 
                ~ replace_na(.x, mean(.x, na.rm = TRUE)))) %>%
  ungroup() %>%
  
  # [Layer 4] 최후의 보루: 전국-월별 평균 대치 (Global Seasonal Mean)
  # 속초나 춘천처럼 특정 시기에 데이터가 전멸한 경우를 완벽 방어
  group_by(month = month(date)) %>%
  mutate(across(c(starts_with("sfc_"), starts_with("mtw_")), 
                ~ replace_na(.x, mean(.x, na.rm = TRUE)))) %>%
  ungroup()

# 5. 파생 변수 생성 (Feature Engineering)
# ---------------------------------------------------------
message(">>> 4단계: 실효습도($H_e$) 및 파생 변수 생성...")

# [Safeguard] He 함수 개선: 초기값이 NA여도 터지지 않게 설계
calc_he_safe <- function(hm, r = 0.7) {
  res <- numeric(length(hm))
  # 초기값이 NA면 보수적으로 65% 설정
  res[1] <- if(is.na(hm[1])) 65 else hm[1]
  for(i in 2:length(hm)) {
    res[i] <- (1-r)*hm[i] + r*res[i-1]
  }
  return(res)
}

df_final_ready <- df_imputed %>%
  group_by(sigungu_cd) %>%
  arrange(date) %>%
  mutate(
    is_dry = if_else(sfc_prcp < 1.0, TRUE, FALSE),
    sfc_he = calc_he_safe(sfc_hm),
    mtw_he = calc_he_safe(mtw_hm),
    sfc_prcp_7d = rollsum(sfc_prcp, k = 7, fill = 0, align = "right"),
    mtw_prcp_7d = rollsum(mtw_prcp, k = 7, fill = 0, align = "right"),
    fire_yn = if_else(fire_cnt > 0, 1, 0)
  ) %>%
  ungroup()

# 6. 최종 무결성 체크 및 저장
# ---------------------------------------------------------
# 보간 후에도 NA가 남았는지 최종 확인 (이제 0이어야 함)
final_na_count <- sum(is.na(df_final_ready %>% dplyr::select(starts_with("sfc_"), starts_with("mtw_"))))

if(final_na_count == 0) {
  message(">>> [성공] 모든 결측치가 복구되었습니다.")
  saveRDS(df_final_ready, "Wildfire_Data/processed_data/df_wildfire_final_cleaned.rds")
  write_csv(df_final_ready, "Wildfire_Data/processed_data/df_wildfire_final_cleaned.csv")
  message(">>> 생성 완료: df_wildfire_final_cleaned.rds")
} else {
  message(">>> [경고] 잔여 결측치 ", final_na_count, "개가 발견되었습니다. 확인이 필요합니다.")
}