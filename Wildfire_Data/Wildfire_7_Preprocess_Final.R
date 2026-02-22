## =========================================================
## Wildfire_7_Preprocess_Final.R
## 목적: 결측치 완전 복구 및 분석용 핵심 변수($H_e$ 등) 생성
## =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(zoo)
})

# 1. 데이터 로드 및 변수명 퇴마(Exorcism) 작업
# ---------------------------------------------------------
df_raw <- read_rds("Wildfire_Data/processed_data/final_wildfire_panel.rds")
colnames(df_raw)
# 2. 초기 스크리닝 (Abnormal Screening)
# ---------------------------------------------------------
message(">>> 1단계: 이상치 스크리닝 및 논리적 결측치 처리...")

df_cleaned <- df_raw %>%
  arrange(sigungu_cd, date) %>%
  mutate(
    # [논리 보정] 강수량 NA는 기상청 원칙에 따라 0으로 채움
    across(c(sfc_prcp, mtw_prcp), ~replace_na(.x, 0)),
    damage_area = if_else(is.na(damage_area) & fire_cnt == 0, 0, damage_area),
    
    # [이상치 제거] 물리적 한계를 벗어나는 센서 오류값 제거
    across(contains("hm"), ~if_else(.x < 0 | .x > 100, NA_real_, .x)),
    across(contains("ta"), ~if_else(.x < -45 | .x > 50, NA_real_, .x)),
    across(contains("ws"), ~if_else(.x < 0 | .x > 60, NA_real_, .x))
  )

# 3. MTW 가용 구간 트리밍 (Availability Trimming)
# ---------------------------------------------------------
message(">>> 2단계: 시군구별 MTW 관측 가용 구간 트리밍...")

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

# 4. 상호 보완적 결측치 복구 (Cross-Imputation Engine)
# ---------------------------------------------------------
message(">>> 3단계: 지상-산악 상관성($r=0.95$)을 이용한 결측치 완전 복구...")

df_imputed <- df_trimmed %>%
  group_by(sigungu_cd) %>%
  mutate(
    # [1] 기온: 지상-산악 편차(Bias) 보정
    temp_diff = mean(mtw_ta - sfc_ta, na.rm = TRUE),
    mtw_ta = if_else(is.na(mtw_ta), sfc_ta + temp_diff, mtw_ta),
    
    # [2] 습도/풍속: 지상-산악 비율(Ratio) 보정
    hm_ratio = mean(mtw_hm / sfc_hm, na.rm = TRUE),
    mtw_hm = if_else(is.na(mtw_hm), sfc_hm * hm_ratio, mtw_hm),
    
    ws_ratio = mean(mtw_ws / sfc_ws, na.rm = TRUE),
    mtw_ws = if_else(is.na(mtw_ws), sfc_ws * ws_ratio, mtw_ws)
  ) %>%
  # [3] 최후의 수단: 시군구별/월별 평균값 (정적 보간)
  group_by(sigungu_cd, month = month(date)) %>%
  mutate(across(c(mtw_ta, mtw_hm, mtw_ws, sfc_ta, sfc_hm, sfc_ws), 
                ~replace_na(.x, mean(.x, na.rm = TRUE)))) %>%
  ungroup() %>%
  dplyr::select(-temp_diff, -hm_ratio, -ws_ratio, -month)

# 5. 파생 변수 생성 (Feature Engineering)
# ---------------------------------------------------------
message(">>> 4단계: 실효습도($H_e$) 및 시계열 변수 생성...")

# 실효습도 가중치 함수 (표준 r=0.7)
calc_he <- function(hm, r = 0.7) {
  res <- numeric(length(hm))
  res[1] <- if(is.na(hm[1])) 50 else hm[1] 
  for(i in 2:length(hm)) {
    res[i] <- if(is.na(hm[i])) res[i-1] else (1-r)*hm[i] + r*res[i-1]
  }
  return(res)
}

df_final_ready <- df_imputed %>%
  group_by(sigungu_cd) %>%
  arrange(date) %>%
  mutate(
    # [Feature 1] 실효습도: 이제 결측치가 없어 끊기지 않음
    sfc_he = calc_he(sfc_hm),
    mtw_he = calc_he(mtw_hm),
    
    # [Feature 2] 7일 누적 강수량: NA가 0으로 채워졌으므로 정확함
    sfc_prcp_7d = rollsum(sfc_prcp, k = 7, fill = 0, align = "right"),
    mtw_prcp_7d = rollsum(mtw_prcp, k = 7, fill = 0, align = "right"),
    
    # [Target] 종속 변수 (Binary)
    fire_yn = if_else(fire_cnt > 0, 1, 0)
  ) %>%
  ungroup()

# 6. 최종 품질 필터링 및 저장
# ---------------------------------------------------------
# 보간 후에도 데이터 품질이 극도로 낮은 곳이 있다면 퇴출
exclude_sgg <- df_final_ready %>%
  group_by(sigungu_cd) %>%
  summarise(na_rate = mean(is.na(mtw_ta))) %>%
  filter(na_rate > 0.05) %>% # 보간 후 결측 5% 초과 지점 (사실상 0개 예상)
  pull(sigungu_cd)

df_model_ready <- df_final_ready %>%
  filter(!(sigungu_cd %in% exclude_sgg))

message(">>> [최종 보고]")
message("- 제외 시군구: ", length(exclude_sgg), "개")
message("- 최종 분석 데이터 규모: ", nrow(df_model_ready), " 행")

# 최종 저장부
saveRDS(df_model_ready, file.path(PROC_DIR, "df_wildfire_final_cleaned.rds"))
write_csv(df_model_ready, file.path(PROC_DIR, "df_wildfire_final_cleaned.csv"))

message(">>> [완료] 전처리 파이프라인 완결!")
message(">>> 생성된 파일: df_wildfire_final_cleaned.rds")