## =========================================================
## Wildfire_9_Annual_Lag_Panel.R 
## 목적: safe_min/max 적용을 통한 수치적 안정성 확보
## =========================================================

df_cleaned <- readRDS("Wildfire_Data/processed_data/df_wildfire_final_cleaned.rds")

# 0. 보조 함수 선언 (근본 함수들)
# ---------------------------------------------------------
safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
}

safe_min <- function(x) {
  valid_x <- x[!is.na(x) & is.finite(x)]
  if (length(valid_x) == 0) return(NA_real_)
  min(valid_x)
}

safe_max <- function(x) {
  valid_x <- x[!is.na(x) & is.finite(x)]
  if (length(valid_x) == 0) return(NA_real_)
  max(valid_x)
}

make_run_max <- function(x) {
  r <- rle(x)
  if (length(r$lengths) == 0) return(0L)
  max(r$lengths[r$values], 0L)
}

# 1. 통합 데이터 생성
# ---------------------------------------------------------
message(">>> [안전 모드] 일일 데이터를 계절별/연간 리스크 지표로 통합 요약 중...")

df_yearly <- df_cleaned %>%
  mutate(
    year = year(date),
    month = month(date)
  ) %>%
  group_by(sigungu_cd, sigungu_nm, year) %>%
  summarise(
    # --- [1] Outcomes (t년 결과) ---
    fire_cnt_year      = sum(fire_cnt, na.rm = TRUE),
    damage_total_year  = sum(damage_area, na.rm = TRUE),
    fire_any_year      = as.integer(any(fire_yn == 1)),
    
    # --- [2] General Weather (Safe 함수로 전면 교체) ---
    ta_sfc_mean    = mean(sfc_ta, na.rm = TRUE),
    hm_sfc_min     = safe_min(sfc_hm),  # min -> safe_min
    ws_sfc_max     = safe_max(sfc_ws),  # max -> safe_max
    prcp_sfc_sum   = sum(sfc_prcp, na.rm = TRUE),
    he_sfc_mean    = mean(sfc_he, na.rm = TRUE),
    
    ta_mtw_mean    = mean(mtw_ta, na.rm = TRUE),
    hm_mtw_min     = safe_min(mtw_hm),  # min -> safe_min
    ws_mtw_max     = safe_max(mtw_ws),  # max -> safe_max
    prcp_mtw_sum   = sum(mtw_prcp, na.rm = TRUE),
    he_mtw_mean    = mean(mtw_he, na.rm = TRUE),
    
    # --- [3] Seasonal Features (계절별 임계치) ---
    spring_dry_days  = sum(is_dry & month %in% 3:5, na.rm = TRUE),
    spring_ws_p95    = safe_quantile(mtw_ws[month %in% 3:5], 0.95),
    spring_he_min    = safe_min(mtw_he[month %in% 3:5]), # 경고 발생 지점 완전 방어
    
    monsoon_prcp_sum = sum(sfc_prcp[month %in% 6:8], na.rm = TRUE),
    winter_dry_days  = sum(is_dry & month %in% c(12,1,2), na.rm = TRUE),
    
    # --- [4] Continuity ---
    max_consecutive_dry_days = make_run_max(is_dry),
    max_dry_run_spring       = make_run_max(is_dry & month %in% 3:5),
    
    .groups = "drop"
  )

# 2. Lag-1 생성 및 최종 무한대 검역
# ---------------------------------------------------------
feature_cols <- setdiff(names(df_yearly), c("sigungu_cd", "sigungu_nm", "year", "fire_cnt_year", "damage_total_year", "fire_any_year"))

df_annual_lagged <- df_yearly %>%
  arrange(sigungu_cd, year) %>%
  group_by(sigungu_cd) %>%
  mutate(across(all_of(feature_cols), ~ lag(.x, 1), .names = "{.col}_lag1")) %>%
  ungroup() %>%
  # [최종 방어선] 시차 적용 후 발생하는 NA/Inf를 물리적 한계값 내의 NA로 정제
  mutate(across(contains("_lag1"), ~ if_else(is.finite(.x), .x, NA_real_))) %>%
  filter(!is.na(ta_sfc_mean_lag1))

saveRDS(df_annual_lagged, "Wildfire_Data/processed_data/df_annual_insurance_panel.rds")
message(">>> [완료] 무결성 패널 생성 성공! (Inf 개수: ", sum(is.infinite(as.matrix(df_annual_lagged))), ")")