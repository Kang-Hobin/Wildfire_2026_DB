## =========================================================
## Wildfire_6_Final_Join.R
## =========================================================

# 1. 산불 데이터 집계 (먼저 수행)
# ---------------------------------------------------------
message(">>> 1. 산불 데이터 집계 중...")
# sgg_lookup이 사전에 정의되어 있어야 합니다 (4번 단계에서 만든 것)
fire_daily_agg <- fire_raw %>%
  mutate(
    date = as.Date(paste(startyear, startmonth, startday, sep="-")),
    damagearea_num = as.numeric(damagearea)
  ) %>%
  left_join(sgg_lookup, by = c("locsi" = "sido_nm", "locgungu" = "sgg_base_nm")) %>% 
  group_by(sigungu_cd, date) %>%
  summarise(
    fire_cnt = n(),
    damage_area = sum(damagearea_num, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(sigungu_cd))

# 2. 메인 패널 초기화 (산불 데이터 결합)
# ---------------------------------------------------------
message(">>> 2. 메인 패널 초기화...")
final_panel <- sgg_base_map %>%
  expand_grid(date = ALL_DATES) %>%
  left_join(fire_daily_agg, by = c("sigungu_cd", "date")) %>%
  mutate(
    fire_cnt = replace_na(fire_cnt, 0),
    damage_area = replace_na(damage_area, 0)
  )

# 3. ASOS/AWS 지상 기상 Waterfall 조인
# ---------------------------------------------------------
message(">>> 3. 지상 기상(ASOS/AWS) 조인 중...")
for (i in 1:3) {
  id_col_a <- paste0("match_ASOS_id_", i)
  final_panel <- final_panel %>%
    left_join(asos_clean %>% rename_with(~paste0(.x, "_A", i), -c(stn_id, date)), 
              by = c("date", setNames("stn_id", id_col_a)))
  
  id_col_w <- paste0("match_AWS_id_", i)
  final_panel <- final_panel %>%
    left_join(aws_clean %>% rename_with(~paste0(.x, "_W", i), -c(stn_id, date)), 
              by = c("date", setNames("stn_id", id_col_w)))
}

# 4. MTW 산악 기상 매핑 및 조인
# ---------------------------------------------------------
message(">>> 4. 산악 기상(MTW) 조인 중...")
mtw_map_fixed <- mtw_annual_map %>%
  mutate(across(starts_with("match_MTW_id_"), ~ {
    val <- as.character(.)
    if_else(str_starts(val, "M"), val, paste0("M", val))
  })) %>%
  dplyr::select(sigungu_cd, map_year, starts_with("match_MTW_"))

final_panel <- final_panel %>%
  mutate(map_year = case_when(
    year(date) < 2018 ~ 2018L, year(date) > 2022 ~ 2022L,
    TRUE ~ as.integer(year(date))
  )) %>%
  left_join(mtw_map_fixed, by = c("sigungu_cd", "map_year"))

for (i in 1:3) {
  id_col <- paste0("match_MTW_id_", i)
  final_panel <- final_panel %>%
    left_join(mtw_clean %>% rename_with(~ paste0(.x, "_M", i), -c(stn_id, date)), 
              by = c("date", setNames("stn_id", id_col)))
}

# 5. [핵심] 변수 최종 확정
# ---------------------------------------------------------
message(">>> 5. 변수 통합 및 최종 패널 생성 ...")

final_panel <- final_panel %>%
  mutate(
    # [지상 기상] AWS 1~3순위 -> ASOS 1~3순위 순차 적용
    sfc_ta = coalesce(ta_avg_W1, ta_avg_W2, ta_avg_W3, ta_avg_A1, ta_avg_A2, ta_avg_A3),
    sfc_hm = coalesce(hm_avg_A1, hm_avg_A2, hm_avg_A3),
    sfc_ws = coalesce(ws_max_W1, ws_max_W2, ws_max_W3, ws_max_A1, ws_max_A2, ws_max_A3),
    
    # [수정] 누락되었던 지상 강수량 변수 추가
    sfc_prcp = coalesce(prcp_day_W1, prcp_day_W2, prcp_day_W3, 
                        prcp_day_A1, prcp_day_A2, prcp_day_A3),
    
    # [산악 기상] MTW 1~3순위
    mtw_ta = coalesce(ta_avg_M1, ta_avg_M2, ta_avg_M3),
    mtw_hm = coalesce(hm_avg_M1, hm_avg_M2, hm_avg_M3),
    mtw_ws = coalesce(ws_max_M1, ws_max_M2, ws_max_M3),
    mtw_prcp = coalesce(prcp_day_M1, prcp_day_M2, prcp_day_M3)
  ) %>%
  # 분석 대상 컬럼 선택 (starts_with를 통해 sfc_prcp도 자동으로 포함됩니다)
  dplyr::select(
    sigungu_cd, sigungu_nm, date, fire_cnt, damage_area,
    starts_with("sfc_"), starts_with("mtw_")
  )

# 6. 저장 및 최종 확인
# ---------------------------------------------------------
saveRDS(final_panel, file.path(PROC_DIR, "final_wildfire_panel.rds"))

message(">>> 현재 패널 컬럼 구성:")
print(colnames(final_panel))
message(">>> [성공] 패널이 저장되었습니다.")