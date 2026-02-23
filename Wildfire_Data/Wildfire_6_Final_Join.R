## =========================================================
## Wildfire_6_Final_Join.R
## =========================================================
library(tidyverse)
library(lubridate)

# 데이터 로드
fire_raw <- read_rds("Wildfire_Data/raw_data/raw_fire_API_2017_2025.rds")
sgg_points <- readRDS("Wildfire_Data/meta_data/sigungu_centroid_master.rds")
asos_clean <- readRDS("Wildfire_Data/processed_data/asos_daily.rds")
aws_clean  <- readRDS("Wildfire_Data/processed_data/aws_daily.rds")
mtw_clean  <- readRDS("Wildfire_Data/processed_data/mtw_raw_standard.rds")
sgg_stn_map <- readRDS("Wildfire_Data/meta_data/sgg_asos_aws_triple_map.rds")
sgg_mtw_map <- readRDS("Wildfire_Data/meta_data/sgg_mtw_annual_map_2018_2022.rds")
sgg_base_map <- sgg_points %>% dplyr::select(sigungu_cd, sigungu_nm)

# 1. 산불 데이터 집계 및 시군구 코드 매핑 (강력한 정규화)
# ---------------------------------------------------------
message(">>> 1. 산불 데이터 집계 및 지옥의 행정구역 매핑 중...")

# [마스터 정비] 매칭 확률을 높이기 위해 여러 버전의 키를 만듭니다.
sgg_points_refined <- sgg_points %>%
  mutate(
    # 정규 키: "강릉시" -> "강릉"
    key_std = str_remove(sgg_base_nm, "시$|군$|구$"),
    # 도시 키: "성남시 분당구" -> "성남" (구 단위가 없는 데이터 구출용)
    key_city = str_extract(sgg_base_nm, "^[^시군구 ]+"),
    # 풀네임 키: "성남 분당" 같은 공백 포함 데이터 대응
    key_full = str_replace_all(sgg_base_nm, "시 |구 |군 ", " ") %>% str_trim()
  )

fire_daily_agg <- fire_raw %>%
  mutate(
    date = as.Date(paste(startyear, startmonth, startday, sep="-")),
    damagearea_num = as.numeric(damagearea),
    locsi = as.character(locsi),
    locgungu = as.character(locgungu) %>% str_trim()
  ) %>%
  # [Step 1] 1차 조인: 정밀 매칭 (강릉 -> 강릉)
  left_join(sgg_points_refined, by = c("locsi" = "sido_nm", "locgungu" = "key_std")) %>%
  
  # [Step 2] 2차 보정: 도시명 기반 매칭 (성남 -> 성남시 수정구 중 하나로)
  # 구 단위가 생략된 대도시 데이터를 구출합니다.
  mutate(sigungu_cd = case_when(
    !is.na(sigungu_cd) ~ sigungu_cd,
    locgungu == "고양" ~ "41281", # 고양시 덕양구로 임시 할당
    locgungu == "용인" ~ "41461", # 용인시 처인구
    locgungu == "성남" ~ "41131", # 성남시 수정구
    locgungu == "수원" ~ "41111", # 수원시 장안구
    locgungu == "안산" ~ "41271", # 안산시 상록구
    locgungu == "안양" ~ "41171", # 안양시 만안구
    locgungu == "포항" ~ "37011", # 포항시 남구
    locgungu == "청주" ~ "33041", # 청주시 상당구
    locgungu == "창원" ~ "38110", # 창원시 의창구
    locgungu == "군위" ~ "27720", # 대구 군위군 (경북/대구 족보 무시 매칭)
    TRUE ~ sigungu_cd
  )) %>%
  
  # [Step 3] 3차 보정: 공백 포함 특수 케이스 (성남 분당 -> 성남시 분당구)
  mutate(locgungu_clean = str_replace_all(locgungu, " ", "")) %>%
  mutate(sigungu_cd = if_else(
    is.na(sigungu_cd),
    sgg_points_refined$sigungu_cd[match(locgungu_clean, str_replace_all(sgg_points_refined$key_std, " ", ""))],
    sigungu_cd
  )) %>%
  
  # 집계
  group_by(sigungu_cd, date) %>%
  summarise(
    fire_cnt = n(),
    damage_area = sum(damagearea_num, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(sigungu_cd))

# 최종 결과 확인
message(">>> [최종 검증] 총 피해 면적: ", sum(fire_daily_agg$damage_area), " ha")

# 2. 메인 패널 생성 (모든 날짜 x 모든 시군구)
ALL_DATES <- seq.Date(as.Date("2017-01-01"), as.Date("2025-12-31"), by="day")

final_panel <- sgg_base_map %>%
  expand_grid(date = ALL_DATES) %>%
  # [추가!] 여기서 매핑 정보를 붙여줘야 match_ASOS_id_1 컬럼들이 생깁니다.
  left_join(sgg_stn_map, by = "sigungu_cd") %>% 
  left_join(fire_daily_agg, by = c("sigungu_cd", "date")) %>%
  mutate(across(c(fire_cnt, damage_area), ~replace_na(.x, 0)))

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

# 4. MTW 산악 기상 매핑 및 조인 (복구 버전)
# ---------------------------------------------------------
message(">>> 4. 산악 기상(MTW) 매핑 정보 결합 중...")

# [핵심] date에서 year를 추출하고, 매핑 데이터가 있는 2018~2022 범위로 맞춥니다.
final_panel <- final_panel %>%
  mutate(map_year = lubridate::year(date)) %>%
  mutate(map_year = case_when(
    map_year < 2018 ~ 2018L,
    map_year > 2022 ~ 2022L,
    TRUE ~ as.integer(map_year)
  ))

# 이제 sgg_mtw_map(map_year 보유)과 조인이 가능해집니다.
final_panel <- final_panel %>%
  left_join(sgg_mtw_map, by = c("sigungu_cd", "map_year"))

# [M 접두사 보정] 만약 ID에 'M'이 없다면 붙여주는 작업 (동료 코드와 정합성)
final_panel <- final_panel %>%
  mutate(across(starts_with("match_MTW_id_"), ~ {
    val <- as.character(.x)
    if_else(str_starts(val, "M"), val, paste0("M", val))
  }))

# MTW Waterfall 조인 루프
# ---------------------------------------------------------
message(">>> 산악 기상 데이터 Waterfall 조인 수행 중...")
for (i in 1:3) {
  id_col <- paste0("match_MTW_id_", i)
  
  final_panel <- final_panel %>%
    left_join(mtw_clean %>% rename_with(~ paste0(.x, "_M", i), -c(stn_id, date)), 
              by = c("date", setNames("stn_id", id_col)))
}

# 5. 변수명 최종 확정 (Prefix 전략 준수)
message(">>> 5. 변수 통합 및 Prefix 적용...")
final_panel <- final_panel %>%
  mutate(
    # 지상 기상 통합
    sfc_ta = coalesce(ta_avg_W1, ta_avg_W2, ta_avg_W3, ta_avg_A1, ta_avg_A2, ta_avg_A3),
    sfc_hm = coalesce(hm_avg_A1, hm_avg_A2, hm_avg_A3),
    sfc_ws = coalesce(ws_max_W1, ws_max_W2, ws_max_W3, ws_max_A1, ws_max_A2, ws_max_A3),
    sfc_prcp = coalesce(prcp_day_W1, prcp_day_W2, prcp_day_W3, prcp_day_A1, prcp_day_A2, prcp_day_A3),
    
    # 산악 기상 통합
    mtw_ta = coalesce(ta_avg_M1, ta_avg_M2, ta_avg_M3),
    mtw_hm = coalesce(hm_avg_M1, hm_avg_M2, hm_avg_M3),
    mtw_ws = coalesce(ws_max_M1, ws_max_M2, ws_max_M3),
    mtw_prcp = coalesce(prcp_day_M1, prcp_day_M2, prcp_day_M3),
    
    # 발생 여부(Binary)
    fire_yn = if_else(fire_cnt > 0, 1, 0)
  ) %>%
  dplyr::select(sigungu_cd, sigungu_nm, date, fire_cnt, fire_yn, damage_area, 
         starts_with("sfc_"), starts_with("mtw_"))

# 최종 저장
saveRDS(final_panel, "Wildfire_Data/processed_data/final_wildfire_panel.rds")