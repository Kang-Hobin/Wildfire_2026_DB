## =========================================================
## Wildfire_4_Meta_Geocoding.R
## 목적: 5자리 시군구 코드 및 접두사(A/W/M) 포함 관측소 좌표 마스터 생성
## =========================================================

library(sf)
library(dplyr)
library(readr)
library(stringr)

# [A] 시군구 중심점 마스터 (5자리 코드 보장)
# ---------------------------------------------------------
message(">>> 시군구 중심점 추출 및 5자리 규격화 중...")
sgg_path <- "Wildfire_Data/meta_data/시군구경계.shp"

if (file.exists(sgg_path)) {
  sgg_master <- st_read(sgg_path, quiet = TRUE) %>%
    st_point_on_surface() %>%
    st_transform(4326) %>%
    mutate(
      lon = st_coordinates(.)[,1],
      lat = st_coordinates(.)[,2]
    ) %>%
    st_drop_geometry() %>%
    transmute(
      # sprintf를 사용하여 앞의 0이 유실되지 않은 5자리 문자열 강제
      sigungu_cd = sprintf("%05d", as.integer(SIGUNGU_CD)),
      sigungu_nm = SIGUNGU_NM,
      sgg_lon = lon,
      sgg_lat = lat
    )
  
  saveRDS(sgg_master, "Wildfire_Data/meta_data/sigungu_centroid_master.rds")
}

# [B] 관측소 위치 마스터 (접두사 전략 적용)
# ---------------------------------------------------------
message(">>> 관측소 위치 마스터 생성 및 접두사(A/W/M) 부여 중...")
meta_path <- "Wildfire_Data/meta_data/Meta_관측지점정보.csv"
mtw_raw_path <- "Wildfire_Data/raw_data/raw_mtw_API_2017_2025.rds"

if (file.exists(meta_path)) {
  meta_raw <- read_csv(meta_path, locale = locale(encoding = "CP949"), 
                       show_col_types = FALSE, name_repair = "unique")
  
  # 공통 정제 로직
  stn_base <- meta_raw %>%
    transmute(
      raw_id = as.integer(지점),
      stn_nm = 지점명,
      stn_lat = as.numeric(위도),
      stn_lon = as.numeric(경도)
    ) %>%
    filter(!is.na(raw_id), !is.na(stn_lat), !is.na(stn_lon))
  
  # 1. ASOS 마스터 (A + 90~296)
  asos_master <- stn_base %>%
    filter(raw_id >= 90 & raw_id <= 296) %>%
    mutate(stn_id = paste0("A", raw_id)) %>%
    dplyr::select(stn_id, stn_nm, stn_lat, stn_lon)
  saveRDS(asos_master, "Wildfire_Data/meta_data/stations_asos_master.rds")
  
  # 2. AWS 마스터 (W + 300~900)
  aws_master <- stn_base %>%
    filter(raw_id >= 300 & raw_id <= 900) %>%
    mutate(stn_id = paste0("W", raw_id)) %>%
    dplyr::select(stn_id, stn_nm, stn_lat, stn_lon)
  saveRDS(aws_master, "Wildfire_Data/meta_data/stations_aws_master.rds")
  
  # 3. MTW 마스터 (M + 5자리 ID)
  if (file.exists(mtw_raw_path)) {
    mtw_ids <- readRDS(mtw_raw_path)$obsid %>% unique() %>% as.integer()
    mtw_master <- stn_base %>%
      filter(raw_id %in% mtw_ids) %>%
      mutate(stn_id = paste0("M", raw_id)) %>%
      dplyr::select(stn_id, stn_nm, stn_lat, stn_lon)
    saveRDS(mtw_master, "Wildfire_Data/meta_data/stations_mtw_master.rds")
  }
}