## =========================================================
## Wildfire_4_Meta_Geocoding.R (Refined)
## =========================================================
library(sf)
library(dplyr)
library(stringr)

# 시도 코드 -> 2글자 약칭 변환 함수 (무결성 확보)
cd_to_sido_short <- function(sigungu_cd) {
  p <- substr(sigungu_cd, 1, 2)
  dplyr::case_when(
    p == "11" ~ "서울", p == "21" ~ "부산", p == "22" ~ "대구", p == "23" ~ "인천",
    p == "24" ~ "광주", p == "25" ~ "대전", p == "26" ~ "울산", p == "29" ~ "세종",
    p == "31" ~ "경기", p == "32" ~ "강원", p == "33" ~ "충북", p == "34" ~ "충남",
    p == "35" ~ "전북", p == "36" ~ "전남", p == "37" ~ "경북", p == "38" ~ "경남",
    p == "39" ~ "제주", TRUE ~ NA_character_
  )
}

# [A] 시군구 중심점 마스터 생성
message(">>> 시군구 중심점 및 조인 키 생성 중...")
sgg_path <- "Wildfire_Data/meta_data/시군구경계.shp"

if (file.exists(sgg_path)) {
  sgg_master <- st_read(sgg_path, quiet = TRUE) %>%
    st_point_on_surface() %>%
    st_transform(4326) %>%
    mutate(
      lon = st_coordinates(.)[,1],
      lat = st_coordinates(.)[,2],
      sigungu_cd = sprintf("%05d", as.integer(SIGUNGU_CD)),
      # 핵심: 조인용 시도 약칭("강원")과 시군구 본명("강릉시") 분리
      sido_nm = cd_to_sido_short(sigungu_cd),
      sgg_base_nm = stringr::word(SIGUNGU_NM, -1) # "강원도 강릉시" -> "강릉시"
    ) %>%
    st_drop_geometry() %>%
    dplyr::select(sigungu_cd, sigungu_nm = SIGUNGU_NM, sido_nm, sgg_base_nm, sgg_lon = lon, sgg_lat = lat)
  
  saveRDS(sgg_master, "Wildfire_Data/meta_data/sigungu_centroid_master.rds")
}