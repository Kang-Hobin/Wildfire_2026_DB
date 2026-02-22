## =========================================================
## Wildfire_5_Spatial_Mapping.R
## 목적: 시군구별 ASOS, AWS, MTW(연도별) 최단거리 매핑 테이블 생성
## =========================================================

library(sf)
library(dplyr)
library(purrr)

# 1. 재료 로드
sgg_points <- st_as_sf(readRDS("Wildfire_Data/meta_data/sigungu_centroid_master.rds"), 
                       coords = c("sgg_lon", "sgg_lat"), crs = 4326)
asos_ref <- readRDS("Wildfire_Data/meta_data/stations_asos_master.rds")
aws_ref  <- readRDS("Wildfire_Data/meta_data/stations_aws_master.rds")
mtw_ref  <- readRDS("Wildfire_Data/meta_data/stations_mtw_master.rds")

# 2. 최단거리 매핑 헬퍼 함수
get_top3_map <- function(target_sf, ref_df, type_label) {
  ref_sf <- st_as_sf(ref_df, coords = c("stn_lon", "stn_lat"), crs = 4326)
  dist_mat <- st_distance(target_sf, ref_sf) # 단위: 미터(m)
  
  top3_idx <- apply(dist_mat, 1, function(x) order(x)[1:3]) %>% t()
  
  map_dfc(1:3, function(j) {
    idx <- top3_idx[, j]
    tibble(
      !!paste0("match_", type_label, "_id_", j) := ref_df$stn_id[idx],
      !!paste0("dist_", type_label, "_m_", j) := as.numeric(dist_mat[cbind(seq_len(nrow(dist_mat)), idx)])
    )
  })
}

# 3. ASOS 및 AWS 매핑 (고정 테이블)
message(">>> ASOS 및 AWS 공간 매칭 수행 중...")
asos_top3 <- get_top3_map(sgg_points, asos_ref, "ASOS")
aws_top3  <- get_top3_map(sgg_points, aws_ref, "AWS")

sgg_asos_aws_map <- readRDS("Wildfire_Data/meta_data/sigungu_centroid_master.rds") %>%
  bind_cols(asos_top3, aws_top3)

saveRDS(sgg_asos_aws_map, "Wildfire_Data/meta_data/sgg_asos_aws_triple_map.rds")

# 4. MTW 연도별 매핑 (2018-2022)
# ---------------------------------------------------------
# MTW는 신규 설치 지점이 많아 연도별 매핑이 중요합니다.
message(">>> MTW 연도별 매핑 테이블 생성 중...")

mtw_annual_maps <- map_dfr(2018:2022, function(yy) {
  # (실제로는 연도별 설치/폐지일 필터링 로직이 추가될 수 있음)
  mtw_top3 <- get_top3_map(sgg_points, mtw_ref, "MTW")
  
  readRDS("Wildfire_Data/meta_data/sigungu_centroid_master.rds") %>%
    bind_cols(mtw_top3) %>%
    mutate(map_year = yy)
})

saveRDS(mtw_annual_maps, "Wildfire_Data/meta_data/sgg_mtw_annual_map_2018_2022.rds")
message(">>> 모든 매핑 테이블 생성 완료!")