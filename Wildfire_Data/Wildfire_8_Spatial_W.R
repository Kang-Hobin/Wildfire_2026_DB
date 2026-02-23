## =========================================================
## Wildfire_8_Spatial_W.R
## 목적: 시군구 공간 가중치 행렬(W) 생성 (252개 시군구 버전)
## 특징: 투영좌표계(EPSG:5179) 적용
## =========================================================

library(sf)
library(dplyr)
library(spdep)
library(Matrix)

# 1. 패널 데이터의 시군구 목록 확인 (252개)
# ---------------------------------------------------------
# 분석 데이터와 W 행렬의 순서를 맞추기 위해 ID 목록을 먼저 확보합니다.
df_clean <- readRDS(file.path("Wildfire_Data/processed_data/df_wildfire_final_cleaned.rds"))
final_sgg_ids <- sort(unique(df_clean$sigungu_cd)) #

# 2. Shapefile 로드 및 정제
# ---------------------------------------------------------
message(">>> 시군구 경계 데이터 로드 및 춘천 포함 정제 중...")

sgg_raw <- st_read(file.path("Wildfire_Data/meta_data/시군구경계.shp"), quiet = TRUE)

sgg_clean <- sgg_raw %>%
  # 춘천(32010)을 제거하던 코드를 삭제하여 252개를 유지합니다.
  filter(SIGUNGU_CD %in% final_sgg_ids) %>% 
  st_make_valid() %>%
  group_by(SIGUNGU_CD) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  # [핵심] 분석 데이터와 동일한 순서로 정렬합니다.
  arrange(SIGUNGU_CD) %>% 
  # 투영좌표계(중부원점 모디파이드 TM)로 변환
  st_transform(5179) %>%
  # 토폴로지 단순화 (인접성 계산 안정화)
  st_simplify(dTolerance = 10)

# 3. 인접 행렬(Neighborhood) 생성 (Queen Contiguity)
# ---------------------------------------------------------
# snap = 1을 주어 아주 미세한 틈이 있어도 인접한 것으로 간주합니다.
nb <- poly2nb(sgg_clean, queen = TRUE, snap = 1)
attr(nb, "region.id") <- sgg_clean$SIGUNGU_CD

# 4. 행-표준화 가중치 행렬($W$) 생성
# ---------------------------------------------------------
# style = "W"는 각 행의 합이 1이 되도록 표준화합니다. $\sum_{j} w_{ij} = 1$
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)
W_sparse <- as_dgRMatrix_listw(listw)

# 행/열 이름 부여 (ID 매핑 보장)
rownames(W_sparse) <- sgg_clean$SIGUNGU_CD
colnames(W_sparse) <- sgg_clean$SIGUNGU_CD

# 5. 최종 객체 구성 및 저장 (Professional Path)
# ---------------------------------------------------------
spatial_weights <- list(
  W = W_sparse,
  sigungu_cd = sgg_clean$SIGUNGU_CD,
  description = "2024 2Q 시군구 252개, Row-standardized matrix"
)

# 파일명을 좀 더 명시적으로 지정
saveRDS(spatial_weights, "Wildfire_Data/meta_data/sgg_spatial_weights_252.rds")

message(">>> [저장 완료] 경로: Wildfire_Data/meta_data/sgg_spatial_weights_252.rds")