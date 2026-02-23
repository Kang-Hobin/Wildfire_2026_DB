## =========================================================
## Wildfire.M1_SAR.R
## 분석 특성: 공간 자기상관(Spatial Autocorrelation)을 반영한 SAR 패널 모형
##            - 고립 지역(Islands) 대응을 위한 급수 전개(Power Series) 방식 적용
##            - 시공간 데이터의 인접 전이 효과(Spatial Spillover) 정량화
## =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(splm)
  library(spdep)
  library(Matrix)
  library(pROC)
})

# ---------------------------------------------------------
# 1. 데이터 및 공간 객체(W) 로드
# ---------------------------------------------------------
message(">>> [1/4] 데이터 및 공간 가중치 행렬 로드 중...")

train_s <- readr::read_rds("Wildfire_Data/processed_data/df_train_final.rds")
test_s  <- readr::read_rds("Wildfire_Data/processed_data/df_test_final.rds")
W_obj   <- readRDS("Wildfire_Data/meta_data/sgg_spatial_weights_252.rds")

W_sparse   <- W_obj$W
region_ids <- as.character(W_obj$sigungu_cd) 

# 공간 행렬 표준화 및 예측용 Matrix 준비
W_mat   <- as.matrix(W_sparse)
W_listw <- spdep::mat2listw(W_mat, style = "W", zero.policy = TRUE)
attr(W_listw$neighbours, "region.id") <- region_ids
W_rs    <- Matrix::Matrix(spdep::listw2mat(W_listw), sparse = TRUE)

# ---------------------------------------------------------
# 2. 필수 함수 정의 (결측치 방어 및 공간 예측)
# ---------------------------------------------------------
message(">>> [2/4] 핵심 함수 정의 중...")

# (A) 패널 균형화 함수: 행렬 연산을 위한 2단계 결측치 방어선
balance_spatial_panel <- function(df, ids, years) {
  df %>%
    complete(year = years, sigungu_cd = ids) %>%
    mutate(across(where(is.numeric), as.numeric)) %>%
    # 산불 지표 초기화
    mutate(
      fire_cnt_year = replace_na(fire_cnt_year, 0),
      damage_capped = replace_na(damage_capped, 0),
      y_freq = log1p(fire_cnt_year),
      y_sev  = log1p(damage_capped / pmax(fire_cnt_year, 1))
    ) %>%
    # 1차 방어: 시군구별 평균 대체
    group_by(sigungu_cd) %>%
    mutate(across(where(is.numeric), ~ replace_na(., mean(., na.rm = TRUE)))) %>%
    ungroup() %>%
    # 2차 방어: 전체 평균 대체 (데이터 부재 지역 대응)
    mutate(across(where(is.numeric), ~ replace_na(., mean(., na.rm = TRUE)))) %>%
    # Factor 레벨 및 정렬
    mutate(sigungu_cd = factor(as.character(sigungu_cd), levels = ids)) %>%
    arrange(year, sigungu_cd)
}

# (B) 공간 계수 추출 함수: splm 객체 내부 구조(arcoef) 대응
get_rho_final <- function(fit) {
  rho <- as.numeric(fit$arcoef)
  if (length(rho) == 0 || is.na(rho)) rho <- as.numeric(fit$errcomp["lambda"])
  if (length(rho) == 0) stop("공간 계수(lambda) 추출 실패.")
  return(rho)
}

# (C) 강건한 공간 예측 함수: Power Series Expansion (섬 지역 대응)
sar_predict_final <- function(fit, newdata, W_matrix, x_vars, order = 2) {
  rho <- get_rho_final(fit)
  b   <- coef(fit)
  X   <- as.matrix(newdata[, x_vars])
  intercept <- if ("(Intercept)" %in% names(b)) b["(Intercept)"] else 0
  xb  <- as.numeric(intercept + X %*% b[x_vars])
  xb[is.na(xb)] <- 0
  
  yhat <- xb
  W_pow <- W_matrix
  curr_rho <- rho
  for (i in 1:order) {
    yhat <- yhat + curr_rho * as.numeric(W_pow %*% xb)
    if (i < order) {
      W_pow <- W_pow %*% W_matrix
      curr_rho <- curr_rho * rho
    }
  }
  return(as.numeric(yhat))
}

# ---------------------------------------------------------
# 3. 모델 적합 및 2025년 예측
# ---------------------------------------------------------
message(">>> [3/4] SAR 모델 적합 및 2025년 지수 산출 중...")

# 독립변수 정의
x_terms <- c("spring_ws_p95_lag1", "spring_he_min_lag1", "max_dry_run_spring_lag1", 
             "monsoon_prcp_sum_lag1", "winter_dry_days_lag1", "ta_sfc_mean_lag1")

# 학습 데이터 균형화 및 모델 적합
all_train_years <- min(train_s$year):max(train_s$year)
train_balanced  <- balance_spatial_panel(train_s, region_ids, all_train_years)

sar_freq_m1 <- spml(y_freq ~ spring_ws_p95_lag1 + spring_he_min_lag1 + max_dry_run_spring_lag1 + 
                      monsoon_prcp_sum_lag1 + winter_dry_days_lag1 + ta_sfc_mean_lag1,
                    data = train_balanced, index = c("sigungu_cd", "year"),
                    listw = W_listw, model = "random", lag = TRUE)

sar_sev_m1  <- spml(y_sev ~ spring_ws_p95_lag1 + spring_he_min_lag1 + max_dry_run_spring_lag1 + 
                      monsoon_prcp_sum_lag1 + winter_dry_days_lag1 + ta_sfc_mean_lag1,
                    data = train_balanced, index = c("sigungu_cd", "year"),
                    listw = W_listw, model = "random", lag = TRUE)

# 테스트 데이터(2025년) 예측
test_2025_balanced <- balance_spatial_panel(test_s, region_ids, 2025)
yhat_freq <- sar_predict_final(sar_freq_m1, test_2025_balanced, W_rs, x_terms)
yhat_sev  <- sar_predict_final(sar_sev_m1, test_2025_balanced, W_rs, x_terms)

test_2025_final <- test_2025_balanced %>%
  mutate(EP_M1_SAR = expm1(yhat_freq) * expm1(yhat_sev))

# ---------------------------------------------------------
# 4. 성능 비교 및 결과 출력
# ---------------------------------------------------------
message(">>> [4/4] 최종 성능 지표 산출 중...")

# ---------------------------------------------------------
# 4. 성능 비교 및 결과 출력 (ROC-AUC 추가)
# ---------------------------------------------------------
message(">>> [4/4] 최종 성능 지표 및 ROC-AUC 산출 중...")

# (1) 실제 발생 여부(Binary) 생성
# 피해액이 0보다 크거나 발생 건수가 0보다 큰 경우를 1로 정의
fire_actual <- ifelse(test_2025_final$damage_capped > 0, 1, 0)

# (2) ROC 객체 생성
# 예측된 빈도 지수(yhat_freq)가 높을수록 발생 확률이 높다고 판단

roc_obj <- pROC::roc(fire_actual, test_2025_final$EP_M1_SAR, quiet = TRUE)
auc_val <- as.numeric(pROC::auc(roc_obj))

# (3) 통합 지표 출력
metrics_m1 <- test_2025_final %>%
  summarise(
    model_type  = "M1_SAR",
    AUC_occ     = auc_val,             # 발생 판별력 추가
    RMSE        = sqrt(mean((damage_capped - EP_M1_SAR)^2, na.rm = TRUE)),
    MAE         = mean(abs(damage_capped - EP_M1_SAR), na.rm = TRUE),
    MedAE       = median(abs(damage_capped - EP_M1_SAR), na.rm = TRUE),
    Spearman    = cor(damage_capped, EP_M1_SAR, method = "spearman"),
    Calib_Ratio = sum(EP_M1_SAR, na.rm = TRUE) / sum(damage_capped, na.rm = TRUE)
  )

print(metrics_m1)

# (4) 시각화: ROC Curve (선택 사항)
plot(roc_obj, main=paste0("ROC Curve (AUC = ", round(auc_val, 3), ")"), col="#2c3e50", lwd=2)
