## =========================================================
## Wildfire.M1_SAR.R
## 분석 특성: 공간 자기상관(Spatial Autocorrelation)을 반영한 SAR 패널 모형
##            - 고립 지역(Islands) 대응을 위한 급수 전개(Power Series) 방식 적용
##            - 시공간 데이터의 인접 전이 효과(Spatial Spillover) 정량화
## =========================================================

# ---------------------------------------------------------
# 0. 데이터 로드 및 전처리
# ---------------------------------------------------------
message(">>> 데이터 로드 및 Factor 설정 중...")
train_s <- readr::read_rds("Wildfire_Data/processed_data/df_train_final.rds")
test_s  <- readr::read_rds("Wildfire_Data/processed_data/df_test_final.rds")

# 시군구 코드의 일관성을 위한 Factor 레벨 동기화
train_s <- train_s %>% mutate(sigungu_cd = factor(sigungu_cd))
test_s  <- test_s  %>% mutate(sigungu_cd = factor(sigungu_cd, levels = levels(train_s$sigungu_cd)))


# ---------------------------------------------------------
# 1. 공간 계수 추출 함수: splm 객체 내부 구조 대응
# ---------------------------------------------------------
# ML 추정 방식(splm_ML)은 rho를 'arcoef'에 저장하며, 
# 경우에 따라 'errcomp' 벡터 내에 존재할 수 있음을 반영함
get_rho_final <- function(fit) {
  # (1) arcoef 주머니 우선 확인 (Atomic vector 에러 방지)
  rho <- as.numeric(fit$arcoef)
  
  # (2) 위 경로가 비어있을 경우 errcomp 내 lambda 검색
  if (length(rho) == 0 || is.na(rho)) {
    rho <- as.numeric(fit$errcomp["lambda"])
  }
  
  if (length(rho) == 0) stop("공간 계수(lambda) 추출 실패. arcoef 또는 errcomp를 확인하세요.")
  return(rho)
}

# ---------------------------------------------------------
# 2. 강건한 공간 예측 함수: Power Series Expansion
# ---------------------------------------------------------
# 고립된 시군구(섬)로 인한 singular matrix 문제를 해결하기 위해 
# 직접적인 역행렬 계산 대신 급수 전개 방식을 사용하여 공간 파급 효과 산출
sar_predict_final <- function(fit, newdata, W_matrix, x_vars, order = 2) {
  rho <- get_rho_final(fit)
  b   <- coef(fit)
  
  # (1) 선형 결합(xb) 계산: 개별 시군구의 기상적 위험(Local Effect)
  X <- as.matrix(newdata[, x_vars])
  intercept <- if ("(Intercept)" %in% names(b)) b["(Intercept)"] else 0
  xb <- as.numeric(intercept + X %*% b[x_vars])
  
  # (2) 결측치 방어: 행렬 연산의 도미노 오염 방지
  xb[is.na(xb)] <- 0
  
  # (3) 급수 전개를 통한 공간 전이(Spillover) 가산
  # y = xb + rho*W*xb + rho^2*W^2*xb (2차수: 이웃의 이웃까지 반영)
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
# 3. 2025년 테스트 데이터 적용 및 지수 산출
# ---------------------------------------------------------
message(">>> M1 SAR 모델 최종 예측치 산출 중...")

# (1) 테스트 데이터 균형화: 252개 시군구 격자 확보
test_2025_balanced <- balance_spatial_panel(test_spatial, ids_in_W, 2025)

# (2) 빈도(N) 및 심도(S) 공간 예측 수행
yhat_freq <- sar_predict_final(sar_freq_m1, test_2025_balanced, W_rs, x_terms)
yhat_sev  <- sar_predict_final(sar_sev_m1, test_2025_balanced, W_rs, x_terms)

# (3) 최종 기대 지급액(EP) 및 성능 지표 결합
# log1p로 학습되었으므로 expm1로 역변환 수행
test_2025_final <- test_2025_balanced %>%
  mutate(
    EP_M1_SAR = expm1(yhat_freq) * expm1(yhat_sev)
  )

# ---------------------------------------------------------
# 4. 최종 모델 비교: M0(전통적) vs M1(공간적)
# ---------------------------------------------------------
# M1 지표 산출
metrics_m1_full <- test_2025_final %>%
  summarise(
    model_type  = "M1_SAR",
    RMSE        = sqrt(mean((damage_capped - EP_M1_SAR)^2, na.rm = TRUE)),
    MAE         = mean(abs(damage_capped - EP_M1_SAR), na.rm = TRUE),
    MedAE       = median(abs(damage_capped - EP_M1_SAR), na.rm = TRUE),
    Spearman    = cor(damage_capped, EP_M1_SAR, method = "spearman"),
    Calib_Ratio = sum(EP_M1_SAR, na.rm = TRUE) / sum(damage_capped, na.rm = TRUE)
  )

# 기존 M0 결과와 병합하여 최종 성능표 생성
comparison_table <- bind_rows(m0_final_perf, metrics_m1_full)

message("--- Final Model Comparison: M0 vs M1 ---")
print(comparison_table)