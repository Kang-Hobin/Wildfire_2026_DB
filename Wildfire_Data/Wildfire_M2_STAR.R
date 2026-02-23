## =========================================================
## Wildfire_M2_STAR.R
##
## 목적:
##  (1) 시군구-연도 패널(252 x T)을 균형화(balancing)하고
##  (2) Frequency(빈도) STAR 계열 5개 모델(M1~M5)을 적합한 뒤
##  (3) Severity(심도) 모델 1개(splm 기반: sar_sev_m1)를 결합하여
##  (4) 2025년 Out-of-Time 예측(EP = N_hat * Sev_hat) 및 성능지표를 산출한다.
##
## 핵심 설계:
##  - Two-part 구조: 빈도 N 모델 5개 + 심도 S 모델 1개 결합
##  - 빈도 모델:
##      M1: Mixed-ST (OLS)
##      M2: Pure-ST  (OLS)
##      M3: Pure-STAR (SAR)
##      M4: Mixed-STAR (SAR + X)
##      M5: Res-STAR (SDM/Durbin)
##  - 심도 모델(sar_sev_m1):
##      splm::spml(ML) 계열 객체라고 가정 (class: splm, splm_ML)
##      공간계수는 sar_sev_m1$arcoef["lambda"] 사용
##
## 안정화 포인트(보고서/공유용):
##  - 스케일링은 train 통계로만 수행하고 test에 동일 적용(누수 방지)
##  - Durbin(lag.X) 항은 newdata에 없을 수 있어 자동 생성(W %*% X)
##  - 심도는 로그 스케일 예측을 역변환할 때 Duan smearing factor로 바이어스 보정
##  - 필요 시(옵션) 심도 예측치를 train 분위수 구간으로 clip (설명 가능)
##
## 출력물:
##  - Wildfire_Data/results/m2_star_5models_perf_final.rds
##  - Wildfire_Data/results/m2_star_5models_perf_final.csv
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(Matrix)
  library(spdep)
  library(spatialreg)
  library(pROC)
})

## ---------------------------------------------------------
## 0. 사용자 설정(경로/옵션)
## ---------------------------------------------------------
PATH_TRAIN <- "Wildfire_Data/processed_data/df_train_final.rds"
PATH_TEST  <- "Wildfire_Data/processed_data/df_test_final.rds"
PATH_W     <- "Wildfire_Data/meta_data/sgg_spatial_weights_252.rds"

## 심도모형(sar_sev_m1) 로딩 옵션
## - 추천: 심도모형을 RDS로 저장해두고 경로로 불러오면 “이 파일만 실행”이 가능해짐
USE_SEV_MODEL_RDS <- FALSE
PATH_SEV_MODEL    <- "Wildfire_Data/models/sar_sev_m1.rds"  # 필요 시 경로 수정

OUT_DIR <- "Wildfire_Data/results"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

OUT_RDS <- file.path(OUT_DIR, "m2_star_5models_perf_final.rds")
OUT_CSV <- file.path(OUT_DIR, "m2_star_5models_perf_final.csv")

## 테스트 연도(Out-of-Time)
TEST_YEAR <- 2025

## Power series 차수: (I - sp W)^(-1) 근사
PS_ORDER <- 2

## 심도 예측 안정화 옵션
CLIP_SEV_PRED <- TRUE    # TRUE: train 분위수 기반으로 y_sev_pred clip
SEV_CLIP_Q    <- c(0.01, 0.99)

## 빈도 모델에 사용할 기상/환경 설명변수(반드시 df에 존재해야 함)
COVARS_LAG1 <- c(
  "spring_ws_p95_lag1",
  "spring_he_min_lag1",
  "max_dry_run_spring_lag1",
  "monsoon_prcp_sum_lag1",
  "winter_dry_days_lag1",
  "ta_sfc_mean_lag1"
)

msg <- function(...) message(sprintf(...))

## ---------------------------------------------------------
## 1. 공통 유틸: W 생성 / 패널 균형화 / 스케일링
## ---------------------------------------------------------

## (1) row-standardized W (고립 지역: zero.policy=TRUE)
build_W_sparse <- function(W_dense, zero_policy = TRUE) {
  lw <- spdep::mat2listw(as.matrix(W_dense), style = "W", zero.policy = zero_policy)
  Wm <- spdep::listw2mat(lw)
  Matrix::Matrix(Wm, sparse = TRUE)
}

## (2) 패널 균형화 + NA 방어
balance_spatial_panel <- function(df, ids, years) {
  df %>%
    mutate(sigungu_cd = as.character(sigungu_cd),
           year = as.integer(year)) %>%
    tidyr::complete(year = years, sigungu_cd = ids) %>%
    mutate(
      fire_cnt_year = tidyr::replace_na(as.numeric(fire_cnt_year), 0),
      damage_capped = tidyr::replace_na(as.numeric(damage_capped), 0),
      fire_any_year = if_else(fire_cnt_year > 0, 1, 0),
      ## 빈도 타겟
      y_freq = log1p(fire_cnt_year),
      ## 심도 타겟(연 평균 “건당” 피해)
      y_sev  = log1p(damage_capped / pmax(fire_cnt_year, 1))
    ) %>%
    ## 1차 방어: 지역 내 평균 대체
    group_by(sigungu_cd) %>%
    mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
    ungroup() %>%
    ## 2차 방어: 전체 평균 대체
    mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
    mutate(sigungu_cd = factor(sigungu_cd, levels = ids)) %>%
    arrange(year, sigungu_cd)
}

## (3) train 기반 scaling (test에 동일 적용, 누수 방지)
scale_with_train <- function(train_df, test_df, cols) {
  mu  <- sapply(train_df[, cols, drop = FALSE], mean, na.rm = TRUE)
  sdv <- sapply(train_df[, cols, drop = FALSE], sd,   na.rm = TRUE)
  sdv[sdv == 0] <- 1
  
  for (cn in cols) {
    train_df[[cn]] <- (train_df[[cn]] - mu[cn]) / sdv[cn]
    test_df[[cn]]  <- (test_df[[cn]]  - mu[cn]) / sdv[cn]
  }
  list(train = train_df, test = test_df, mu = mu, sd = sdv)
}

## ---------------------------------------------------------
## 2. 핵심: "lm / sarlm / splm" 모두 예측 가능한 robust predictor
## ---------------------------------------------------------

## (A) 모델별 공간계수(sp) 추출
## - sarlm: rho/lambda
## - splm_ML: arcoef["lambda"] (네가 확인한 구조)
get_spatial_param <- function(fit) {
  if (!is.null(fit$rho))    return(as.numeric(fit$rho))
  if (!is.null(fit$lambda)) return(as.numeric(fit$lambda))
  
  if (!is.null(fit$arcoef)) {
    if ("lambda" %in% names(fit$arcoef)) return(as.numeric(fit$arcoef["lambda"]))
    if ("rho" %in% names(fit$arcoef))    return(as.numeric(fit$arcoef["rho"]))
  }
  return(0)
}

## (B) 계수 추출(모델 클래스 혼재 대응)
get_coef_safe <- function(fit) {
  b <- tryCatch(coef(fit), error = function(e) NULL)
  if (!is.null(b)) return(b)
  if (!is.null(fit$coefficients)) return(fit$coefficients)  # splm fallback
  stop("Cannot extract coefficients from fit.")
}

## (C) Durbin 계수(lag.X)가 있을 때 newdata에 자동 생성
## - lagsarlm(Durbin=...)은 coef 이름이 "lag.xxx" 형태로 들어올 수 있음
## - 예측 시 newdata에 lag.xxx 열이 없으면 W %*% xxx로 생성해준다.
ensure_lag_terms <- function(newdata, coef_names, W_matrix) {
  lag_terms <- coef_names[startsWith(coef_names, "lag.")]
  if (length(lag_terms) == 0) return(newdata)
  
  for (lt in lag_terms) {
    if (lt %in% names(newdata)) next
    v <- sub("^lag\\.", "", lt)
    if (v %in% names(newdata)) {
      newdata[[lt]] <- as.numeric(W_matrix %*% as.numeric(newdata[[v]]))
    } else {
      newdata[[lt]] <- 0
    }
  }
  newdata
}

## (D) 예측 함수(모델 행렬에 의존하지 않고 "계수명-컬럼명" 매칭으로 xb 구성)
## - 장점: splm/부분적 formula 문제 회피, “Some coefficients not found” 최소화
predict_star_robust <- function(fit, newdata, W_matrix, order = 2) {
  b <- get_coef_safe(fit)
  sp <- get_spatial_param(fit)
  
  coef_names <- names(b)
  newdata <- ensure_lag_terms(newdata, coef_names, W_matrix)
  
  ## Intercept 제외 변수명
  vars <- setdiff(coef_names, "(Intercept)")
  
  ## newdata에 없는 변수는 0으로 생성(최후 방어)
  missing_vars <- setdiff(vars, names(newdata))
  if (length(missing_vars) > 0) {
    for (v in missing_vars) newdata[[v]] <- 0
  }
  
  X <- as.matrix(newdata[, vars, drop = FALSE])
  
  xb <- if ("(Intercept)" %in% coef_names) as.numeric(b["(Intercept)"]) else 0
  xb <- xb + as.numeric(X %*% b[vars])
  xb[is.na(xb)] <- 0
  
  ## 공간 파급 효과(근사)
  yhat <- xb
  if (!is.na(sp) && sp != 0) {
    W_pow <- W_matrix
    sp_k <- sp
    for (k in 1:order) {
      yhat <- yhat + sp_k * as.numeric(W_pow %*% xb)
      if (k < order) {
        W_pow <- W_pow %*% W_matrix
        sp_k  <- sp_k * sp
      }
    }
  }
  as.numeric(yhat)
}

## ---------------------------------------------------------
## 3. 평가 지표(robust)
## ---------------------------------------------------------
evaluate_models <- function(actual_cnt, actual_loss, n_hat, ep_hat, model_name) {
  actual_occ <- as.numeric(actual_cnt > 0)
  
  n_hat  <- pmax(as.numeric(n_hat), 0)
  ep_hat <- pmax(as.numeric(ep_hat), 0)
  actual_loss <- as.numeric(actual_loss)
  
  ## AUC: occurrence는 N_hat로 평가
  auc_val <- NA_real_
  if (length(unique(actual_occ)) >= 2) {
    auc_val <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(actual_occ, n_hat, quiet = TRUE))),
      error = function(e) NA_real_
    )
  }
  
  rmse  <- sqrt(mean((actual_loss - ep_hat)^2, na.rm = TRUE))
  mae   <- mean(abs(actual_loss - ep_hat), na.rm = TRUE)
  medae <- median(abs(actual_loss - ep_hat), na.rm = TRUE)
  
  ## Spearman: 정보가 있는 관측치만(실제>0 또는 예측>0)
  idx_inf <- which((actual_loss > 0) | (ep_hat > 0))
  spearman <- NA_real_
  if (length(idx_inf) >= 3) {
    a <- actual_loss[idx_inf]
    p <- ep_hat[idx_inf]
    if (sd(a, na.rm = TRUE) > 0 && sd(p, na.rm = TRUE) > 0) {
      spearman <- suppressWarnings(cor(a, p, method = "spearman", use = "complete.obs"))
    } else {
      spearman <- 0
    }
  }
  
  denom_all <- sum(actual_loss, na.rm = TRUE)
  calib_all <- ifelse(denom_all > 0, sum(ep_hat, na.rm = TRUE) / denom_all, NA_real_)
  
  idx_pos <- which(actual_loss > 0)
  denom_pos <- sum(actual_loss[idx_pos], na.rm = TRUE)
  calib_pos <- ifelse(denom_pos > 0, sum(ep_hat[idx_pos], na.rm = TRUE) / denom_pos, NA_real_)
  
  mean_pos_actual <- mean(actual_loss[idx_pos], na.rm = TRUE)
  mean_pos_pred   <- mean(ep_hat[idx_pos], na.rm = TRUE)
  calib_mean_pos  <- ifelse(is.finite(mean_pos_actual) && mean_pos_actual > 0,
                            mean_pos_pred / mean_pos_actual, NA_real_)
  
  tibble(
    model = model_name,
    AUC_occ = auc_val,
    RMSE = rmse,
    MAE = mae,
    MedAE = medae,
    Spearman = spearman,
    Calib = calib_all,
    Calib_Pos = calib_pos,
    Calib_Mean_Pos = calib_mean_pos,
    Pred_ZeroShare = mean(ep_hat == 0, na.rm = TRUE)
  )
}

## ---------------------------------------------------------
## 4. 메인 파이프라인: Load -> Panel -> Fit -> Predict -> Eval -> Save
## ---------------------------------------------------------
msg(">>> [1/6] Loading data...")

df_train_raw <- readRDS(PATH_TRAIN)
df_test_raw  <- readRDS(PATH_TEST)
df_all <- bind_rows(df_train_raw, df_test_raw)

W_obj <- readRDS(PATH_W)
region_ids <- as.character(W_obj$sigungu_cd)
years_all  <- sort(unique(as.integer(df_all$year)))

W_rs <- build_W_sparse(W_obj$W, zero_policy = TRUE)
dimnames(W_rs) <- list(region_ids, region_ids)

msg(">>> Regions: %d | Years: %d (%d-%d)",
    length(region_ids), length(years_all), min(years_all), max(years_all))

msg(">>> [2/6] Building balanced panel + temporal/spatial lags...")

panel_bal <- balance_spatial_panel(df_all, region_ids, years_all)

## 시차(1년) 생성
panel_bal <- panel_bal %>%
  group_by(sigungu_cd) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    y_freq_lag1 = lag(y_freq, 1),
    y_sev_lag1  = lag(y_sev, 1)
  ) %>%
  ungroup() %>%
  filter(!is.na(y_freq_lag1), !is.na(y_sev_lag1))

## 연도별 공간시차 Wy_freq_lag1
panel_bal <- panel_bal %>%
  group_by(year) %>%
  mutate(Wy_freq_lag1 = as.numeric(W_rs %*% y_freq_lag1)) %>%
  ungroup()

## 공변량 체크
missing_covars <- setdiff(COVARS_LAG1, names(panel_bal))
if (length(missing_covars) > 0) {
  stop("Missing covariates in data: ", paste(missing_covars, collapse = ", "))
}

msg(">>> [3/6] Train/Test split + scaling (no leakage)...")

train_st <- panel_bal %>% filter(year < TEST_YEAR)
test_st  <- panel_bal %>% filter(year == TEST_YEAR)

stopifnot(nrow(test_st) == length(region_ids))

## 빈도 모델용 scaling (빈도에서만 표준화 적용)
scale_cols_freq <- c("y_freq_lag1", "Wy_freq_lag1", COVARS_LAG1)
scaled_freq <- scale_with_train(train_st, test_st, scale_cols_freq)
train_st <- scaled_freq$train
test_st  <- scaled_freq$test

msg(">>> [4/6] Fit frequency models M1~M5...")

## M1: Mixed-ST (OLS)
f_m1 <- as.formula(paste(
  "y_freq ~ y_freq_lag1 + Wy_freq_lag1 +",
  paste(COVARS_LAG1, collapse = " + ")
))
m1_st <- lm(f_m1, data = train_st)

## M2: Pure-ST (OLS)
m2_st <- lm(y_freq ~ y_freq_lag1 + Wy_freq_lag1, data = train_st)

## STAR 모델용 block W(listw) 구성
train_years <- sort(unique(train_st$year))
T_train <- length(train_years)

## (매우 중요) block W는 (year, sigungu_cd) 정렬을 전제로 함
train_st <- train_st %>% arrange(year, sigungu_cd)
test_st  <- test_st  %>% arrange(year, sigungu_cd)

W_block_train <- Matrix::kronecker(Matrix::Diagonal(T_train), W_rs)
W_time_train  <- spdep::mat2listw(as.matrix(W_block_train), style = "W", zero.policy = TRUE)

## M3: Pure-STAR (SAR)
m3_star <- spatialreg::lagsarlm(
  y_freq ~ y_freq_lag1,
  data = train_st,
  listw = W_time_train,
  method = "LU",
  zero.policy = TRUE
)

## M4: Mixed-STAR (SAR + X)
f_m4 <- as.formula(paste(
  "y_freq ~ y_freq_lag1 +",
  paste(COVARS_LAG1, collapse = " + ")
))
m4_star <- spatialreg::lagsarlm(
  f_m4,
  data = train_st,
  listw = W_time_train,
  method = "LU",
  zero.policy = TRUE
)

## M5: Res-STAR (SDM/Durbin)
m5_star <- spatialreg::lagsarlm(
  f_m4,
  data = train_st,
  listw = W_time_train,
  method = "LU",
  zero.policy = TRUE,
  Durbin = as.formula(paste("~", paste(COVARS_LAG1, collapse = " + ")))
)

msg(">>> [5/6] Load/Fit severity model + predict Sev_hat (with smearing)...")

## 심도 모형 로드 옵션(이 파일 단독 실행 목적)
if (USE_SEV_MODEL_RDS) {
  if (!file.exists(PATH_SEV_MODEL)) stop("Severity model RDS not found: ", PATH_SEV_MODEL)
  sar_sev_m1 <- readRDS(PATH_SEV_MODEL)
}

## 심도 모형 존재 확인
if (!exists("sar_sev_m1")) {
  stop("Object sar_sev_m1 not found. Load it in advance or set USE_SEV_MODEL_RDS=TRUE.")
}

## (A) 로그 스케일 심도 예측
y_sev_pred <- predict_star_robust(sar_sev_m1, test_st, W_rs, order = PS_ORDER)

## (B) 선택 옵션: 심도 예측치 clip (train 분위수 기반)
## - 하드코딩(-5, 0.001 등) 대신, 데이터 기반이라 설명이 쉬움
if (CLIP_SEV_PRED) {
  q_lo <- quantile(train_st$y_sev, SEV_CLIP_Q[1], na.rm = TRUE)
  q_hi <- quantile(train_st$y_sev, SEV_CLIP_Q[2], na.rm = TRUE)
  y_sev_pred <- pmin(pmax(y_sev_pred, q_lo), q_hi)
}

## (C) Duan smearing factor (train residual로 자동 추정)
## - 로그 변환 모델을 역변환할 때 발생하는 바이어스 보정
## - 가능하면 fitted.values가 존재해야 함(대부분 splm_ML에 존재)
smear <- 1.0
if (!is.null(sar_sev_m1$fitted.values)) {
  resid_train <- as.numeric(train_st$y_sev) - as.numeric(sar_sev_m1$fitted.values)
  smear <- mean(exp(resid_train), na.rm = TRUE)
} else {
  msg(">>> Warning: sar_sev_m1$fitted.values not found. Smearing factor set to 1.0")
}

## ---------------------------------------------------------
## (D) 역변환: Sev_hat = E[sev | X] 
## ---------------------------------------------------------
# expm1()은 x가 음수일 때 -1 ~ 0 사이를 반환함. 
# pmax(..., 0)를 하기 전에, exp(y_sev_pred) 방식으로 접근하여 0 증발을 원천 차단.

# 1. exp(y_sev_pred)는 이론적으로는 양수여야 함
sev_hat_raw <- exp(y_sev_pred) * smear

# 2. log1p의 역변환은 exp(y)-1 이지만, 0으로 수렴하는 것을 막기 위해 
#    -1을 하기 전의 기댓값이 1보다 작으면 아주 작은 양수(epsilon)를 부여함.
sev_hat <- pmax(sev_hat_raw - 1, 1e-8) 

## ---------------------------------------------------------
## (E) 진단 및 최후의 수단 (Fallback)
## ---------------------------------------------------------
msg(">>> Severity Hat Summary: Min=%.6f, Mean=%.6f, Max=%.6f", 
    min(sev_hat), mean(sev_hat), max(sev_hat))

if (all(sev_hat <= 1e-6)) {
  msg(">>> [Emergency] sev_hat still too low. Applying Train-Mean Smearing...")
  # 모델 예측이 완전히 망가졌을 경우를 대비해, 훈련 데이터의 평균 심도를 하한선으로 강제 주입
  train_avg_sev <- mean(expm1(train_st$y_sev[train_st$y_sev > 0]), na.rm = TRUE)
  sev_hat <- pmax(sev_hat, train_avg_sev * 0.1) # 훈련 평균의 10%라도 보장
}

# 이제 절대로 0이 나올 수 없음
stopifnot(any(sev_hat > 0))

msg(">>> [6/6] Predict N_hat, combine EP, evaluate, save...")

models_n <- list(m1_st, m2_st, m3_star, m4_star, m5_star)

## 빈도 예측: y_freq_hat -> N_hat = expm1(y_freq_hat)
n_preds <- lapply(models_n, function(m) {
  pmax(expm1(predict_star_robust(m, test_st, W_rs, order = PS_ORDER)), 0)
})

## EP 결합
ep_preds <- lapply(n_preds, function(n) pmax(n * sev_hat, 0))

## 평가
actual_cnt  <- test_st$fire_cnt_year
actual_loss <- test_st$damage_capped

res <- bind_rows(
  evaluate_models(actual_cnt, actual_loss, n_preds[[1]], ep_preds[[1]], "M1: Mixed-ST (OLS)"),
  evaluate_models(actual_cnt, actual_loss, n_preds[[2]], ep_preds[[2]], "M2: Pure-ST (OLS)"),
  evaluate_models(actual_cnt, actual_loss, n_preds[[3]], ep_preds[[3]], "M3: Pure-STAR (SAR)"),
  evaluate_models(actual_cnt, actual_loss, n_preds[[4]], ep_preds[[4]], "M4: Mixed-STAR (SAR+X)"),
  evaluate_models(actual_cnt, actual_loss, n_preds[[5]], ep_preds[[5]], "M5: Res-STAR (SDM/Durbin)")
)

print(res %>% arrange(desc(Spearman)))

saveRDS(res, OUT_RDS)
write.csv(res, OUT_CSV, row.names = FALSE)

msg(">>> Done. Saved: %s", OUT_RDS)
msg(">>> Done. Saved: %s", OUT_CSV)

## (선택) 추가 진단 출력: 모델별 총 예측 EP 합/실제 합
## - Calib 해석에 도움
diag_tbl <- tibble(
  model = res$model,
  sum_ep = sapply(ep_preds, sum, na.rm = TRUE),
  sum_actual = sum(actual_loss, na.rm = TRUE),
  smear = smear
)
print(diag_tbl)


