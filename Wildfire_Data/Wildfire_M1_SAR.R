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

# ---------------------------------------------------------
# 5. 결과 저장: M0, M2와의## =========================================================
## Wildfire_M1_SAR_2part.R
##
## 목적:
##  - SAR 기반 2-part (Occurrence–Total Severity) 모델로
##    2025 Out-of-Time EP(=P_hat * E[Loss|positive]) 예측 및 성능 산출
##
## 핵심:
##  - splm::spml(lag=TRUE) 기반 SAR 패널 모형 2개
##    Part1: fire_any ~ X  (SAR-LPM; 예측값을 [0,1]로 clip)
##    Part2: y_loss=log1p(damage_capped) ~ fire_any + X (SAR; 예측 시 fire_any=1로 조건부화)
##  - islands(고립 8개) + disjoint subgraphs 대응:
##    component-wise solve + 실패 시 power series fallback
##  - 패널 균형화 + 결측치 대체는 train 통계로만 수행(누수 방지)
##
## 출력:
##  - Wildfire_Data/results/m1_sar_perf.rds
##  - Wildfire_Data/results/m1_sar_test_pred_2025.rds
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(Matrix)
  library(spdep)
  library(splm)
  library(pROC)
})

# ---------------------------------------------------------
# 0) User settings
# ---------------------------------------------------------
PATH_TRAIN <- "Wildfire_Data/processed_data/df_train_final.rds"
PATH_TEST  <- "Wildfire_Data/processed_data/df_test_final.rds"
PATH_W     <- "Wildfire_Data/meta_data/sgg_spatial_weights_252.rds"

OUT_DIR <- "Wildfire_Data/results"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

OUT_PERF_RDS <- file.path(OUT_DIR, "m1_sar_perf.rds")
OUT_PRED_RDS <- file.path(OUT_DIR, "m1_sar_test_pred_2025.rds")

TEST_YEAR <- 2025
PS_FALLBACK_ORDER <- 20

COVARS_LAG1 <- c(
  "spring_ws_p95_lag1",
  "spring_he_min_lag1",
  "max_dry_run_spring_lag1",
  "monsoon_prcp_sum_lag1",
  "winter_dry_days_lag1",
  "ta_sfc_mean_lag1"
)

msg <- function(...) message(sprintf(...))

# ---------------------------------------------------------
# 1) Build W (row-standardized; islands allowed)
# ---------------------------------------------------------
build_W <- function(W_dense, region_ids) {
  lw <- spdep::mat2listw(as.matrix(W_dense), style = "W", zero.policy = TRUE)
  attr(lw$neighbours, "region.id") <- region_ids
  W_rs <- Matrix::Matrix(spdep::listw2mat(lw), sparse = TRUE)
  dimnames(W_rs) <- list(region_ids, region_ids)
  list(listw = lw, W = W_rs)
}

# ---------------------------------------------------------
# 2) Balanced template + train-only imputation (no leakage)
# ---------------------------------------------------------
make_balanced_template <- function(df_all, ids, years) {
  df_all %>%
    mutate(sigungu_cd = as.character(sigungu_cd),
           year = as.integer(year)) %>%
    tidyr::complete(year = years, sigungu_cd = ids) %>%
    mutate(sigungu_cd = factor(sigungu_cd, levels = ids)) %>%
    arrange(year, sigungu_cd)
}

impute_with_train_stats <- function(panel_df, train_mask) {
  df <- panel_df
  num_cols <- names(df)[sapply(df, is.numeric)]
  
  train_df <- df[train_mask, , drop = FALSE]
  
  by_id_means <- train_df %>%
    group_by(sigungu_cd) %>%
    summarise(across(all_of(num_cols), ~ mean(.x, na.rm = TRUE)), .groups="drop")
  
  global_means <- sapply(train_df[, num_cols, drop = FALSE], mean, na.rm = TRUE)
  
  df <- df %>%
    left_join(by_id_means, by="sigungu_cd", suffix=c("", ".idmean"))
  
  for (cn in num_cols) {
    idm <- paste0(cn, ".idmean")
    df[[cn]] <- ifelse(is.na(df[[cn]]), df[[idm]], df[[cn]])
    df[[cn]] <- ifelse(is.na(df[[cn]]), global_means[[cn]], df[[cn]])
    df[[idm]] <- NULL
  }
  df
}

# ---------------------------------------------------------
# 3) Component-wise safe spatial filter: y = (I - rho W)^(-1) xb
# ---------------------------------------------------------
get_components_from_W <- function(W) {
  A <- (W != 0) | (Matrix::t(W) != 0)
  A <- Matrix::drop0(A)
  n <- nrow(A)
  
  comp_id <- integer(n)
  cur <- 0L
  
  neighbors <- vector("list", n)
  for (i in 1:n) {
    idx <- Matrix::which(A[i, ], arr.ind = TRUE)
    if (length(idx) == 0) {
      neighbors[[i]] <- integer(0)
    } else if (is.matrix(idx)) {
      neighbors[[i]] <- as.integer(idx[, 2])
    } else {
      neighbors[[i]] <- as.integer(idx)
    }
  }
  
  for (i in 1:n) {
    if (comp_id[i] != 0L) next
    cur <- cur + 1L
    
    stack <- c(i)
    comp_id[i] <- cur
    
    while (length(stack) > 0) {
      v <- stack[[length(stack)]]
      stack <- stack[-length(stack)]
      
      nb <- neighbors[[v]]
      if (length(nb) == 0) next
      
      new_nodes <- nb[comp_id[nb] == 0L]
      if (length(new_nodes) > 0) {
        comp_id[new_nodes] <- cur
        stack <- c(stack, new_nodes)
      }
    }
  }
  
  comp_id
}

spatial_filter_safe <- function(xb, W, rho, ps_order = 20) {
  n <- length(xb)
  if (is.na(rho) || rho == 0) return(as.numeric(xb))
  
  comp_id <- get_components_from_W(W)
  y <- numeric(n)
  
  for (cid in sort(unique(comp_id))) {
    idx <- which(comp_id == cid)
    if (length(idx) == 1) {
      y[idx] <- xb[idx]  # isolate
      next
    }
    
    Wc <- W[idx, idx, drop=FALSE]
    Ic <- Matrix::Diagonal(n=length(idx), x=1)
    
    ok <- TRUE
    yc <- tryCatch(
      as.numeric(Matrix::solve(Ic - rho * Wc, xb[idx])),
      error = function(e) { ok <<- FALSE; rep(NA_real_, length(idx)) }
    )
    
    if (!ok || any(!is.finite(yc))) {
      yk <- xb[idx]
      Wpow_x <- xb[idx]
      rho_k <- rho
      for (k in 1:ps_order) {
        Wpow_x <- as.numeric(Wc %*% Wpow_x)
        yk <- yk + rho_k * Wpow_x
        rho_k <- rho_k * rho
      }
      yc <- yk
    }
    
    y[idx] <- yc
  }
  
  as.numeric(y)
}

# ---------------------------------------------------------
# 4) spml extract + robust predict
# ---------------------------------------------------------
get_spml_rho <- function(fit) {
  if (!is.null(fit$arcoef)) {
    if ("lambda" %in% names(fit$arcoef)) return(as.numeric(fit$arcoef["lambda"]))
    return(as.numeric(fit$arcoef[1]))
  }
  if (!is.null(fit$errcomp) && "lambda" %in% names(fit$errcomp)) {
    return(as.numeric(fit$errcomp["lambda"]))
  }
  stop("Cannot extract SAR parameter (lambda/rho) from spml fit.")
}

get_coef_safe <- function(fit) {
  b <- tryCatch(coef(fit), error=function(e) NULL)
  if (!is.null(b)) return(b)
  if (!is.null(fit$coefficients)) return(fit$coefficients)
  stop("Cannot extract coefficients.")
}

predict_spml_sar <- function(fit, newdata, W, x_vars, ps_order = 20) {
  rho <- get_spml_rho(fit)
  b <- get_coef_safe(fit)
  
  intercept <- if ("(Intercept)" %in% names(b)) as.numeric(b["(Intercept)"]) else 0
  
  for (v in x_vars) if (!v %in% names(newdata)) newdata[[v]] <- 0
  X <- as.matrix(newdata[, x_vars, drop=FALSE])
  xb <- as.numeric(intercept + X %*% b[x_vars])
  xb[!is.finite(xb)] <- 0
  
  spatial_filter_safe(xb, W, rho, ps_order = ps_order)
}

# ---------------------------------------------------------
# 5) Metrics (M1 format)
# ---------------------------------------------------------
evaluate_m1 <- function(actual_loss, occ_actual, occ_score, ep_hat) {
  actual_loss <- as.numeric(actual_loss)
  occ_actual  <- as.numeric(occ_actual)
  occ_score   <- as.numeric(occ_score)
  ep_hat      <- pmax(as.numeric(ep_hat), 0)
  
  auc_val <- NA_real_
  if (length(unique(occ_actual)) >= 2) {
    auc_val <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(occ_actual, occ_score, quiet = TRUE))),
      error = function(e) NA_real_
    )
  }
  
  rmse  <- sqrt(mean((actual_loss - ep_hat)^2, na.rm = TRUE))
  mae   <- mean(abs(actual_loss - ep_hat), na.rm = TRUE)
  medae <- median(abs(actual_loss - ep_hat), na.rm = TRUE)
  
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
  
  calib <- sum(ep_hat, na.rm = TRUE) / sum(actual_loss, na.rm = TRUE)
  
  tibble(
    model_type = "M1_SAR_2PART",
    AUC_occ = auc_val,
    RMSE = rmse,
    MAE = mae,
    MedAE = medae,
    Spearman = spearman,
    Calib_Ratio = calib
  )
}

# ---------------------------------------------------------
# 6) Main
# ---------------------------------------------------------
msg(">>> [1/5] Loading data & W...")
df_train_raw <- readRDS(PATH_TRAIN)
df_test_raw  <- readRDS(PATH_TEST)
df_all_raw   <- bind_rows(df_train_raw, df_test_raw)

W_obj <- readRDS(PATH_W)
region_ids <- as.character(W_obj$sigungu_cd)
years_all  <- sort(unique(as.integer(df_all_raw$year)))

W_pack <- build_W(W_obj$W, region_ids)
W_listw <- W_pack$listw
W_rs <- W_pack$W

msg(">>> Regions: %d | Years: %d (%d-%d)",
    length(region_ids), length(years_all), min(years_all), max(years_all))

msg(">>> [2/5] Balanced panel + targets (no leakage impute)...")
panel <- make_balanced_template(df_all_raw, region_ids, years_all)

panel <- panel %>%
  mutate(
    fire_cnt_year = as.numeric(fire_cnt_year),
    damage_capped = as.numeric(damage_capped),
    fire_cnt_year = tidyr::replace_na(fire_cnt_year, 0),
    damage_capped = tidyr::replace_na(damage_capped, 0),
    fire_any = if_else(damage_capped > 0, 1, 0),
    y_loss = log1p(damage_capped)
  )

train_mask <- panel$year < TEST_YEAR
panel <- impute_with_train_stats(panel, train_mask)

missing_cov <- setdiff(COVARS_LAG1, names(panel))
if (length(missing_cov) > 0) stop("Missing covariates: ", paste(missing_cov, collapse=", "))

train_df <- panel %>% filter(year < TEST_YEAR) %>% arrange(sigungu_cd, year)
test_df  <- panel %>% filter(year == TEST_YEAR) %>% arrange(sigungu_cd, year)

stopifnot(nrow(test_df) == length(region_ids))

msg(">>> [3/5] Fit SAR 2-part (Occurrence + Total Loss)...")

# Part1: Occurrence (SAR-LPM)
f_occ <- as.formula(paste("fire_any ~", paste(COVARS_LAG1, collapse=" + ")))
m_occ_sar <- splm::spml(
  f_occ, data = train_df,
  index = c("sigungu_cd", "year"),
  listw = W_listw,
  model = "random",
  lag = TRUE
)

# Part2: Total loss on full panel (match W), then condition-on-positive by setting fire_any=1 at prediction time
f_loss <- as.formula(paste("y_loss ~ fire_any +", paste(COVARS_LAG1, collapse=" + ")))
m_loss_sar <- splm::spml(
  f_loss, data = train_df,
  index = c("sigungu_cd", "year"),
  listw = W_listw,
  model = "random",
  lag = TRUE
)

msg(">>> [4/5] Predict 2025...")
# Part1 prediction (clip to [0,1])
occ_score <- predict_spml_sar(m_occ_sar, test_df, W_rs, COVARS_LAG1, ps_order = PS_FALLBACK_ORDER)
P_hat_sar <- pmin(pmax(occ_score, 0), 1)

# Part2 prediction: conditional positive loss
test_loss_nd <- test_df
test_loss_nd$fire_any <- 1
x_vars_loss <- c("fire_any", COVARS_LAG1)

y_loss_hat_pos <- predict_spml_sar(m_loss_sar, test_loss_nd, W_rs, x_vars_loss, ps_order = PS_FALLBACK_ORDER)

# Smearing for log1p inverse
smear <- 1.0
resid_loss <- tryCatch(as.numeric(residuals(m_loss_sar)), error=function(e) NULL)
if (!is.null(resid_loss)) smear <- mean(exp(resid_loss), na.rm = TRUE)

Loss_hat_pos_sar <- pmax(exp(y_loss_hat_pos) * smear - 1, 1e-8)

# EP
EP_hat_sar <- pmax(P_hat_sar * Loss_hat_pos_sar, 0)

msg(">>> [5/5] Evaluate & save...")

actual_loss <- as.numeric(test_df$damage_capped)
actual_occ  <- as.numeric(test_df$fire_any)

metrics_m1 <- evaluate_m1(
  actual_loss = actual_loss,
  occ_actual  = actual_occ,
  occ_score   = P_hat_sar,
  ep_hat      = EP_hat_sar
)

print(metrics_m1)

# Save performance table (keep legacy file name)
saveRDS(metrics_m1, OUT_PERF_RDS)

# Save 2025 predictions (for mapping/diagnostics)
test_2025_final <- test_df %>%
  mutate(
    P_hat_sar = P_hat_sar,
    Loss_hat_pos_sar = Loss_hat_pos_sar,
    EP_M1_SAR = EP_hat_sar
  )

saveRDS(test_2025_final, OUT_PRED_RDS)

msg("Saved: %s", OUT_PERF_RDS)
msg("Saved: %s", OUT_PRED_RDS) 통합 비교를 위한 RDS 출력
# ---------------------------------------------------------
message(">>> [5/5] M1_SAR 성능 지표 및 예측 결과 저장 중...")

# (1) 성능 지표 테이블 저장
# M0와 컬럼명을 일치시켜 저장합니다. (metrics_m1_full 객체 기준)
saveRDS(metrics_m1, "Wildfire_Data/results/m1_sar_perf.rds")

# (2) (선택 사항) 2025년 예측치가 포함된 전체 데이터프레임 저장
# 나중에 지도(Map)를 그리거나 잔차 분석을 할 때 요긴하게 쓰입니다.
saveRDS(test_2025_final, "Wildfire_Data/results/m1_sar_test_pred_2025.rds")

cat("\n--- M1_SAR 결과 저장 완료 (Wildfire_Data/results/) ---\n")