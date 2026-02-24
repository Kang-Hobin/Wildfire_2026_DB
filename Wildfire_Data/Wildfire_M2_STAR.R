## =========================================================
## Wildfire_M2_STAR_v2.R
##
## 목적:
##  (1) 시군구-연도 패널 균형화 + train 통계로 결측대체(누수 방지)
##  (2) 빈도(N) STAR 계열 5개 모델(M1~M5) 적합
##  (3) 심도(S: 건당 피해) SAR 모델 1개 적합
##  (4) 2025 OOT 예측: EP = N_hat * Sev_hat
##  (5) 성능지표 저장
##
## 핵심 안정성:
##  - islands / disjoint subgraphs: component-wise solve + fallback(power series)
##  - AUC_occ는 EP가 아니라 N_hat로 계산 (정의 정상화)
##  - 역변환 bias correction: Duan smearing(모델별 잔차 기반) 사용
##    (분산보정 exp(mu+s2/2)은 객체/가정 불일치가 잦아 v2에선 배제)
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

# ---------------------------------------------------------
# 0) Settings
# ---------------------------------------------------------
PATH_TRAIN <- "Wildfire_Data/processed_data/df_train_final.rds"
PATH_TEST  <- "Wildfire_Data/processed_data/df_test_final.rds"
PATH_W     <- "Wildfire_Data/meta_data/sgg_spatial_weights_252.rds"

OUT_DIR <- "Wildfire_Data/results"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

OUT_RDS <- file.path(OUT_DIR, "m2_star_5models_perf_final.rds")
OUT_CSV <- file.path(OUT_DIR, "m2_star_5models_perf_final.csv")

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
# 1) W helpers (row-standardized; islands allowed)
# ---------------------------------------------------------
build_W <- function(W_dense, region_ids) {
  lw <- spdep::mat2listw(as.matrix(W_dense), style="W", zero.policy=TRUE)
  attr(lw$neighbours, "region.id") <- region_ids
  W_rs <- Matrix::Matrix(spdep::listw2mat(lw), sparse=TRUE)
  dimnames(W_rs) <- list(region_ids, region_ids)
  list(listw = lw, W = W_rs)
}

# block-diagonal W for pooled SAR over years (training only)
build_block_listw <- function(W_rs, T_train, zero_policy = TRUE) {
  W_block <- Matrix::kronecker(Matrix::Diagonal(T_train), W_rs)
  spdep::mat2listw(as.matrix(W_block), style="W", zero.policy = zero_policy)
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
  
  train_df <- df[train_mask, , drop=FALSE]
  
  by_id_means <- train_df %>%
    group_by(sigungu_cd) %>%
    summarise(across(all_of(num_cols), ~ mean(.x, na.rm=TRUE)), .groups="drop")
  
  global_means <- sapply(train_df[, num_cols, drop=FALSE], mean, na.rm=TRUE)
  
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

scale_with_train <- function(train_df, test_df, cols) {
  mu  <- sapply(train_df[, cols, drop=FALSE], mean, na.rm=TRUE)
  sdv <- sapply(train_df[, cols, drop=FALSE], sd, na.rm=TRUE)
  sdv[sdv == 0] <- 1
  
  for (cn in cols) {
    train_df[[cn]] <- (train_df[[cn]] - mu[cn]) / sdv[cn]
    test_df[[cn]]  <- (test_df[[cn]]  - mu[cn]) / sdv[cn]
  }
  list(train=train_df, test=test_df, mu=mu, sd=sdv)
}

# ---------------------------------------------------------
# 3) Component-wise safe spatial filter: y = (I - rho W)^(-1) xb
#    (islands + disjoint subgraphs safe)
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
      y[idx] <- xb[idx]
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
# 4) Robust prediction for lm / sarlm (lagsarlm)
#    - build xb by matching coefficient names to newdata columns
#    - Durbin lag.X terms auto-generated: "lag.var" = W %*% var
#    - spatial reduced form applied using single-year W_rs (consistent for OOT year)
# ---------------------------------------------------------
get_spatial_param <- function(fit) {
  if (!is.null(fit$rho)) return(as.numeric(fit$rho))
  if (!is.null(fit$lambda)) return(as.numeric(fit$lambda))
  return(0)
}

get_coef_safe <- function(fit) {
  b <- tryCatch(coef(fit), error=function(e) NULL)
  if (!is.null(b)) return(b)
  if (!is.null(fit$coefficients)) return(fit$coefficients)
  stop("Cannot extract coefficients.")
}

ensure_lag_terms <- function(newdata, coef_names, W_matrix) {
  lag_terms <- coef_names[startsWith(coef_names, "lag.")]
  if (length(lag_terms) == 0) return(newdata)
  
  for (lt in lag_terms) {
    if (lt %in% names(newdata)) next
    
    v <- sub("^lag\\.", "", lt)
    
    if (!v %in% names(newdata)) {
      newdata[[lt]] <- 0
      next
    }
    
    x <- as.numeric(newdata[[v]])
    
    # Case 1) single-year slice: length matches W dimension
    if (length(x) == nrow(W_matrix)) {
      newdata[[lt]] <- as.numeric(W_matrix %*% x)
      
    } else {
      # Case 2) panel stacked: compute year-wise W * x_year
      if (!("year" %in% names(newdata))) {
        stop("ensure_lag_terms: panel-length detected but 'year' column not found.")
      }
      newdata[[lt]] <- apply_W_by_year(x, newdata$year, W_matrix)
    }
  }
  
  newdata
}

predict_mu_star <- function(fit, newdata, W_matrix, ps_order = 20) {
  b <- get_coef_safe(fit)
  rho <- get_spatial_param(fit)
  
  coef_names <- names(b)
  newdata <- ensure_lag_terms(newdata, coef_names, W_matrix)
  
  vars <- setdiff(coef_names, "(Intercept)")
  missing_vars <- setdiff(vars, names(newdata))
  if (length(missing_vars) > 0) {
    for (v in missing_vars) newdata[[v]] <- 0
  }
  
  X <- as.matrix(newdata[, vars, drop=FALSE])
  xb <- if ("(Intercept)" %in% coef_names) as.numeric(b["(Intercept)"]) else 0
  xb <- xb + as.numeric(X %*% b[vars])
  xb[!is.finite(xb)] <- 0
  
  # If panel-stacked, apply reduced form within each year slice (block diagonal logic)
  if (length(xb) == nrow(W_matrix)) {
    return(spatial_filter_safe(xb, W_matrix, rho, ps_order = ps_order))
  } else {
    if (!("year" %in% names(newdata))) stop("predict_mu_star: panel-length detected but 'year' column not found.")
    y <- numeric(length(xb))
    yrs <- unique(as.integer(newdata$year))
    for (yy in yrs) {
      idx <- which(as.integer(newdata$year) == yy)
      if (length(idx) != nrow(W_matrix)) {
        stop(sprintf("predict_mu_star: year=%s has %d rows, expected %d (check balancing/sorting).",
                     yy, length(idx), nrow(W_matrix)))
      }
      y[idx] <- spatial_filter_safe(xb[idx], W_matrix, rho, ps_order = ps_order)
    }
    return(as.numeric(y))
  }
}

# Duan smearing (log1p scale -> original)
smearing_factor <- function(fit) {
  r <- tryCatch(as.numeric(residuals(fit)), error=function(e) NULL)
  if (is.null(r) || length(r) == 0) return(1.0)
  mean(exp(r), na.rm = TRUE)
}

inv_log1p_smear <- function(mu_hat, smear) {
  # E[Z] approx = smear * exp(mu_hat) - 1   (with Z = expm1(Y))
  pmax(smear * exp(mu_hat) - 1, 0)
}

# ---------------------------------------------------------
# 5) Metrics (AUC uses N_hat, NOT EP)
# ---------------------------------------------------------
evaluate_models <- function(actual_cnt, actual_loss, n_hat, ep_hat, model_name) {
  actual_cnt <- as.numeric(actual_cnt)
  actual_loss <- as.numeric(actual_loss)
  
  occ_actual <- as.numeric(actual_cnt > 0)
  n_hat  <- pmax(as.numeric(n_hat), 0)
  ep_hat <- pmax(as.numeric(ep_hat), 0)
  
  # AUC uses N_hat (NOT EP)
  auc_val <- NA_real_
  if (length(unique(occ_actual)) >= 2) {
    auc_val <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(occ_actual, n_hat, quiet = TRUE))),
      error = function(e) NA_real_
    )
  }
  
  rmse  <- sqrt(mean((actual_loss - ep_hat)^2, na.rm=TRUE))
  mae   <- mean(abs(actual_loss - ep_hat), na.rm=TRUE)
  medae <- median(abs(actual_loss - ep_hat), na.rm=TRUE)
  
  idx_inf <- which((actual_loss > 0) | (ep_hat > 0))
  spearman <- NA_real_
  if (length(idx_inf) >= 3) {
    a <- actual_loss[idx_inf]
    p <- ep_hat[idx_inf]
    if (sd(a, na.rm=TRUE) > 0 && sd(p, na.rm=TRUE) > 0) {
      spearman <- suppressWarnings(cor(a, p, method="spearman", use="complete.obs"))
    } else {
      spearman <- 0
    }
  }
  
  calib <- sum(ep_hat, na.rm=TRUE) / sum(actual_loss, na.rm=TRUE)
  
  tibble(
    model = model_name,
    AUC_occ = auc_val,
    RMSE = rmse,
    MAE = mae,
    MedAE = medae,
    Spearman = spearman,
    Calib = calib,
    Pred_ZeroShare = mean(ep_hat == 0, na.rm=TRUE)
  )
}

# ---------------------------------------------------------
# 6) Main pipeline
# ---------------------------------------------------------
msg(">>> [1/7] Loading data & W...")
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

msg(">>> [2/7] Balancing panel + targets (no leakage impute)...")
panel <- make_balanced_template(df_all_raw, region_ids, years_all) %>%
  mutate(
    fire_cnt_year = as.numeric(fire_cnt_year),
    damage_capped = as.numeric(damage_capped),
    fire_cnt_year = tidyr::replace_na(fire_cnt_year, 0),
    damage_capped = tidyr::replace_na(damage_capped, 0),
    
    # targets
    y_freq = log1p(fire_cnt_year),
    y_sev  = log1p(damage_capped / pmax(fire_cnt_year, 1))
  )

train_mask <- panel$year < TEST_YEAR
panel <- impute_with_train_stats(panel, train_mask)

missing_cov <- setdiff(COVARS_LAG1, names(panel))
if (length(missing_cov) > 0) stop("Missing covariates: ", paste(missing_cov, collapse=", "))

# temporal lags
panel <- panel %>%
  group_by(sigungu_cd) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    y_freq_lag1 = lag(y_freq, 1),
    y_sev_lag1  = lag(y_sev, 1)
  ) %>%
  ungroup() %>%
  filter(!is.na(y_freq_lag1), !is.na(y_sev_lag1))

# spatial lag Wy_freq_lag1 by year
panel <- panel %>%
  group_by(year) %>%
  mutate(Wy_freq_lag1 = as.numeric(W_rs %*% y_freq_lag1)) %>%
  ungroup()

train_df <- panel %>% filter(year < TEST_YEAR) %>% arrange(year, sigungu_cd)
test_df  <- panel %>% filter(year == TEST_YEAR) %>% arrange(year, sigungu_cd)

stopifnot(nrow(test_df) == length(region_ids))

msg(">>> [3/7] Scaling for frequency models (train stats only)...")
scale_cols <- c("y_freq_lag1", "Wy_freq_lag1", COVARS_LAG1)
scaled <- scale_with_train(train_df, test_df, scale_cols)
train_df <- scaled$train
test_df  <- scaled$test

msg(">>> [4/7] Fit Frequency models (M1~M5)...")

# M1: Mixed-ST (OLS)
f_m1 <- as.formula(paste(
  "y_freq ~ y_freq_lag1 + Wy_freq_lag1 +",
  paste(COVARS_LAG1, collapse=" + ")
))
m1_st <- lm(f_m1, data=train_df)

# M2: Pure-ST (OLS)
m2_st <- lm(y_freq ~ y_freq_lag1 + Wy_freq_lag1, data=train_df)

# STAR pooled SAR over years: block-diagonal W for training
train_years <- sort(unique(train_df$year))
T_train <- length(train_years)
W_time_train <- build_block_listw(W_rs, T_train, zero_policy = TRUE)

# M3: Pure-STAR (SAR)
m3_star <- spatialreg::lagsarlm(
  y_freq ~ y_freq_lag1,
  data = train_df,
  listw = W_time_train,
  method = "LU",
  zero.policy = TRUE
)

# M4: Mixed-STAR (SAR + X)
f_m4 <- as.formula(paste("y_freq ~ y_freq_lag1 +", paste(COVARS_LAG1, collapse=" + ")))
m4_star <- spatialreg::lagsarlm(
  f_m4,
  data = train_df,
  listw = W_time_train,
  method = "LU",
  zero.policy = TRUE
)

# M5: Res-STAR (SDM/Durbin)
m5_star <- spatialreg::lagsarlm(
  f_m4,
  data = train_df,
  listw = W_time_train,
  method = "LU",
  zero.policy = TRUE,
  Durbin = as.formula(paste("~", paste(COVARS_LAG1, collapse=" + ")))
)

msg(">>> [5/7] Fit Severity model (SAR on y_sev)...")
# Severity pooled SAR over years (same block W idea)
f_sev <- as.formula(paste("y_sev ~ y_sev_lag1 +", paste(COVARS_LAG1, collapse=" + ")))
sev_star <- spatialreg::lagsarlm(
  f_sev,
  data = train_df,
  listw = W_time_train,
  method = "LU",
  zero.policy = TRUE
)

# Smearing factors (model-specific)
smear_m1 <- smearing_factor(m1_st)
smear_m2 <- smearing_factor(m2_st)
smear_m3 <- smearing_factor(m3_star)
smear_m4 <- smearing_factor(m4_star)
smear_m5 <- smearing_factor(m5_star)
smear_sev <- smearing_factor(sev_star)

msg(">>> [6a/7] Learn calibration factors from TRAIN (no leakage)...")

## Addon : Calibration Adjustment ---
# Panel-safe helper: apply single-year W within each year slice
# Assumes newdata is sorted by (year, sigungu_cd) OR at least grouped contiguous by year.
apply_W_by_year <- function(x, year_vec, W_rs) {
  year_vec <- as.integer(year_vec)
  out <- numeric(length(x))
  yrs <- unique(year_vec)
  
  for (yy in yrs) {
    idx <- which(year_vec == yy)
    # Expect idx length equals nrow(W_rs)=252
    if (length(idx) != nrow(W_rs)) {
      stop(sprintf("apply_W_by_year: year=%s has %d rows, expected %d (check balancing/sorting).",
                   yy, length(idx), nrow(W_rs)))
    }
    out[idx] <- as.numeric(W_rs %*% as.numeric(x[idx]))
  }
  out
}

# --- Predict severity on TRAIN ---
mu_sev_train <- predict_mu_star(sev_star, train_df, W_rs, ps_order = PS_FALLBACK_ORDER)
sev_hat_train <- inv_log1p_smear(mu_sev_train, smear_sev)
sev_hat_train <- pmax(sev_hat_train, 1e-8)

# --- Predict frequency on TRAIN for each model ---
mu1_tr <- predict_mu_star(m1_st,   train_df, W_rs, ps_order = PS_FALLBACK_ORDER)
mu2_tr <- predict_mu_star(m2_st,   train_df, W_rs, ps_order = PS_FALLBACK_ORDER)
mu3_tr <- predict_mu_star(m3_star, train_df, W_rs, ps_order = PS_FALLBACK_ORDER)
mu4_tr <- predict_mu_star(m4_star, train_df, W_rs, ps_order = PS_FALLBACK_ORDER)
mu5_tr <- predict_mu_star(m5_star, train_df, W_rs, ps_order = PS_FALLBACK_ORDER)

n1_tr <- inv_log1p_smear(mu1_tr, smear_m1)
n2_tr <- inv_log1p_smear(mu2_tr, smear_m2)
n3_tr <- inv_log1p_smear(mu3_tr, smear_m3)
n4_tr <- inv_log1p_smear(mu4_tr, smear_m4)
n5_tr <- inv_log1p_smear(mu5_tr, smear_m5)

ep1_tr <- pmax(n1_tr * sev_hat_train, 0)
ep2_tr <- pmax(n2_tr * sev_hat_train, 0)
ep3_tr <- pmax(n3_tr * sev_hat_train, 0)
ep4_tr <- pmax(n4_tr * sev_hat_train, 0)
ep5_tr <- pmax(n5_tr * sev_hat_train, 0)

actual_loss_train <- as.numeric(train_df$damage_capped)

sum_actual_tr <- sum(actual_loss_train, na.rm = TRUE)

safe_factor <- function(sum_actual, sum_pred) {
  if (!is.finite(sum_pred) || sum_pred <= 0) return(1.0)
  as.numeric(sum_actual / sum_pred)
}

c1 <- safe_factor(sum_actual_tr, sum(ep1_tr, na.rm=TRUE))
c2 <- safe_factor(sum_actual_tr, sum(ep2_tr, na.rm=TRUE))
c3 <- safe_factor(sum_actual_tr, sum(ep3_tr, na.rm=TRUE))
c4 <- safe_factor(sum_actual_tr, sum(ep4_tr, na.rm=TRUE))
c5 <- safe_factor(sum_actual_tr, sum(ep5_tr, na.rm=TRUE))

calib_factors <- tibble(
  model = c("M1: Mixed-ST (OLS)", "M2: Pure-ST (OLS)", "M3: Pure-STAR (SAR)", "M4: Mixed-STAR (SAR+X)", "M5: Res-STAR (SDM/Durbin)"),
  calib_factor = c(c1, c2, c3, c4, c5)
)

print(calib_factors)

msg(">>> [6/7] Predict 2025 (mu -> inverse log1p with smearing)...")

# Severity prediction: mu on log1p scale, then inverse with smearing
mu_sev_hat <- predict_mu_star(sev_star, test_df, W_rs, ps_order = PS_FALLBACK_ORDER)
sev_hat <- inv_log1p_smear(mu_sev_hat, smear_sev)
# extra guard
sev_hat <- pmax(sev_hat, 1e-8)

# Frequency predictions (N_hat) for each model
mu1 <- predict_mu_star(m1_st,   test_df, W_rs, ps_order = PS_FALLBACK_ORDER)
mu2 <- predict_mu_star(m2_st,   test_df, W_rs, ps_order = PS_FALLBACK_ORDER)
mu3 <- predict_mu_star(m3_star, test_df, W_rs, ps_order = PS_FALLBACK_ORDER)
mu4 <- predict_mu_star(m4_star, test_df, W_rs, ps_order = PS_FALLBACK_ORDER)
mu5 <- predict_mu_star(m5_star, test_df, W_rs, ps_order = PS_FALLBACK_ORDER)

n1 <- inv_log1p_smear(mu1, smear_m1)
n2 <- inv_log1p_smear(mu2, smear_m2)
n3 <- inv_log1p_smear(mu3, smear_m3)
n4 <- inv_log1p_smear(mu4, smear_m4)
n5 <- inv_log1p_smear(mu5, smear_m5)

# Combine EP
ep1 <- pmax(n1 * sev_hat, 0)
ep2 <- pmax(n2 * sev_hat, 0)
ep3 <- pmax(n3 * sev_hat, 0)
ep4 <- pmax(n4 * sev_hat, 0)
ep5 <- pmax(n5 * sev_hat, 0)

msg(">>> [7/7] Evaluate (raw vs calibrated) & save...")

actual_cnt  <- test_df$fire_cnt_year
actual_loss <- test_df$damage_capped

# --- Apply calibration to TEST EP (EP only) ---
ep1_adj <- ep1 * c1
ep2_adj <- ep2 * c2
ep3_adj <- ep3 * c3
ep4_adj <- ep4 * c4
ep5_adj <- ep5 * c5

res_raw <- bind_rows(
  evaluate_models(actual_cnt, actual_loss, n1, ep1, "M1: Mixed-ST (OLS) [RAW]"),
  evaluate_models(actual_cnt, actual_loss, n2, ep2, "M2: Pure-ST (OLS) [RAW]"),
  evaluate_models(actual_cnt, actual_loss, n3, ep3, "M3: Pure-STAR (SAR) [RAW]"),
  evaluate_models(actual_cnt, actual_loss, n4, ep4, "M4: Mixed-STAR (SAR+X) [RAW]"),
  evaluate_models(actual_cnt, actual_loss, n5, ep5, "M5: Res-STAR (SDM/Durbin) [RAW]")
)

res_adj <- bind_rows(
  evaluate_models(actual_cnt, actual_loss, n1, ep1_adj, "M1: Mixed-ST (OLS) [CAL]"),
  evaluate_models(actual_cnt, actual_loss, n2, ep2_adj, "M2: Pure-ST (OLS) [CAL]"),
  evaluate_models(actual_cnt, actual_loss, n3, ep3_adj, "M3: Pure-STAR (SAR) [CAL]"),
  evaluate_models(actual_cnt, actual_loss, n4, ep4_adj, "M4: Mixed-STAR (SAR+X) [CAL]"),
  evaluate_models(actual_cnt, actual_loss, n5, ep5_adj, "M5: Res-STAR (SDM/Durbin) [CAL]")
)

strip_suffix <- function(x) gsub("\\s*\\[(RAW|CAL)\\]\\s*$", "", x)

res <- bind_rows(res_raw, res_adj)

strip_suffix <- function(x) gsub("\\s*\\[(RAW|CAL)\\]\\s*$", "", x)

res <- res %>%
  mutate(
    version = ifelse(grepl("\\[CAL\\]$", model), "CAL", "RAW"),
    model_base = strip_suffix(model)
  ) %>%
  left_join(
    calib_factors %>% rename(model_base = model),
    by = "model_base"
  ) %>%
  mutate(
    calib_factor = ifelse(version == "CAL", calib_factor, NA_real_)
  ) %>%
  dplyr::select(-model_base)

print(res %>% arrange(desc(Spearman)))

saveRDS(res, OUT_RDS)
write.csv(res, OUT_CSV, row.names = FALSE)

msg("Saved: %s", OUT_RDS)
msg("Saved: %s", OUT_CSV)

