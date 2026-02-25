## =========================================================
## Wildfire_M2_STAR_v2_occ5.R
##
## Purpose:
##  - Keep original M2 (Freq×Sev) code intact; this is an alternative M2' pipeline
##  - Build Occurrence-based 2-part model under the SAME 5 model family (M1~M5)
##    Part1: Occurrence probability P_hat via 5 models (ST/STAR variants)
##    Part2: Conditional positive total loss L^+ via STAR pooled SAR on log1p(loss)
##  - 2025 OOT: EP = P_hat * Lhat_pos
##  - Metrics + optional train-only calibration
##
## Notes:
##  - Part1 uses linear probability-style models; P_hat is clipped to [0,1]
##  - Part2 uses Duan smearing for inverse log1p
##  - AUC_occ uses P_hat (NOT EP)
##  - Islands/disjoint graphs: component-wise solve + fallback(power series)
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

OUT_RDS <- file.path(OUT_DIR, "m2_occ5_perf_final.rds")
OUT_CSV <- file.path(OUT_DIR, "m2_occ5_perf_final.csv")

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

# Panel-safe helper: apply single-year W within each year slice
apply_W_by_year <- function(x, year_vec, W_rs) {
  year_vec <- as.integer(year_vec)
  out <- numeric(length(x))
  yrs <- unique(year_vec)
  
  for (yy in yrs) {
    idx <- which(year_vec == yy)
    if (length(idx) != nrow(W_rs)) {
      stop(sprintf("apply_W_by_year: year=%s has %d rows, expected %d (check balancing/sorting).",
                   yy, length(idx), nrow(W_rs)))
    }
    out[idx] <- as.numeric(W_rs %*% as.numeric(x[idx]))
  }
  out
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
    
    if (length(x) == nrow(W_matrix)) {
      newdata[[lt]] <- as.numeric(W_matrix %*% x)
    } else {
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
  # inverse for log1p(Z): E[Z] approx = smear * exp(mu_hat) - 1
  pmax(smear * exp(mu_hat) - 1, 0)
}

# ---------------------------------------------------------
# 5) Metrics (AUC uses P_hat, NOT EP)
# ---------------------------------------------------------
evaluate_models_occ <- function(actual_cnt, actual_loss, p_hat, ep_hat, model_name,
                                tol = 1e-8, eps = 1e-8) {
  actual_cnt  <- as.numeric(actual_cnt)
  L <- pmax(as.numeric(actual_loss), 0)
  
  occ_actual <- as.numeric(actual_cnt > 0)
  
  p_hat <- pmin(pmax(as.numeric(p_hat), 0), 1)
  P <- pmax(as.numeric(ep_hat), 0)
  
  auc_val <- NA_real_
  if (length(unique(occ_actual)) >= 2) {
    auc_val <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(occ_actual, p_hat, quiet = TRUE))),
      error = function(e) NA_real_
    )
  }
  
  rmse   <- sqrt(mean((L - P)^2, na.rm = TRUE))
  mae    <- mean(abs(L - P), na.rm = TRUE)
  medae  <- median(abs(L - P), na.rm = TRUE)
  rmsle  <- sqrt(mean((log1p(L) - log1p(P))^2, na.rm = TRUE))
  spear  <- suppressWarnings(cor(L, P, method = "spearman", use = "complete.obs"))
  
  P0 <- P <= tol
  br_rate_down <- mean((L > 0) & P0, na.rm = TRUE)
  br_rate_up   <- mean((L == 0) & (!P0), na.rm = TRUE)
  
  br_sev_down <- sum(pmax(L - P, 0), na.rm = TRUE) / (sum(L, na.rm = TRUE) + eps)
  br_sev_up   <- sum(pmax(P - L, 0), na.rm = TRUE) / (sum(P, na.rm = TRUE) + eps)
  
  tibble(
    model_type   = model_name,
    AUC_occ      = auc_val,
    RMSE         = rmse,
    MAE          = mae,
    MedianAE     = medae,
    RMSLE        = rmsle,
    Spearman     = spear,
    BR_rate_down = br_rate_down,
    BR_rate_up   = br_rate_up,
    BR_sev_down  = br_sev_down,
    BR_sev_up    = br_sev_up
  )
}

safe_factor <- function(sum_actual, sum_pred) {
  if (!is.finite(sum_pred) || sum_pred <= 0) return(1.0)
  as.numeric(sum_actual / sum_pred)
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
    
    # targets for 2-part
    fire_any = as.numeric(fire_cnt_year > 0),
    y_loss   = log1p(damage_capped)
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
    fire_any_lag1 = lag(fire_any, 1),
    y_loss_lag1   = lag(y_loss, 1)
  ) %>%
  ungroup() %>%
  filter(!is.na(fire_any_lag1), !is.na(y_loss_lag1))

# spatial lag of occurrence-lag by year: Wy_fire_any_lag1
panel <- panel %>%
  group_by(year) %>%
  mutate(Wy_fire_any_lag1 = as.numeric(W_rs %*% fire_any_lag1)) %>%
  ungroup()

train_df <- panel %>% filter(year < TEST_YEAR) %>% arrange(year, sigungu_cd)
test_df  <- panel %>% filter(year == TEST_YEAR) %>% arrange(year, sigungu_cd)

stopifnot(nrow(test_df) == length(region_ids))

msg(">>> [3/7] Scaling for occurrence models (train stats only)...")
# Occurrence is binary, but we scale lag terms/covars (same style as original M2)
scale_cols <- c("fire_any_lag1", "Wy_fire_any_lag1", COVARS_LAG1)
scaled <- scale_with_train(train_df, test_df, scale_cols)
train_df <- scaled$train
test_df  <- scaled$test

msg(">>> [4/7] Fit Occurrence models (M1~M5, ST/STAR family)...")

# M1: Mixed-ST (OLS) on fire_any
f_m1 <- as.formula(paste(
  "fire_any ~ fire_any_lag1 + Wy_fire_any_lag1 +",
  paste(COVARS_LAG1, collapse=" + ")
))
m1_occ <- lm(f_m1, data=train_df)

# M2: Pure-ST (OLS) on fire_any
m2_occ <- lm(fire_any ~ fire_any_lag1 + Wy_fire_any_lag1, data=train_df)

# STAR pooled SAR over years: block-diagonal W for training
train_years <- sort(unique(train_df$year))
T_train <- length(train_years)
W_time_train <- build_block_listw(W_rs, T_train, zero_policy = TRUE)

# M3: Pure-STAR (SAR) on fire_any
m3_occ <- spatialreg::lagsarlm(
  fire_any ~ fire_any_lag1,
  data = train_df,
  listw = W_time_train,
  method = "LU",
  zero.policy = TRUE
)

# M4: Mixed-STAR (SAR + X) on fire_any
f_m4 <- as.formula(paste("fire_any ~ fire_any_lag1 +", paste(COVARS_LAG1, collapse=" + ")))
m4_occ <- spatialreg::lagsarlm(
  f_m4,
  data = train_df,
  listw = W_time_train,
  method = "LU",
  zero.policy = TRUE
)

# M5: Res-STAR (SDM/Durbin) on fire_any
m5_occ <- spatialreg::lagsarlm(
  f_m4,
  data = train_df,
  listw = W_time_train,
  method = "LU",
  zero.policy = TRUE,
  Durbin = as.formula(paste("~", paste(COVARS_LAG1, collapse=" + ")))
)

msg(">>> [5/7] Fit Conditional Loss model (STAR pooled SAR on y_loss)...")
# Part 2: log1p(total loss) with lag + fire_any + covars
f_loss <- as.formula(paste("y_loss ~ y_loss_lag1 + fire_any +", paste(COVARS_LAG1, collapse=" + ")))
m_loss <- spatialreg::lagsarlm(
  f_loss,
  data = train_df,
  listw = W_time_train,
  method = "LU",
  zero.policy = TRUE
)

smear_loss <- smearing_factor(m_loss)

msg(">>> [6a/7] Learn calibration factors from TRAIN (no leakage)...")
# Train: conditional positive loss prediction (set fire_any=1)
train_loss_nd <- train_df
train_loss_nd$fire_any <- 1

x_vars_loss <- c("y_loss_lag1", "fire_any", COVARS_LAG1)
mu_loss_tr_pos <- predict_mu_star(m_loss, train_loss_nd, W_rs, ps_order = PS_FALLBACK_ORDER)
Loss_tr_pos <- pmax(exp(mu_loss_tr_pos) * smear_loss - 1, 1e-8)

# Train: occurrence probability per model (clip)
P1_tr <- pmin(pmax(as.numeric(predict_mu_star(m1_occ, train_df, W_rs, ps_order = PS_FALLBACK_ORDER)), 0), 1)
P2_tr <- pmin(pmax(as.numeric(predict_mu_star(m2_occ, train_df, W_rs, ps_order = PS_FALLBACK_ORDER)), 0), 1)
P3_tr <- pmin(pmax(as.numeric(predict_mu_star(m3_occ, train_df, W_rs, ps_order = PS_FALLBACK_ORDER)), 0), 1)
P4_tr <- pmin(pmax(as.numeric(predict_mu_star(m4_occ, train_df, W_rs, ps_order = PS_FALLBACK_ORDER)), 0), 1)
P5_tr <- pmin(pmax(as.numeric(predict_mu_star(m5_occ, train_df, W_rs, ps_order = PS_FALLBACK_ORDER)), 0), 1)

EP1_tr <- pmax(P1_tr * Loss_tr_pos, 0)
EP2_tr <- pmax(P2_tr * Loss_tr_pos, 0)
EP3_tr <- pmax(P3_tr * Loss_tr_pos, 0)
EP4_tr <- pmax(P4_tr * Loss_tr_pos, 0)
EP5_tr <- pmax(P5_tr * Loss_tr_pos, 0)

actual_loss_train <- as.numeric(train_df$damage_capped)
sum_actual_tr <- sum(actual_loss_train, na.rm = TRUE)

c1 <- safe_factor(sum_actual_tr, sum(EP1_tr, na.rm=TRUE))
c2 <- safe_factor(sum_actual_tr, sum(EP2_tr, na.rm=TRUE))
c3 <- safe_factor(sum_actual_tr, sum(EP3_tr, na.rm=TRUE))
c4 <- safe_factor(sum_actual_tr, sum(EP4_tr, na.rm=TRUE))
c5 <- safe_factor(sum_actual_tr, sum(EP5_tr, na.rm=TRUE))

calib_factors <- tibble(
  model = c("M1: Mixed-ST (OLS)", "M2: Pure-ST (OLS)", "M3: Pure-STAR (SAR)",
            "M4: Mixed-STAR (SAR+X)", "M5: Res-STAR (SDM/Durbin)"),
  calib_factor = c(c1, c2, c3, c4, c5)
)

print(calib_factors)

msg(">>> [6/7] Predict 2025 (P_hat + conditional Loss_hat_pos)...")

# Test: conditional positive loss (set fire_any=1)
test_loss_nd <- test_df
test_loss_nd$fire_any <- 1

mu_loss_hat_pos <- predict_mu_star(m_loss, test_loss_nd, W_rs, ps_order = PS_FALLBACK_ORDER)
Loss_hat_pos <- pmax(exp(mu_loss_hat_pos) * smear_loss - 1, 1e-8)

# Test: occurrence probability per model (clip)
P1 <- pmin(pmax(as.numeric(predict_mu_star(m1_occ, test_df, W_rs, ps_order = PS_FALLBACK_ORDER)), 0), 1)
P2 <- pmin(pmax(as.numeric(predict_mu_star(m2_occ, test_df, W_rs, ps_order = PS_FALLBACK_ORDER)), 0), 1)
P3 <- pmin(pmax(as.numeric(predict_mu_star(m3_occ, test_df, W_rs, ps_order = PS_FALLBACK_ORDER)), 0), 1)
P4 <- pmin(pmax(as.numeric(predict_mu_star(m4_occ, test_df, W_rs, ps_order = PS_FALLBACK_ORDER)), 0), 1)
P5 <- pmin(pmax(as.numeric(predict_mu_star(m5_occ, test_df, W_rs, ps_order = PS_FALLBACK_ORDER)), 0), 1)

# EP
EP1 <- pmax(P1 * Loss_hat_pos, 0)
EP2 <- pmax(P2 * Loss_hat_pos, 0)
EP3 <- pmax(P3 * Loss_hat_pos, 0)
EP4 <- pmax(P4 * Loss_hat_pos, 0)
EP5 <- pmax(P5 * Loss_hat_pos, 0)

msg(">>> [7/7] Evaluate (RAW vs CAL) & save...")

actual_cnt  <- test_df$fire_cnt_year
actual_loss <- test_df$damage_capped

tol <- 1e-8
eps <- 1e-8

# Apply calibration to TEST EP (EP only)
EP1_adj <- EP1 * c1
EP2_adj <- EP2 * c2
EP3_adj <- EP3 * c3
EP4_adj <- EP4 * c4
EP5_adj <- EP5 * c5

res_raw <- bind_rows(
  evaluate_models_occ(actual_cnt, actual_loss, P1, EP1, "M1: Mixed-ST (OLS) [RAW]", tol=tol, eps=eps),
  evaluate_models_occ(actual_cnt, actual_loss, P2, EP2, "M2: Pure-ST (OLS) [RAW]", tol=tol, eps=eps),
  evaluate_models_occ(actual_cnt, actual_loss, P3, EP3, "M3: Pure-STAR (SAR) [RAW]", tol=tol, eps=eps),
  evaluate_models_occ(actual_cnt, actual_loss, P4, EP4, "M4: Mixed-STAR (SAR+X) [RAW]", tol=tol, eps=eps),
  evaluate_models_occ(actual_cnt, actual_loss, P5, EP5, "M5: Res-STAR (SDM/Durbin) [RAW]", tol=tol, eps=eps)
)

res_cal <- bind_rows(
  evaluate_models_occ(actual_cnt, actual_loss, P1, EP1_adj, "M1: Mixed-ST (OLS) [CAL]", tol=tol, eps=eps),
  evaluate_models_occ(actual_cnt, actual_loss, P2, EP2_adj, "M2: Pure-ST (OLS) [CAL]", tol=tol, eps=eps),
  evaluate_models_occ(actual_cnt, actual_loss, P3, EP3_adj, "M3: Pure-STAR (SAR) [CAL]", tol=tol, eps=eps),
  evaluate_models_occ(actual_cnt, actual_loss, P4, EP4_adj, "M4: Mixed-STAR (SAR+X) [CAL]", tol=tol, eps=eps),
  evaluate_models_occ(actual_cnt, actual_loss, P5, EP5_adj, "M5: Res-STAR (SDM/Durbin) [CAL]", tol=tol, eps=eps)
)

strip_suffix <- function(x) gsub("\\s*\\[(RAW|CAL)\\]\\s*$", "", x)

res <- bind_rows(res_raw, res_cal) %>%
  mutate(
    version    = if_else(grepl("\\[CAL\\]$", model_type), "CAL", "RAW"),
    model_base = strip_suffix(model_type)
  ) %>%
  left_join(
    calib_factors %>% rename(model_base = model),
    by = "model_base"
  ) %>%
  mutate(
    calib_factor = if_else(version == "CAL", calib_factor, NA_real_)
  ) %>%
  dplyr::select(-model_base)

print(res %>% arrange(desc(Spearman)))

saveRDS(res, OUT_RDS)
readr::write_csv(res, OUT_CSV)

msg("Saved: %s", OUT_RDS)
msg("Saved: %s", OUT_CSV)