## =========================================================
## Wildfire_M1_SAR_additive_tau_grid.R
## 목적: SAR 2-part에서 Loss 파트를 additive(small+excess)로 분해하여
##      tau grid 성능 비교
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

PATH_TRAIN <- "Wildfire_Data/processed_data/df_train_final.rds"
PATH_TEST  <- "Wildfire_Data/processed_data/df_test_final.rds"
PATH_W     <- "Wildfire_Data/meta_data/sgg_spatial_weights_252.rds"

OUT_DIR <- "Wildfire_Data/results"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
OUT_RDS <- file.path(OUT_DIR, "m1_sar_additive_tau_grid_perf.rds")
OUT_CSV <- file.path(OUT_DIR, "m1_sar_additive_tau_grid_perf.csv")

TEST_YEAR <- 2025
PS_FALLBACK_ORDER <- 20
TAU_GRID <- seq(5, 50, by = 5)

COVARS_LAG1 <- c(
  "spring_ws_p95_lag1","spring_he_min_lag1","max_dry_run_spring_lag1",
  "monsoon_prcp_sum_lag1","winter_dry_days_lag1","ta_sfc_mean_lag1"
)

msg <- function(...) message(sprintf(...))

# ---------- W ----------
build_W <- function(W_dense, region_ids) {
  lw <- spdep::mat2listw(as.matrix(W_dense), style="W", zero.policy=TRUE)
  attr(lw$neighbours, "region.id") <- region_ids
  W_rs <- Matrix::Matrix(spdep::listw2mat(lw), sparse=TRUE)
  dimnames(W_rs) <- list(region_ids, region_ids)
  list(listw=lw, W=W_rs)
}

# ---------- panel template + train-only impute ----------
make_balanced_template <- function(df_all, ids, years) {
  df_all %>%
    mutate(sigungu_cd=as.character(sigungu_cd), year=as.integer(year)) %>%
    tidyr::complete(year=years, sigungu_cd=ids) %>%
    mutate(sigungu_cd=factor(sigungu_cd, levels=ids)) %>%
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
  
  df <- df %>% left_join(by_id_means, by="sigungu_cd", suffix=c("", ".idmean"))
  for (cn in num_cols) {
    idm <- paste0(cn, ".idmean")
    df[[cn]] <- ifelse(is.na(df[[cn]]), df[[idm]], df[[cn]])
    df[[cn]] <- ifelse(is.na(df[[cn]]), global_means[[cn]], df[[cn]])
    df[[idm]] <- NULL
  }
  df
}

# ---------- component-wise filter (same as your M1_SAR_2part) ----------
get_components_from_W <- function(W) {
  A <- (W != 0) | (Matrix::t(W) != 0)
  A <- Matrix::drop0(A)
  n <- nrow(A)
  
  comp_id <- integer(n)
  cur <- 0L
  
  neighbors <- vector("list", n)
  for (i in 1:n) {
    idx <- Matrix::which(A[i, ], arr.ind=TRUE)
    if (length(idx)==0) neighbors[[i]] <- integer(0)
    else if (is.matrix(idx)) neighbors[[i]] <- as.integer(idx[,2])
    else neighbors[[i]] <- as.integer(idx)
  }
  
  for (i in 1:n) {
    if (comp_id[i]!=0L) next
    cur <- cur + 1L
    stack <- c(i); comp_id[i] <- cur
    while (length(stack)>0) {
      v <- stack[[length(stack)]]
      stack <- stack[-length(stack)]
      nb <- neighbors[[v]]
      if (length(nb)==0) next
      new_nodes <- nb[comp_id[nb]==0L]
      if (length(new_nodes)>0) {
        comp_id[new_nodes] <- cur
        stack <- c(stack, new_nodes)
      }
    }
  }
  comp_id
}

spatial_filter_safe <- function(xb, W, rho, ps_order=20) {
  n <- length(xb)
  if (is.na(rho) || rho==0) return(as.numeric(xb))
  comp_id <- get_components_from_W(W)
  y <- numeric(n)
  
  for (cid in sort(unique(comp_id))) {
    idx <- which(comp_id==cid)
    if (length(idx)==1) { y[idx] <- xb[idx]; next }
    
    Wc <- W[idx, idx, drop=FALSE]
    Ic <- Matrix::Diagonal(n=length(idx), x=1)
    
    ok <- TRUE
    yc <- tryCatch(
      as.numeric(Matrix::solve(Ic - rho*Wc, xb[idx])),
      error=function(e){ ok<<-FALSE; rep(NA_real_, length(idx)) }
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

get_spml_rho <- function(fit) {
  if (!is.null(fit$arcoef)) {
    if ("lambda" %in% names(fit$arcoef)) return(as.numeric(fit$arcoef["lambda"]))
    return(as.numeric(fit$arcoef[1]))
  }
  if (!is.null(fit$errcomp) && "lambda" %in% names(fit$errcomp)) {
    return(as.numeric(fit$errcomp["lambda"]))
  }
  stop("Cannot extract SAR parameter from spml.")
}

get_coef_safe <- function(fit) {
  b <- tryCatch(coef(fit), error=function(e) NULL)
  if (!is.null(b)) return(b)
  if (!is.null(fit$coefficients)) return(fit$coefficients)
  stop("Cannot extract coefficients.")
}

predict_spml_sar <- function(fit, newdata, W, x_vars, ps_order=20) {
  rho <- get_spml_rho(fit)
  b <- get_coef_safe(fit)
  intercept <- if ("(Intercept)" %in% names(b)) as.numeric(b["(Intercept)"]) else 0
  for (v in x_vars) if (!v %in% names(newdata)) newdata[[v]] <- 0
  X <- as.matrix(newdata[, x_vars, drop=FALSE])
  xb <- as.numeric(intercept + X %*% b[x_vars])
  xb[!is.finite(xb)] <- 0
  spatial_filter_safe(xb, W, rho, ps_order=ps_order)
}

smearing_factor <- function(fit) {
  r <- tryCatch(as.numeric(residuals(fit)), error=function(e) NULL)
  if (is.null(r) || length(r)==0) return(1.0)
  mean(exp(r), na.rm=TRUE)
}

inv_log1p_smear <- function(mu_hat, smear) pmax(exp(mu_hat)*smear - 1, 0)

evaluate_m1_diag <- function(actual_loss, actual_occ,
                             P_hat, loss_hat_pos, EP_hat,
                             model_name, tau) {
  
  actual_loss <- as.numeric(actual_loss)
  actual_occ  <- as.numeric(actual_occ)
  P_hat       <- as.numeric(P_hat)
  loss_hat_pos <- as.numeric(loss_hat_pos)
  EP_hat      <- as.numeric(EP_hat)
  
  # --- 기본 지표 ---
  auc_val <- NA_real_
  if (length(unique(actual_occ)) >= 2) {
    auc_val <- as.numeric(
      pROC::auc(pROC::roc(actual_occ, P_hat, quiet = TRUE))
    )
  }
  
  rmse  <- sqrt(mean((actual_loss - EP_hat)^2, na.rm=TRUE))
  mae   <- mean(abs(actual_loss - EP_hat), na.rm=TRUE)
  medae <- median(abs(actual_loss - EP_hat), na.rm=TRUE)
  
  calib <- sum(EP_hat, na.rm=TRUE) / sum(actual_loss, na.rm=TRUE)
  
  # --- Spearman 분해 ---
  # 1) 전체 EP
  sp_ep_all <- suppressWarnings(
    cor(actual_loss, EP_hat, method="spearman", use="complete.obs")
  )
  
  # 2) 양성 구간 EP
  idx_pos <- which(actual_loss > 0)
  sp_ep_pos <- NA_real_
  sp_loss_pos <- NA_real_
  if (length(idx_pos) >= 5) {
    sp_ep_pos <- suppressWarnings(
      cor(actual_loss[idx_pos],
          EP_hat[idx_pos],
          method="spearman",
          use="complete.obs")
    )
    sp_loss_pos <- suppressWarnings(
      cor(actual_loss[idx_pos],
          loss_hat_pos[idx_pos],
          method="spearman",
          use="complete.obs")
    )
  }
  
  # 3) 발생확률 순위
  sp_p_only <- suppressWarnings(
    cor(actual_occ, P_hat, method="spearman", use="complete.obs")
  )
  
  tibble(
    tau = tau,
    model = model_name,
    AUC_occ = auc_val,
    RMSE = rmse,
    MAE = mae,
    MedAE = medae,
    Calib = calib,
    Spearman_EP_all = sp_ep_all,
    Spearman_EP_pos = sp_ep_pos,
    Spearman_Loss_pos = sp_loss_pos,
    Spearman_P_only = sp_p_only
  )
}

# ---------- main ----------
msg(">>> Load data & W")
df_train_raw <- readRDS(PATH_TRAIN)
df_test_raw  <- readRDS(PATH_TEST)
df_all_raw <- bind_rows(df_train_raw, df_test_raw)

W_obj <- readRDS(PATH_W)
region_ids <- as.character(W_obj$sigungu_cd)
years_all  <- sort(unique(as.integer(df_all_raw$year)))

W_pack <- build_W(W_obj$W, region_ids)
W_listw <- W_pack$listw
W_rs <- W_pack$W

panel <- make_balanced_template(df_all_raw, region_ids, years_all) %>%
  mutate(
    fire_cnt_year = as.numeric(fire_cnt_year),
    damage_capped = as.numeric(damage_capped),
    fire_cnt_year = tidyr::replace_na(fire_cnt_year, 0),
    damage_capped = tidyr::replace_na(damage_capped, 0),
    fire_any = if_else(damage_capped>0, 1, 0)
  )

train_mask <- panel$year < TEST_YEAR
panel <- impute_with_train_stats(panel, train_mask)

train_df <- panel %>% filter(year < TEST_YEAR) %>% arrange(sigungu_cd, year)
test_df  <- panel %>% filter(year == TEST_YEAR) %>% arrange(sigungu_cd, year)
stopifnot(nrow(test_df) == length(region_ids))

# Occurrence model fixed (once)
msg(">>> Fit Occurrence SAR-LPM (fixed)")
f_occ <- as.formula(paste("fire_any ~", paste(COVARS_LAG1, collapse=" + ")))
m_occ_sar <- splm::spml(
  f_occ, data=train_df,
  index=c("sigungu_cd","year"), listw=W_listw,
  model="random", lag=TRUE
)

occ_score <- predict_spml_sar(m_occ_sar, test_df, W_rs, COVARS_LAG1, ps_order=PS_FALLBACK_ORDER)
P_hat <- pmin(pmax(occ_score, 0), 1)

actual_loss <- as.numeric(test_df$damage_capped)
actual_occ  <- as.numeric(test_df$fire_any)

# tau grid on loss decomposition
res_list <- lapply(TAU_GRID, function(tau){
  
  msg("... tau=%d", tau)
  
  # define targets on FULL panel (to keep W conformable)
  train_df2 <- train_df %>%
    mutate(
      y_small = log1p(pmin(damage_capped, tau)),
      ex_any  = as.integer(damage_capped > tau),
      y_excess = log1p(pmax(damage_capped - tau, 0))
    )
  
  test_df2 <- test_df %>%
    mutate(
      ex_any = as.integer(damage_capped > tau)
    )
  
  # small SAR
  f_small <- as.formula(paste("y_small ~", paste(COVARS_LAG1, collapse=" + ")))
  m_small <- splm::spml(
    f_small, data=train_df2,
    index=c("sigungu_cd","year"), listw=W_listw,
    model="random", lag=TRUE
  )
  
  # excess-any SAR-LPM (rare: still full panel target)
  f_ex_occ <- as.formula(paste("ex_any ~", paste(COVARS_LAG1, collapse=" + ")))
  m_ex_occ <- splm::spml(
    f_ex_occ, data=train_df2,
    index=c("sigungu_cd","year"), listw=W_listw,
    model="random", lag=TRUE
  )
  
  # excess-size SAR on y_excess (full panel, zeros included)
  f_ex_size <- as.formula(paste("y_excess ~", paste(COVARS_LAG1, collapse=" + ")))
  m_ex_size <- splm::spml(
    f_ex_size, data=train_df2,
    index=c("sigungu_cd","year"), listw=W_listw,
    model="random", lag=TRUE
  )
  
  # predict components
  mu_small <- predict_spml_sar(m_small, test_df2, W_rs, COVARS_LAG1, ps_order=PS_FALLBACK_ORDER)
  smear_small <- smearing_factor(m_small)
  small_hat <- pmin(inv_log1p_smear(mu_small, smear_small), tau)
  
  ex_score <- predict_spml_sar(m_ex_occ, test_df2, W_rs, COVARS_LAG1, ps_order=PS_FALLBACK_ORDER)
  p_ex_hat <- pmin(pmax(ex_score, 0), 1)
  
  mu_ex <- predict_spml_sar(m_ex_size, test_df2, W_rs, COVARS_LAG1, ps_order=PS_FALLBACK_ORDER)
  smear_ex <- smearing_factor(m_ex_size)
  ex_hat <- inv_log1p_smear(mu_ex, smear_ex)
  
  loss_hat_pos <- small_hat + p_ex_hat * ex_hat
  EP_hat <- P_hat * loss_hat_pos
  
  evaluate_m1_diag(
    actual_loss,
    actual_occ,
    P_hat,
    loss_hat_pos,
    EP_hat,
    model_name="M1_SAR_ADD_LOSS",
    tau=tau
  )
})

res_grid <- bind_rows(res_list)
print(
  res_grid %>%
    dplyr::select(tau,
           Spearman_EP_all,
           Spearman_EP_pos,
           Spearman_Loss_pos,
           Spearman_P_only,
           Calib) %>%
    arrange(desc(Spearman_EP_pos))
)

saveRDS(res_grid, OUT_RDS)
write.csv(res_grid, OUT_CSV, row.names=FALSE)

msg("Saved: %s", OUT_RDS)
msg("Saved: %s", OUT_CSV)