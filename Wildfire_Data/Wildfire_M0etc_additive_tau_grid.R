## =========================================================
## Wildfire_M0_additive_tau_grid.R
## 목적:
##  - M0 additive(severity: small+excess) 구조에서
##    tau grid를 돌려 breakpoint 후보를 탐색
##
## NOTE:
##  - Occurrence는 기존처럼 NB(prob_nb) 고정
##  - Severity는 fire_any_year==1에서만 학습
##  - small: min(damage, tau)
##  - excess: max(damage-tau, 0) with 2-part (logit + size)
## =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(MASS)
  library(quantreg)
})

# ---------------------------------------------------------
# Settings
# ---------------------------------------------------------
PATH_TRAIN <- "Wildfire_Data/processed_data/df_train_final.rds"
PATH_TEST  <- "Wildfire_Data/processed_data/df_test_final.rds"

OUT_DIR <- "Wildfire_Data/results"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

OUT_RDS <- file.path(OUT_DIR, "m0_additive_tau_grid_perf.rds")
OUT_CSV <- file.path(OUT_DIR, "m0_additive_tau_grid_perf.csv")

# grid (추천: 5~50에서 5단위, 필요시 조정)
TAU_GRID <- seq(5, 50, by = 5)

# excess-positive 최소 샘플 수 (너무 적으면 intercept-only)
MIN_EXCESS_POS <- 20

# covariates
x_var_sev <- c("spring_ws_p95_lag1", "max_dry_run_spring_lag1", "monsoon_prcp_sum_lag1",
               "ta_mtw_mean_lag1", "hm_mtw_min_lag1")

msg <- function(...) message(sprintf(...))

# ---------------------------------------------------------
# helpers (English code as requested)
# ---------------------------------------------------------
fit_additive_for_tau <- function(train_s, test_s, tau, min_excess_pos = 20) {
  
  # 1) Occurrence: NB
  model_cnt_nb <- MASS::glm.nb(
    fire_cnt_year ~ ta_sfc_mean_lag1 + spring_he_min_lag1 + spring_ws_p95_lag1 +
      spring_dry_days_lag1 + monsoon_prcp_sum_lag1 + sigungu_cd,
    data = train_s
  )
  
  test2 <- test_s %>%
    mutate(
      mu_nb   = predict(model_cnt_nb, newdata = ., type = "response"),
      prob_nb = 1 - dnbinom(0, mu = mu_nb, size = model_cnt_nb$theta)
    )
  
  # AUC on occurrence
  auc_nb <- NA_real_
  if (length(unique(test2$fire_any_year)) >= 2) {
    auc_nb <- as.numeric(pROC::auc(pROC::roc(test2$fire_any_year, test2$prob_nb, quiet = TRUE)))
  }
  
  # 2) Severity training only on fire==1
  train_pos <- train_s %>%
    filter(fire_any_year == 1) %>%
    mutate(
      damage_capped = as.numeric(damage_capped),
      y_small = pmin(damage_capped, tau),
      y_excess = pmax(damage_capped - tau, 0),
      ex_any = as.integer(y_excess > 0),
      log1p_small  = log1p(y_small),
      log1p_excess = log1p(y_excess)
    )
  
  test_pos <- test2 %>%
    mutate(
      damage_capped = as.numeric(damage_capped),
      y_small = pmin(damage_capped, tau),
      y_excess = pmax(damage_capped - tau, 0),
      ex_any = as.integer(y_excess > 0)
    )
  
  # small models
  f_small_logn <- as.formula(paste("log1p_small ~", paste(x_var_sev, collapse = " + ")))
  m_small_logn <- lm(f_small_logn, data = train_pos)
  smear_small  <- mean(exp(residuals(m_small_logn)), na.rm = TRUE)
  
  f_small_qr <- as.formula(paste("y_small ~", paste(x_var_sev, collapse = " + ")))
  m_small_qr <- quantreg::rq(f_small_qr, tau = 0.5, data = train_pos)
  
  pred_small_mu <- as.numeric(predict(m_small_logn, newdata = test_pos))
  small_hat_ols <- pmax(exp(pred_small_mu) * smear_small - 1, 0)
  small_hat_qr  <- pmax(as.numeric(predict(m_small_qr, newdata = test_pos)), 0)
  
  small_hat_ols <- pmin(small_hat_ols, tau)
  small_hat_qr  <- pmin(small_hat_qr,  tau)
  
  # excess 2-part
  f_ex_occ <- as.formula(paste("ex_any ~", paste(x_var_sev, collapse = " + ")))
  m_ex_occ <- glm(f_ex_occ, data = train_pos, family = binomial(link = "logit"))
  
  train_ex_pos <- train_pos %>% filter(ex_any == 1)
  
  if (nrow(train_ex_pos) < min_excess_pos) {
    # intercept-only fallback
    m_ex_size_logn <- lm(log1p_excess ~ 1, data = train_ex_pos)
    m_ex_size_qr   <- quantreg::rq(log1p_excess ~ 1, tau = 0.5, data = train_ex_pos)
  } else {
    f_ex_size_logn <- as.formula(paste("log1p_excess ~", paste(x_var_sev, collapse = " + ")))
    m_ex_size_logn <- lm(f_ex_size_logn, data = train_ex_pos)
    
    f_ex_size_qr <- as.formula(paste("log1p_excess ~", paste(x_var_sev, collapse = " + ")))
    m_ex_size_qr <- quantreg::rq(f_ex_size_qr, tau = 0.5, data = train_ex_pos)
  }
  
  smear_excess <- mean(exp(residuals(m_ex_size_logn)), na.rm = TRUE)
  
  p_ex_hat <- as.numeric(predict(m_ex_occ, newdata = test_pos, type = "response"))
  
  mu_ex_hat <- as.numeric(predict(m_ex_size_logn, newdata = test_pos))
  excess_hat_ols <- pmax(exp(mu_ex_hat) * smear_excess - 1, 0)
  
  mu_ex_hat_qr <- as.numeric(predict(m_ex_size_qr, newdata = test_pos))
  excess_hat_qr <- pmax(exp(mu_ex_hat_qr) - 1, 0)  # QR on log1p scale
  
  sev_add_ols <- small_hat_ols + p_ex_hat * excess_hat_ols
  sev_add_qr  <- small_hat_qr  + p_ex_hat * excess_hat_qr
  
  # EP
  ep_ols <- test2$prob_nb * sev_add_ols
  ep_qr  <- test2$prob_nb * sev_add_qr
  
  # metrics
  actual <- as.numeric(test2$damage_capped)
  
  metric_row <- function(ep, name) {
    tibble(
      tau = tau,
      model_type = name,
      AUC_occ = auc_nb,
      RMSE = sqrt(mean((actual - ep)^2, na.rm = TRUE)),
      MAE  = mean(abs(actual - ep), na.rm = TRUE),
      MedAE = median(abs(actual - ep), na.rm = TRUE),
      Spearman = suppressWarnings(cor(actual, ep, method = "spearman", use = "complete.obs")),
      Calib_Ratio = sum(ep, na.rm = TRUE) / sum(actual, na.rm = TRUE),
      ex_pos_n = nrow(train_ex_pos),
      ex_pos_fallback = nrow(train_ex_pos) < min_excess_pos
    )
  }
  
  bind_rows(
    metric_row(ep_ols, "EP_NB_ADD_OLS"),
    metric_row(ep_qr,  "EP_NB_ADD_QR")
  )
}

# ---------------------------------------------------------
# main
# ---------------------------------------------------------
msg(">>> Loading train/test...")
train_s <- readRDS(PATH_TRAIN)
test_s  <- readRDS(PATH_TEST)

train_s <- train_s %>% mutate(sigungu_cd = factor(sigungu_cd))
test_s  <- test_s  %>% mutate(sigungu_cd = factor(sigungu_cd, levels = levels(train_s$sigungu_cd)))

missing_cov <- setdiff(x_var_sev, names(train_s))
if (length(missing_cov) > 0) stop("Missing covariates in train: ", paste(missing_cov, collapse = ", "))

msg(">>> Running tau grid: %s", paste(TAU_GRID, collapse = ", "))

res_list <- lapply(TAU_GRID, function(tau) {
  msg("... tau = %d", tau)
  fit_additive_for_tau(train_s, test_s, tau, min_excess_pos = MIN_EXCESS_POS)
})

res_grid <- bind_rows(res_list)

print(res_grid %>% arrange(model_type, RMSE))

saveRDS(res_grid, OUT_RDS)
write.csv(res_grid, OUT_CSV, row.names = FALSE)

msg("Saved: %s", OUT_RDS)
msg("Saved: %s", OUT_CSV)

# quick topline (by Spearman within model_type)
topline <- res_grid %>%
  group_by(model_type) %>%
  slice_max(order_by = Spearman, n = 1, with_ties = FALSE) %>%
  ungroup()

msg(">>> Topline by Spearman:")
print(topline)