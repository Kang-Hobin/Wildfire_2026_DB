## =========================================================
## Wildfire_M0.R
## 분석 특성: 시공간 상관 미반영 전통 모형 (Benchmark)
## 목적: Hurdle(Two-part) 모형 성능 산출 (M1과 동일 지표 적용)
## =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(car)
  library(pROC)
  library(mgcv)
  library(MASS)
  library(quantreg)
})

# ---------------------------------------------------------
# 1. 데이터 로드 및 전처리
# ---------------------------------------------------------
message(">>> [1/5] 데이터 로드 및 Factor 설정 중...")
train_s <- readr::read_rds("Wildfire_Data/processed_data/df_train_final.rds")
test_s  <- readr::read_rds("Wildfire_Data/processed_data/df_test_final.rds")

# 시군구 코드 동기화
train_s <- train_s %>% mutate(sigungu_cd = factor(sigungu_cd))
test_s  <- test_s  %>% mutate(sigungu_cd = factor(sigungu_cd, levels = levels(train_s$sigungu_cd)))

# ---------------------------------------------------------
# 2. [Step 1] Occurrence Model (발생 판별)
# ---------------------------------------------------------
message(">>> [2/5] Step 1: 발생(Occurrence) 모델 구축 중...")

x_var_occ <- c("spring_dry_days_lag1", "max_dry_run_spring_lag1", "spring_ws_p95_lag1", 
               "spring_he_min_lag1", "monsoon_prcp_sum_lag1", "winter_dry_days_lag1", "ta_sfc_mean_lag1")

# (A) 로지스틱 회귀
model_occ_logit <- glm(as.formula(paste("fire_any_year ~", paste(c(x_var_occ, "sigungu_cd"), collapse = " + "))), 
                       data = train_s, family = binomial(link = "logit"))

# (B) 일반화 가산 모형 (GAM, 시군구 랜덤효과 처리)
model_occ_gam <- mgcv::gam(
  fire_any_year ~ s(ta_sfc_mean_lag1, k = 5) + s(spring_he_min_lag1, k = 5) + 
    s(prcp_sfc_sum_lag1, k = 5) + spring_ws_p95_lag1 + 
    spring_dry_days_lag1 + monsoon_prcp_sum_lag1 + 
    s(sigungu_cd, bs = "re"),
  data = train_s, family = binomial(link = "logit"), method = "REML"
)

# (C) 음이항 회귀 (Count 기반 발생 확률)
model_cnt_nb <- MASS::glm.nb(fire_cnt_year ~ ta_sfc_mean_lag1 + spring_he_min_lag1 + spring_ws_p95_lag1 + 
                               spring_dry_days_lag1 + monsoon_prcp_sum_lag1 + sigungu_cd, data = train_s)

# 테스트 데이터 예측
test_s <- test_s %>%
  mutate(
    prob_logit = predict(model_occ_logit, newdata = ., type = "response"),
    # NB 확률: P(N >= 1)
    mu_nb    = predict(model_cnt_nb, newdata = ., type = "response"),
    prob_gam   = predict(model_occ_gam, newdata = ., type = "response"),
    prob_nb  = 1 - dnbinom(0, mu = mu_nb, size = model_cnt_nb$theta)
  )

# Step 1 평가: AUC
auc_res <- list(
  Logit = as.numeric(pROC::auc(pROC::roc(test_s$fire_any_year, test_s$prob_logit, quiet=TRUE))),
  GAM   = as.numeric(pROC::auc(pROC::roc(test_s$fire_any_year, test_s$prob_gam, quiet=TRUE))),
  NB    = as.numeric(pROC::auc(pROC::roc(test_s$fire_any_year, test_s$prob_nb, quiet=TRUE)))
)

# ---------------------------------------------------------
# 3. [Step 2] Severity Model (피해 심도)
# ---------------------------------------------------------
message(">>> [3/5] Step 2: 심도(Severity) 모델 구축 중...")

train_sev_s <- train_s %>% filter(fire_any_year == 1)
x_var_sev   <- c("spring_ws_p95_lag1", "max_dry_run_spring_lag1", "monsoon_prcp_sum_lag1", "ta_mtw_mean_lag1", "hm_mtw_min_lag1")

# (M1) Lognormal (OLS)
model_sev_logn <- lm(as.formula(paste("log_damage_capped ~", paste(x_var_sev, collapse = " + "))), data = train_sev_s)
smear_factor   <- mean(exp(residuals(model_sev_logn)), na.rm = TRUE)

# (M2) Quantile Regression (Median)
model_sev_qr50 <- quantreg::rq(as.formula(paste("damage_capped ~", paste(x_var_sev, collapse = " + "))), 
                               tau = 0.5, data = train_sev_s)

# 심도 예측 (전체 테스트셋 대상 - EP 산출용)
test_s <- test_s %>%
  mutate(
    pred_log_m1 = predict(model_sev_logn, newdata = .),
    sev_hat_m1  = exp(pred_log_m1) * smear_factor - 1,
    sev_hat_m2  = pmax(as.numeric(predict(model_sev_qr50, newdata = .)), 0)
  )

# ---------------------------------------------------------
# 4. 최종 통합 성능 평가: Expected Payout (EP)
# ---------------------------------------------------------
message(">>> [4/5] 최종 통합 성능 지표 산출 중...")

test_s <- test_s %>%
  mutate(
    EP_NB_OLS = prob_nb * sev_hat_m1,
    EP_NB_QR  = prob_nb * sev_hat_m2
  )

# 통합 5대 지표 산출 (M1과 동일 규격)
m0_final_perf <- test_s %>%
  pivot_longer(cols = c(EP_NB_OLS, EP_NB_QR), names_to = "model_type", values_to = "EP_val") %>%
  group_by(model_type) %>%
  summarise(
    AUC_occ     = auc_res$NB, # NB 발생 확률 기준 AUC
    RMSE        = sqrt(mean((damage_capped - EP_val)^2, na.rm = TRUE)),
    MAE         = mean(abs(damage_capped - EP_val), na.rm = TRUE),
    MedAE       = median(abs(damage_capped - EP_val), na.rm = TRUE),
    Spearman    = cor(damage_capped, EP_val, method = "spearman"),
    Calib_Ratio = sum(EP_val, na.rm = TRUE) / sum(damage_capped, na.rm = TRUE)
  )

# ---------------------------------------------------------
# 5. 결과 출력
# ---------------------------------------------------------
message("\n--- [Result] M0 Hurdle Model Performance ---")
print(m0_final_perf)
print(auc_res) # M0 occ_AUC

# M1과의 비교를 위해 객체 저장
saveRDS(m0_final_perf, "Wildfire_Data/results/m0_final_perf.rds")

