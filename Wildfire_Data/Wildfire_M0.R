## =========================================================
## Wildfire_M0.R
## 분석 특성: 시공간 상관 미반영 전통 모형의 성능 확인
## 목적: Hurdle(Two-part) 모형을 통한 기초 기대 손실(Expected Loss) 산출
## =========================================================

library(tidyverse)
library(car)
library(pROC)
library(mgcv)
library(MASS)
library(quantreg)

# ---------------------------------------------------------
# 1. 데이터 로드 및 전처리 (R10 전처리 데이터 활용)
# ---------------------------------------------------------
message(">>> 데이터 로드 및 Factor 설정 중...")
train_s <- readr::read_rds("Wildfire_Data/processed_data/df_train_final.rds")
test_s  <- readr::read_rds("Wildfire_Data/processed_data/df_test_final.rds")

# 시군구 코드의 일관성을 위한 Factor 레벨 동기화
train_s <- train_s %>% mutate(sigungu_cd = factor(sigungu_cd))
test_s  <- test_s  %>% mutate(sigungu_cd = factor(sigungu_cd, levels = levels(train_s$sigungu_cd)))


# ---------------------------------------------------------
# 2. [Step 1] Occurrence Model: 산불 발생 여부 및 확률 예측
# ---------------------------------------------------------
message(">>> Step 1: 발생(Occurrence) 모델 구축 중...")

# 2.1 다중공선성(VIF) 검동 및 변수 정제
# 초기 후보 변수군 정의
x_var_occ_init <- c(
  "spring_dry_days_lag1", "spring_ws_p95_lag1", "spring_he_min_lag1",
  "max_dry_run_spring_lag1", "monsoon_prcp_sum_lag1", 
  "winter_dry_days_lag1", "ta_sfc_mean_lag1", "prcp_sfc_sum_lag1"
)

# VIF 체크를 통한 변수 선별 (prcp_sfc_sum_lag1 제거 결정)
x_var_occ_final <- c(
  "spring_dry_days_lag1", "max_dry_run_spring_lag1", "spring_ws_p95_lag1", 
  "spring_he_min_lag1", "monsoon_prcp_sum_lag1", "winter_dry_days_lag1", "ta_sfc_mean_lag1"
)

# 2.2 모델 적합 (Logit, GAM, NB)
# (A) 로지스틱 회귀 (고정효과 포함)
formula_occ <- as.formula(paste("fire_any_year ~", paste(c(x_var_occ_final, "sigungu_cd"), collapse = " + ")))
model_occ_refined <- glm(formula_occ, data = train_s, family = binomial(link = "logit"))

# (B) 일반화 가산 모형 (GAM, 시군구 랜덤효과 처리)
model_occ_gam <- mgcv::gam(
  fire_any_year ~ s(ta_sfc_mean_lag1, k = 5) + s(spring_he_min_lag1, k = 5) + 
    s(prcp_sfc_sum_lag1, k = 5) + spring_ws_p95_lag1 + 
    spring_dry_days_lag1 + monsoon_prcp_sum_lag1 + 
    s(sigungu_cd, bs = "re"),
  data = train_s, family = binomial(link = "logit"), method = "REML"
)

# (C) 음이항 회귀 (NB, 발생 횟수 기반 확률 도출용)
model_cnt_nb <- MASS::glm.nb(
  fire_cnt_year ~ ta_sfc_mean_lag1 + spring_he_min_lag1 + prcp_sfc_sum_lag1 +
    spring_ws_p95_lag1 + spring_dry_days_lag1 + monsoon_prcp_sum_lag1 +
    sigungu_cd,
  data = train_s
)

# 2.3 발생 모형 성능 평가 (AUC & Brier Score)
test_s <- test_s %>%
  mutate(
    prob_logit = predict(model_occ_refined, newdata = ., type = "response"),
    prob_gam   = predict(model_occ_gam, newdata = ., type = "response")
  )

# NB 발생 확률 변환: P(N >= 1) = 1 - P(N = 0)
mu_nb <- predict(model_cnt_nb, newdata = test_s, type = "response")
theta_nb <- model_cnt_nb$theta
test_s$prob_nb <- 1 - dnbinom(0, mu = mu_nb, size = theta_nb)

# 성능 지표 산출
brier_score <- function(y, p) mean((y - p)^2, na.rm = TRUE)
auc_res <- list(
  Logit = pROC::auc(pROC::roc(test_s$fire_any_year, test_s$prob_logit, quiet=TRUE)),
  GAM   = pROC::auc(pROC::roc(test_s$fire_any_year, test_s$prob_gam, quiet=TRUE)),
  NB    = pROC::auc(pROC::roc(test_s$fire_any_year, test_s$prob_nb, quiet=TRUE))
)
brier_res <- list(
  Logit = brier_score(test_s$fire_any_year, test_s$prob_logit),
  GAM   = brier_score(test_s$fire_any_year, test_s$prob_gam),
  NB    = brier_score(test_s$fire_any_year, test_s$prob_nb)
)


# ---------------------------------------------------------
# 3. [Step 2] Severity Model: 산불 피해 심도 예측 (100ha Cap 적용)
# ---------------------------------------------------------
message(">>> Step 2: 심도(Severity) 모델 구축 중...")

# 3.1 심도 모델용 데이터셋 분리 (발생 건수만 추출)
train_sev_s <- train_s %>% filter(fire_any_year == 1)
test_sev_s  <- test_s  %>% filter(fire_any_year == 1)

# 심도용 변수 정의 (확산 중심)
x_var_sev <- c("spring_ws_p95_lag1", "max_dry_run_spring_lag1", 
               "monsoon_prcp_sum_lag1", "ta_mtw_mean_lag1", "hm_mtw_min_lag1")

# 3.2 모델 적합 (Lognormal & Quantile Regression)
formula_sev <- as.formula(paste("log_damage_capped ~", paste(x_var_sev, collapse = " + ")))

# (M1) Lognormal 모델 (OLS)
model_sev_logn <- lm(formula_sev, data = train_sev_s)
smear_factor   <- mean(exp(residuals(model_sev_logn)), na.rm = TRUE)

# (M2) 분위수 회귀 (Quantile Regression, 중앙값)
model_sev_qr50 <- quantreg::rq(update(formula_sev, damage_capped ~ .), tau = 0.5, data = train_sev_s)

# 3.3 심도 예측 및 성능 비교
test_s <- test_sev_s %>%
  mutate(
    # M1 예측: 로그 역변환 및 스미어링 보정 적용
    pred_log_m1 = predict(model_sev_logn, newdata = .),
    sev_hat_m1  = (exp(pred_log_m1) * smear_factor) - 1,
    
    # M2 예측: 원본 스케일 중앙값 예측
    sev_hat_m2  = pmax(as.numeric(predict(model_sev_qr50, newdata = .)), 0)
  )

perf_comp <- test_s %>%
  pivot_longer(cols = c(sev_hat_m1, sev_hat_m2), names_to = "model", values_to = "pred") %>%
  group_by(model) %>%
  summarise(
    RMSE = sqrt(mean((damage_capped - pred)^2, na.rm = TRUE)),
    MAE  = mean(abs(damage_capped - pred), na.rm = TRUE),
    MedAE = median(abs(damage_capped - pred), na.rm = TRUE),
    Spearman = cor(damage_capped, pred, method = "spearman")
  )


# ---------------------------------------------------------
# 4. 최종 통합 성능 평가: Expected Payout (EP) 산출 및 비교
# ---------------------------------------------------------
message(">>> 최종 통합 성능 지표 산출 중...")

test_s <- test_s %>%
  mutate(
    # [M0-OLS 조합] NB(발생) * Lognormal(심도)
    EP_OLS = prob_nb * sev_hat_m1,
    
    # [M0-QR 조합] NB(발생) * Quantreg(심도)
    EP_QR  = prob_nb * sev_hat_m2
  )

m0_final_perf <- test_s %>%
  pivot_longer(cols = c(EP_OLS, EP_QR), names_to = "model_type", values_to = "EP_val") %>%
  group_by(model_type) %>%
  summarise(
    RMSE = sqrt(mean((damage_capped - EP_val)^2, na.rm = TRUE)),
    MAE  = mean(abs(damage_capped - EP_val), na.rm = TRUE),
    MedAE = median(abs(damage_capped - EP_val), na.rm = TRUE),
    Spearman = cor(damage_capped, EP_val, method = "spearman"),
    Calib_Ratio = sum(EP_val, na.rm = TRUE) / sum(damage_capped, na.rm = TRUE)
  )

# ---------------------------------------------------------
# 5. 결과 출력
# ---------------------------------------------------------
cat("\n--- Step 1: Occurrence Performance ---\n")
print(auc_res)
print(brier_res)

cat("\n--- Step 2: Severity Performance ---\n")
print(perf_comp)

cat("\n--- Final M0 Hurdle Model Performance ---\n")
print(m0_final_perf)