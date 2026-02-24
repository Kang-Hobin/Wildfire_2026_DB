## =========================================================
## Wildfire_10_EDA_and_Train_Test_Split.R
## 목적: Y 변수 무결성 확보 및 Out-of-Time 데이터 분리
## =========================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# 0. 환경 설정
if (!dir.exists("Wildfire_Data/figures")) dir.create("Wildfire_Data/figures", recursive = TRUE)

# 1. 데이터 로드 (9번에서 생성된 최신 패널)
df <- readRDS("Wildfire_Data/processed_data/df_annual_insurance_panel.rds")

# 2. Y 변수 생성 및 100ha 캡핑 (보험 수리적 전처리)
# ---------------------------------------------------------
message(">>> Y 변수 생성 및 수치 검증 중...")

df_final <- df %>%
  mutate(
    # [Y1] Hurdle 1: 발생 여부 (0 or 1)
    fire_any_year = as.integer(fire_any_year),
    fire_any_year_f = factor(fire_any_year, levels = c(0, 1), labels = c("No Fire", "Fire")),
    
    # [Y2] Hurdle 2: 피해 규모 (Severity) - 100ha 캡핑
    # 극단적인 대형 산불(Outlier)이 전체 회귀선을 왜곡하는 것을 방지
    damage_capped = pmin(damage_total_year, 100),
    log_damage_capped = log1p(damage_capped),
    
    # [Y3] Count: 발생 횟수 (Frequency)
    fire_cnt_year = as.integer(fire_cnt_year)
  ) %>%
  # Y 변수에 결측치가 있으면 모델이 터지므로 최종 확인
  filter(!is.na(fire_any_year), !is.na(damage_total_year))

# 3. 데이터 분포 시각화 (Hurdle Logic 검증)
# ---------------------------------------------------------
# (시각화 코드는 기존과 동일하므로 중략하지만, df_final을 참조하도록 유지)
p1 <- ggplot(df_final, aes(x = fire_any_year_f, fill = fire_any_year_f)) +
  geom_bar(alpha = 0.8) + scale_fill_manual(values = c("gray70", "firebrick")) +
  labs(title = "(A) Occurrence (Y1)") + theme_minimal()

p4 <- ggplot(df_final %>% filter(damage_capped > 0), aes(x = damage_capped)) +
  geom_histogram(fill = "darkorange", bins = 50) +
  labs(title = "(D) Capped Severity (Y2)") + theme_minimal()

combined_plot <- (p1 + p4) # 필요한 그래프 위주로 병합
ggsave("Wildfire_Data/figures/EDA_Hurdle_Logic.png", combined_plot, width = 10, height = 5)

# 4. [핵심 수정] Train-Test Split (2024년 기준 고정)
# ---------------------------------------------------------
# 2017-2024: 모델 학습 (충분한 시계열 확보)
# 2025: 최종 미래 예측력 테스트 (가장 최신 데이터)
train_df <- df_final %>% filter(year <= 2024)
test_df  <- df_final %>% filter(year >= 2025)

# Y 데이터 누락 여부 최종 체크
message(">>> [검증] Train Set Y 변수 확인: ", 
        all(c("fire_any_year", "damage_capped", "fire_cnt_year") %in% colnames(train_df)))

# 이질적 분포에 대한 확인
tau <- 10 # 산불 대응 단계 확산 대응 초기 임계점

small_log_cap100 <- log1p(df_final$damage_capped[df_final$damage_capped < tau])
large_log_cap100 <- log1p(df_final$damage_capped[df_final$damage_capped >= tau])
# truncated cases
capped_log_cap100 <- log1p(df_final$damage_total_year[df_final$damage_capped >= 100])


par(mfrow=c(1,2))## =========================================================
## Wildfire_10_EDA_and_Train_Test_Split.R
## 목적: Y 변수 무결성 확보 및 Out-of-Time 데이터 분리
## =========================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# 0. 환경 설정
if (!dir.exists("Wildfire_Data/figures")) dir.create("Wildfire_Data/figures", recursive = TRUE)

# 1. 데이터 로드 (9번에서 생성된 최신 패널)
df <- readRDS("Wildfire_Data/processed_data/df_annual_insurance_panel.rds")

# 2. Y 변수 생성 및 100ha 캡핑 (보험 수리적 전처리)
# ---------------------------------------------------------
message(">>> Y 변수 생성 및 수치 검증 중...")

df_final <- df %>%
  mutate(
    # [Y1] Hurdle 1: 발생 여부 (0 or 1)
    fire_any_year = as.integer(fire_any_year),
    fire_any_year_f = factor(fire_any_year, levels = c(0, 1), labels = c("No Fire", "Fire")),
    
    # [Y2] Hurdle 2: 피해 규모 (Severity) - 100ha 캡핑
    # 극단적인 대형 산불(Outlier)이 전체 회귀선을 왜곡하는 것을 방지
    damage_capped = pmin(damage_total_year, 100),
    log_damage_capped = log1p(damage_capped),
    
    # [Y3] Count: 발생 횟수 (Frequency)
    fire_cnt_year = as.integer(fire_cnt_year)
  ) %>%
  # Y 변수에 결측치가 있으면 모델이 터지므로 최종 확인
  filter(!is.na(fire_any_year), !is.na(damage_total_year))

# 3. 데이터 분포 시각화 (Hurdle Logic 검증)
# ---------------------------------------------------------
# (시각화 코드는 기존과 동일하므로 중략하지만, df_final을 참조하도록 유지)
p1 <- ggplot(df_final, aes(x = fire_any_year_f, fill = fire_any_year_f)) +
  geom_bar(alpha = 0.8) + scale_fill_manual(values = c("gray70", "firebrick")) +
  labs(title = "(A) Occurrence (Y1)") + theme_minimal()

p4 <- ggplot(df_final %>% filter(damage_capped > 0), aes(x = damage_capped)) +
  geom_histogram(fill = "darkorange", bins = 50) +
  labs(title = "(D) Capped Severity (Y2)") + theme_minimal()

combined_plot <- (p1 + p4) # 필요한 그래프 위주로 병합
ggsave("Wildfire_Data/figures/EDA_Hurdle_Logic.png", combined_plot, width = 10, height = 5)

# 4. [핵심 수정] Train-Test Split (2024년 기준 고정)
# ---------------------------------------------------------
# 2017-2024: 모델 학습 (충분한 시계열 확보)
# 2025: 최종 미래 예측력 테스트 (가장 최신 데이터)
train_df <- df_final %>% filter(year <= 2024)
test_df  <- df_final %>% filter(year >= 2025)

# Y 데이터 누락 여부 최종 체크
message(">>> [검증] Train Set Y 변수 확인: ", 
        all(c("fire_any_year", "damage_capped", "fire_cnt_year") %in% colnames(train_df)))

# 이질적 분포에 대한 확인
tau <- 10 # 산불 대응 단계 확산 대응 초기 임계점

small_log_cap100 <- log1p(df_final$damage_capped[df_final$damage_capped < tau])
large_log_cap100 <- log1p(df_final$damage_capped[df_final$damage_capped >= tau])
# truncated cases
capped_log_cap100 <- log1p(df_final$damage_total_year[df_final$damage_capped >= 100])


par(mfrow=c(1,2))
hist(small_log_cap100, breaks=30, main="log payout of Small fires")
hist(large_log_cap100, breaks=30, main="log payout of 10ha+ fires")
par(mfrow=c(1,1))

# 5. 결과 저장
# ---------------------------------------------------------
saveRDS(train_df, "Wildfire_Data/processed_data/df_train_final.rds")
saveRDS(test_df, "Wildfire_Data/processed_data/df_test_final.rds")

cat("\n--- Split Summary ---\n")
cat("Train Set (2017-2024):", nrow(train_df), "obs\n")
cat("Test Set  (2025):", nrow(test_df), "obs\n")
cat("Avg Damage in Train:", mean(train_df$damage_capped, na.rm=TRUE), "ha\n")
cat("Avg Damage in Test:", mean(test_df$damage_capped, na.rm=TRUE), "ha\n")

hist(small_log_cap100, breaks=30, main="log Small")
hist(large_log_cap100, breaks=30, main="log Large")

# 5. 결과 저장
# ---------------------------------------------------------
saveRDS(train_df, "Wildfire_Data/processed_data/df_train_final.rds")
saveRDS(test_df, "Wildfire_Data/processed_data/df_test_final.rds")

cat("\n--- Split Summary ---\n")
cat("Train Set (2017-2024):", nrow(train_df), "obs\n")
cat("Test Set  (2025):", nrow(test_df), "obs\n")
cat("Avg Damage in Train:", mean(train_df$damage_capped, na.rm=TRUE), "ha\n")
cat("Avg Damage in Test:", mean(test_df$damage_capped, na.rm=TRUE), "ha\n")
