## =========================================================
## Wildfire_res_summary.R  (Fixed: explicit 3-file version)
##
## 목적:
##  - 결과 RDS 3개(M0, M1, M2)를 명시적으로 로드 (자동 탐색 X)
##  - 공통 컬럼 정규화 후 통합 비교표 생성
##  - M2(STAR)는 RAW/CAL 비교표(wide)까지 생성
##  - family별 topline 1개(의도 반영) 생성
##
## 입력(수정해서 사용):
##  - PATH_M0: benchmark 결과 (예: m0_final_perf.rds)
##  - PATH_M1: SAR 2-part 결과 (예: m1_sar_2part_perf_2025.rds)
##  - PATH_M2: STAR 5-model 결과 (예: m2_star_5models_perf_final.rds)
##
## 출력:
##  - Wildfire_Data/results/res_summary_all.csv
##  - Wildfire_Data/results/res_summary_m2_raw_cal.csv (가능할 때)
##  - Wildfire_Data/results/res_summary_topline.csv
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

RESULT_DIR <- "Wildfire_Data/results"

# ---------------------------------------------------------
# 0) USER: set your 3 result files here
# ---------------------------------------------------------
PATH_M0 <- file.path(RESULT_DIR, "m0_final_perf.rds")

# 아래는 예시. 너의 실제 저장 파일명으로 맞춰줘.
# (이전 대화 기준으로는 sar2part_vs_ols2part_perf_2025.rds가 있었는데,
#  지금은 SAR만 저장한다고 했으니 새 파일명으로 바꾸는 걸 추천)
PATH_M1 <- file.path(RESULT_DIR, "m1_sar_perf.rds")

PATH_M2 <- file.path(RESULT_DIR, "m2_star_5models_perf_final.rds")

OUT_CSV_ALL <- file.path(RESULT_DIR, "res_summary_all.csv")
OUT_CSV_M2  <- file.path(RESULT_DIR, "res_summary_m2_raw_cal.csv")
OUT_CSV_TOP <- file.path(RESULT_DIR, "res_summary_topline.csv")

msg <- function(...) message(sprintf(...))

# ---------------------------------------------------------
# helpers
# ---------------------------------------------------------
safe_read_tbl <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  obj <- readRDS(path)
  if (!inherits(obj, "data.frame")) stop("RDS is not a data.frame/tibble: ", path)
  as_tibble(obj)
}

normalize_cols <- function(df) {
  df <- as_tibble(df)
  
  # model name
  if (!"model" %in% names(df)) {
    if ("model_type" %in% names(df)) df <- df %>% rename(model = model_type)
  }
  
  # AUC column
  if (!"AUC_occ" %in% names(df)) {
    if ("AUC" %in% names(df)) df <- df %>% rename(AUC_occ = AUC)
    if ("auc" %in% names(df)) df <- df %>% rename(AUC_occ = auc)
  }
  
  # calibration column
  if (!"Calib" %in% names(df)) {
    if ("Calib_Ratio" %in% names(df)) df <- df %>% rename(Calib = Calib_Ratio)
  }
  
  keep <- intersect(
    c("model","AUC_occ","RMSE","MAE","MedAE","Spearman","Calib","version","calib_factor"),
    names(df)
  )
  
  df %>% dplyr::select(all_of(keep))
}

strip_suffix <- function(x) gsub("\\s*\\[(RAW|CAL)\\]\\s*$", "", x)

add_family <- function(df, family, source_file) {
  df %>%
    mutate(
      family = family,
      source_file = source_file
    )
}

infer_version <- function(df) {
  df %>%
    mutate(
      version = case_when(
        !is.na(version) ~ version,
        str_detect(model, "\\[CAL\\]$") ~ "CAL",
        str_detect(model, "\\[RAW\\]$") ~ "RAW",
        TRUE ~ NA_character_
      )
    )
}

# ---------------------------------------------------------
# 1) load 3 tables
# ---------------------------------------------------------
msg(">>> Loading 3 result tables...")

m0 <- safe_read_tbl(PATH_M0) %>% normalize_cols() %>% add_family("M0 (Benchmark)", basename(PATH_M0))
m1 <- safe_read_tbl(PATH_M1) %>% normalize_cols() %>% add_family("M1 (SAR 2-part)", basename(PATH_M1))
m2 <- safe_read_tbl(PATH_M2) %>% normalize_cols() %>% add_family("M2 (STAR)", basename(PATH_M2))

res_all <- bind_rows(m0, m1, m2) %>% infer_version()

msg(">>> Combined table:")
print(res_all)

write.csv(res_all, OUT_CSV_ALL, row.names = FALSE)
msg("Saved: %s", OUT_CSV_ALL)

# ---------------------------------------------------------
# 2) M2 RAW vs CAL wide summary (if present)
# ---------------------------------------------------------
m2_only <- res_all %>% filter(family == "M2 (STAR)")

if (nrow(m2_only) > 0 && any(m2_only$version %in% c("RAW","CAL"))) {
  # prefer making wide only when both exist
  if (all(c("RAW","CAL") %in% unique(na.omit(m2_only$version)))) {
    m2_wide <- m2_only %>%
      mutate(model_base = strip_suffix(model)) %>%
      dplyr::select(model_base, version, AUC_occ, RMSE, MAE, MedAE, Spearman, Calib, calib_factor) %>%
      distinct() %>%
      pivot_wider(
        names_from = version,
        values_from = c(RMSE, MAE, MedAE, Spearman, Calib, calib_factor),
        names_glue = "{.value}_{version}"
      ) %>%
      arrange(desc(Spearman_RAW))
    
    msg(">>> M2 RAW vs CAL (wide):")
    print(m2_wide)
    
    write.csv(m2_wide, OUT_CSV_M2, row.names = FALSE)
    msg("Saved: %s", OUT_CSV_M2)
  } else {
    msg(">>> M2 wide summary skipped: RAW and CAL not both present.")
  }
} else {
  msg(">>> M2 not found or no version labels present.")
}

# ---------------------------------------------------------
# 3) Topline: 1 row per family (rules)
#   - M2: prefer CAL if exists else RAW else max Spearman
#   - M1: choose max Spearman (assumes file already SAR-only)
#   - M0: choose max Spearman
# ---------------------------------------------------------
top_m0 <- res_all %>% filter(family == "M0 (Benchmark)") %>%
  slice_max(order_by = Spearman, n = 1, with_ties = FALSE)

top_m1 <- res_all %>% filter(family == "M1 (SAR 2-part)") %>%
  slice_max(order_by = Spearman, n = 1, with_ties = FALSE)

m2_pool <- res_all %>% filter(family == "M2 (STAR)")
top_m2 <- m2_pool
if (nrow(m2_pool) > 0) {
  if (any(m2_pool$version == "CAL", na.rm = TRUE)) {
    top_m2 <- m2_pool %>% filter(version == "CAL")
  } else if (any(m2_pool$version == "RAW", na.rm = TRUE)) {
    top_m2 <- m2_pool %>% filter(version == "RAW")
  }
  top_m2 <- top_m2 %>% slice_max(order_by = Spearman, n = 1, with_ties = FALSE)
}

topline <- bind_rows(top_m0, top_m1, top_m2) %>%
  arrange(family)

msg(">>> Topline (1 per family):")
print(topline)

write.csv(topline, OUT_CSV_TOP, row.names = FALSE)
msg("Saved: %s", OUT_CSV_TOP)

msg(">>> Done.")
