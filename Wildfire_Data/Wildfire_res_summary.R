## =========================================================
## Wildfire_res_summary.R
##
## 목적:
##  - Wildfire_Data/results/에 저장된 결과표 3종(M0/M1/M2)을 자동 로드
##  - 공통 컬럼로 정리 후 한 눈에 비교 가능한 요약 테이블 생성
##  - (특히) M2 STAR 결과는 RAW vs CAL을 wide 형식으로 요약
##
## 출력:
##  - 콘솔: 통합표 / M2 RAW-CAL 요약 / 베스트 모델 요약
##  - 파일(옵션): Wildfire_Data/results/res_summary_all.csv
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

RESULT_DIR <- "Wildfire_Data/results"
OUT_CSV_ALL <- file.path(RESULT_DIR, "res_summary_all.csv")
OUT_CSV_M2  <- file.path(RESULT_DIR, "res_summary_m2_raw_cal.csv")
OUT_CSV_TOP <- file.path(RESULT_DIR, "res_summary_topline.csv")

msg <- function(...) message(sprintf(...))

# ---------------------------------------------------------
# 0) helpers
# ---------------------------------------------------------
safe_read_rds <- function(path) {
  tryCatch(readRDS(path), error = function(e) NULL)
}

as_tibble_safe <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "data.frame")) return(as_tibble(x))
  # sometimes objects are lists with a data.frame inside
  if (is.list(x)) {
    for (nm in names(x)) {
      if (inherits(x[[nm]], "data.frame")) return(as_tibble(x[[nm]]))
    }
  }
  NULL
}

normalize_cols <- function(df) {
  # normalize common column names across scripts
  # expected candidates:
  # model / model_type / model_type.x
  # AUC_occ / AUC / auc
  # Calib / Calib_Ratio
  df <- as_tibble(df)
  
  # model column
  if (!"model" %in% names(df)) {
    if ("model_type" %in% names(df)) df <- df %>% rename(model = model_type)
  }
  
  # AUC column
  if (!"AUC_occ" %in% names(df)) {
    if ("AUC" %in% names(df)) df <- df %>% rename(AUC_occ = AUC)
    if ("auc" %in% names(df)) df <- df %>% rename(AUC_occ = auc)
  }
  
  # Calib column
  if (!"Calib" %in% names(df)) {
    if ("Calib_Ratio" %in% names(df)) df <- df %>% rename(Calib = Calib_Ratio)
  }
  
  # keep only columns we care about if present
  keep <- intersect(
    c("model","AUC_occ","RMSE","MAE","MedAE","Spearman","Calib","version","calib_factor"),
    names(df)
  )
  df <- df %>% dplyr::select(all_of(keep))
  
  df
}

infer_family <- function(path, df) {
  fn <- basename(path)
  
  # rough family tagging
  family <- case_when(
    str_detect(fn, regex("^m0_|m0", ignore_case = TRUE)) ~ "M0 (Benchmark)",
    str_detect(fn, regex("^m1_|sar2part|m1", ignore_case = TRUE)) ~ "M1 (SAR 2-part)",
    str_detect(fn, regex("^m2_|star|5models", ignore_case = TRUE)) ~ "M2 (STAR)",
    TRUE ~ "Other"
  )
  
  # if model names include RAW/CAL, mark as M2 star-ish
  if ("model" %in% names(df) && any(str_detect(df$model, "\\[CAL\\]|\\[RAW\\]"))) {
    family <- "M2 (STAR)"
  }
  
  family
}

strip_suffix <- function(x) gsub("\\s*\\[(RAW|CAL)\\]\\s*$", "", x)

# ---------------------------------------------------------
# 1) auto-discover candidate rds files
# ---------------------------------------------------------
msg(">>> Scanning result directory: %s", RESULT_DIR)
rds_files <- list.files(RESULT_DIR, pattern = "\\.rds$", full.names = TRUE)

# prefer likely files but keep all
candidates <- rds_files[str_detect(basename(rds_files),
                                   regex("m0|m1|m2|sar2part|star|perf", ignore_case=TRUE))]
if (length(candidates) == 0) {
  stop("No candidate .rds files found under: ", RESULT_DIR)
}

msg(">>> Found %d candidate RDS files.", length(candidates))
print(basename(candidates))

# ---------------------------------------------------------
# 2) read + normalize + bind
# ---------------------------------------------------------
all_tbls <- list()

for (p in candidates) {
  obj <- safe_read_rds(p)
  df  <- as_tibble_safe(obj)
  if (is.null(df)) next
  
  df <- normalize_cols(df)
  fam <- infer_family(p, df)
  
  df <- df %>%
    mutate(
      source_file = basename(p),
      family = fam
    )
  
  all_tbls[[length(all_tbls) + 1]] <- df
}

if (length(all_tbls) == 0) stop("No readable data.frame/tibble found in candidate RDS files.")

res_all <- bind_rows(all_tbls)

msg(">>> Combined rows: %d", nrow(res_all))
print(res_all)

# save combined csv
write.csv(res_all, OUT_CSV_ALL, row.names = FALSE)
msg("Saved combined table: %s", OUT_CSV_ALL)

# ---------------------------------------------------------
# 3) M2 STAR: RAW vs CAL wide summary (if exists)
# ---------------------------------------------------------
res_m2 <- res_all %>% filter(family == "M2 (STAR)")

if (nrow(res_m2) > 0 && any(str_detect(res_m2$model, "\\[(RAW|CAL)\\]"))) {
  res_m2_wide <- res_m2 %>%
    mutate(
      version = ifelse(!is.na(version), version,
                       ifelse(str_detect(model, "\\[CAL\\]$"), "CAL", "RAW")),
      model_base = strip_suffix(model)
    ) %>%
    dplyr::select(model_base, version, AUC_occ, RMSE, MAE, MedAE, Spearman, Calib, calib_factor, source_file) %>%
    distinct() %>%
    pivot_wider(
      names_from = version,
      values_from = c(RMSE, MAE, MedAE, Spearman, Calib, calib_factor),
      names_glue = "{.value}_{version}"
    ) %>%
    arrange(desc(Spearman_RAW))
  
  msg(">>> [M2] RAW vs CAL summary (wide):")
  print(res_m2_wide)
  
  write.csv(res_m2_wide, OUT_CSV_M2, row.names = FALSE)
  msg("Saved M2 RAW/CAL summary: %s", OUT_CSV_M2)
} else {
  msg(">>> [M2] RAW/CAL rows not detected. Skipping wide summary.")
}

# ---------------------------------------------------------
# 4) Topline: 1 best per family (with sensible preference rules)
#   - M2: prefer CAL if exists, else RAW
#   - M1: prefer SAR models (name contains "SAR")
#   - Exclude obviously broken/old outputs if needed
# ---------------------------------------------------------

strip_suffix <- function(x) gsub("\\s*\\[(RAW|CAL)\\]\\s*$", "", x)

res_all2 <- res_all %>%
  mutate(
    version2 = case_when(
      !is.na(version) ~ version,
      str_detect(model, "\\[CAL\\]$") ~ "CAL",
      str_detect(model, "\\[RAW\\]$") ~ "RAW",
      TRUE ~ NA_character_
    ),
    model_base = strip_suffix(model)
  ) %>%
  # (선택) 구버전/이상치 파일 걸러내기: 파일명에 "perf"만 있고 final이 아닌 등
  # 필요하면 아래 조건을 더 빡세게 바꿔도 됨
  filter(!str_detect(source_file, regex("m2_star_perf", ignore_case = TRUE))) %>%
  # Calib이 음수거나 Spearman이 극단적으로 이상한 경우는 보통 계산/정의 오류였던 산출물
  filter(is.na(Spearman) | Spearman > -0.5)

# ---- M2: CAL 우선, 없으면 RAW ----
m2_pref <- res_all2 %>% 
  filter(family == "M2 (STAR)") %>%
  {
    if (any(.$version2 == "CAL", na.rm = TRUE)) filter(., version2 == "CAL")
    else if (any(.$version2 == "RAW", na.rm = TRUE)) filter(., version2 == "RAW")
    else .
  } %>%
  slice_max(order_by = Spearman, n = 1, with_ties = FALSE)

# ---- M1: SAR 포함 모델 우선 ----
m1_pref <- res_all2 %>%
  filter(family == "M1 (SAR 2-part)") %>%
  {
    if (any(str_detect(.$model, regex("SAR", ignore_case = TRUE)))) {
      filter(., str_detect(model, regex("SAR", ignore_case = TRUE)))
    } else .
  } %>%
  slice_max(order_by = Spearman, n = 1, with_ties = FALSE)

# ---- M0: 그냥 Spearman 최고 ----
m0_pref <- res_all2 %>%
  filter(family == "M0 (Benchmark)") %>%
  slice_max(order_by = Spearman, n = 1, with_ties = FALSE)

# ---- Other families (if any) ----
other_pref <- res_all2 %>%
  filter(!family %in% c("M0 (Benchmark)", "M1 (SAR 2-part)", "M2 (STAR)")) %>%
  group_by(family) %>%
  slice_max(order_by = Spearman, n = 1, with_ties = FALSE) %>%
  ungroup()

topline <- bind_rows(m0_pref, m1_pref, m2_pref, other_pref) %>%
  dplyr::select(-model_base) %>%
  arrange(family)

msg(">>> Topline (preferred rules applied):")
print(topline)

write.csv(topline, OUT_CSV_TOP, row.names = FALSE)
msg("Saved topline: %s", OUT_CSV_TOP)