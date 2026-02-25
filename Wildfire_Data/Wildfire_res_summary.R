# =========================================================
# Wildfire_Final_Compare_Table.R
# - Load M0/M1/M2 performance RDS
# - Harmonize columns to the standard metric set
# - Create final comparison table + rankings
# - Save combined RDS/CSV
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

# ---------------------------------------------------------
# 0) Paths
# ---------------------------------------------------------
IN_M0 <- "Wildfire_Data/results/m0_final_perf.rds"
IN_M1 <- "Wildfire_Data/results/m1_sar_perf.rds"
IN_M2 <- "Wildfire_Data/results/m2_star_5models_perf_final.rds"

OUT_DIR <- "Wildfire_Data/results"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

OUT_RDS <- file.path(OUT_DIR, "final_perf_all_models.rds")
OUT_CSV <- file.path(OUT_DIR, "final_perf_all_models.csv")

OUT_RANK_CSV <- file.path(OUT_DIR, "final_perf_all_models_ranked.csv")

# Standard metric columns (must match M0/M1/M2)
STD_COLS <- c(
  "model_type",
  "AUC_occ",
  "RMSE",
  "MAE",
  "MedianAE",
  "RMSLE",
  "Spearman",
  "BR_rate_down",
  "BR_rate_up",
  "BR_sev_down",
  "BR_sev_up"
)

# ---------------------------------------------------------
# 1) Coerce helper
# ---------------------------------------------------------
coerce_perf <- function(x) {
  x <- as.data.frame(x)
  
  # Legacy name fixes
  if ("model" %in% names(x) && !"model_type" %in% names(x)) {
    x$model_type <- x$model
    x$model <- NULL
  }
  if ("MedAE" %in% names(x) && !"MedianAE" %in% names(x)) {
    x$MedianAE <- x$MedAE
    x$MedAE <- NULL
  }
  
  # If any missing standard cols, create as NA
  missing <- setdiff(STD_COLS, names(x))
  for (cn in missing) x[[cn]] <- NA_real_
  
  # Keep only standard cols in order
  x <- x[, STD_COLS, drop = FALSE]
  
  # Coerce numeric cols
  num_cols <- setdiff(STD_COLS, "model_type")
  for (cn in num_cols) x[[cn]] <- as.numeric(x[[cn]])
  
  as_tibble(x)
}

# ---------------------------------------------------------
# 2) Load & combine
# ---------------------------------------------------------
m0 <- coerce_perf(readRDS(IN_M0))
m1 <- coerce_perf(readRDS(IN_M1))
m2 <- coerce_perf(readRDS(IN_M2))

final <- bind_rows(m0, m1, m2)

# Optional: add group labels for readability (M0/M1/M2)
final <- final %>%
  mutate(
    model_group = case_when(
      str_detect(model_type, "^EP_") ~ "M0",
      str_detect(model_type, "^M1_") ~ "M1",
      str_detect(model_type, "^M[1-5]:") ~ "M2",
      TRUE ~ "Other"
    )
  )

# Show final table
print(final)

# Save combined
saveRDS(final, OUT_RDS)
readr::write_csv(final, OUT_CSV)

message(sprintf("Saved combined RDS: %s", OUT_RDS))
message(sprintf("Saved combined CSV: %s", OUT_CSV))

# ---------------------------------------------------------
# 3) Ranking table (for quick final comparison)
#    - Lower is better: RMSE, MAE, MedianAE, RMSLE, BR_* (all)
#    - Higher is better: AUC_occ, Spearman
# ---------------------------------------------------------
ranked <- final %>%
  mutate(
    rank_AUC      = dense_rank(desc(AUC_occ)),
    rank_Spear    = dense_rank(desc(Spearman)),
    rank_RMSE     = dense_rank(RMSE),
    rank_MAE      = dense_rank(MAE),
    rank_MedianAE = dense_rank(MedianAE),
    rank_RMSLE    = dense_rank(RMSLE),
    rank_BRrd     = dense_rank(BR_rate_down),
    rank_BRru     = dense_rank(BR_rate_up),
    rank_BRsd     = dense_rank(BR_sev_down),
    rank_BRsu     = dense_rank(BR_sev_up)
  ) %>%
  # Simple composite (equal weights): lower is better
  mutate(
    rank_total = rank_AUC + rank_Spear +
      rank_RMSE + rank_MAE + rank_MedianAE + rank_RMSLE +
      rank_BRrd + rank_BRru + rank_BRsd + rank_BRsu
  ) %>%
  arrange(rank_total, RMSE, MAE)

print(ranked %>% dplyr::select(model_group, model_type, rank_total, everything()))

readr::write_csv(ranked, OUT_RANK_CSV)
message(sprintf("Saved ranked CSV: %s", OUT_RANK_CSV))
