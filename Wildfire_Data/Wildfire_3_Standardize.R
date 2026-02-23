## =========================================================
## Wildfire_3_Standardize.R
## 목적: ASOS/AWS/MTW ID 중복 방지를 위한 접두사 부여
## =========================================================

# [A] ASOS 표준화 (ID 앞에 'A' 추가)
asos_std <- read_csv(file.path("wildfire_Data/raw_data/ASOS_raw.csv"), 
                     locale = locale(encoding = "CP949"), show_col_types = FALSE) %>%
  transmute(
    stn_id   = paste0("A", as.integer(지점)), # 예: A175
    date     = as.Date(일시),
    ta_avg   = `평균기온(°C)`,
    hm_avg   = `평균 상대습도(%)`,
    ws_max   = `최대 풍속(m/s)`,
    prcp_day = `일강수량(mm)`
  )
saveRDS(asos_std, file.path("wildfire_Data/processed_data/asos_daily.rds"))

# [B] AWS 표준화 (ID 앞에 'W' 추가)
raw_bytes <- readBin(file.path("wildfire_Data/raw_data/AWS_raw.csv"), "raw", file.info(file.path("wildfire_Data/raw_data/AWS_raw.csv"))$size)
aws_txt <- iconv(list(raw_bytes), from = "EUC-KR", to = "UTF-8")[[1]]

aws_std <- data.table::fread(text = aws_txt, data.table = FALSE) %>%
  transmute(
    stn_id   = paste0("W", as.integer(지점)), # 예: W175
    date     = as.Date(일시),
    ta_avg   = `평균기온(°C)`,
    ws_max   = `최대 순간 풍속(m/s)`,
    prcp_day = `일강수량(mm)`
  )
saveRDS(aws_std, file.path("wildfire_Data/processed_data/ws_daily.rds"))

# [C] MTW 표준화 (ID 앞에 'M' 추가)
mtw_std <- readRDS(file.path("wildfire_Data/raw_data/raw_mtw_API_2017_2025.rds")) %>%
  transmute(
    stn_id = paste0("M", as.integer(obsid)), # 예: M47192
    date   = as.Date(tm),
    ta_avg = as.numeric(tm2m),
    hm_avg = as.numeric(hm2m),
    ws_max = as.numeric(ws10m),
    prcp_day = as.numeric(rn)
  )
saveRDS(mtw_std, file.path("wildfire_Data/processed_data/mtw_raw_standard.rds"))