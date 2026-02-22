## =========================================================
## Wildfire_1_Setup.R
## 목적 : 초기 폴더 구조 생성 및 패키지 준비
## =========================================================
pkgs <- c("dplyr","purrr","tibble","stringr","httr","xml2","lubridate","readr", "data.table")
install.packages(setdiff(pkgs, rownames(installed.packages())))
lapply(pkgs, library, character.only = TRUE)

# 폴더 자동 생성 (이미 만드셨지만 코드에 넣어두면 재현성이 높아집니다)
dirs <- c("Wildfire_Data/raw_data", "Wildfire_Data/meta_data", "Wildfire_Data/processed_data")
sapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# API 키 설정 (환경 변수 권장)
# Sys.setenv(DATA_GO_KR_KEY = "본인의_서비스키_입력")
API_KEY <- Sys.getenv("DATA_GO_KR_KEY")