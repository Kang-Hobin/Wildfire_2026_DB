## =========================================================
## Wildfire_1_Setup.R
## =========================================================

# 1. 필수 패키지 목록 (lwgeom, units 추가)
pkgs <- c("dplyr", "purrr", "tibble", "stringr", "httr", "xml2", 
          "lubridate", "readr", "data.table", "sf", "lwgeom", "units")

# 2. 미설치 패키지 확인 및 설치
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) {
  install.packages(to_install)
  message(paste("다음 패키지를 설치했습니다:", paste(to_install, collapse = ", ")))
}

# 3. 라이브러리 로드
lapply(pkgs, library, character.only = TRUE)

message("--- 환경 설정 완료: 모든 라이브러리가 준비되었습니다. ---")