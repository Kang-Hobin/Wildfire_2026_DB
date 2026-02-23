## =========================================================
## Wildfire_2_API_Harvest.R
## 목적: 산림청(산불) 및 기상청(산림기상) API를 통해 원천 데이터 수집
## 생성 파일: 
##   1) Wildfire_Data/raw_data/raw_fire_API_2017_2025.rds
##   2) Wildfire_Data/raw_data/raw_mtw_API_2017_2025.rds
##   3) Wildfire_Data/raw_data/mtw_1200_failed.csv (오류 로그)
## =========================================================

# 1. 환경 설정 및 API 키 체크
library(dplyr)
library(httr)
library(xml2)
library(purrr)
library(readr)
library(lubridate)

API_KEY <- Sys.getenv("DATA_GO_KR_KEY")
if (API_KEY == "") {
  stop("API 키가 설정되지 않았습니다. Sys.setenv(DATA_GO_KR_KEY = '본인키')를 먼저 실행하세요.")
}

# 2. 공통 헬퍼 함수 (API 통신 및 XML 파싱)
api_get_text <- function(url, query, timeout_sec = 120){ # 타임아웃 120초로 연장
  res <- tryCatch({
    httr::GET(url, query = query, httr::timeout(timeout_sec))
  }, error = function(e) {
    message("연결 에러 발생: ", e$message)
    return(NULL)
  })
  
  if (is.null(res)) return(NULL)
  
  httr::stop_for_status(res)
  httr::content(res, as = "text", encoding = "UTF-8")
}

parse_openapi_items_safe <- function(xml_txt){
  doc <- tryCatch(xml2::read_xml(xml_txt), error = function(e) return(NULL))
  if (is.null(doc)) return(list(ok = FALSE, data = tibble(), meta = list(resultMsg = "Invalid XML")))
  
  rc <- xml2::xml_text(xml2::xml_find_first(doc, ".//resultCode"))
  rm <- xml2::xml_text(xml2::xml_find_first(doc, ".//resultMsg"))
  
  if (!identical(rc, "00")) {
    return(list(ok = FALSE, data = tibble(), meta = list(resultCode = rc, resultMsg = rm)))
  }
  
  items <- xml2::xml_find_all(doc, ".//items/item")
  df <- purrr::map_dfr(items, function(it){
    ch <- xml2::xml_children(it)
    tibble::as_tibble(stats::setNames(as.list(xml2::xml_text(ch)), xml2::xml_name(ch)))
  })
  
  totalCount <- suppressWarnings(as.integer(xml2::xml_text(xml2::xml_find_first(doc, ".//totalCount"))))
  list(ok = TRUE, data = df, meta = list(resultCode = rc, resultMsg = rm, totalCount = totalCount))
}

# 3. [A] 산불 발생 데이터 수집 함수 (페이징 처리 완결판)
# ---------------------------------------------------------
get_fire_stats_all <- function(st_date, ed_date) {
  message(">>> [Phase 1] 전체 데이터 건수 확인 중...")
  
  # 1. 우선 1건만 요청해서 전체 개수(totalCount) 파악
  init_q <- list(ServiceKey = API_KEY, searchStDt = st_date, searchEdDt = ed_date, 
                 numOfRows = 1, pageNo = 1, `_type` = "xml")
  init_txt <- api_get_text(FIRE_URL, init_q)
  init_parsed <- parse_openapi_items_safe(init_txt)
  
  if (!init_parsed$ok) stop("초기 연결 실패: ", init_parsed$meta$resultMsg)
  
  total_cnt <- init_parsed$meta$totalCount
  rows_per_page <- 1000
  total_pages <- ceiling(total_cnt / rows_per_page) # 총 페이지 수 계산
  
  message(paste(">>> 총", total_cnt, "건 발견 (약", total_pages, "회 추가 호출 필요)"))
  
  # 2. 모든 페이지를 순회하며 데이터 수집
  all_data <- tibble()
  
  for (p in 1:total_pages) {
    message(paste(">>> 산불 데이터 수집 중... (", p, "/", total_pages, "페이지)"))
    
    q <- list(ServiceKey = API_KEY, searchStDt = st_date, searchEdDt = ed_date, 
              numOfRows = rows_per_page, pageNo = p, `_type` = "xml")
    
    txt <- api_get_text(FIRE_URL, q)
    parsed <- parse_openapi_items_safe(txt)
    
    if (parsed$ok) {
      all_data <- bind_rows(all_data, parsed$data)
    } else {
      message(paste("!!! ", p, "페이지 수집 실패. 건너뜁니다."))
    }
    
    Sys.sleep(0.2) # API 서버 매너 타임
  }
  
  return(all_data)
}

# 4. [B] 산림 기상(MTW) 데이터 수집 함수 (일별 반복 및 캐시 기능)
MTW_URL <- "http://apis.data.go.kr/1400377/mtweather/mountListSearch"

fetch_mtw_robust <- function(start_date, end_date, cache_path = "Wildfire_Data/raw_data/mtw_cache.rds") {
  dates <- seq.Date(as.Date(start_date), as.Date(end_date), by = "day")
  
  # 캐시 불러오기 (중단된 부분부터 다시 시작 가능)
  if (file.exists(cache_path)) {
    acc <- readRDS(cache_path)
    done_dates <- as.Date(unique(substr(acc$tm, 1, 8)), format = "%Y%m%d")
    dates <- dates[!(dates %in% done_dates)]
    message(paste("캐시 확인: 이미 수집된", length(done_dates), "일치를 제외하고 진행합니다."))
  } else {
    acc <- tibble()
  }
  
  failed_log <- tibble()
  
  for (d in as.character(dates)) {
    d_str <- format(as.Date(d), "%Y%m%d")
    tm_str <- paste0(d_str, "1200") # 매일 낮 12시 데이터 기준
    
    message(paste("MTW 수집 중:", d_str))
    
    res <- tryCatch({
      txt <- api_get_text(MTW_URL, list(ServiceKey = API_KEY, tm = tm_str, numOfRows = 1000, `_type` = "xml"))
      parse_openapi_items_safe(txt)
    }, error = function(e) list(ok = FALSE, meta = list(resultMsg = "Request Error")))
    
    if (res$ok) {
      acc <- bind_rows(acc, res$data)
      saveRDS(acc, cache_path) # 매일 단위로 임시 저장
    } else {
      failed_log <- bind_rows(failed_log, tibble(date = d, msg = res$meta$resultMsg))
      write_csv(failed_log, "Wildfire_Data/raw_data/mtw_1200_failed.csv")
    }
    Sys.sleep(0.1) # 서버 부하 방지용 짧은 휴식
  }
  return(acc)
}

## =========================================================
## 5. 실행 및 파일 저장
## =========================================================

# 기간 설정 (2017-01-01 ~ 2025-12-31)
ST <- "20170101"
ED <- "20251231"

# [실행 1] 산불 데이터
fire_raw_final <- get_fire_stats_all(ST, ED)
saveRDS(fire_raw_final, "Wildfire_Data/raw_data/raw_fire_API_2017_2025.rds")

# [실행 2] 산림 기상 데이터
mtw_raw_final <- fetch_mtw_robust("2017-01-01", "2025-12-31")
saveRDS(mtw_raw_final, "Wildfire_Data/raw_data/raw_mtw_API_2017_2025.rds")

message("=========================================================")
message("수집 완료! Wildfire_Data/raw_data 폴더를 확인하세요.")
message("=========================================================")