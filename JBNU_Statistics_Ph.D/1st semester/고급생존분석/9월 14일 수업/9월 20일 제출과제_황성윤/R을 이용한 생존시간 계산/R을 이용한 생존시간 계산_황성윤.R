install.packages("tidyverse")
install.packages("lubridate")
library(tidyverse)
library(lubridate)

## base R을 이용한 생존시간 계산
### 태어난 날과 금일 날짜 입력
date_ex <- 
  tibble(
    birth_date = c("1989-04-18"), 
    last_fup_date = c("2020-09-14")
  )
glimpse(date_ex)
### date 형식으로 데이터 변환
(date_r <- date_ex %>% 
      mutate(
      birth_date = as.Date(birth_date, format = "%Y-%m-%d"), 
      last_fup_date = as.Date(last_fup_date, format = "%Y-%m-%d") 
    ))
### Calculating survival times
date_r %>% 
  mutate(
    os_yrs = 
      as.numeric(
        difftime(last_fup_date, 
                 birth_date, 
                 units = "days"))
  )

## lubridate 패키지를 이용한 생존시간 계산
### Formatting dates
(date_l <- date_ex %>% 
    mutate(
      birth_date = ymd(birth_date), 
      last_fup_date = ymd(last_fup_date)
    ))
### Calculating survival times
date_l %>% 
  mutate(
    os_yrs = 
      (as.duration(birth_date %--% last_fup_date) / dyears(1)) * 365.25
  )
