install.packages("survival")
library(survival)

#### retinopathy data
# 목적 : 당뇨병 환자의 시력 손실을 지연시키기 위한 레이저치료법의 효과를 확인하기 위해 수집
# 총 197명의 환자에 대해서 한쪽 눈에는 레이저치료를 실시하고 다른 눈에는 치료를 실시하지 않은 상태에서 
# 시력상실이 발생할 때까지의 시간(futime)을 관찰함.
# 그룹 : 변수 trt (0 = no treatment, 1 = laser)
# 중도절단 여부 : 변수 status (0 = censored, 1 = visual loss)
data(retinopathy, packages="survival")
head(retinopathy, 20)

### Q1 : 그룹별로 관측값과 중도절단(censored) 개수의 비율 산출
retinopathy_notreat <- subset(retinopathy, trt==0)
head(retinopathy_notreat, 10)
# 레이저치료를 실시하지 않은 그룹에서 관측값 개수 산출
sum(retinopathy_notreat$status) 
# 레이저치료를 실시하지 않은 그룹에서 중도절단 개수 산출
nrow(retinopathy_notreat) - sum(retinopathy_notreat$status) 
# 비율 = 101:96

retinopathy_laser <- subset(retinopathy, trt==1)
head(retinopathy_laser, 10)
# 레이저치료를 실시한 그룹에서 관측값 개수 산출
sum(retinopathy_laser$status) 
# 레이저치료를 실시한 그룹에서 중도절단 개수 산출
nrow(retinopathy_laser) - sum(retinopathy_laser$status)
# 비율 = 54:143

### Q2 : 그룹별로 관측값 데이터(status==1)만을 고려한 생존시간 평균 산출
# 레이저치료를 실시하지 않은 그룹
aggregate(futime~status, retinopathy_notreat, mean) # answer = 18.949
# 레이저치료를 실시한 그룹
aggregate(futime~status, retinopathy_laser, mean) # answer = 18.222

### Q3 : 그룹별로 관측값과 중도절단 데이터를 모두 포함한 생존시간 평균 산출
tapply(retinopathy$futime, retinopathy$trt, mean) # answer = 32.288(no treat), 38.871(laser)


#### 생존분석 사례 조사
## url : http://kostat.go.kr/file_total/2014_pp_03.pdf
## 논문제목 : 생존분석을 이용한 보육교사 이직률 결정요인 분석
## 데이터수집 : 고용노동부와 한국고용정보원의 '대졸자 직업이동 경로조사 : GOMS'
## 데이터 : 매년 2~3년제 대학 이상 고등교육 과정을 이수한 졸업자를 모집단으로 하여 
# 졸업년도 다음 해 9월에서 11월경 조사하고, 그로부터 2년마다 1회 추적조사 추가 실시
## 변수
# 종속변수 : 보육교사 재직기간
# 설명변수 : 시간당 임금, 학력, 시설규모, 시설 소재 지역, 종사상 지위, 혼인상태, 성별, 연령 
## 분석방법 : Kaplan-Meier estimation, Cox propotional hazard model, Competing risk model
## 분석문제와 결과 : 시간당 임금, 혼인상태가 보육교사 이직률에 결정적인 역할을 하는 것으로 나타남.
# 다만, 학력의 효과는 2005년 자료와 2008년 자료에서 다르게 나타났음.