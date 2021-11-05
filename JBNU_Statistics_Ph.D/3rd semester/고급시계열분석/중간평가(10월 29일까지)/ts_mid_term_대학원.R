######################################
###### Mid-Term
###### Time series analysis
###### 2021.10.
######################################


######################################
######################################
library(lmtest)
library(data.table)
library(ggplot2)

setwd("C:\\R-Project\\DAT\\Time Series Data")
######################################
######################################


######################################
## 1. data1.csv
######################################
# ‘data1.csv’는 모의 실험에 의해 생성된 시계열자료이다.
# 1) 적절한 추세 모형을 적합시킨 후 잔차분석을 하여라
# 2) 1)에서 적합한 추세모형에서 1~10 시차 후의 예측값을 구하여라.
# 3) 적절한 평활법을 적용한 후 잔차분석을 하여라.
# 4) 3)의 결과를 이용하여 1~10시차 후의 예측값을 구하여라.
# 5) 실제 1~10차 후의 관측값이 ‘data1_new.csv’일 때, 
#     2), 4)의 결과를 이용하여 어느 모형이 더 적합했는지에 대해 비교하여라.

dt1 <- read.csv('data1.csv')[,2:3]
dt1$t2 <- (1:100)^2

plot(dt1$t, dt1$z, type='l')  ## 선형추세

### 선형 추세 모형
m1 <- lm(z~t, data=dt1)
summary(m1)  ## 추세가 유의하

### 잔차분석
resid_m1 <- resid(m1)
plot(resid_m1, type='l')
abline(h=0, lty=2)

dwtest(m1)
dwtest(m1, alternative="less")

t.test(resid_m1)

shapiro.test(resid_m1)  ##H0 : normal distribution
qqnorm(resid(m1))
qqline(resid_m1, col = 2)

acf(resid_m1)
pacf(resid_m1)
sum(resid_m1^2) ##SSE


### 예측값 
pred_y1 <- predict(m1, newdata = data.frame(t=101:110))


### 평활법 : 2중지수 평활법
m2 <- HoltWinters(dt1$z,  gamma=F)
m2

plot(m2)


### 평활법 : 단순지수 평활법
m21 <- HoltWinters(dt1$z, beta=F,  gamma=F)
m21

plot(m21)


m2$SSE
m21$SSE


### 평활법 : 계절 평활법
a <- ts(dt1$z, frequency=12)

m22 <- HoltWinters(a)
m22

plot(m22)

m22$SSE




### 잔차분석
resid_m2 <- resid(m2)
plot(resid_m2, type='l')
abline(h=0, lty=2)

dwtest(resid_m2~1)
dwtest(resid_m2~1, alternative="less")

shapiro.test(resid_m2)  ##H0 : normal distribution
qqnorm(resid_m2)
qqline(resid_m2, col = 2)

acf(resid_m2)
sum(resid_m2^2) ##SSE

### 예측값
pred_y2 <- predict(m2, n.ahead = 10)
plot(m2, pred_y2)

### 실제데이터와 비교
newdt <- read.csv("data1_new.csv")[,2:3]

plot(newdt$z, type='l', ylim=c(5,15))
lines(pred_y1, col='red', lty=2)
lines(as.numeric(pred_y2), col='blue', lty=2)

MSE_trend <- sum((newdt$z-pred_y1)^2)/10
MSE_trend
MSE_smooth <- sum((newdt$z-as.numeric(pred_y2))^2)/10
MSE_smooth

library(forecast)
a <- forecast(m2, h=10)
sum((a$mean - pred_y1)^2)



######################################
## 2. usapass.txt
######################################
# ‘usapass.txt’는 미국 월별 비행기 승객 수(단위 : 천 명)의  시계열자료이다. log 변환 후 아래의 분석을 수행하시오.
# 1) 왜 log 변환이 필요한지에 대해  간단히 설명하여라.
# 2) 적절한 추세모형 적합 후 잔차분석을 하여라.
# 3) 적절한 평활법을 적용한 후 잔차분석을 하여라.
# 4) 적절한 분해법에 의해 각 성분을 분해해여 시계열 그림을 그려라.
# 5) 4)에서 추정된 불규칙성분을 통해 적용된 분해법이 적절했는지 논하여라.

z <- scan("usapass.txt") 
t <- 1:length(z)


plot(z, type='l')
lnz <- log(z)
plot(lnz, type='l')

####### 추세+계절성분

#### 1. 선형추세 +지시함수

dt1 <- data.frame( t = 1:length(z),
                   s = as.factor(rep(1:12,12)),
                   lnz=lnz)

m1 <- lm(lnz ~ t+s, data=dt1 )   
summary(m1)  ## 모두 유의

plot(lnz, type='l')
lines(fitted(m1), col='red', lty=2)

### 1. 선형추세 + 지시함수 : 잔차분석
resid_m1 <- resid(m1)
plot(resid_m1, type='l')
abline(h=0, lty=2)

dwtest(m1)  ##독립아닙

shapiro.test(resid_m1)  ##H0 : normal distribution
qqnorm(resid_m1)
qqline(resid_m1, col = 2)

acf(resid_m1) ## WN 아님
sum(resid_m1^2) ##SSE


#### 2. 선형추세 + 삼각함수
dt2 <-as.data.table(dt1[,c(1,3)]) 

dt2 <- cbind(dt2,
             dt2[, lapply(as.list(1:5), 
                          function(i) sin(2*pi*i/12*t))])

names(dt2)[-(1:2)]<-paste("sin", c(12,6,4,3,2.4), sep="_")

dt2 <- cbind(dt2,
             dt2[, lapply(as.list(1:5), 
                          function(i) cos(2*pi*i/12*t))])

names(dt2)[-(1:7)]<-paste("cos", c(12,6,4,3,2.4), sep="_")

m2 <- lm(lnz ~ ., data=dt2 )   
summary(m2)  ## 모두 유의

plot(lnz, type='l')
lines(fitted(m2), col='red', lty=2)

### 2. 선형추세 + 삼각함수 : 잔차분석
resid_m2 <- resid(m2)
plot(resid_m2, type='l')
abline(h=0, lty=2)

dwtest(m2)  ##독립아닙

shapiro.test(resid_m2)  ##H0 : normal distribution
qqnorm(resid_m2)
qqline(resid_m2, col = 2)

acf(resid_m2) ## WN 아님
sum(resid_m2^2) ##SSE



### 3. 평활법 : 계절 평활법
lnz <- ts(lnz, frequency=12)

m3 <- HoltWinters(lnz, seasonal='additive')
m3

plot(m3)

### 잔차분석
resid_m3 <- resid(m3)
plot(resid_m3, type='l')
abline(h=0, lty=2)

dwtest(resid_m3~1)

shapiro.test(resid_m3)  ##H0 : normal distribution
qqnorm(resid_m3)
qqline(resid_m3, col = 2)

acf(resid_m3)
sum(resid_m3^2) ##SSE


### 4. 분해법 : 이동평균법
m4 <- stl(lnz, s.window=12)
plot(m4)


### 4, 분해법 : 이동평균법 - 잔차분석
resid_m4 <- m4$time.series[,3]
plot(resid_m4, type='l')
abline(h=0, lty=2)

dwtest(resid_m4~1)

shapiro.test(resid_m4)  ##H0 : normal distribution
qqnorm(resid_m4)
qqline(resid_m4, col = 2)

acf(resid_m4)
sum(resid_m4^2) ##SSE


### 4. 분해법 : 추세법
m_t <- lm(lnz ~ t, data=dt1 ) 
dt1$trend <- fitted(m_t)

plot(lnz, type='l')
lines(dt1$trend, lty=2, col='red')


dif <- dt1$lnz - dt1$trend
m_s <- lm(dif~s, data=dt1)
dt1$season <- fitted(m_s)

plot(fitted(m_s), type='l')

plot(lnz, type='l')
lines(dt1$trend+dt1$season, lty=2, col='red')


### 4, 분해법 : 추세법 - 잔차분석
resid_m4 <- dt1$lnz - dt1$trend - dt1$season
plot(resid_m4, type='l')
abline(h=0, lty=2)

dwtest(resid_m4~1)

shapiro.test(resid_m4)  ##H0 : normal distribution
qqnorm(resid_m4)
qqline(resid_m4, col = 2)

acf(resid_m4)
pacf(resid_m4)
sum(resid_m4^2) ##SSE








######################################
## 3. 
######################################

## Zt = 1 + 0.9Z(t-1) + et, t=1,2,...,100

set.seed(1234)

sim_ar <- function(n, mu, phi){
  ### n : sample size
  ### mu : mean
  ### phi : p-dim coefficients
  p <- length(phi)
  z <- mu + rnorm(n+100)  #epsilon ~ WN(0, sigma^2), iid N(0,1)
  # delta <- (1-sum(phi))*mu
  delta <- 1
  for (k in (length(phi)+1):(n+100)){
    z[k] <- delta + sum(z[(k-1):(k-p)]*phi) + rnorm(1)
  }
  # return(z[-(1:100)])
  return(z)
}

plot(sim_ar(100, 10, 0.9), type='l')

par(mfrow=c(2,1))

z1 <- 10+arima.sim(n=100, list(order=c(1,0,0), ar=0.9))
plot.ts(z1)

z2 <- arima.sim(n=100, list(order=c(1,0,0), ar=0.9))
plot.ts(z1)

a <- sim_ar(100, 10, 0.9)
plot(a, type='l')


#############################
##
################################
arima.sim(n=100, list(order=c(0,0,1), ma=-0.5))

## zt = et - theta*e(t-1)  ~ MA(1), theta
## zt= et - 0.5 * e(t-1), theta=0.5











