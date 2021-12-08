######################################
###### Time series analysis
###### CH 10 Seasonal ARIMA
###### 202102
######################################

##############  package
library(sarima)
library(ggplot2)
library(forecast)
library(data.table)
library(gridExtra)
library(fUnitRoots)
library(portes)
library(lmtest)


setwd("C:\\R-Project\\DAT\\Time Series Data")


#################################################
##### EX 10.2
#################################################

tour <- scan("tourist.txt")

dt <- data.table(t = 1:length(tour),
                 z = tour)

## 시계열그림 / ACF PACF

p3 <- ggplot(dt, aes(t, z)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ ggtitle('국내 입국 관광객 자료의 시계열그림')+
  theme_bw()

p1 <- ggAcf(dt$z) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_blank())

p2 <- ggPacf(dt$z) + 
  theme_bw() +ylim(-1,1) +
  theme(plot.title = element_blank())

grid.arrange(p3, p1, p2, nrow = 2,
             layout_matrix = rbind(c(1,1),
                                   c(2,3)))

##### 정상화 : 로그변환

dt[, lnz := log(z)]
dt[, sqrtz := sqrt(z)]
dt[, boxcoxz := BoxCox(z,lambda= BoxCox.lambda(z))]

# 세 변환 비교
melt.dt <- melt(dt, id=1)
ggplot(melt.dt, aes(t, value)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ 
  facet_wrap(variable~.,nrow=2, scales = "free_y")+
  theme_bw()

bptest(lm(z~t2, dt))
bptest(lm(lnz~t, dt))
bptest(lm(sqrtz~t, dt))
bptest(lm(boxcoxz~t, dt))

p3 <- ggplot(dt, aes(t, lnz)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ ggtitle('로그변환된 자료의 시계열그림')+
  theme_bw()

p1 <- ggAcf(dt$lnz) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_blank())

p2 <- ggPacf(dt$lnz) + 
  theme_bw() +ylim(-1,1) +
  theme(plot.title = element_blank())

grid.arrange(p3, p1, p2, nrow = 2,
             layout_matrix = rbind(c(1,1),
                                   c(2,3)))

##### 정상화 : 계절차분

dt[, df12_lnz := c(rep(0,12), diff(lnz,12))]
sdt <- dt[-(1:12)]

p3 <- ggplot(sdt, aes(t, df12_lnz)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ ggtitle('로그변환과 계절차분된 자료의 시계열그림')+
  theme_bw()

p1 <- ggAcf(sdt$df12_lnz) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_blank())

p2 <- ggPacf(sdt$df12_lnz) + 
  theme_bw() +ylim(-1,1) +
  theme(plot.title = element_blank())

grid.arrange(p3, p1, p2, nrow = 2,
             layout_matrix = rbind(c(1,1),
                                   c(2,3)))

## 단위근검정 H0 : phi=1
adfTest(sdt$df12_lnz, lags = 0, type = "c")
adfTest(sdt$df12_lnz, lags = 1, type = "c")
adfTest(sdt$df12_lnz, lags = 2, type = "c")


##### 정상화 : 차분
sdt[, df1_df12_lnz := c(0, diff(df12_lnz))]
ssdt <- sdt[-1,]

p3 <- ggplot(ssdt, aes(t, df1_df12_lnz)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ ggtitle('로그변환과 계절차분 및 1차 차분된 자료의 시계열그림')+
  theme_bw()

p1 <- ggAcf(ssdt$df1_df12_lnz, lag.max = 40) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_blank())

p2 <- ggPacf(ssdt$df1_df12_lnz, lag.max = 40) + 
  theme_bw() +ylim(-1,1) +
  theme(plot.title = element_blank())

grid.arrange(p3, p1, p2, nrow = 2,
             layout_matrix = rbind(c(1,1),
                                   c(2,3)))

## 단위근검정 H0 : phi=1
adfTest(sdt$df1_df12_lnz, lags = 0, type = "nc")
adfTest(sdt$df1_df12_lnz, lags = 1, type = "nc")
adfTest(sdt$df1_df12_lnz, lags = 2, type = "nc")


## 여러가지 acf/PACF 
## sma, sar : 12시차 24시차만 유의
## arima(1,0,0)(1,0,0)_12 : 12시차 양쪽으로 유의

### 잔차검정


## df1_df12_lnz 계절차분, 차분을 시행한 데이터 => ma의 차수 q=1, sma 차수 Q=1
## lnz  => ARIMA(0,1,1)(0,1,1)_12 : 선택
## lnz  => ARIMA(1,1,0)(1,1,0)_12
## lnz  => ARIMA(1,1,0)(0,1,1)_12

### 모형적합 

fit1 = arima(dt$lnz, order = c(0,1,1), seasonal = list(order = c(0,1,1), 
                                                      period = 12))
fit1

## 잔차분석
dt[, res := as.numeric(resid(fit1))]

p3 <- ggplot(dt, aes(t, res)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ ggtitle('잔차의 시계열그림')+
  theme_bw()

p1 <- ggAcf(dt$res) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_blank())

p2 <- ggPacf(dt$res) + 
  theme_bw() +ylim(-1,1) +
  theme(plot.title = element_blank())

grid.arrange(p3, p1, p2, nrow = 2,
             layout_matrix = rbind(c(1,1),
                                   c(2,3)))

LjungBox(fit1, lags=seq(6,24,6))# 잔차의 포트맨토검정  

### 모형적합  ARIMA(0,1,1)(1,1,0)_12

fit2 = arima(dt$lnz, order = c(0,1,1), 
             seasonal = list(order = c(1,1,0),  period = 12));  fit2

### 잔차검정

dt[, res2 := as.numeric(resid(fit2))]

p3 <- ggplot(dt, aes(t, res2)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ ggtitle('잔차의 시계열그림')+
  theme_bw()

p1 <- ggAcf(dt$res2) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_blank())

p2 <- ggPacf(dt$res2) + 
  theme_bw() +ylim(-1,1) +
  theme(plot.title = element_blank())

grid.arrange(p3, p1, p2, nrow = 2,
             layout_matrix = rbind(c(1,1),
                                   c(2,3)))
LjungBox(fit2, lags=seq(6,24,6))# 잔차의 포트맨토검정 



###################################################
##########
auto.arima(ssdt$df1_df12_lnz)  ### ARIMA(2,1,0)(0,1,0)_12
auto.arima(ssdt$df1_df12_lnz, trace=T, ic='aic')


fit3 = arima(ssdt$df1_df12_lnz, order = c(0,0,1), 
             seasonal = list(order = c(0,0,1),  period = 12))
fit3

##################################################
##########################  simulation 
##################################################

### ARIMA(p,d,q)(P,D,Q)_s
##sma : ARIMA(0,0,0)(0,0,1)_12
x <- sim_sarima(n=100,
                model=list(sma=-0.6, nseasons=12, sigma2 = 1)) # SMA(1)
y <- sim_sarima(n=100,
                model=list(sma=0.6, nseasons=12, sigma2 = 1)) # SMA(1)

p1 <- ggAcf(x, lag.max=60) + ylab("") +  
  theme_bw() + ylim(-1,1) + ggtitle('ACF : THETA=0.6') + 
  theme(plot.title = element_text(size=10))

p2 <- ggPacf(x, lag.max=60) + ylab("") +
  theme_bw() +ylim(-1,1)  + ggtitle('PACF : THETA=0.6')+ 
  theme(plot.title = element_text(size=10))

p3 <- ggAcf(y, lag.max=60) + ylab("") +
  theme_bw() + ylim(-1,1) + ggtitle('ACF : THETA=-0.6')+ 
  theme(plot.title = element_text(size=10))

p4 <- ggPacf(y, lag.max=60) + ylab("") +
  theme_bw() +ylim(-1,1) + ggtitle('PACF : THETA=-0.6')+ 
  theme(plot.title = element_text(size=10))


grid.arrange(p1, p2, p3, p4, nrow=2)


#######################################
#######################################
### train-data , test_data
### 

