######################################
###### Time series analysis
###### CH 8 ARMA Fit
###### 202102
######################################


##############  package
library(data.table)
library(ggplot2)
library(forecast)
library(gridExtra)
library(lmtest)  ##coeftest
library(portes)  ##LjungBox
library(astsa)
library(fUnitRoots)  # library for function adfTest
library(tseries)  ##JB test
setwd("C:\\R-Project\\DAT\\Time Series Data")

#################################################
##### EX 8.6 가스로 자료 분석
#################################################

z <- scan("gas.txt", what=list(0,0))

dt <- data.table( t = 1:length(z[[1]]),
                  rate = z[[1]],
                  co2 = z[[2]])
## 시계열그림 / ACF PACF
p3 <- ggplot(dt, aes(t, rate)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ ggtitle('Time series plot of the input Gas feed rate')+
  theme_bw()

p1 <- ggAcf(dt$rate) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_blank())

p2 <- ggPacf(dt$rate) + 
  theme_bw() +ylim(-1,1) +
  theme(plot.title = element_blank())



grid.arrange(p3, p1, p2, nrow = 2,
             layout_matrix = rbind(c(1,1),
                                   c(2,3)))

LjungBox(dt$rate, lags=seq(6,24,6))#포트맨토 검정 

## 단위근 검정 
adfTest(dt$rate, lags = 0, type = "nc")
adfTest(dt$rate, lags = 1, type = "nc")
adfTest(dt$rate, lags = 2, type = "nc")

fit1 <- arima(dt$rate, order=c(3,0,0))
summary(fit1)
coeftest(fit1)

fit2 <- arima(dt$rate, order=c(3,0,0), include.mean = F)
summary(fit2)
coeftest(fit2)




### 잔차분석 

dt[, resid := as.numeric(resid(fit2))]

p3 <- ggplot(dt, aes(t, resid)) + 
  geom_line(col='steelblue') +
  geom_hline(yintercept = 0, col='grey', lty=2)+
  xlab("")+ylab('')+ ggtitle('Time series plot of Residual')+
  theme_bw()

p1 <- ggAcf(dt$resid) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_blank())

p2 <- ggPacf(dt$resid) + 
  theme_bw() +ylim(-1,1) +
  theme(plot.title = element_blank())

grid.arrange(p3, p1, p2, nrow = 2,
             layout_matrix = rbind(c(1,1),
                                   c(2,3)))


## 정규성 검정
shapiro.test(dt$resid)  ##H0 : normal distribution

ggplot(dt, aes(t, resid)) + 
  geom_line(col='steelblue') +
  geom_hline(yintercept = 0, col='grey', lty=2)+
  xlab("")+ylab('')+ ggtitle('Time serise plot of Residual')+
  theme_bw()

qqnorm(dt$resid)
qqline(dt$resid, col = 2)

LjungBox(fit2, lags=seq(6,24,6))# 잔차의 포트맨토 검정 

jarque.bera.test(dt$resid)  ##JB test H0: normal

sarima(dt$rate, 3, 0, 0) 
sarima.for(dt$rate, 12, 3, 0, 0)


#################################################
##### EX 8.7 : 과대적합 
#################################################

z <- scan("eg8_7.txt")

dt <- data.table( t = 1:length(z),
                  z=z)

p3 <- ggplot(dt, aes(t, z)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ ggtitle('Time serise plot')+
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

## 모형적합 AR(1)
fit <- arima(dt$z,order=c(1,0,0))
summary(fit)
coeftest(fit)


### 잔차분석 
dt[, resid := as.numeric(resid(fit))]

p3 <- ggplot(dt, aes(t, resid)) + 
  geom_line(col='steelblue') +
  geom_hline(yintercept = 0, col='grey', lty=2)+
  xlab("")+ylab('')+ ggtitle('Time series plot of Residual')+
  theme_bw()

p1 <- ggAcf(dt$resid) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_blank())

p2 <- ggPacf(dt$resid) + 
  theme_bw() +ylim(-1,1) +
  theme(plot.title = element_blank())

grid.arrange(p3, p1, p2, nrow = 2,
             layout_matrix = rbind(c(1,1),
                                   c(2,3)))

## 정규성 검정
shapiro.test(dt$resid)  ##H0 : normal distribution
jarque.bera.test(dt$resid)  ##JB test H0: normal

ggplot(dt, aes(t, resid)) + 
  geom_line(col='steelblue') +
  geom_hline(yintercept = 0, col='grey', lty=2)+
  xlab("")+ylab('')+ ggtitle('Time serise plot of Residual')+
  theme_bw()

qqnorm(dt$resid)
qqline(dt$resid, col = 2)

LjungBox(fit, lags=seq(6,24,3))# 잔차의 포트맨토 검정 

sarima(dt$z, 1,0,0) 
sarima.for(dt$z, 25, 1,0,0) 

## 모형적합 AR(2), ARMA(1,1)

fit1 <- sarima(dt$z, 1,0,0)    
fit2 <- sarima(dt$z, 2,0,0)    
fit3 <- sarima(dt$z, 1,0,1) 

fit1$ttable
fit2$ttable
fit3$ttable

fit1$AIC
fit2$AIC
fit3$AIC

fit1$fit$sigma2
fit2$fit$sigma2
fit3$fit$sigma2


#################################################
##### EX 8.8 : 단위근 검정
#################################################r

z <- scan("elecstock.txt")

dt <- data.table( t = 1:length(z),
                  z=z)

p3 <- ggplot(dt, aes(t, z)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ ggtitle('Time series plot')+
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

LjungBox(dt$z, lags=seq(6,24,6)) 


## 단위근 검정 H0 : phi=1
adfTest(dt$z, lags = 0, type = "c")
adfTest(dt$z, lags = 1, type = "c")
adfTest(dt$z, lags = 2, type = "c") 

# function adf.test를 이용할 수도 있음
adf.test(dt$z)    # ADF 검정
pp.test(dt$z)     # PP 검정

## 차분
dt[, diff_z := c(0, diff(z))]

p3 <- ggplot(dt, aes(t, diff_z)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ ggtitle('Time series plot')+
  theme_bw()

p1 <- ggAcf(dt$diff_z) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_blank())

p2 <- ggPacf(dt$diff_z) + 
  theme_bw() +ylim(-1,1) +
  theme(plot.title = element_blank())

grid.arrange(p3, p1, p2, nrow = 2,
             layout_matrix = rbind(c(1,1),
                                   c(2,3)))

LjungBox(dt$diff_z, lags=seq(6,24,6))

t.test(dt$diff_z)


## 모형적합 ARIMA(0,1,0)
fit <- arima(dt$z,order=c(0,1,0))
summary(fit)
coeftest(fit)

summary(fit)


### 잔차분석 
dt[, resid := as.numeric(resid(fit))]

p3 <- ggplot(dt, aes(t, resid)) + 
  geom_line(col='steelblue') +
  geom_hline(yintercept = 0, col='grey', lty=2)+
  xlab("")+ylab('')+ ggtitle('Time series plot of Residual')+
  theme_bw()

p1 <- ggAcf(dt$resid) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_blank())

p2 <- ggPacf(dt$resid) + 
  theme_bw() +ylim(-1,1) +
  theme(plot.title = element_blank())

grid.arrange(p3, p1, p2, nrow = 2,
             layout_matrix = rbind(c(1,1),
                                   c(2,3)))


## 정규성 검정
shapiro.test(dt$resid)  ##H0 : normal distribution
jarque.bera.test(dt$resid)  ##JB test H0: normal

ggplot(dt, aes(t, resid)) + 
  geom_line(col='steelblue') +
  geom_hline(yintercept = 0, col='grey', lty=2)+
  xlab("")+ylab('')+ ggtitle('Time series plot of Residual')+
  theme_bw()

qqnorm(dt$resid)
qqline(dt$resid, col = 2)

LjungBox(fit, lags=seq(6,24,3))# 잔차의 포트맨토 검정 




#################################################
##### EX 8.9 : female
#################################################

z <- scan("female.txt")

dt <- data.table( t = seq.Date(as.Date("1983-01-01"), 
                               by='month', 
                               length.out=length(z)),
                  z=z)

p3 <- ggplot(dt, aes(t, z)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ ggtitle('Time series plot')+
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

LjungBox(dt$z, lags=seq(6,24,6)) 


## 단위근 검정 
adfTest(dt$z, lags = 0, type = "ct")
adfTest(dt$z, lags = 1, type = "ct")
adfTest(dt$z, lags = 2, type = "ct") 

# function adf.test를 이용할 수도 있음
library(tseries)   # library for function adf.test & pp.test
adf.test(dt$z)    # ADF 검정
pp.test(dt$z)     # PP 검정


## 선형모형적합 
fit1 <- lm(z~t, data=dt)  # 선형모형 적합
summary(fit1)
coefficients(fit1)

dt[, resid_lm := residuals(fit1)]

p3 <- ggplot(dt, aes(t, resid_lm)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ ggtitle('Time series plot')+
  theme_bw()

p1 <- ggAcf(dt$resid_lm) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_blank())

p2 <- ggPacf(dt$resid_lm) + 
  theme_bw() +ylim(-1,1) +
  theme(plot.title = element_blank())

grid.arrange(p3, p1, p2, nrow = 2,
             layout_matrix = rbind(c(1,1),
                                   c(2,3)))

LjungBox(dt$resid_lm, lags=seq(6,24,6)) 

adfTest(dt$resid_lm, lags = 0, type = "nc")
adfTest(dt$resid_lm, lags = 1, type = "nc")
adfTest(dt$resid_lm, lags = 2, type = "nc") 

jarque.bera.test(dt$resid_lm)

###### lm (선형추세 모형) -> 잔차검정
###### et ~ AR(1)

fit_resid <- arima(dt$resid_lm, order=c(1,0,0), include.mean = F)
fit_resid
coeftest(fit_resid)

#### Zt = alpha + beta*t + et
#### et = phi*e(t-1)+ut
#### ut ~ WN(0, sigma^2)

#### Zt = T + I
#### Zt- T = It ~ stationary process

#######################################
## 차분
dt[, diff_z := c(0, diff(z))]

p3 <- ggplot(dt, aes(t, diff_z)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ ggtitle('Time series plot')+
  theme_bw()

p1 <- ggAcf(dt$diff_z) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_blank())

p2 <- ggPacf(dt$diff_z) + 
  theme_bw() +ylim(-1,1) +
  theme(plot.title = element_blank())

grid.arrange(p3, p1, p2, nrow = 2,
             layout_matrix = rbind(c(1,1),
                                   c(2,3)))

LjungBox(dt$diff_z, lags=seq(6,24,6))
t.test(dt$diff_z)

## 모형적합 ARIMA(0,1,0)
fit <- arima(dt$z,order=c(0,1,0),
             include.mean = TRUE)
summary(fit)



### 잔차분석 
dt[, resid := as.numeric(resid(fit))]

p3 <- ggplot(dt, aes(t, resid)) + 
  geom_line(col='steelblue') +
  geom_hline(yintercept = 0, col='grey', lty=2)+
  xlab("")+ylab('')+ ggtitle('Time series plot of Residual')+
  theme_bw()

p1 <- ggAcf(dt$resid) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_blank())

p2 <- ggPacf(dt$resid) + 
  theme_bw() +ylim(-1,1) +
  theme(plot.title = element_blank())

grid.arrange(p3, p1, p2, nrow = 2,
             layout_matrix = rbind(c(1,1),
                                   c(2,3)))


## 정규성 검정
shapiro.test(dt$resid)  ##H0 : normal distribution
jarque.bera.test(dt$resid)  ##JB test H0: normal

ggplot(dt, aes(t, resid)) + 
  geom_line(col='steelblue') +
  geom_hline(yintercept = 0, col='grey', lty=2)+
  xlab("")+ylab('')+ ggtitle('Time series plot of Residual')+
  theme_bw()

qqnorm(dt$resid)
qqline(dt$resid, col = 2)

LjungBox(fit, lags=seq(6,24,3))# 잔차의 포트맨토 검정 
