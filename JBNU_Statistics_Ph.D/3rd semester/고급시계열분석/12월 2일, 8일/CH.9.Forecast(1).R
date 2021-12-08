######################################
###### Time series analysis
###### CH 9 Forecast
###### 202102
######################################

##############  package

library(astsa)
library(data.table)
library(ggplot2)
library(gridExtra)
library(forecast)
library(portes)  ##LjungBox
library(fUnitRoots)  # library for function adfTest
library(tseries)  ##JB test
setwd("C:\\R-Project\\DAT\\Time Series Data")

#################################################
##### EX 9.1 AR(1)
#################################################

z <- scan("eg8_7.txt")

dt <- data.table( t = 1:length(z),z = z)
## 시계열그림 / ACF PACF

p3 <- ggplot(dt, aes(t, z)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ ggtitle('Simulated AR(1) Process')+
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
fit <- arima(dt$z, order=c(1,0,0), method='CSS')

forecast_fit <- forecast(fit, 25)

head(forecast_fit)
head(forecast_fit$mean)
head(forecast_fit$lower)
head(forecast_fit$upper)

dt_tmp <- data.table(l = 1:25,
                     hat_Zn_l = as.numeric(forecast_fit$mean),
                     upper_95 = as.numeric(forecast_fit$upper[,2]),
                     lower_95 = as.numeric(forecast_fit$lower[,2]))

ggplot(dt_tmp, aes(l, hat_Zn_l)) + 
  geom_line(col='steelblue') +
  geom_line(aes(l, upper_95), col='red')+
  geom_line(aes(l, lower_95), col='red')+
  xlab("")+ylab('')+
  theme_bw()


sarima_fit <- sarima.for(z, 25,1,0,0)

sarima_fit$pred
sarima_fit$se

head(sarima_fit$pred + qnorm(0.975)*sarima_fit$se)  ##95% 신뢰구간 상한
head(sarima_fit$pred - qnorm(0.975)*sarima_fit$se)  ##95% 신뢰구간 하한



#################################################
##### EX 9.5 IMA(1,1)
#################################################
z <- scan("eg9_5.txt")

dt <- data.table( t = 1:length(z),
                  z = z)
## 시계열그림 / ACF PACF

p3 <- ggplot(dt, aes(t, z)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+ ggtitle('')+
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

adfTest(dt$z, lags = 0, type = "c")
adfTest(dt$z, lags = 1, type = "c")
adfTest(dt$z, lags = 2, type = "c")

dt[, diff_z := c(0, diff(z))]


p3 <- ggplot(dt, aes(t, diff_z)) + 
  geom_line(col='steelblue') +
  xlab("")+ylab('')+
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

adfTest(dt$diff_z, lags = 0, type = "nc")
adfTest(dt$diff_z, lags = 1, type = "nc")
adfTest(dt$diff_z, lags = 2, type = "nc")

fit1 <- arima(dt$diff_z, order=c(0,0,1), include.mean = F)
fit1
##### Wt=(1-B)Zt, Wt=et-theta*e(t-1)  
#####      => (1-B)Zt=et-theta*e(t-1) ## ARIMA(0,1,1)


fit <- arima(dt$z, order=c(0,1,1), include.mean = F)
fit
##### (1-B)Zt=et-theta*e(t-1) ## ARIMA(0,1,1)

### 잔차분석   

dt[, resid := as.numeric(resid(fit))]

p3 <- ggplot(dt, aes(t, resid)) + 
  geom_line(col='steelblue') +
  geom_hline(yintercept = 0, col='grey', lty=2)+
  xlab("")+ylab('')+ ggtitle('Time serise plot of Residual')+
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

qqnorm(dt$resid)
qqline(dt$resid, col = 2)

LjungBox(fit, lags=seq(6,24,6))

jarque.bera.test(dt$resid)  ##JB test H0: normal

## 예측

sarima.for(z, 25, 0,1,1)
sarima.for(dt$diff_z, 25, 0,0,1)
