##############  package
library(data.table)
library(ggplot2)
library(forecast)
library(gridExtra)

##############  경로지정
setwd("C:\\R-Project\\DAT\\Time Series Data")

### 분산이 일정하지 않은 경우

data("AirPassengers")
z <- as.numeric(AirPassengers)

plot(z, type='l')


tmp.dat <- data.table(
  t=seq.Date(as.Date("1949-01-01"), 
             by='month', 
             length.out=144),
  original=z,
  log = log(z),
  sqrt = sqrt(z),
  Boxcox_0.15 = BoxCox(z,lambda= BoxCox.lambda(z, method='loglik')) 
)

res <- BoxCox.lambda(z, method='loglik')
res

melt.tmp <- melt(tmp.dat, id=1)
melt.tmp

ggplot(melt.tmp, aes(t,value)) + 
  geom_line(col='steelblue', lwd=1) +
  xlab("")+ylab("")+
  ggtitle("Time Series")+
  theme_bw()+
  facet_wrap(.~variable, nrow = 2, scales = "free")


### 확률보행과정 정상화
## Zt = Z(t-1)+et = et + e(t-1) + e(t-2) + ... + e1

e <- rnorm(100)
plot(e, type='l')
plot(cumsum(e), type='l')

tmp.dat <- data.table(
  t=1:100,
  z=cumsum(e)
)
tmp.dat[, dz := c(0, diff(z))]

######### 확률보행과정의 ACF/PACF

## AR(1), phi=1 => non-stationary ===> ARIMA(0,1,0)

p1 <- ggplot(tmp.dat, aes(t, z)) + 
  geom_line(col='steelblue', lwd=1) +
  xlab("")+ylab("")+
  ggtitle("Time Series plot : Random Walk")+
  # scale_x_continuous(breaks = seq(0,30, by = 6))+
  # geom_hline(yintercept = 0, lty=2, col='grey')+
  theme_bw()

p2 <- ggAcf(tmp.dat$z) + 
  theme_bw() + ylim(-1,1) +
  scale_x_continuous(breaks = seq(1,20, by = 2))+
  ggtitle("SACF")

p3 <- ggPacf(tmp.dat$z) + 
  theme_bw() +ylim(-1,1) +
  scale_x_continuous(breaks = seq(1,20, by = 2))+
  ggtitle("SPACF")

grid.arrange(p1, p2, p3, nrow = 2,
             layout_matrix = rbind(c(1, 1),
                                   c(2, 3)))


######### 차분한 확률보행과정의 ACF/PACF
p1 <- ggplot(tmp.dat, aes(t, dz)) + 
  geom_line(col='steelblue', lwd=1) +
  xlab("")+ylab("")+
  ggtitle("Time Series plot : diff(Random Walk)")+
  # scale_x_continuous(breaks = seq(0,30, by = 6))+
  geom_hline(yintercept = 0, lty=2, col='grey')+
  theme_bw()

p2 <- ggAcf(tmp.dat$dz) + 
  theme_bw() + ylim(-1,1) +
  scale_x_continuous(breaks = seq(1,20, by = 2))+
  ggtitle("SACF")

p3 <- ggPacf(tmp.dat$dz) + 
  theme_bw() +ylim(-1,1) +
  scale_x_continuous(breaks = seq(1,20, by = 2))+
  ggtitle("SPACF")

grid.arrange(p1, p2, p3, nrow = 2,
             layout_matrix = rbind(c(1, 1),
                                   c(2, 3)))





####################################
z <-scan("D:/Workplace/GSS_MBA_for_Big_Data/JBNU_Statistics_Ph.D/3rd semester/고급시계열분석/프로그램자료모음/제5판 시계열분석 data/depart.txt") 

tmp.dat <- data.table(
  t=seq.Date(as.Date("1984-01-01"), 
             by='month', 
             length.out=60),
  z=z
)

tmp.dat[, lnz:= log(z)]
tmp.dat[, diff_z := c(0, diff(lnz))]

melt_dt <- melt(tmp.dat, id=1)

ggplot(melt_dt, aes(t, value)) + 
  geom_line(col='steelblue', lwd=1) +
  xlab("")+ylab("")+
  # ggtitle("Time Series plot : diff(Random Walk)")+
  # scale_x_continuous(breaks = seq(0,30, by = 6))+
  # geom_hline(yintercept = 0, lty=2, col='grey')+
  theme_bw()+
  facet_wrap(~variable, scales = 'free_y')

acf(tmp.dat$diff_z, lag.max = 36)
pacf(tmp.dat$diff_z, lag.max = 36)

####### 차분한 데이터의 ACF/PACF
p2 <- ggAcf(tmp.dat[,diff_z]) + 
  theme_bw() + ylim(-1,1) +
  ggtitle("SACF")

p3 <- ggPacf(tmp.dat[,diff_z ]) + 
  theme_bw() +ylim(-1,1) +
  ggtitle("SPACF")

grid.arrange(p2, p3, nrow = 2)

###### 계절차분 
tmp.dat[, diff_s_z := c(rep(0,13), diff(diff(lnz), lag=12))]

p1 <- ggplot(tmp.dat[-(1:13)], aes(t, diff_s_z)) + 
  geom_line(col='steelblue', lwd=1) +
  xlab("")+ylab("")+
  # ggtitle("Time Series plot : diff(Random Walk)")+
  # scale_x_continuous(breaks = seq(0,30, by = 6))+
  # geom_hline(yintercept = 0, lty=2, col='grey')+
  theme_bw()


p2 <- ggAcf(tmp.dat[-(1:13),diff_s_z ]) + 
  theme_bw() + ylim(-1,1) +
  ggtitle("SACF")

p3 <- ggPacf(tmp.dat[-(1:13),diff_s_z ]) + 
  theme_bw() +ylim(-1,1) +
  ggtitle("SPACF")

grid.arrange(p1, p2, p3, nrow = 2,
             layout_matrix = rbind(c(1, 1),
                                   c(2, 3)))

acf(tmp.dat[-(1:13),diff_s_z ], lag.max = 36)



####### ARIMA(1,1,0)
n <- 10000

dt <- data.table(t = 1:n,
                 x = as.numeric(arima.sim(n=n, list(order=c(1,0,0), ar=0.8))),
                 z = as.numeric(arima.sim(n=n-1, list(order=c(1,1,0), ar=0.8))))

p1 <- ggplot(dt, aes(t, x)) + 
  geom_line(col='steelblue', lwd=1) +
  xlab("")+ylab("")+
  ggtitle("AR(1)")+
  theme_bw()


p2 <- ggplot(dt, aes(t, z)) + 
  geom_line(col='steelblue', lwd=1) +
  xlab("")+ylab("")+
  ggtitle("ARIMA(1,1,0)")+
  theme(plot.title = element_text(size = 10 ))+
  theme_bw()

p3 <- ggAcf(dt[,x]) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_text(size = 10 ))+
  theme(text = element_text(size=10))+
  ggtitle("SACF of AR(1)")

p4 <- ggAcf(dt[,z]) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_text(size = 10 ))+
  theme(text = element_text(size=10))+
  ggtitle("SACF of ARIMA(1,1,0)")

p5 <- ggPacf(dt[,x]) + 
  theme_bw() + ylim(-1,1) +  theme(plot.title = element_text(size = 10 ))+
  theme(text = element_text(size=10))+
  ggtitle("SPACF of AR(1)")


p6 <- ggPacf(dt[,z]) + 
  theme_bw() + ylim(-1,1) +  theme(plot.title = element_text(size = 10 ))+
  theme(text = element_text(size=10))+
  ggtitle("SPACF of ARIMA(1,1,0)")


grid.arrange(p1, p3, p5, p2, p4, p6, nrow = 2)


####### ARIMA(0,1,1)


dt <- data.table(t = 1:10000,
                 x = as.numeric(arima.sim(n=10000, list(order=c(0,0,1), ma=-0.8))),
                 z = as.numeric(arima.sim(n=9999, list(order=c(0,1,1), ma=-0.8))))

p1 <- ggplot(dt, aes(t, x)) + 
  geom_line(col='steelblue', lwd=1) +
  xlab("")+ylab("")+
  
  ggtitle("MA(1)")+
  theme(plot.title = element_text(size = 10 ))+
  theme_bw()


p2 <- ggplot(dt, aes(t, z)) + 
  geom_line(col='steelblue', lwd=1) +
  xlab("")+ylab("")+
  ggtitle("ARIMA(0,1,1)")+
  theme(plot.title = element_text(size = 10 ))+
  theme_bw()

p3 <- ggAcf(dt[,x]) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_text(size = 10 ))+
  theme(text = element_text(size=10))+
  ggtitle("SACF of MA(1)")

p4 <- ggAcf(dt[,z]) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_text(size = 10 ))+
  theme(text = element_text(size=10))+
  ggtitle("SACF of ARIMA(0,1,1)")

p5 <- ggPacf(dt[,x]) + 
  theme_bw() + ylim(-1,1) +  theme(plot.title = element_text(size = 10 ))+
  theme(text = element_text(size=10))+
  ggtitle("SPACF of MA(1)")


p6 <- ggPacf(dt[,z]) + 
  theme_bw() + ylim(-1,1) +  theme(plot.title = element_text(size = 10 ))+
  theme(text = element_text(size=10))+
  ggtitle("SPACF of ARIMA(0,1,1)")


grid.arrange(p1, p3, p5, p2,p4, p6, nrow = 2)



###################################################################
#######   과대차분
###################################################################

## ARIMA(0,1,1)
dt <- data.table(t = 1:100,
                 z = as.numeric(arima.sim(n=99, list(order=c(0,1,1), ma=0.5))))

p1 <- ggAcf(dt[,z]) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_text(size = 10 ))+
  theme(text = element_text(size=10))+
  ggtitle("SACF of ARIMA(0,1,1)")


p2 <- ggAcf(dt[,diff(z)]) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_text(size = 10 ))+
  theme(text = element_text(size=10))+
  ggtitle("SACF of MA(1)")


p3 <- ggAcf(dt[,diff(diff(z))]) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_text(size = 10 ))+
  theme(text = element_text(size=10))+
  ggtitle("SACF of diff(MA(1))")

grid.arrange(p1, p2, p3, nrow = 1)



