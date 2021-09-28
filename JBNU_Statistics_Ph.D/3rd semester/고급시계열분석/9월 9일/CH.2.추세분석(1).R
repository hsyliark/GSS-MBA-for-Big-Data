
##############  package
library(lmtest)
library(data.table)
library(ggplot2)


##############  경로지정
setwd("C:\\R-Project\\DAT\\Time Series Data")



####################################################
##### 국내 총인구 
####################################################

z <- scan("D:/Workplace/GSS_MBA_for_Big_Data/JBNU_Statistics_Ph.D/3rd semester/고급시계열분석/프로그램자료모음/제5판 시계열분석 data/population.txt") 
pop <- round(z/10000)
plot(1:length(z), pop, type='l')


tmp.data <- data.table(
  day = seq.Date(as.Date("1960-01-01"), 
                 by='year', length.out=length(z)),
  pop = round(z/10000),
  t = 1:length(z),
  t2 = (1:length(z))^2
)

ggplot(tmp.data, aes(day, pop)) + 
  geom_line(col='skyblue') +
  geom_point(col='steelblue')+
  xlab("")+ylab("")+
  theme_bw()


######## 1차 선형 추세

m1 <- lm(pop~t, data=tmp.data)
summary(m1)
dwtest(m1)

tmp.data[, fitted_pop := fitted(m1)]
# tmp.data$fitted_pop <- fitted(m1)

ggplot(tmp.data, aes(day, pop)) + 
  geom_line( col='skyblue') +
  geom_point(col='steelblue')+
  geom_line(aes(day,fitted_pop ), lty=2, col='red', size=1)+
  xlab("")+ylab("")+
  scale_x_date(date_labels = "%Y") +
  theme_bw()


tmp.data[, res := resid(m1)]

ggplot(tmp.data, aes(day, res)) + geom_line(col='skyblue') +
  geom_point(col='steelblue')+
  xlab("")+ylab("")+
  # scale_x_date(date_breaks = "1 year", date_labels = "%Y-%m") +
  geom_hline(yintercept = 0, col='grey', lty=2) +
  theme_bw()


######## 2차 선형 추세

m2 <- lm(pop~t+t2, data=tmp.data)
summary(m2)
dwtest(m2)

tmp.data[, fitted_pop_2 := as.numeric(fitted(m2))]

ggplot(tmp.data, aes(day, pop)) + 
  geom_line( col='skyblue') +
  geom_point(col='steelblue')+
  geom_line(aes(day,fitted_pop_2 ), lty=2, col='red', size=1)+
  xlab("")+ylab("")+
  scale_x_date(date_labels = "%Y") +
  theme_bw()


tmp.data[, res2 := resid(m2)]

ggplot(tmp.data, aes(day, res2)) + geom_line(col='skyblue') +
  geom_point(col='steelblue')+
  xlab("")+ylab("")+
  # scale_x_date(date_breaks = "1 year", date_labels = "%Y-%m") +
  geom_hline(yintercept = 0, col='grey', lty=2) +
  theme_bw()


########## 로그변환/2차 추세
######## 2차 선형 추세

m3 <- lm(log(pop)~t+t2, data=tmp.data)
summary(m3)
dwtest(m3)

tmp.data[, fitted_pop_3 := as.numeric(fitted(m3))]

ggplot(tmp.data, aes(day, log(pop))) + 
  geom_line( col='skyblue') +
  geom_point(col='steelblue')+
  geom_line(aes(day,fitted_pop_3 ), lty=2, col='red', size=1)+
  xlab("")+ylab("")+
  scale_x_date(date_labels = "%Y") +
  theme_bw()


tmp.data[, res3 := resid(m3)]

ggplot(tmp.data, aes(day, res3)) + geom_line(col='skyblue') +
  geom_point(col='steelblue')+
  xlab("")+ylab("")+
  # scale_x_date(date_breaks = "1 year", date_labels = "%Y-%m") +
  geom_hline(yintercept = 0, col='grey', lty=2) +
  theme_bw()





####################################################
##### 주기성분을 갖는 시계열
####################################################


##########  삼각함수

x <- seq(0,24,0.01 )

s<-12

par(mfrow=c(3,1))
plot(x, sin(2*pi*x/s), type='l', col='steelblue', lwd=2,
     xlab="", ylab="", main=paste0('sin :: ', "frequency=", s))
abline(h=0, lty=2)
abline(v= seq(0, 24, by=s), lty=2)

plot(x, cos(2*pi*x/s), type='l', col='steelblue', lwd=2,
     xlab="", ylab="", main=paste0('cos :: ', "frequency=", s))
abline(h=0, lty=2)
abline(v= seq(0, 24, by=s), lty=2)

plot(x, sin(2*pi*x/s)+cos(2*pi*x/s), type='l', col='steelblue', lwd=2,
     xlab="", ylab="", main=paste0('sin+cos :: ', "frequency=", s))
abline(h=0, lty=2)
abline(v= seq(0, 24, by=s), lty=2)
abline(v= seq(1.5, 24, by=s), lty=2)

plot(x, 1.5* sin(2*pi*x/s)+0.7* cos(2*pi*x/s), type='l', col='steelblue', lwd=2,
     xlab="", ylab="", main=paste0('sin+cos :: ', "frequency=", s))
abline(h=0, lty=2)
abline(v= seq(0, 24, by=s), lty=2)
abline(v= seq(1.5, 24, by=s), lty=2)

par(mfrow=c(3,1))


s<-12
plot(x, sin(2*pi*x/s), type='l', col='steelblue', lwd=2,
     xlab="", ylab="", main=paste0('sin :: ', "frequency=", s))
abline(h=0, lty=2)
abline(v= seq(0, 24, by=s), lty=2)

s<-6
plot(x, sin(2*pi*x/s), type='l', col='steelblue', lwd=2,
     xlab="", ylab="", main=paste0('sin :: ', "frequency=", s))
abline(h=0, lty=2)
abline(v= seq(0, 24, by=s), lty=2)
 
plot(x, sin(2*pi*x/12)+sin(2*pi*x/6), type='l', col='steelblue', lwd=2,
     xlab="", ylab="", main=paste0('sin :: ', "frequency=", s))
abline(h=0, lty=2)
abline(v= seq(0, 24, by=12), lty=2)

plot(x, 2*sin(2*pi*x/12)+0.8*sin(2*pi*x/6), type='l', col='steelblue', lwd=2,
     xlab="", ylab="", main=paste0('sin :: ', "frequency=", s))
abline(h=0, lty=2)
abline(v= seq(0, 24, by=12), lty=2)


par(mfrow=c(4,1))
s<-12
plot(x, sin(2*pi*x/s), type='l', col='steelblue', lwd=2,
     xlab="", ylab="", main=paste0('sin :: ', "frequency=", s))
abline(h=0, lty=2)
abline(v= seq(0, 24, by=s), lty=2)

s<-6
plot(x, sin(2*pi*x/s), type='l', col='steelblue', lwd=2,
     xlab="", ylab="", main=paste0('sin :: ', "frequency=", s))
abline(h=0, lty=2)
abline(v= seq(0, 24, by=s), lty=2)

s<-3
plot(x, sin(2*pi*x/s), type='l', col='steelblue', lwd=2,
     xlab="", ylab="", main=paste0('sin :: ', "frequency=", s))
abline(h=0, lty=2)
abline(v= seq(0, 24, by=s), lty=2)

plot(x, sin(2*pi*x/12)+sin(2*pi*x/6)+sin(2*pi*x/3), type='l', col='steelblue', lwd=2,
     xlab="", ylab="", main=paste0('sin + sin + sin :: ', "frequency=", 12))
abline(h=0, lty=2)
abline(v= seq(0, 24, by=12), lty=2)



y <- sin(2*pi*x/12)+sin(2*pi*x/6)+sin(2*pi*x/3)+
  cos(2*pi*x/12)+cos(2*pi*x/6)+cos(2*pi*x/3)

par(mfrow=c(2,1))

plot(x, sin(2*pi*x/12)+sin(2*pi*x/6)+sin(2*pi*x/3), type='l', col='steelblue', lwd=2,
     xlab="", ylab="", main=paste0('sin + sin + sin :: ', "frequency=", 12))
abline(h=0, lty=2)
abline(v= seq(0, 24, by=12), lty=2)

plot(x,y, type='l', col='steelblue', lwd=2,
     xlab="", ylab="", main="frequency=12")
abline(h=0, lty=2)
abline(v= seq(0, 24, by=12), lty=2)

par(mfrow=c(1,1))

###############################

n <- 100;  
t <- 1:n
a1 <- -0.8;  
a2 <- 1.4
phi1 <- pi/3;  
phi2 <- 3*pi/4
first <- a1*sin(pi*t/6+phi1)    # 첫 번째 주기성분
second <- a2*sin(pi*t/3+phi2)  # 두 번째 주기성분

dt <- data.table(t=t,
                 first = first,
                 second=second,
                 z = first+second)
# install.packages('gridExtra')
library(gridExtra)
p1 <- ggplot(dt, aes(t,first)) + geom_line(col='skyblue', size=1) +
  geom_point(col='steelblue', size=1)+
  ylim(-2.5,2)+xlab("")+
  scale_x_continuous(breaks = seq(1,100, by = 12))+
  theme_bw()

p2 <- ggplot(dt, aes(t,second)) + geom_line(col='skyblue', size=1) +
  geom_point(col='steelblue', size=1)+
  ylim(-2.5,2)+xlab("")+
  scale_x_continuous(breaks = seq(1,100, by = 12))+
  theme_bw()

p3 <- ggplot(dt, aes(t,z)) + geom_line(col='skyblue', size=1) +
  geom_point(col='steelblue', size=1)+
  scale_x_continuous(breaks = seq(1,100, by = 12))+
  ylim(-2.5,2)+xlab("")+
  theme_bw()

grid.arrange(p1, p2, p3, nrow = 3)

####################################################
##### 백화점 매출액 - 지시함수
####################################################

z <-scan("D:/Workplace/GSS_MBA_for_Big_Data/JBNU_Statistics_Ph.D/3rd semester/고급시계열분석/프로그램자료모음/제5판 시계열분석 data/depart.txt") 

dep <- ts(z, frequency=12, start=c(1984,1))

tmp.data <- data.table(
  day = seq.Date(as.Date("1984-01-01"), 
                 by='month', length.out=length(z)),
  z=z  
)
tmp.data[, lndep := log(z)]
tmp.data[, y := as.factor(as.integer(cycle(dep)))]
tmp.data[, trend := 1:length(z)]

ggplot(tmp.data, aes(day, z)) + geom_line(col='skyblue') +
  geom_point(col='steelblue')+
  xlab("")+ylab("")+
  scale_x_date(date_breaks = "1 year", date_labels = "%Y-%m") +
  theme_bw()

ggplot(tmp.data, aes(day, lndep)) + geom_line(col='skyblue') +
  geom_point(col='steelblue')+
  xlab("")+ylab("")+
  scale_x_date(date_breaks = "1 year", date_labels = "%Y-%m") +
  theme_bw()


reg <- lm(lndep ~ 0+trend+y, data=tmp.data )   
dwtest(reg)
summary(reg)

tmp.data[,res := resid(reg)]
tmp.data[, fitted_lndep := fitted(reg)]

ggplot(tmp.data, aes(day, lndep)) + geom_line(col='skyblue', lwd=1) +
  geom_line(aes(day, fitted_lndep), col='orange', lwd=1) +
  xlab("")+ylab("")+
  scale_x_date(date_breaks = "1 year", date_labels = "%Y-%m") +
  theme_bw()

ggplot(tmp.data, aes(day, res)) + geom_line(col='skyblue') +
  geom_point(col='steelblue')+
  xlab("")+ylab("")+
  scale_x_date(date_breaks = "1 year", date_labels = "%Y-%m") +
  geom_hline(yintercept = 0, col='grey', lty=2) +
  theme_bw()



####################################################
##### 백화점 매출액 - 삼각함수
####################################################

tmp.data_sub <- tmp.data[,.(lndep, trend)]

tmp.data_sub <- cbind(tmp.data_sub,
                      tmp.data_sub[, lapply(as.list(1:5), 
                                            function(i) sin(2*pi*i/12*trend))])

names(tmp.data_sub)[-(1:2)]<-paste("sin", c(12,6,4,3,2.4), sep="_")

tmp.data_sub <- cbind(tmp.data_sub,
                      tmp.data_sub[, lapply(as.list(1:5), 
                                            function(i) cos(2*pi*i/12*trend))])

names(tmp.data_sub)[-(1:7)]<-paste("cos", c(12,6,4,3,2.4), sep="_")

reg_2 <- lm(lndep ~., data=tmp.data_sub)
summary(reg_2)
dwtest(reg_2)
dwtest(reg_2, alternative = "less")

tmp.data_sub[, res := resid(reg_2)]
tmp.data_sub[, day:=tmp.data$day]
tmp.data_sub[, fitted_lndep := fitted(reg_2)]

ggplot(tmp.data_sub, aes(day, lndep)) + geom_line(col='skyblue', lwd=1) +
  geom_line(aes(day, fitted_lndep), col='orange', lwd=1) +
  xlab("")+ylab("")+
  scale_x_date(date_breaks = "1 year", date_labels = "%Y-%m") +
  theme_bw()

ggplot(tmp.data_sub, aes(day, res)) + geom_line(col='skyblue') +
  geom_point(col='steelblue')+
  xlab("")+ylab("")+
  scale_x_date(date_breaks = "1 year", date_labels = "%Y-%m") +
  geom_hline(yintercept = 0, col='grey', lty=2) +
  theme_bw()




####################################################
##### 국내 총인구 
####################################################
z<- scan("catv.txt")
k = 70000000
t <- 1:length(z)

tmp.data <- data.table(
  day = seq.Date(as.Date("1970-01-01"), by='year', length.out=length(z)),
  t = 1: length(z),
  catv=z,
  lncatv = log(k/z -1)
)

ggplot(tmp.data, aes(day, catv)) + geom_line(col='steelblue')+
  geom_point(col='steelblue')+
  xlab("")+ylab("catv")+
  # scale_x_date(date_breaks = "5 year", date_labels = "%Y") +
  theme_bw()
ggplot(tmp.data, aes(day, lncatv)) + geom_line(col='steelblue')+
  geom_point(col='steelblue')+
  xlab("")+ylab("ln(k/catv-1)")+
  # scale_x_date(date_breaks = "5 year", date_labels = "%Y") +
  theme_bw()

fit <- lm(lncatv ~t, tmp.data)  
summary(fit)

dwtest(fit)

tmp.data[, fiited_catv := k/(exp(fitted(fit))+1)]
tmp.data[, resid := catv-fiited_catv]

ggplot(tmp.data, aes(day, catv)) + geom_line(col='steelblue', lwd=1)+
  geom_line( aes(day, fiited_catv), col='orange', lty=2, lwd=1)+
  xlab("")+ylab("catv")+
  theme_bw()


ggplot(tmp.data, aes(t, resid)) + 
  geom_line(col='steelblue', lwd=1) + 
  geom_point(col='steelblue')+
  geom_hline(yintercept = 0, lty=2)+
  xlab("")+ylab("residual")+
  theme_bw()

