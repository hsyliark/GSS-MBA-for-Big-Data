#=========================================================================
# title           : 고급시계열
# subtitle        : 전력 수요 예측 자료의 시계열 자료 분석
# 
#=========================================================================

#=========================================================================
# 1. Pacakge load 
#=========================================================================
#install.packages("forecast")
#install.packages("tseries")
#install.packages("imputeTS")
#install.packages("zoo")
#install.packages("dplyr")

#패키지 호출
library(forecast);library(tseries);library(imputeTS);library(zoo);library(dplyr);


#=========================================================================
# 2. 자료 읽기(load) 및 전처리 (Pre-processing for analysis) 
#=========================================================================

setwd("C:\\R-Project\\DAT\\etc")

#자료 읽기
data <- read.delim("./전력수요_시계열자료.txt",header=T,stringsAsFactors = F)

#문자형 시간에 시간 속성 부여
data$time<-as.POSIXct(data$time)

#시계열 자료에서는 시간 간격 및 순서가 중요
#시간순 정렬을 하면서 빈 시간을 채움
f.time<-seq(min(data$time),max(data$time),by='15 mins')
summary(f.time)
i.dat<-data[match(f.time,data$time),]

#13개의 Time missing이 있는 것을 확인
length(which(is.na(i.dat$time)))
i.dat$time<-f.time

#요일 정보 추가
i.dat$week<-weekdays(i.dat$time)


#전체 값의 범위을 살펴봄
boxplot(i.dat$value,main="S사 전력 사용량 범위",ylab="Kwh")

#전체 기간의 패턴을 살펴봄
plot(x=as.POSIXct(i.dat$time),y=as.numeric(i.dat$value),main="S사 전력사용량 추이",xlab="시간",ylab="Kwh")
#전체 패턴을 살펴봄 Lower bound에 약간의 out-lier로 생각 됨.
#특별히 추세와 이분산은 보이지 않음

#이상치 제거
out.seq<-quantile(i.dat$value,probs=c(0.001),na.rm=T)
i.dat$value[which(i.dat$value<=out.seq[1] )]<-NA

#이상치 제거가 적절히 되었는지 확인함.
plot(x=as.POSIXct(i.dat$time),y=as.numeric(i.dat$value),main="S사 전력사용량 추이",xlab="시간",ylab="Kwh")

#NA interpolation
length(which(is.na(i.dat$value))) #31개의 missing이 있는 것을 확인
i.dat$value<-na.interpolation(i.dat$value,option="linear") #linear interpolation 수행
length(which(is.na(i.dat$value))) #모두 interploation 된 것을 확인

plot(x=as.POSIXct(i.dat$time),y=as.numeric(i.dat$value),main="S사 전력사용량 추이",xlab="시간",ylab="Kwh")

#전체로 보니 패턴 파악이 어려워 일주일 단위로 살펴봄
seqs = seq(min(f.time),length.out=28, by="7 days") # 전체 시간을 7일 간격으로 분리함

#=========================================================================
# 3. 적절한 분석을 위한 탐색적 자료 분석 EDA (Explanatory Data Analysis)
#=========================================================================

#1주일의 시작과 끝 시점을 생성
week.info<-rollapply(seqs,2,function(x){
  w.info<-data.frame(Start=as.POSIXct(x[1]),End=as.POSIXct(x[2]))
}) 

#자료를 1주일 단위로 분리함
week.list<-list()
for( i in 1:nrow(week.info)){
  week.list[[i]]<-which(i.dat[,1]>=week.info[i,1] & i.dat[,1]<week.info[i,2])  
}
week.dat<-lapply(week.list,function(y){return(i.dat[y,])})  

length(week.dat)
#총 26주 + 2일의 자료임

#26개의 plot의 상하한을 고정하기위해 전체값을 기준으로 상하한 생성
y.bound<-c(round(min(i.dat$value,na.rm=T)*0.9,-2)
           ,round(max(i.dat$value,na.rm=T)*1.1,-2))

#플랏 작성
lapply(week.dat,function(week){
  x=week[,1]
  y=week[,2]
  plot(x,y,xlab="Time",ylab="Kwh",axes=F,type='l',ylim=y.bound)
  title("S사 전력소비",line=1,cex=2)
  title(paste0(substr(min(x),1,10)," ~ ",substr(max(x),1,10)),line=0,cex=0.5)
  axis(1,as.POSIXct(paste0(seq(as.POSIXct(substr(min(x),1,10)),length.out=8,by="1 days")," 00:00:00")),c("토","일","월","화","수","목","금","토"))
  axis(2,label=round(seq(y.bound[1],y.bound[2],length.out=4),0),at=round(seq(y.bound[1],y.bound[2],length.out=4),0))
})

#자기 상관 계수 및 부분 자기 상관 계수 확인 
acf(data$value,lag.max = 200)
pacf(data$value,lag.max = 200)
#AR(1~2)*SAR(1)*SMA(1) K(주기)=96으로 생각됨

#27개의 주간 Plot을 살펴보니 휴일 정보가 중요한 요인으로 작용하는 것으로 생각 됨
#1일 안에 각 시점별로 값의 패턴이 있음
#평일/휴일 , 각시점별  moving average, exponentail smoothing을 수행하기 위한 전처리

#요일정보 추가함
i.dat$week<-weekdays(i.dat$time)

#공휴일 정보 추가
#광복절
idx1<-which(i.dat$time>=as.POSIXct("2017-08-15 00:00:00") & i.dat$time<as.POSIXct("2017-08-16 00:00:00"))
length(idx1) # 1일치 96개가 선택 되었는지 확인
#임시공휴일,개천절,추석
idx2<-which(i.dat$time>=as.POSIXct("2017-10-02 00:00:00") & i.dat$time<as.POSIXct("2017-10-07 00:00:00"))
length(idx2) # 5일치 5*96개가 선택 되었는지 확인
#한글날
idx3<-which(i.dat$time>=as.POSIXct("2017-10-09 00:00:00") & i.dat$time<as.POSIXct("2017-10-10 00:00:00"))
length(idx3) # 1일치 96개가 선택 되었는지 확인
#크리스마스
idx4<-which(i.dat$time>=as.POSIXct("2017-12-25 00:00:00") & i.dat$time<as.POSIXct("2017-12-26 00:00:00"))
length(idx4) # 1일치 96개가 선택 되었는지 확인

#공휴일 변수 생성
holiday<-rep(0,nrow(i.dat))
holiday[c(idx1,idx2,idx3,idx4)]<-1
i.dat$holiday<-holiday
head(i.dat)

#요일 정보의 평일 휴일 변환
i.dat$week[-which(i.dat$week=="일요일" | i.dat$week=="토요일")]<-0
i.dat$week[which(i.dat$week=="일요일" | i.dat$week=="토요일")]<-1

#통합 휴일 변수 생성
i.dat$holiday<-as.numeric(i.dat$week)+as.numeric(i.dat$holiday)
i.dat<-i.dat[,-3]
head(i.dat,200)


#============================================================
# 4 시계열 자료의 예측 
#============================================================
#평일과 휴일을 분리한 데이터에서 각 시각별 예측값을 산출하기 위한 자료 전처리 함수
Pre.pro<-function(dat,lookback,day,today,holi){
  #dat : [time][value][holiday] 로 구성된 데이터
  #lookback : 과거 몇 일 까지 예측에 사용할 것인가 ex) 5 (과거 평일 5일 혹은 과거 휴일 5일)
  #day :  몇 일 예측값을 산출 할 것인가 ex) 1 (내일 하루)
  #today : 오늘 날짜 (기준시점) 
  #holiday : 내일이 휴일인가? 0-평일 1-휴일
  f.tr.dat<-dat%>%filter(holiday==holi)
  t.end<-as.POSIXct(paste0(today," 23:59:59"))
  tr.dat<-f.tr.dat%>%filter(time<=t.end)
  tr.dat<-tr.dat[rev(seq(nrow(tr.dat),by=-1,length.out=96*lookback)),]
  week.dat<-matrix(tr.dat$value,ncol=lookback,nrow=96,byrow=F)
  
  #미래 자료
  p.start<-as.POSIXct(paste0(today," 23:59:59")) # 예측 분석 시작 시점
  te.dat<-dat%>%filter((as.POSIXct(p.start)<time))
  te.dat<-te.dat[c(1:(96*day)),]
  
  #tr.dat lookback 만큼의 훈련자료셋
  #tr.dat day 만큼의 테스트셋
  #week.dat row : 시각 column : lookback
  return(list(tr.dat=tr.dat,te.dat=te.dat,week.dat=week.dat))
}

#함수에 today로 사용할 벡터 생성
pred.period<-seq(as.Date('2017-12-22 00:00:00'),as.Date('2017-12-30 00:00:00'),by='1 days')

#함수에 holiday로 사용할 벡터 생성
pred.holi<-c(1,1,1,0,0,0,0,1,1)

#pred.today와 pred.holiday가 동시에 변하면서 함수실행
f.dat<-mapply(function(today,holi){Pre.pro(i.dat,5,1,today,holi)},today=pred.period,holi=pred.holi,SIMPLIFY = F)
length(f.dat) #9일치 자료 예측을 위한 자료 생성

#============================================================
# 4.1 이동평균 방법 수행 
#============================================================

pred.val<-lapply(f.dat,function(x){return(apply(x[[3]],1,mean,na.rm=T))})

#============================================================
# 4.2 지수 평활 방법 수행
#============================================================

pred.val.exp<-lapply(f.dat,function(x){return(apply(x[[3]],1,function(y){
  return(as.numeric(forecast(ses(y,order=lookback,alpha=0.7),h=1)$mean))
}))})


#============================================================
# 4.3 ARIMA*SARIMA를 이용한 예측
#============================================================

pred.val.arma<-lapply(f.dat,function(x){
  tseries <- ts(x[[1]][,2], frequency = 96)
  #adf.test(tseries, alternative="stationary")
  arma_fit<-auto.arima(tseries,trace=T)
  #tsdiag(par.auto)
  pred<- forecast(arma_fit, h = 96)
  return(pred)
})

#============================================================
# 4.4 예측값의 평가 
#============================================================

#총 9일치의 예측값을 합침
#이동평균
pred.val<-do.call(cbind,pred.val)
pred.val<-as.numeric(pred.val)

#지수평활
pred.val.exp<-do.call(cbind,pred.val.exp)
pred.val.exp<-as.numeric(pred.val.exp)

#ARIMA*SARIMA
pred.val.arma<-do.call(cbind,lapply(pred.val.arma,function(x){return(x$mean)}))
pred.val.arma<-as.numeric(pred.val.arma)

#총 9일치의 예측치와 비교 할 실제 값을 하나의 자료로 불러옴
true.val<-lapply(f.dat,function(x){return(x[[2]])})
true.val<-do.call(rbind,true.val)


#얼마나 예측이 잘되었는가를 MAPE를 통하여 측정
mape <- function(actual,pred){
  mape <- mean(abs((actual - pred)/actual))*100
  return (mape)
}
#결과는 Plot에 함께 작성

#============================================================
# 4.5 예측값과 실제값을 비교할 Plot 작성
#============================================================

par(mfrow=c(2,2),oma=c(1,1,4,1),font=4)
#이동평균 결과
plot(y=true.val$value,x=true.val$time,type='l',ylim=c(100,1000),col='gray90',lwd=2.5,ylab="Kwh",xlab="Time",main="Simple Moving Average")
lines(y=pred.val,x=true.val$time,lty=2,col='darkblue',lwd=1.5)
legend(x=c(true.val$time[1],true.val$time[200]),y=c(800,1000),c("True","Pred"),lty=c(1,2),lwd=2,col=c("gray90",'darkblue'))
title(paste0("MAPE : ",round(mape(true.val$value,pred.val),2)),cex=0.1,adj=0,line=0.2,font=1)

#지수평활 결과
plot(y=true.val$value,x=true.val$time,type='l',ylim=c(100,1000),col='gray90',lwd=2.5,ylab="Kwh",xlab="Time",main="Simple Exponential Smoothing(alpha=0.7)")
lines(y=pred.val.exp,x=true.val$time,lty=2,col='red',lwd=1.5)
legend(x=c(true.val$time[1],true.val$time[200]),y=c(800,1000),c("True","Pred"),lty=c(1,2),lwd=2,col=c("gray90",'red'))
title(paste0("MAPE : ",round(mape(true.val$value,pred.val.exp),2)),cex=0.1,adj=0,line=0.2,font=1)

#ARIMA * SARIMA 결과
plot(y=true.val$value,x=true.val$time,type='l',ylim=c(100,1000),col='gray90',lwd=2.5,ylab="Kwh",xlab="Time",main="ARIMA / SARIMA")  
lines(y=pred.val.arma,x=true.val$time,lty=2,col='green',lwd=1.5)  
legend(x=c(true.val$time[1],true.val$time[200]),y=c(800,1000),c("True","Pred"),lty=c(1,2),lwd=2,col=c("gray90",'green'))
title(paste0("MAPE : ",round(mape(true.val$value,pred.val.arma),2)),cex=0.1,adj=0,line=0.2,font=1)
  
#Graphic device 종료
#dev.off()






