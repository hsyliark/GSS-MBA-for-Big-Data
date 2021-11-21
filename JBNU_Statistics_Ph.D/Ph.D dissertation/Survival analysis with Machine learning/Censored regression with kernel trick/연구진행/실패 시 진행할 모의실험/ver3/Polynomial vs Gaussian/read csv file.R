dat.res5 <- read.csv("C:/Users/Hi/Desktop/ver3/Polynomial vs Gaussian/설명변수 5/중도절단 0/p5cen0.csv",sep=",",header=T)
dat.res5 <- dat.res4[,-1]
dat.res5_1 <-dat.res5[dat.res5$number=="50",]
dat.res5_2 <-dat.res5[dat.res5$number=="100",]
dat.res5_3 <-dat.res5[dat.res5$number=="200",]
