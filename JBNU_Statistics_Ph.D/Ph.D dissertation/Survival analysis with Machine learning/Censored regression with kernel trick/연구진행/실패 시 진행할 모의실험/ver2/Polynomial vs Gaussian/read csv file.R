dat.res1 <- read.csv("C:/Users/Hi/Desktop/ver2/Polynomial vs Gaussian/설명변수 3/중도절단 0/p3cen0.csv",sep=",",header=T)
dat.res1 <- dat.res1[,-1]
dat.res1_1 <-dat.res1[dat.res1$number=="50",]
dat.res1_2 <-dat.res1[dat.res1$number=="100",]
dat.res1_3 <-dat.res1[dat.res1$number=="200",]
