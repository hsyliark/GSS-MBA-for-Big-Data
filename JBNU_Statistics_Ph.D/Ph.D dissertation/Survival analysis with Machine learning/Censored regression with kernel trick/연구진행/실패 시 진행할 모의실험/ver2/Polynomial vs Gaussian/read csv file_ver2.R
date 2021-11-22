dat.res16 <- read.csv("C:/Users/Hi/Desktop/ver2/Polynomial vs Gaussian/설명변수 9/중도절단 50/p9cen50.csv",sep=",",header=T)
dat.res16 <- dat.res16[,-1]
dat.res16_1 <-dat.res16[dat.res16$number=="50",]
dat.res16_2 <-dat.res16[dat.res16$number=="100",]
dat.res16_3 <-dat.res16[dat.res16$number=="200",]
