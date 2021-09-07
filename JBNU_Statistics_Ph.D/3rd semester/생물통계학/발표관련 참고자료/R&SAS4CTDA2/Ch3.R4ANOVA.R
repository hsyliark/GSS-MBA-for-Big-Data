
###################################################
### code chunk number 2: R4ANOVA.rnw:77-78
###################################################
dat   = read.csv("DBP.csv", header=T)


###################################################
### code chunk number 3: R4ANOVA.rnw:81-85
###################################################
library(xtable)
print(xtable(dat, digits=0, align=rep("c", 1+ncol(dat)),
caption="Diastolic Blood Pressure Trial Data.",label = "ANOVA.data.DBP"),
table.placement = "htbp",caption.placement = "top",include.rownames=FALSE)


###################################################
### code chunk number 4: data
###################################################
# read the data into R
dat  = read.csv("DBP.csv",header=T)
# create the difference
dat$diff = dat$DBP5-dat$DBP1
# print the first a few observation
head(dat)


###################################################
### code chunk number 5: ANOVA.fig4DBP1
###################################################
# call boxplot
boxplot(diff~TRT, dat, xlab="Treatment", 
	ylab="DBP Changes", las=1)


###################################################
### code chunk number 6: R4ANOVA.rnw:330-331
###################################################
# call boxplot
boxplot(diff~TRT, dat, xlab="Treatment", 
	ylab="DBP Changes", las=1)


###################################################
### code chunk number 7: R4ANOVA.rnw:341-343
###################################################
# call t-test with equal variance 
t.test(diff~TRT, dat, var.equal=T)


###################################################
### code chunk number 8: R4ANOVA.rnw:350-351
###################################################
t.test(diff~TRT, dat, var.equal=F)


###################################################
### code chunk number 9: R4ANOVA.rnw:355-356
###################################################
var.test(diff~TRT, dat)


###################################################
### code chunk number 10: R4ANOVA.rnw:363-364
###################################################
wilcox.test(diff~TRT, dat)


###################################################
### code chunk number 11: R4ANOVA.rnw:369-375
###################################################
# data from treatment A
diff.A = dat[dat$TRT=="A",]$diff
# data from treatment B
diff.B = dat[dat$TRT=="B",]$diff
# call t.test for one-sided test
t.test(diff.A, diff.B,alternative="less")


###################################################
### code chunk number 12: ANOVA.boot
###################################################
# load the library "bootstrap"
library(bootstrap)
# define a function to calculate the mean difference 
#	between treatment groups A to B:
mean.diff = function(bn,dat) 
  diff(tapply(dat[bn,]$diff, dat[bn,]$TRT,mean))


###################################################
### code chunk number 13: R4ANOVA.rnw:396-400
###################################################
# number of bootstrap
nboot     = 1000
# call "bootstrap" function
boot.mean = bootstrap(1:dim(dat)[1], nboot, mean.diff,dat)


###################################################
### code chunk number 14: ANOVA.fig4DBP1.boot
###################################################
# extract the mean differences  
x = boot.mean$thetastar
# calcualte the bootstrap quantiles
x.quantile = quantile(x, c(0.025,0.5, 0.975))
# show the quantiles
print(x.quantile)
# make a histogram
hist(boot.mean$thetastar, xlab="Mean Differences", main="")
# add the vertical lines for the quantiles 
abline(v=x.quantile,lwd=2, lty=c(4,1,4)) 


###################################################
### code chunk number 15: R4ANOVA.rnw:419-420
###################################################
# extract the mean differences  
x = boot.mean$thetastar
# calcualte the bootstrap quantiles
x.quantile = quantile(x, c(0.025,0.5, 0.975))
# show the quantiles
print(x.quantile)
# make a histogram
hist(boot.mean$thetastar, xlab="Mean Differences", main="")
# add the vertical lines for the quantiles 
abline(v=x.quantile,lwd=2, lty=c(4,1,4)) 


###################################################
### code chunk number 16: R4ANOVA.rnw:433-434
###################################################
aggregate(dat[,3:7], list(TRT=dat$TRT), mean)


###################################################
### code chunk number 17: R4ANOVA.rnw:440-448
###################################################
# call reshpae
Dat = reshape(dat, direction="long", 
varying=c("DBP1","DBP2","DBP3","DBP4","DBP5"),
 idvar = c("Subject","TRT","Age","Sex","diff"),sep="")
colnames(Dat) = c("Subject","TRT","Age","Sex","diff","Time","DBP")
Dat$Time = as.factor(Dat$Time)
# show the first 6 observations
head(Dat)


###################################################
### code chunk number 18: aov
###################################################
# test treatment "A"
datA   = Dat[Dat$TRT=="A",]
test.A = aov(DBP~Time, datA)
summary(test.A)
# test treatment "B"
datB   = Dat[Dat$TRT=="B",]
test.B = aov(DBP~Time, datB)
summary(test.B)


###################################################
### code chunk number 19: tukey
###################################################
TukeyHSD(test.A)
TukeyHSD(test.B)


###################################################
### code chunk number 20: R4ANOVA.rnw:493-495
###################################################
mod2 = aov(DBP~ TRT*Time, Dat)
summary(mod2)


###################################################
### code chunk number 21: ANOVA.fig4DBP1.interaction
###################################################
par(mfrow=c(2,1),mar=c(5,3,1,1))
with(Dat,interaction.plot(Time,TRT,DBP,las=1,legend=T))
with(Dat,interaction.plot(TRT,Time,DBP,las=1,legend=T))


###################################################
### code chunk number 22: R4ANOVA.rnw:507-508
###################################################
par(mfrow=c(2,1),mar=c(5,3,1,1))
with(Dat,interaction.plot(Time,TRT,DBP,las=1,legend=T))
with(Dat,interaction.plot(TRT,Time,DBP,las=1,legend=T))


###################################################
### code chunk number 23: R4ANOVA.rnw:519-520
###################################################
TukeyHSD(aov(DBP ~ TRT*Time,Dat)) 


###################################################
### code chunk number 24: diff
###################################################
# attached the data into this R session
attach(dat)
# create the changes from baseline
diff2to1 = DBP2-DBP1
diff3to1 = DBP3-DBP1
diff4to1 = DBP4-DBP1
diff5to1 = DBP5-DBP1


###################################################
### code chunk number 25: R4ANOVA.rnw:555-557
###################################################
# calculate the correlations
cor(cbind(diff2to1,diff3to1,diff4to1,diff5to1))


###################################################
### code chunk number 26: R4ANOVA.rnw:563-568
###################################################
# calculate the mean changes
MCh =aggregate(cbind(diff2to1,diff3to1,diff4to1,diff5to1),
         list(TRT=TRT), mean)
# print the chanhe
print(MCh)


###################################################
### code chunk number 27: R4ANOVA.rnw:576-586
###################################################
# call "manova" to fit a manova
maov1=manova(cbind(diff2to1,diff3to1,diff4to1,diff5to1)~TRT,dat)
# then F-test with Pillai (default in R)
summary(maov1)
# F-test with Hotelling-Lawley
summary(maov1, test="Hotelling-Lawley")
# F-test with Wilks
summary(maov1, test="Wilks")
# F-test with Roy
summary(maov1, test="Roy")


###################################################
### code chunk number 28: R4ANOVA.rnw:600-603
###################################################
n  = c(168, 182, 165,188) 
p4 = c(.41, .62, .73, .77)
x4 = c(69, 113, 120, 145)


###################################################
### code chunk number 29: R4ANOVA.rnw:608-609
###################################################
prop.test(x4, n)


###################################################
### code chunk number 30: R4ANOVA.rnw:615-616
###################################################
prop.test(x4[c(1,3)], n[c(1,3)])


###################################################
### code chunk number 31: R4ANOVA.rnw:619-620
###################################################
prop.test(x4[c(2,3)], n[c(2,3)])


###################################################
### code chunk number 32: R4ANOVA.rnw:623-624
###################################################
prop.test(x4[c(3,4)], n[c(3,4)])


###################################################
### code chunk number 33: R4ANOVA.rnw:631-640
###################################################
# create a dataframe for the Ulcer trial
Ulcer = data.frame(
# use ``factor" to create the treatment factor
trt   = factor(rep(c("0 mg C","400 mg C","800 mg C","1600 mg C"), 
each=2),levels=c("0 mg C","400 mg C","800 mg C","1600 mg C") ), 
Heal  = c("Yes","No","Yes","No","Yes","No","Yes","No"), 
y     = c(x4[1],n[1]-x4[1],x4[2],n[2]-x4[2],x4[3],
	 n[3]-x4[3],x4[4],n[4]-x4[4]))
Ulcer


###################################################
### code chunk number 34: R4ANOVA.rnw:644-646
###################################################
tab.Ulcer = xtabs(y~trt+Heal,Ulcer)
tab.Ulcer


###################################################
### code chunk number 35: ANOVA.fig4ulcer1
###################################################
# layout for the plot
par(mfrow=c(1,2), mar=c(4,2,1,1))
# call ``dotchart"
dotchart(tab.Ulcer)
# call ``mosaicplot"
mosaicplot(tab.Ulcer,color=T,las=1, main=" ", 
	xlab="Treatment",ylab="Heal Status" )


###################################################
### code chunk number 36: R4ANOVA.rnw:661-662
###################################################
# layout for the plot
par(mfrow=c(1,2), mar=c(4,2,1,1))
# call ``dotchart"
dotchart(tab.Ulcer)
# call ``mosaicplot"
mosaicplot(tab.Ulcer,color=T,las=1, main=" ", 
	xlab="Treatment",ylab="Heal Status" )


###################################################
### code chunk number 37: R4ANOVA.rnw:669-670
###################################################
margin.table(tab.Ulcer,1)


###################################################
### code chunk number 38: R4ANOVA.rnw:674-675
###################################################
summary(tab.Ulcer)


