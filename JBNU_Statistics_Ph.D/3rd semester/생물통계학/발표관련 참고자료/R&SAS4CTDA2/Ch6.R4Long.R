### R code from vignette source 'R4Long.rnw'

###################################################
### code chunk number 1: R4Long.rnw:10-27
###################################################
options(continue=" ") # from Ross Ihaka to remove the "+" continuation
rm(list=ls())
badfonts = FALSE
.detach.stuff = function() {
s1 = grep("package|Autoloads",search())
nattach = length(search())
xattach = search()[c(-1,-s1)]
for (i in xattach)
eval(substitute(detach(i),list(i=i)))
}
.detach.stuff()
## ps.options(family="NimbusSan")
ps.options(colormodel="cmyk")
options(width=53, digits=3)
## palette(gray((0:8)/8))
#library(lattice)
#trellis.par.set(canonical.theme(color=FALSE))


###################################################
### code chunk number 2: Long.dat1
###################################################
dat  = read.csv("DBP.csv",header=T)
# show first a few observations
head(dat)


###################################################
### code chunk number 3: R4Long.rnw:198-199
###################################################
summary(dat)


###################################################
### code chunk number 4: R4Long.rnw:207-215
###################################################
# reshape the data into "long" direction
Dat = reshape(dat, direction="long", 
 varying=c("DBP1","DBP2","DBP3","DBP4","DBP5"),
 idvar = c("Subject","TRT","Age","Sex"),sep="")
# rename the variables
colnames(Dat) = c("Subject","TRT","Age","Sex","Time","DBP")
# print the first 6 observations
head(Dat)


###################################################
### code chunk number 5: Long.figdat1.1
###################################################
library(lattice)
print(xyplot(DBP~Time|as.factor(Subject),type="l",
   groups=TRT,lty=c(1,8),lwd=2,layout=c(10,4), Dat))


###################################################
### code chunk number 6: R4Long.rnw:227-228
###################################################
library(lattice)
print(xyplot(DBP~Time|as.factor(Subject),type="l",
   groups=TRT,lty=c(1,8),lwd=2,layout=c(10,4), Dat))


###################################################
### code chunk number 7: Long.figdat1.2
###################################################
print(xyplot(DBP~Time|TRT,type="l",Dat, groups=as.factor(Subject)))


###################################################
### code chunk number 8: R4Long.rnw:239-240
###################################################
print(xyplot(DBP~Time|TRT,type="l",Dat, groups=as.factor(Subject)))


###################################################
### code chunk number 9: Long.figdat1.3
###################################################
print(bwplot(DBP~as.factor(Time)|TRT,Dat, xlab="Time"))


###################################################
### code chunk number 10: R4Long.rnw:255-256
###################################################
print(bwplot(DBP~as.factor(Time)|TRT,Dat, xlab="Time"))


###################################################
### code chunk number 11: Long.intslope
###################################################
num.Subj = 40
# initiate the intercept and slope
intercept = slope = numeric(num.Subj)
# loop-over
for(i in 1:num.Subj){
# fit regression model
mod          = lm(DBP~Time, Dat[Dat$Subject==i,])
# extract the intercept and slope
intercept[i] = coef(mod)[1]
slope[i]     = coef(mod)[2]
}
# make a dataframe "dat.coef" 
dat.coef = data.frame(Subject=dat$Subject,TRT=dat$TRT,
Intercept = intercept, Slope=slope)
# print it out
dat.coef


###################################################
### code chunk number 12: Long.figdat1.interceptslope
###################################################
# Make histogram for both intercept and slope
int.hist   = hist(intercept,plot=F)
slope.hist = hist(slope,plot=F)
# make layout for plotting       
top        = max(c(int.hist$counts, slope.hist$counts))
nf         = layout(matrix(c(2,0,1,3),2,2,byrow=T),
		c(3,1), c(1,3),T)
par(mar=c(5,4,1,1))
# plot the intercept and slope
plot(Slope~Intercept,las=1,dat.coef,xlab="Intercept",
  ylab="Slope",pch=as.character(TRT))
par(mar=c(0,4,1,1))
# add the intercept histogram
barplot(int.hist$counts, axes=FALSE, 
	ylim=c(0, top), space=0)
par(mar=c(5,0,1,1))
# add the slope histogram
barplot(slope.hist$counts, axes=FALSE, 
	xlim=c(0, top), space=0, horiz=TRUE)


###################################################
### code chunk number 13: R4Long.rnw:308-309
###################################################
# Make histogram for both intercept and slope
int.hist   = hist(intercept,plot=F)
slope.hist = hist(slope,plot=F)
# make layout for plotting       
top        = max(c(int.hist$counts, slope.hist$counts))
nf         = layout(matrix(c(2,0,1,3),2,2,byrow=T),
		c(3,1), c(1,3),T)
par(mar=c(5,4,1,1))
# plot the intercept and slope
plot(Slope~Intercept,las=1,dat.coef,xlab="Intercept",
  ylab="Slope",pch=as.character(TRT))
par(mar=c(0,4,1,1))
# add the intercept histogram
barplot(int.hist$counts, axes=FALSE, 
	ylim=c(0, top), space=0)
par(mar=c(5,0,1,1))
# add the slope histogram
barplot(slope.hist$counts, axes=FALSE, 
	xlim=c(0, top), space=0, horiz=TRUE)


###################################################
### code chunk number 14: R4Long.rnw:315-321
###################################################
# fit model 1 with interation
mod1.coef = lm(Slope~Intercept*TRT, dat.coef)
summary(mod1.coef)
# fit model 2 without interaction
mod2.coef = lm(Slope~Intercept+TRT, dat.coef)
summary(mod2.coef)


###################################################
### code chunk number 15: R4Long.rnw:326-330
###################################################
# test slope difference
t.test(Slope~TRT, dat.coef)
# test intercept difference
t.test(Intercept~TRT, dat.coef)


###################################################
### code chunk number 16: R4Long.rnw:339-347
###################################################
# load the library lme4
library(lmerTest)
# Fit Model 1
mod1DBP =  lmer(DBP~TRT*Time+(Time|Subject), Dat)
# Fit Model 2
mod2DBP =  lmer(DBP~TRT*Time+(1|Subject), Dat)
# model comparison   
anova(mod1DBP, mod2DBP)


###################################################
### code chunk number 17: R4Long.rnw:351-352
###################################################
summary(mod2DBP)


###################################################
### code chunk number 18: R4Long.rnw:357-363
###################################################
# fit Model 3 
mod3DBP =  lmer(DBP~TRT+Time+(Time|Subject), Dat)
# fit Model 4
mod4DBP =  lmer(DBP~TRT+Time+(1|Subject), Dat)
# model comparison
anova(mod3DBP, mod4DBP)


###################################################
### code chunk number 19: R4Long.rnw:368-369
###################################################
summary(mod3DBP)


###################################################
### code chunk number 20: Long.figdat1.qqmath
###################################################
print(qqmath(~resid(mod3DBP)|TRT,Dat, 
    xlab="Theoretical Normal", ylab="Residuals"))


###################################################
### code chunk number 21: R4Long.rnw:388-389
###################################################
print(qqmath(~resid(mod3DBP)|TRT,Dat, 
    xlab="Theoretical Normal", ylab="Residuals"))


###################################################
### code chunk number 22: Long.figdat1.resid2
###################################################
print(bwplot(resid(mod3DBP)~as.factor(Time)|TRT,Dat, 
	xlab="Time",ylim=c(-5,5), ylab="Residuals"))


###################################################
### code chunk number 23: R4Long.rnw:403-404
###################################################
print(bwplot(resid(mod3DBP)~as.factor(Time)|TRT,Dat, 
	xlab="Time",ylim=c(-5,5), ylab="Residuals"))


###################################################
### code chunk number 24: R4Long.rnw:411-415
###################################################
# fit Model 3 include ``Age" effect
mod5DBP =  lmer(DBP~TRT+Time+Age+(Time|Subject), Dat)
# call anova to test ``Age" effect
anova(mod3DBP, mod5DBP)


###################################################
### code chunk number 25: R4Long.rnw:419-423
###################################################
# fit Model 6 including ``Age" and ``Sex"
mod6DBP =  lmer(DBP~TRT+Time+Age+Sex+(Time|Subject), Dat)
# test the ``Sex" effect
anova(mod5DBP, mod6DBP)


###################################################
### code chunk number 26: Long.dat2
###################################################
dat  = read.csv("Ulcer.csv", header=T)
# print the first 6 obs
head(dat)


###################################################
### code chunk number 27: R4Long.rnw:442-453
###################################################
# total n for each TRT
n  = tapply(rep(1, dim(dat)[1]),dat$TRT,sum)
# number for time 1
n1 = tapply(dat$Time1,dat$TRT,sum)
# number for time 2
n2 = tapply(dat$Time2,dat$TRT,sum)
# number for time 4
n4 = tapply(dat$Time4,dat$TRT,sum)
print(rbind(n,n1,n2,n4))
# proportions
print( round(rbind(n1/n,n2/n,n4/n),2))


###################################################
### code chunk number 28: R4Long.rnw:458-471
###################################################
Dat     = reshape(dat, direction="long", 
varying = c("Time0","Time1","Time2","Time4"),
 idvar  = c("Subject","TRT","WeekH"),sep="")
colnames(Dat) = c("Subject","TRT","WeekH","Time","Heal")
# sort the data by Subject: very important for gee
Dat      = Dat[order(Dat$Subject),] 
# Remove the baseline for model fitting
Dat      = Dat[Dat$Time > 0,]
# make the TRT and Time as factors
Dat$TRT  = as.factor(Dat$TRT)
Dat$Time = as.factor(Dat$Time) 
# show the first 6 observations
head(Dat)


###################################################
### code chunk number 29: R4Long.rnw:475-481
###################################################
# fit Model 1: with interaction
mod1glm = glm(Heal~TRT*Time, family=binomial, Dat)
# fit Model 2: without interaction
mod2glm = glm(Heal~TRT+Time, family=binomial, data=Dat)
# test these two model using Chi-Square test
anova(mod1glm,mod2glm, test="Chi")


###################################################
### code chunk number 30: R4Long.rnw:485-486
###################################################
summary(mod2glm)


###################################################
### code chunk number 31: R4Long.rnw:492-497
###################################################
# load the ``multcomp" library
library(multcomp)
# multiple comparisons
glht.mod2glm = glht(mod2glm, mcp(TRT="Tukey", Time="Tukey"))
summary(glht.mod2glm)


###################################################
### code chunk number 32: R4Long.rnw:505-512
###################################################
# load MASS library
library(MASS)
# fit the Model 3
mod3glm = glmmPQL(Heal~TRT, random=~1|Subject, 
	family=binomial, Dat)
# print the summary
summary(mod3glm)


###################################################
### code chunk number 33: R4Long.rnw:519-522
###################################################
# fit Model 4
mod4glm = glm(Heal~TRT, family=binomial, Dat)
summary(mod4glm)


###################################################
### code chunk number 34: R4Long.rnw:530-537
###################################################
# load the ``gee" library
library(gee)
# fit the gee model with independent patient effect
fit.gee1 = gee(Heal~TRT,id=Subject,family=binomial, 
data=Dat,corstr="independence", scale.fix=T)
# print the summary
summary(fit.gee1)


###################################################
### code chunk number 35: R4Long.rnw:543-546
###################################################
fit.gee2 = gee(Heal~TRT,id=Subject,family=binomial, 
data=Dat,corstr="exchangeable", scale.fix=T)
summary(fit.gee2)


