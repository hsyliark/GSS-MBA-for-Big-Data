### R code from vignette source 'R4ANCOVA.rnw'

###################################################
### code chunk number 1: R4ANCOVA.rnw:11-34
###################################################
options(continue=" ") # from Ross Ihaka to remove the "+" continution
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
# lattice setting
library(lattice)
#trellis.par.set(canonical.theme(color=FALSE))
ltheme = canonical.theme(color = FALSE) # in-built B&W theme
ltheme$strip.background$col = "transparent " # change strip bg
lattice.options(default.theme=ltheme) # set as default
# random seed
set.seed(12345678)


###################################################
### code chunk number 2: R4ANCOVA.rnw:334-337
###################################################
dat  = read.csv("DBP.csv",header=T)
# creat the difference
dat$diff = dat$DBP5-dat$DBP1


###################################################
### code chunk number 3: ANCOVA.BaseDBP
###################################################
boxplot(DBP1~TRT, dat, las=1, 
	xlab="Treatment", ylab="DBP at Baseline")


###################################################
### code chunk number 4: R4ANCOVA.rnw:350-351
###################################################
boxplot(DBP1~TRT, dat, las=1, 
	xlab="Treatment", ylab="DBP at Baseline")


###################################################
### code chunk number 5: R4ANCOVA.rnw:360-361
###################################################
t.test(DBP1~TRT, dat)


###################################################
### code chunk number 6: R4ANCOVA.rnw:368-374
###################################################
# call function table to make the 2 by 2 table 
SexbyTRT = table(dat$TRT,dat$Sex)
# print it
SexbyTRT
# call prop.test to test the difference
prop.test(SexbyTRT)


###################################################
### code chunk number 7: R4ANCOVA.rnw:381-385
###################################################
# Fit the main effect model on "Sex" and "Age"
bm1=lm(DBP1~Sex+Age, dat)
# Show the result
summary(bm1)


###################################################
### code chunk number 8: ANCOVA.fig4Baseline
###################################################
# plot the ``Age" to ``DBP1"
plot(DBP1~Age,las=1,pch=as.character(Sex), dat, 
	xlab="Age", ylab="Baseline DBP")
# add the regression lines using ``abline"
abline(bm1$coef[1], bm1$coef[3],lwd=2, lty=1)
abline(bm1$coef[1]+bm1$coef[2], bm1$coef[3],lwd=2, lty=4)


###################################################
### code chunk number 9: R4ANCOVA.rnw:402-403
###################################################
# plot the ``Age" to ``DBP1"
plot(DBP1~Age,las=1,pch=as.character(Sex), dat, 
	xlab="Age", ylab="Baseline DBP")
# add the regression lines using ``abline"
abline(bm1$coef[1], bm1$coef[3],lwd=2, lty=1)
abline(bm1$coef[1]+bm1$coef[2], bm1$coef[3],lwd=2, lty=4)


###################################################
### code chunk number 10: R4ANCOVA.rnw:415-421
###################################################
# start with full model
m0 = lm(diff~TRT*Age*Sex, dat)
# stepwise model selection
m1 = step(m0)
# output the ANOVA
anova(m1)


###################################################
### code chunk number 11: R4ANCOVA.rnw:425-431
###################################################
# fit the reduced model
m2 = lm(diff~TRT+Age, dat)
# output the anova
anova(m2)
# output the model fit
summary(m2)


###################################################
### code chunk number 12: ANCOVA.fig4diff
###################################################
plot(diff~Age,las=1,pch=as.character(TRT), dat, 
xlab="Age", ylab="DBP Change")
abline(m2$coef[1], m2$coef[3],lwd=2, lty=1)
abline(m2$coef[1]+m2$coef[2], m2$coef[3],lwd=2, lty=4)


###################################################
### code chunk number 13: R4ANCOVA.rnw:445-446
###################################################
plot(diff~Age,las=1,pch=as.character(TRT), dat, 
xlab="Age", ylab="DBP Change")
abline(m2$coef[1], m2$coef[3],lwd=2, lty=1)
abline(m2$coef[1]+m2$coef[2], m2$coef[3],lwd=2, lty=4)


###################################################
### code chunk number 14: R4ANCOVA.rnw:466-473
###################################################
# attached the data into this R session
attach(dat)
# create the changes from baseline
diff2to1 = DBP2-DBP1
diff3to1 = DBP3-DBP1
diff4to1 = DBP4-DBP1
diff5to1 = DBP5-DBP1


###################################################
### code chunk number 15: R4ANCOVA.rnw:479-490
###################################################
# call "manova" to fit a manova adjusting for "Age"
macov1=manova(cbind(diff2to1,diff3to1,diff4to1,
			diff5to1)~TRT+Age,dat)
# then F-test with Pillai (default in R)
summary(macov1)
# F-test with Hotelling-Lawley
summary(macov1, test="Hotelling-Lawley")
# F-test with Wilks
summary(macov1, test="Wilks")
# F-test with Roy
summary(macov1, test="Roy")


###################################################
### code chunk number 16: R4ANCOVA.rnw:502-506
###################################################
library(flexmix)
data("betablocker")  
betablocker$Center = as.factor(betablocker$Center)
colnames(betablocker) = c("Deaths","Total","Center","TRT")


###################################################
### code chunk number 17: R4ANCOVA.rnw:508-509
###################################################
betablocker


###################################################
### code chunk number 18: beta.glm
###################################################
# fit a logistic regression using glm
beta.glm = glm(cbind(Deaths,Total-Deaths)~TRT+Center,
		family=binomial,data=betablocker)
# print the model fitting
anova(beta.glm)


###################################################
### code chunk number 19: R4ANCOVA.rnw:522-523
###################################################
summary(beta.glm)


###################################################
### code chunk number 20: R4ANCOVA.rnw:529-531
###################################################
est.dp = sum(resid(beta.glm, type="pearson")^2)/beta.glm$df.res
est.dp


###################################################
### code chunk number 21: R4ANCOVA.rnw:535-536
###################################################
summary(beta.glm, dispersion=est.dp)


###################################################
### code chunk number 22: R4ANCOVA.rnw:540-545
###################################################
# fit quasi-likelihood for binomial data
beta.glm2 = glm(cbind(Deaths,Total- Deaths)~TRT+Center,
	family=quasibinomial,data=betablocker)
# print the model fit
summary(beta.glm2)


###################################################
### code chunk number 23: R4ANCOVA.rnw:554-556
###################################################
library(HSAUR)
data(polyps)


###################################################
### code chunk number 24: R4ANCOVA.rnw:561-562
###################################################
polyps


###################################################
### code chunk number 25: R4ANCOVA.rnw:566-570
###################################################
# Poisson Regression
m0.polyps = glm(number~treat*age, polyps, family=poisson())
# print the model fit
summary(m0.polyps)


###################################################
### code chunk number 26: R4ANCOVA.rnw:574-576
###################################################
est.dp = sum(resid(m0.polyps, type="pearson")^2)/m0.polyps$df.res
est.dp


###################################################
### code chunk number 27: R4ANCOVA.rnw:579-580
###################################################
summary(m0.polyps, dispersion=est.dp)


###################################################
### code chunk number 28: R4ANCOVA.rnw:584-592
###################################################
# refit the model without interaction
m1.polyps = glm(number~treat+age, polyps, family=poisson())
# estimate the dispersion parameter
est.dp = sum(resid(m1.polyps, type="pearson")^2)/m1.polyps$df.res
# print the estimated dispersion parameter
est.dp
# print the model fit adjusting the over dispersion
summary(m1.polyps, dispersion=est.dp)


###################################################
### code chunk number 29: R4ANCOVA.rnw:597-601
###################################################
# fit the quasi Poisson
m2.polyps = glm(number~treat+age, polyps, family=quasipoisson())
# print the model fit
summary(m2.polyps)


###################################################
### code chunk number 30: R4ANCOVA.rnw:607-613
###################################################
# load the MASS library
library(MASS)
# fit the negative binomial model
m3.polyps = glm.nb(number~treat+age, polyps)
# print the model fit
summary(m3.polyps)


