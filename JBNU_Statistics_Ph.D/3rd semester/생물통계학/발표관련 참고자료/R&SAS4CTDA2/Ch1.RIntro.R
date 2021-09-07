### R code from vignette source 'RIntro.rnw'


###################################################
### code chunk number 9: RI-simuNorm
###################################################
help(rnorm)


###################################################
### code chunk number 10: RI-input
###################################################
# simulated input values
n      = 100
mu     = 100
sd     = 10
mu.d   = 20
age.mu = 50
age.sd = 10


###################################################
### code chunk number 11: RI-placebo
###################################################
# fix the seed for random number generation 
set.seed(123)
# use "rnorm" to generate random normal
age         = rnorm(n, age.mu, age.sd)
bp.base     = rnorm(n,mu,sd)
bp.end      = rnorm(n,mu,sd)
# take the difference between endpoint and baseline
bp.diff     = bp.end-bp.base
# put the data together using "cbind" to column-bind
dat4placebo = round(cbind(age,bp.base,bp.end,bp.diff))


###################################################
### code chunk number 12: RI-headplacebo
###################################################
head(dat4placebo)


###################################################
### code chunk number 13: RI-drug
###################################################
age      = rnorm(n, age.mu, age.sd)
bp.base  = rnorm(n,mu,sd)
bp.end   = rnorm(n,mu-mu.d,sd)
bp.diff  = bp.end-bp.base
dat4drug = round(cbind(age,bp.base,bp.end,bp.diff))


###################################################
### code chunk number 14: RI-dat
###################################################
# make a dataframe to hold all data
dat     = data.frame(rbind(dat4placebo,dat4drug))
# make "trt" as a factor for treatment.
dat$trt = as.factor(rep(c("Placebo", "Drug"), each=n))


###################################################
### code chunk number 15: RIntro.rnw:353-357
###################################################
# check the data dimension
dim(dat)
# print the first 6 obervations to see the variable names
head(dat)


###################################################
### code chunk number 16: RI.boxplot4placebo
###################################################
# call boxplot
boxplot(dat4placebo, las=1, main="Placebo")


###################################################
### code chunk number 17: RIntro.rnw:373-374
###################################################
# call boxplot
boxplot(dat4placebo, las=1, main="Placebo")


###################################################
### code chunk number 18: RI.boxplot4drug
###################################################
boxplot(dat4drug, las=1, main="Drug")


###################################################
### code chunk number 19: RIntro.rnw:387-388
###################################################
boxplot(dat4drug, las=1, main="Drug")


###################################################
### code chunk number 20: RI.xyplot.dat
###################################################
#load the lattice library
library(lattice)
# call xyplot function and print it
print(xyplot(bp.diff~age|trt, data=dat,xlab="Age", 
strip=strip.custom(bg="white"), 
ylab="Blood Pressure Difference",lwd=3,cex=1.3,pch=20,
type=c("p", "r")))


###################################################
### code chunk number 21: RIntro.rnw:421-422
###################################################
#load the lattice library
library(lattice)
# call xyplot function and print it
print(xyplot(bp.diff~age|trt, data=dat,xlab="Age", 
strip=strip.custom(bg="white"), 
ylab="Blood Pressure Difference",lwd=3,cex=1.3,pch=20,
type=c("p", "r")))


###################################################
### code chunk number 22: RI-lm
###################################################
lm1 = lm(bp.diff~trt*age, data=dat)
summary(lm1)


###################################################
### code chunk number 23: tablm
###################################################
# load the xtable library
library(xtable)
# call xtable to make the table
print(xtable(lm1, caption="ANOVA Table for Simulated 
Clinical Trial Data", label = "tab4RI.coef"),
table.placement = "htbp",caption.placement = "top")


###################################################
### code chunk number 24: RI.fig4lm
###################################################
layout(matrix(1:4, nrow=2))
plot(lm1)


###################################################
### code chunk number 25: RIntro.rnw:469-470
###################################################
layout(matrix(1:4, nrow=2))
plot(lm1)


