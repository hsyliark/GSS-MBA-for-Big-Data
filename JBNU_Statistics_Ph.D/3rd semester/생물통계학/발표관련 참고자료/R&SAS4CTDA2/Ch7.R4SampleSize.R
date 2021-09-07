#

###################################################
### chunk number 2: SS.pow4num
###################################################
#line 152 "R4SampleSize.rnw"
# Type-I error
alpha = 0.05
# Type-II error
beta  = c(0.05,0.1,0.15,0.2,0.25,0.3)
# power
pow   = 1-beta
# numerator in the sample size
num = 2*(qnorm(1-alpha/2)+qnorm(1-beta))^2
# plot the power to the numerator
plot(pow, num, xlab="Power",las=1, ylab="Numerator")
# add the line to it
lines(pow, num)
# use arrows to show the values of numerator 
for(i in 1:length(pow)){
arrows(pow[i],0, pow[i], num[i], length=0.13)
arrows(pow[i],num[i], pow[length(beta)],num[i], length=0.13)
}


###################################################
### chunk number 3: 
###################################################
#line 174 "R4SampleSize.rnw"
#line 152 "R4SampleSize.rnw#from line#174#"
# Type-I error
alpha = 0.05
# Type-II error
beta  = c(0.05,0.1,0.15,0.2,0.25,0.3)
# power
pow   = 1-beta
# numerator in the sample size
num = 2*(qnorm(1-alpha/2)+qnorm(1-beta))^2
# plot the power to the numerator
plot(pow, num, xlab="Power",las=1, ylab="Numerator")
# add the line to it
lines(pow, num)
# use arrows to show the values of numerator 
for(i in 1:length(pow)){
arrows(pow[i],0, pow[i], num[i], length=0.13)
arrows(pow[i],num[i], pow[length(beta)],num[i], length=0.13)
}
#line 175 "R4SampleSize.rnw"


###################################################
### chunk number 4: 
###################################################
#line 196 "R4SampleSize.rnw"
power.t.test(delta=0.5, sd=1, power=0.8)


###################################################
### chunk number 5: 
###################################################
#line 201 "R4SampleSize.rnw"
power.t.test(delta=0.5, sd=1, power=0.8,
	 alternative = c("one.sided"))


###################################################
### chunk number 6: 
###################################################
#line 209 "R4SampleSize.rnw"
power.t.test(n=64,delta=0.5, sd=1)


###################################################
### chunk number 7: 
###################################################
#line 214 "R4SampleSize.rnw"
power.t.test(n=64,sd=1,power=0.8)


###################################################
### chunk number 8: SS.pow2size
###################################################
#line 220 "R4SampleSize.rnw"
# use pow from 0.2 to 0.9 by 0.05
pow  = seq(0.2, 0.9, by=0.05)
# keep track of the size using for-loop
size = NULL
for(i in 1:length(pow))
size[i] = power.t.test(delta=0.5, sd=1, power=pow[i])$n
# plot the size to power
plot(pow, size, las=1,type="b", xlab="Power", 
	ylab="Sample Size Required")


###################################################
### chunk number 9: 
###################################################
#line 234 "R4SampleSize.rnw"
#line 220 "R4SampleSize.rnw#from line#234#"
# use pow from 0.2 to 0.9 by 0.05
pow  = seq(0.2, 0.9, by=0.05)
# keep track of the size using for-loop
size = NULL
for(i in 1:length(pow))
size[i] = power.t.test(delta=0.5, sd=1, power=pow[i])$n
# plot the size to power
plot(pow, size, las=1,type="b", xlab="Power", 
	ylab="Sample Size Required")
#line 235 "R4SampleSize.rnw"


###################################################
### chunk number 10: 
###################################################
#line 259 "R4SampleSize.rnw"
# load the library
library(samplesize)
# sample size calculation
n.indep.t.test.eq(power = 0.8, alpha = 0.95,
	 mean.diff = 0.8, sd.est = 0.83)


###################################################
### chunk number 11: 
###################################################
#line 269 "R4SampleSize.rnw"
n.indep.t.test.neq(power = 0.8, alpha = 0.95, 
	mean.diff = 0.8, sd.est = 0.83, k=0.5)


###################################################
### chunk number 12: 
###################################################
#line 282 "R4SampleSize.rnw"
n.welch.test(power = 0.8, alpha = 0.95,
	 mean.diff = 2, sd.est1 = 1, sd.est2 = 2)


###################################################
### chunk number 13: 
###################################################
#line 304 "R4SampleSize.rnw"
n.wilcox.ord(beta = 0.2, alpha = 0.05, t = 0.5, 
	p = c(0.66, 0.15, 0.19), q = c(0.61, 0.23, 0.16))


###################################################
### chunk number 14: 
###################################################
#line 349 "R4SampleSize.rnw"
power.prop.test(p1 = .75, p2 = .50, power = .80)


###################################################
### chunk number 15: 
###################################################
#line 355 "R4SampleSize.rnw"
power.prop.test(n = 60, p1 = .75, p2 = .5)


###################################################
### chunk number 16: SS.pow2alpha
###################################################
#line 364 "R4SampleSize.rnw"
# set up the power range
pow   = seq(0.5, 0.9, by=0.05)
# a for-loop to calculate alpha
alpha = NULL
for(i in 1:length(pow)){
alpha[i] = power.prop.test(n=60, p1=0.75, p2=0.5, 
	power=pow[i], sig.level=NULL)$sig.level
}
# make the plot
plot(pow, alpha, las=1,type="b", lwd=2, xlab="Power", 
		ylab="Significance Level")
# add a segment for alpha=0.05
segments(pow[1], 0.05, 0.816,0.05, lwd=2)
# point to the power=0.816 for alpha=0.05
arrows(0.816,0.05, 0.816,0, lwd=2)


###################################################
### chunk number 17: 
###################################################
#line 384 "R4SampleSize.rnw"
#line 364 "R4SampleSize.rnw#from line#384#"
# set up the power range
pow   = seq(0.5, 0.9, by=0.05)
# a for-loop to calculate alpha
alpha = NULL
for(i in 1:length(pow)){
alpha[i] = power.prop.test(n=60, p1=0.75, p2=0.5, 
	power=pow[i], sig.level=NULL)$sig.level
}
# make the plot
plot(pow, alpha, las=1,type="b", lwd=2, xlab="Power", 
		ylab="Significance Level")
# add a segment for alpha=0.05
segments(pow[1], 0.05, 0.816,0.05, lwd=2)
# point to the power=0.816 for alpha=0.05
arrows(0.816,0.05, 0.816,0, lwd=2)
#line 385 "R4SampleSize.rnw"


###################################################
### chunk number 18: 
###################################################
#line 410 "R4SampleSize.rnw"
# load the library into R
library(pwr)
# display the help menu 
library(help=pwr)


###################################################
### chunk number 19: h
###################################################
#line 418 "R4SampleSize.rnw"
h = ES.h(0.75,0.5)
print(h)


###################################################
### chunk number 20: 
###################################################
#line 428 "R4SampleSize.rnw"
pwr.2p.test(h=h,power=0.8,sig.level=0.05)


###################################################
### chunk number 21: prop
###################################################
#line 435 "R4SampleSize.rnw"
power.prop.test(p1=0.75, p2=0.5, power=0.8)


###################################################
### chunk number 22: 
###################################################
#line 448 "R4SampleSize.rnw"
library(gsDesign)


###################################################
### chunk number 23: 
###################################################
#line 454 "R4SampleSize.rnw"
help(nBinomial)


###################################################
### chunk number 24:  eval=FALSE
###################################################
## #line 488 "R4SampleSize.rnw"
## power.prop.test(p1=0.75,p2=0.5, power=0.8)


###################################################
### chunk number 25: 
###################################################
#line 492 "R4SampleSize.rnw"
nBinomial(p1=0.75, p2=0.5, alpha=.05, beta=0.2, sided=2) 


###################################################
### chunk number 26: 
###################################################
#line 502 "R4SampleSize.rnw"
nBinomial(p1=.35, p2=.2, beta=.2, ratio=0.5,outtype=2, 
	alpha=.05, sided=1)


###################################################
### chunk number 27: 
###################################################
#line 509 "R4SampleSize.rnw"
nBinomial(p1=.35, p2=.2, beta=.2, alpha=.05, sided=1)


###################################################
### chunk number 28: SS.riskreduction
###################################################
#line 518 "R4SampleSize.rnw"
# sequence of control event rate
p1  =  seq(.1, .3, .01)
# reduce by 30%  to calculate the sample size required
p2 <- p1 *.7
y1 <- nBinomial(p1, p2, beta=.2, outtype=1, alpha=.025, sided=1)
# reduce by 40%  to calculate the sample size required
p2 <- p1 * .6
y2 <- nBinomial(p1, p2, beta=.2, outtype=1, alpha=.025, sided=1)
# reduce by 50%  to calculate the sample size required
p2 <- p1 * .5
y3 <- nBinomial(p1, p2, beta=.2, outtype=1, alpha=.025, sided=1)
# reduce by 60%  to calculate the sample size required
p2 <- p1 * .4
y4 <- nBinomial(p1, p2, beta=.2, outtype=1, alpha=.025, sided=1)
# make the plot for 30% reduction 
plot(p1, y1, type="l", las=1,ylab="Sample Size",
     xlab="Control Group Event Rate", ylim=c(0, 3000), lwd=2)
# add a line for 40% reduction
lines(p1, y2, lty=2, lwd=2)
# add a line for 50% reduction
lines(p1, y3, lty=3, lwd=2)
# add a line for 60% reduction
lines(p1, y4, lty=4, lwd=2)
# add a legend
legend("topright",lty=c(1,2, 3, 4), lwd=2,
 legend=c("30 pct reduction", "40 pct reduction",
                "50 pct reduction", "60 pct reduction"))


###################################################
### chunk number 29: 
###################################################
#line 550 "R4SampleSize.rnw"
#line 518 "R4SampleSize.rnw#from line#550#"
# sequence of control event rate
p1  =  seq(.1, .3, .01)
# reduce by 30%  to calculate the sample size required
p2 <- p1 *.7
y1 <- nBinomial(p1, p2, beta=.2, outtype=1, alpha=.025, sided=1)
# reduce by 40%  to calculate the sample size required
p2 <- p1 * .6
y2 <- nBinomial(p1, p2, beta=.2, outtype=1, alpha=.025, sided=1)
# reduce by 50%  to calculate the sample size required
p2 <- p1 * .5
y3 <- nBinomial(p1, p2, beta=.2, outtype=1, alpha=.025, sided=1)
# reduce by 60%  to calculate the sample size required
p2 <- p1 * .4
y4 <- nBinomial(p1, p2, beta=.2, outtype=1, alpha=.025, sided=1)
# make the plot for 30% reduction 
plot(p1, y1, type="l", las=1,ylab="Sample Size",
     xlab="Control Group Event Rate", ylim=c(0, 3000), lwd=2)
# add a line for 40% reduction
lines(p1, y2, lty=2, lwd=2)
# add a line for 50% reduction
lines(p1, y3, lty=3, lwd=2)
# add a line for 60% reduction
lines(p1, y4, lty=4, lwd=2)
# add a legend
legend("topright",lty=c(1,2, 3, 4), lwd=2,
 legend=c("30 pct reduction", "40 pct reduction",
                "50 pct reduction", "60 pct reduction"))
#line 551 "R4SampleSize.rnw"


###################################################
### chunk number 30: SS.nSurv1
###################################################
#line 629 "R4SampleSize.rnw"
# calculate sample size (denoted by ss1)
ss1 = nSurvival(lambda1=.3,lambda2=.2,
	eta =.1,Ts=3,Tr=1,sided=1,alpha=.025,beta=0.2)
# print the required ss
print(ss1)


###################################################
### chunk number 31: 
###################################################
#line 653 "R4SampleSize.rnw"
# call nSurvival to calculate sample size
ss2 =  nSurvival(lambda1=-log(.5) / 6, 
	lambda2=-log(.5) / 6 *.7, eta=-log(.95)/12, 
	Tr=30 ,Ts=36, type="rr", entry="unif")
# print the required sample size
print(ss2$Sample.size)
# print the required number of events
print(ss2$Num.events)


###################################################
### chunk number 32: 
###################################################
#line 677 "R4SampleSize.rnw"
library(gsDesign)


###################################################
### chunk number 33: 
###################################################
#line 715 "R4SampleSize.rnw"
# load the default gsDesign
x = gsDesign()
# print the output
x


###################################################
### chunk number 34: 
###################################################
#line 724 "R4SampleSize.rnw"
x$n.I


###################################################
### chunk number 35: 
###################################################
#line 730 "R4SampleSize.rnw"
# calculate the sample size for non-group sequential design
n.fix = nBinomial(p1=0.3, p2=0.15)
# print the fixed sample size
n.fix


###################################################
### chunk number 36: 
###################################################
#line 737 "R4SampleSize.rnw"
n.GS = n.fix*x$n.I
n.GS


###################################################
### chunk number 37: 
###################################################
#line 743 "R4SampleSize.rnw"
x.new = gsDesign(n.fix=n.fix)
x.new


###################################################
### chunk number 38: SS.fig4gsDesign
###################################################
#line 752 "R4SampleSize.rnw"
print(plot(x, plottype=1))


###################################################
### chunk number 39: 
###################################################
#line 758 "R4SampleSize.rnw"
#line 752 "R4SampleSize.rnw#from line#758#"
print(plot(x, plottype=1))
#line 759 "R4SampleSize.rnw"


###################################################
### chunk number 40: SS.fig4gsDesignspending
###################################################
#line 766 "R4SampleSize.rnw"
print(plot(x, plottype=5))


###################################################
### chunk number 41: 
###################################################
#line 772 "R4SampleSize.rnw"
#line 766 "R4SampleSize.rnw#from line#772#"
print(plot(x, plottype=5))
#line 773 "R4SampleSize.rnw"


###################################################
### chunk number 42: 
###################################################
#line 873 "R4SampleSize.rnw"
# minimum treatment effect
delta       = 1.2
# number of repeated measurements
n           = 5
# trial duration
tau         = 2
# within-subject variability
sig2within  = 7
# between-subject variability
sig2between = 2
# significance level
alpha       = 0.05
# desired power
pow         = 0.9
# variance of slopes
sig2        = 12*(n-1)*sig2within/(tau^2*n*(n+1))+sig2between
print(sig2)
# calculate sample size
N           = (qnorm(1-alpha/2)+qnorm(pow))^2*2*sig2/delta^2
cat("Sample size needed=",N, sep=" ","\n\n")


###################################################
### chunk number 43: 
###################################################
#line 900 "R4SampleSize.rnw"
# sample size for each treatment
N           = 40
delta       = 1.2
n           = 2
tau         = 2
sig2within  = 7
sig2between = 2
alpha       = 0.05
sig2        = 12*(n-1)*sig2within/(tau^2*n*(n+1))+sig2between
pow.get     = pnorm(sqrt( N*delta^2/(2*sig2))-qnorm(1-alpha/2))
cat("Power obtained=", pow.get,sep=" ", "\n\n")


###################################################
### chunk number 44: 
###################################################
#line 915 "R4SampleSize.rnw"
pow.Long =  
  function(N, n, delta, tau, sig2within, sig2between, alpha){
sig2     = 12*(n-1)*sig2within/(tau^2*n*(n+1))+sig2between
pow.get  = pnorm(sqrt( N*delta^2/(2*sig2))-qnorm(1-alpha/2))
pow.get
}


###################################################
### chunk number 45: 
###################################################
#line 925 "R4SampleSize.rnw"
# the sample size inputs
N = seq(20, 100, by=20)
# the number of repeated measurements
n = seq(2,10, by=2)
# power matrix
pow = matrix(0, ncol=length(N), nrow=length(n))
colnames(pow) = n
rownames(pow) = N
# loop to calculate the power
for (i in 1:length(N)){
	for(j in 1:length(n)){
	pow[i,j] = pow.Long(N[i], n[j], delta, tau, 
			sig2within, sig2between, alpha)
				}	
			}
# print the power matrix
pow


###################################################
### chunk number 46: SS.fig4longpow
###################################################
#line 947 "R4SampleSize.rnw"
persp(N, n, pow, theta = 30, phi = 30, expand = 0.5, 
col = "lightblue", xlab="Sample Size (N)",
ylab="# of Measurements (n)",
zlab="Power") 


###################################################
### chunk number 47: 
###################################################
#line 956 "R4SampleSize.rnw"
#line 947 "R4SampleSize.rnw#from line#956#"
persp(N, n, pow, theta = 30, phi = 30, expand = 0.5, 
col = "lightblue", xlab="Sample Size (N)",
ylab="# of Measurements (n)",
zlab="Power") 
#line 957 "R4SampleSize.rnw"


###################################################
### chunk number 48: 
###################################################
#line 1016 "R4SampleSize.rnw"
# the treatment effect
delta  = 0.075
# number of repeated measurements
n       = 5
# study duration
tau    = 2
# correlation
rho    =  0.5
# significance level
alpha = 0.05
# Type-II error
beta  = 0.1
# associated power
pow  = 1-beta
# sigma calculation
sig2 = 3*(n-1)*(1-rho)/(tau^2*n*(n+1))
sig2
# sample size calculation
N = (qnorm(1-alpha/2)+qnorm(pow))^2*2*sig2/delta^2
cat("Sample size needed=",N, sep=" ","\n\n")


###################################################
### chunk number 49: 
###################################################
#line 1082 "R4SampleSize.rnw"
# CV range
CV =seq(0.1, 0.7, by=0.1)
CV
# PC range
PC =  seq(0.1, 0.5, by=0.05)
PC
# Sample size calculation
SS = matrix(0, ncol=length(PC), nrow=length(CV))
colnames(SS) = PC
rownames(SS)= CV
# for-loop to calculate the sample size for each PC and CV combination
for (i in 1:length(CV)){
	for(j in 1:length(PC)){
	SS[i,j] = ceiling(8*CV[i]^2/PC[j]^2*(1+(1-PC[j])^2))
					}
				}
# print out the table
SS


###################################################
### chunk number 50: SS.fig4PC2CV
###################################################
#line 1105 "R4SampleSize.rnw"
persp(CV, PC, SS, theta = 30, phi = 30, expand = 0.5, 
col = "lightblue", xlab="CV", ylab="PC", zlab="Sample Size") 


###################################################
### chunk number 51: 
###################################################
#line 1113 "R4SampleSize.rnw"
#line 1105 "R4SampleSize.rnw#from line#1113#"
persp(CV, PC, SS, theta = 30, phi = 30, expand = 0.5, 
col = "lightblue", xlab="CV", ylab="PC", zlab="Sample Size") 
#line 1114 "R4SampleSize.rnw"


