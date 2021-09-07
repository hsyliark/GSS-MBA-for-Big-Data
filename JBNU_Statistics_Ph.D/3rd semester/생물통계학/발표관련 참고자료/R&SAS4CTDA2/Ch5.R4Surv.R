### R code from vignette source 'R4Surv2.rnw'


###################################################
### code chunk number 2: Surv.dat1
###################################################
dat        = read.csv("CTCarcinoma.csv", header=T)
# show first 6 observations using R function "head"
head(dat)


###################################################
### code chunk number 3: Surv.km
###################################################
# load the R library
library(survival)
# fit Kaplan-Meier 
fit.km = survfit(Surv(Time,Status==0)~TRT,
		type=c("kaplan-meier"),dat)
# print the model fitting
fit.km


###################################################
### code chunk number 4: Surv.fig4dat1
###################################################
plot(fit.km, lty=c(1,4,8),lwd=2,xlab="Time in Weeks",ylab="S(t)")
legend("bottomleft",title="Line Types",
c("S+CT","S+CT+IT","S+IT"),lty=c(1,4,8), lwd=2) 


###################################################
### code chunk number 5: R4Surv2.rnw:551-552
###################################################
plot(fit.km, lty=c(1,4,8),lwd=2,xlab="Time in Weeks",ylab="S(t)")
legend("bottomleft",title="Line Types",
c("S+CT","S+CT+IT","S+IT"),lty=c(1,4,8), lwd=2) 


###################################################
### code chunk number 6: Surv.fig4hazard1
###################################################
fit.fleming =survfit(Surv(Time,Status==0)~TRT,
		dat,type='fleming') 
plot(fit.fleming,lty=c(1,4,8),lwd=2,fun="cumhaz", 
xlab="Time in Weeks", ylab="Cumulative Hazard") 
legend("topleft",title="Line Types",
c("S+CT","S+CT+IT","S+IT"),lty=c(1,4,8),lwd=2) 


###################################################
### code chunk number 7: R4Surv2.rnw:613-614
###################################################
fit.fleming =survfit(Surv(Time,Status==0)~TRT,
		dat,type='fleming') 
plot(fit.fleming,lty=c(1,4,8),lwd=2,fun="cumhaz", 
xlab="Time in Weeks", ylab="Cumulative Hazard") 
legend("topleft",title="Line Types",
c("S+CT","S+CT+IT","S+IT"),lty=c(1,4,8),lwd=2) 


###################################################
### code chunk number 8: R4Surv2.rnw:620-623
###################################################
# use "survdiff" to test difference
fit.diff = survdiff(Surv(Time, Status==0)~TRT,data=dat)
fit.diff


###################################################
### code chunk number 9: R4Surv2.rnw:634-641
###################################################
# fit exponential model
fit.exp =survreg(Surv(Time, Status==0)~TRT,dat,
		dist="exponential")
summary(fit.exp)
# fit Weibull model
fit.Weibull =survreg(Surv(Time, Status==0)~TRT,dat)
summary(fit.Weibull)


###################################################
### code chunk number 10: R4Surv2.rnw:669-676
###################################################
# fit exponential model +Age
fit.exp.age = survreg(Surv(Time, Status==0)~TRT+Age,
		dat,dist="exponential")
summary(fit.exp.age)
# fit Weibull model+Age
fit.Weibull.age = survreg(Surv(Time,Status==0)~TRT+Age,dat)
summary(fit.Weibull.age)


###################################################
### code chunk number 11: R4Surv2.rnw:684-690
###################################################
# fit Cox
fit.Cox = coxph(Surv(Time, Status==0)~TRT,dat)
summary(fit.Cox)
# fit Cox +Age
fit.Cox.age = coxph(Surv(Time, Status==0)~TRT+Age,dat)
summary(fit.Cox.age)


###################################################
### code chunk number 12: R4Surv2.rnw:705-707
###################################################
dat        = read.csv("BreastCancer.csv", header=T)
head(dat)


###################################################
### code chunk number 13: R4Surv2.rnw:720-726
###################################################
cria.tau = function(data){
l   = data$tL;r = data$tU
# sort all the time points
tau = sort(unique(c(l,r[is.finite(r)])))
return(tau)
}


###################################################
### code chunk number 14: R4Surv2.rnw:730-741
###################################################
S.ini = function(tau){
# take the ordered time
m    = length(tau)
# fit the Kaplan-Meier
ekm  = survfit(Surv(tau[1:m-1],rep(1,m-1))~1)
# Output the estimated Survival
So   = c(1,ekm$surv)
# estimate the step
p    = -diff(So)
return(p)
}


###################################################
### code chunk number 15: R4Surv2.rnw:745-778
###################################################
cria.A = function(data,tau){
tau12  = cbind(tau[-length(tau)],tau[-1])
interv = function(x,inf,sup) 
		ifelse(x[1]>=inf & x[2]<=sup,1,0)
A      = apply(tau12,1,interv,inf=data$tL,sup=data$tU)
id.lin.zero = which(apply(A==0, 1, all))
if(length(id.lin.zero)>0) A = A[-id.lin.zero, ]
return(A)
}
# Turnbull function
Turnbull = function(p, A, data, eps=1e-3, 
	iter.max=200, verbose=FALSE){
n =nrow(A);m=ncol(A);Q=matrix(1,m)
iter    = 0
repeat {
iter = iter + 1; diff = (Q-p)
maxdiff = max(abs(as.vector(diff)))
if (verbose) print(maxdiff)
if (maxdiff<eps | iter>=iter.max) break
Q  = p; C  = A%*%p; p=p*((t(A)%*%(1/C))/n)
}
cat("Iterations = ", iter,"\n")
cat("Max difference = ", maxdiff,"\n")
cat("Convergence criteria: Max difference < 1e-3","\n")
dimnames(p) = list(NULL,c("P Estimate"))
surv        = round(c(1,1-cumsum(p)),digits=5)
right       = data$tU
if(any(!(is.finite(right)))){
t = max(right[is.finite(right)])
return(list(time=tau[tau<t],surv=surv[tau<t]))
}
else return(list(time=tau,surv=surv))
}


###################################################
### code chunk number 16: R4Surv2.rnw:782-794
###################################################
# get the data for TRT=1
dat1 = dat[dat$TRT==1,]
dat1$tU[is.na(dat1$tU)] = Inf
# sort the time points
tau  = cria.tau(dat1)
# Estimate the initial Survival
p    = S.ini(tau=tau)
# run Turnbull and name it as "mb1"
A    = cria.A(data=dat1,tau=tau)
mb1  = Turnbull(p,A,dat1)
# print the estimates
mb1


###################################################
### code chunk number 17: R4Surv2.rnw:799-806
###################################################
dat0  = dat[dat$TRT==0,]
dat0$tU[is.na(dat0$tU)] = Inf
tau   = cria.tau(dat0)
p     = S.ini(tau=tau)
A     = cria.A(data=dat0,tau=tau)
mb0   = Turnbull(p,A,dat0)
mb0


###################################################
### code chunk number 18: figSurv.Turnbull1
###################################################
# plot the TRT=1
plot(mb1$time,mb1$surv,lty=1,lwd=2,type="s", ylim=c(0,1),
xlim=range(c(0,60)),xlab="Time in Months",ylab="S(t)")
# add a line for TRT=0
lines(mb0$time,mb0$surv,lty=4,lwd=2,type="s")
# put a legend
legend("topright",title="Line Types",lty=c(1,4),lwd=2,
c("Radiation Only","Radiation+Chemotherapy"))


###################################################
### code chunk number 19: R4Surv2.rnw:826-827
###################################################
# plot the TRT=1
plot(mb1$time,mb1$surv,lty=1,lwd=2,type="s", ylim=c(0,1),
xlim=range(c(0,60)),xlab="Time in Months",ylab="S(t)")
# add a line for TRT=0
lines(mb0$time,mb0$surv,lty=4,lwd=2,type="s")
# put a legend
legend("topright",title="Line Types",lty=c(1,4),lwd=2,
c("Radiation Only","Radiation+Chemotherapy"))


###################################################
### code chunk number 20: Surv.midpoint
###################################################
# get the midpoint
time   = dat$tL+((dat$tU-dat$tL)/2)
# get the censorship
Status = ifelse(is.finite(time),1,0)
# replace the NA with left-time
time   = ifelse(is.finite(time),time,dat$tL)
# fit Kaplan-Meier model
ekm    = survfit(Surv(time, Status)~TRT,
		type=c("kaplan-meier"),dat)
# print the output
ekm


###################################################
### code chunk number 21: figSurv.midpoint
###################################################
plot(mb1$time,mb1$surv,lty=1,lwd=3,type="s",ylim=c(0,1),
xlim=range(c(0,50)), xlab="Time in Months",ylab="S(t)")
legend("bottomleft",title="Line Types",lty=c(1,4),lwd=3,
c("Radiotherapy Only","Radiotherapy+Chemotherapy"))
lines(mb0$time,mb0$surv,lty=4,lwd=3,type="s")
lines(ekm[1]$time,ekm[1]$surv,type="s",lty=4,lwd=1)
lines(ekm[2]$time,ekm[2]$surv,type="s",lty=1,lwd=1)


###################################################
### code chunk number 22: R4Surv2.rnw:866-867
###################################################
plot(mb1$time,mb1$surv,lty=1,lwd=3,type="s",ylim=c(0,1),
xlim=range(c(0,50)), xlab="Time in Months",ylab="S(t)")
legend("bottomleft",title="Line Types",lty=c(1,4),lwd=3,
c("Radiotherapy Only","Radiotherapy+Chemotherapy"))
lines(mb0$time,mb0$surv,lty=4,lwd=3,type="s")
lines(ekm[1]$time,ekm[1]$surv,type="s",lty=4,lwd=1)
lines(ekm[2]$time,ekm[2]$surv,type="s",lty=1,lwd=1)


###################################################
### code chunk number 23: R4Surv2.rnw:883-889
###################################################
# load the library
library(interval)
# fit the NPMLE by calling "icfit" and name it as "fit
fit.Int = icfit(Surv(tL,tU,type="interval2")~TRT,data=dat)
# print the summary
summary(fit.Int)


###################################################
### code chunk number 24: figSurv.Int
###################################################
plot(fit.Int, XLAB="Time in Months", YLAB="Survival Probability",
 col=c("lightblue","pink"), LEGEND=F,estpar=list(col=c("blue","red"),lwd=3,lty=1))
legend("bottomleft", legend=c("Radiotherapy Only","Radiotherapy+Chemotherapy"),
   col=c("red","blue"),lwd=2)


###################################################
### code chunk number 25: R4Surv2.rnw:907-908
###################################################
plot(fit.Int, XLAB="Time in Months", YLAB="Survival Probability",
 col=c("lightblue","pink"), LEGEND=F,estpar=list(col=c("blue","red"),lwd=3,lty=1))
legend("bottomleft", legend=c("Radiotherapy Only","Radiotherapy+Chemotherapy"),
   col=c("red","blue"),lwd=2)


###################################################
### code chunk number 26: R4Surv2.rnw:921-927
###################################################
# create a new dataset "breast" and keep "dat" for future use
breast = dat
# replace 0's with NA as left-censored
breast$tL[breast$tL==0]= NA
# print the first a few observations
head(breast)


###################################################
### code chunk number 27: R4Surv2.rnw:932-939
###################################################
# fit exponential
fit.exp=survreg(Surv(tL,tU,type="interval2")~TRT,
		breast,dist="exponential")
summary(fit.exp)
# fit Weibull
fit.Weibull=survreg(Surv(tL,tU,type="interval2")~TRT,data=breast)
summary(fit.Weibull)


###################################################
### code chunk number 28: R4Surv2.rnw:959-963
###################################################
# call "ictest" to test "TRT" difference
test.Int = ictest(Surv(tL,tU,type="interval2")~TRT,data=dat)
# print the summary
test.Int


