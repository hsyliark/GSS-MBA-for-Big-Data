### R code from vignette source 'R4Bayes.rnw'


###################################################
### code chunk number 2: R4Bayes.rnw:292-293
###################################################
library(R2WinBUGS)


###################################################
### code chunk number 3: R4Bayes.rnw:297-298
###################################################
library(help=R2WinBUGS)


###################################################
### code chunk number 4: rbugs
###################################################
library(rbugs)


###################################################
### code chunk number 5: pack
###################################################
library(MCMCpack)


###################################################
### code chunk number 6: R4Bayes.rnw:365-379
###################################################
# set seed to 123
set.seed(123)
# n=30 patients
n        = 30
# known sigma
sigma  = 2
# population mean
mu      = 3
# simulate data
y         =  rnorm(n, mu,sigma)
# print the data
# the mean and variance of the simulated data
mean(y)
var(y)


###################################################
### code chunk number 7: R4Bayes.rnw:383-393
###################################################
# the prior parameters
mu0 = 2; tau0 = .5 
# the weight
w           =  tau0^2/(tau0^2+sigma^2/n)
# the posterior mean
muP       =  w*mean(y) + (1-w)*mu0
# the posterior standard deviation
sigmaP = sqrt(1/(1/tau0^2+n/sigma^2))
# direct simulation of posterior normal
Bayes1.norm2norm = rnorm(10000, muP,sigmaP)


###################################################
### code chunk number 8: R4Bayes.rnw:397-398
###################################################
quantile(Bayes1.norm2norm, c(0.025,0.25,0.5,0.75,0.975))


###################################################
### code chunk number 9: R4Bayes.rnw:402-406
###################################################
# call the function
Bayes2.norm2norm = MCnormalnormal(y, sigma^2, mu0, tau0^2, 10000)
# print the summary
summary(Bayes2.norm2norm)


###################################################
### code chunk number 10: Bayes.fig4norm2norm
###################################################
x     = seq(0,6,0.01)
plot(x, dnorm(x, mu0,tau0), type="l", lwd=1,las=1, ylim=c(0,1.4),
   xlab="mu", ylab="density")
lines(x, dnorm(x, mean(y), sigma/sqrt(n)), lty=8, lwd=1)
lines(density(Bayes1.norm2norm), lty=8, lwd=3)
legend("topright", c("Prior","Likelihood", "Posterior"), lwd=c(1,1,3), 
lty=c(1,8,4))


###################################################
### code chunk number 11: R4Bayes.rnw:421-422
###################################################
x     = seq(0,6,0.01)
plot(x, dnorm(x, mu0,tau0), type="l", lwd=1,las=1, ylim=c(0,1.4),
   xlab="mu", ylab="density")
lines(x, dnorm(x, mean(y), sigma/sqrt(n)), lty=8, lwd=1)
lines(density(Bayes1.norm2norm), lty=8, lwd=3)
legend("topright", c("Prior","Likelihood", "Posterior"), lwd=c(1,1,3), 
lty=c(1,8,4))


###################################################
### code chunk number 12: R4Bayes.rnw:434-441
###################################################
# total patients for each treatment
n  = c(168, 182, 165,188)
# number healed
x = c(69, 113, 120, 145)
# the observed proportion
p = x/n
p


###################################################
### code chunk number 13: R4Bayes.rnw:454-459
###################################################
# the objective function
obj = function(parm){
a = parm[1]; b = parm[2]
( pbeta(0.50,a,b) -0.75)^2 +( pbeta(0.95,a,b)- 0.85)^2
}


###################################################
### code chunk number 14: R4Bayes.rnw:463-466
###################################################
# call optim to search the root with initial values at (3,3)
out = optim(c(3,3), obj)
print(out)


###################################################
### code chunk number 15: R4Bayes.rnw:471-473
###################################################
pbeta(0.5,out$par[1], out$par[2])
pbeta(0.95,out$par[1], out$par[2])


###################################################
### code chunk number 16: R4Bayes.rnw:484-488
###################################################
# direct simulation
Bayes1.betabin = rbeta(10000, 120.062, 45.183)
# print the quantiles
quantile(Bayes1.betabin, c(0.025,0.25,0.5,0.75,0.975))


###################################################
### code chunk number 17: R4Bayes.rnw:494-500
###################################################
# keep the parameters
x3 = 120; n3 =165; a = 0.062; b=0.183
# call the MCbinomialbeta function for 10000 simulation 
Bayes2.betabin = MCbinomialbeta(x3, n3, a, b, mc=10000)
# print the summary
summary(Bayes2.betabin)


###################################################
### code chunk number 18: Bayes.fig4betabin
###################################################
p = seq(0,1,0.01)
plot(p, dbeta(p,a, b),  lwd=3, type="l",ylim=c(0,13),
  xlab="Healing Rate", ylab="density")
lines(density(Bayes1.betabin), lty=4, lwd=3)
lines(density(Bayes2.betabin), lty=8, lwd=3)
legend("topleft", c("Prior", "Direct Simulation","MCbinomialbeta"), 
lwd=3, lty=c(1,4,8))


###################################################
### code chunk number 19: R4Bayes.rnw:518-519
###################################################
p = seq(0,1,0.01)
plot(p, dbeta(p,a, b),  lwd=3, type="l",ylim=c(0,13),
  xlab="Healing Rate", ylab="density")
lines(density(Bayes1.betabin), lty=4, lwd=3)
lines(density(Bayes2.betabin), lty=8, lwd=3)
legend("topleft", c("Prior", "Direct Simulation","MCbinomialbeta"), 
lwd=3, lty=c(1,4,8))


###################################################
### code chunk number 20: R4Bayes.rnw:538-541
###################################################
dat = read.csv("DBP.csv",header=T)
# we are interested in the blood pressure change
dat$diff = dat$DBP5-dat$DBP1


###################################################
### code chunk number 21: bayes.lm
###################################################
# fit the Bayes regression model with 1000 burn-in
BayesMod  = MCMCregress(diff~TRT+Age, dat)
# print the MCMC result
summary(BayesMod)


###################################################
### code chunk number 22: Bayes.fig4lm
###################################################
# make the margin
par(mar = c(3,2,1.5,1))
# make the MCMC plot
plot(BayesMod)


###################################################
### code chunk number 23: R4Bayes.rnw:580-581
###################################################
# make the margin
par(mar = c(3,2,1.5,1))
# make the MCMC plot
plot(BayesMod)


###################################################
### code chunk number 24: R4Bayes.rnw:589-592
###################################################
library(flexmix)
data("betablocker")  
betablocker$Center = as.factor(betablocker$Center)


###################################################
### code chunk number 25: center1
###################################################
# extract center 1
beta1 = betablocker[betablocker$Center == 1,
			c("Deaths","Total","Treatment")]
# print the center 1 data
beta1 


###################################################
### code chunk number 26: R4Bayes.rnw:606-611
###################################################
# make a dataframe
beta1 = data.frame(trt = c(rep("TRT", 38),rep("Cont",39)),
	death = c(rep(1,3), rep(0,38-3), rep(1,3), rep(0,39-3))) 
# print the first 6 observations
head(beta1)


###################################################
### code chunk number 27: R4Bayes.rnw:617-621
###################################################
# fit logistic regression
glm.beta = glm(death ~trt,binomial,beta1)
# print the result
summary(glm.beta)


###################################################
### code chunk number 28: R4Bayes.rnw:633-637
###################################################
 ## Call MCMClogit with default 
Bayes1.beta = MCMClogit(death~trt, data=beta1) 
# print the summary for MCMC
summary(Bayes1.beta)


###################################################
### code chunk number 29: Bayes.fig4glm
###################################################
plot(Bayes1.beta)


###################################################
### code chunk number 30: R4Bayes.rnw:648-649
###################################################
plot(Bayes1.beta)


###################################################
### code chunk number 31: R4Bayes.rnw:657-661
###################################################
# fit the Bayesian logistic regression with multivariate normal prior
Bayes2.beta = MCMClogit(death~trt, B0=.001,data=beta1) 
# print the fit
summary(Bayes2.beta)


###################################################
### code chunk number 32: R4Bayes.rnw:669-671
###################################################
library(HSAUR)
data(polyps)


###################################################
### code chunk number 33: R4Bayes.rnw:683-687
###################################################
 ## Call MCMCpoissont with default 
Bayes.polyps <- MCMCpoisson(number ~ treat+age, polyps)
# print the summary for MCMC
summary(Bayes.polyps)


###################################################
### code chunk number 34: Bayes.fig4poisson
###################################################
# set a beta margin for plotting
par(mar = c(3,2,1.5,1))
# plot the MCMC
plot(Bayes.polyps)


###################################################
### code chunk number 35: R4Bayes.rnw:703-704
###################################################
# set a beta margin for plotting
par(mar = c(3,2,1.5,1))
# plot the MCMC
plot(Bayes.polyps)


###################################################
### code chunk number 36: R4Bayes.rnw:735-744
###################################################
# call the hypergeo library
library(hypergeo)
# make a function call 
pXgtY = function(sx,tx,sy,ty){
tmp1 = beta(sx+sy,tx+ty)/(beta(sx,tx)*beta(sy,ty)*sy)
tmp2 = genhypergeo(U=c(sx+sy,sy+ty,1),L=c(sy+1,sx+tx+sy+ty), 
	check_mod=FALSE,z=1)
tmp1*tmp2
}


###################################################
### code chunk number 37: R4Bayes.rnw:749-761
###################################################
# compare 800 mg C to 400 mg C
p800to400 = pXgtY(x[3], n[3]-x[3], x[2],n[2]-x[2])
p800to400
# compare 800 mg C to 0 mg C
p800to0 = pXgtY(x[3], n[3]-x[3], x[1],n[1]-x[1])
p800to0
# compare 1600 mg C to 800 mg C
p1600to800 = pXgtY(x[4], n[4]-x[4], x[3],n[3]-x[3])
p1600to800
# compare 1600 mg C to 400 mg C
p1600to400 = pXgtY(x[4], n[4]-x[4], x[2],n[2]-x[2])
p1600to400


###################################################
### code chunk number 38: Bayes.fig4betabin2
###################################################
#  make p from 0.3 to 0.9 by 100 points 
n.pts =100
p      = seq(0.3,0.9,length=n.pts)
# the prior and the distributions from 4 treatments
pb0  = dbeta(p,1,1)
pb1  = dbeta(p,x[1],n[1]-x[1])
pb2  = dbeta(p,x[2],n[2]-x[2])
pb3  = dbeta(p,x[3],n[3]-x[3])
pb4  = dbeta(p,x[4],n[4]-x[4])

# the maximum to set the yaxis limit
ymax= max(pb1,pb2,pb3,pb4)
# plot the prior and posterior
plot(p,pb0, lwd=1, lty=8, las=1,type="l",xlab="Healing Probability",
 ylab="Density",ylim=c(0,ymax))
lines(p,pb1, lwd=3, lty=1)
lines(p,pb2, lwd=3, lty=2)
lines(p,pb3, lwd=3, lty=3)
lines(p,pb4, lwd=3, lty=4)
legend("topleft", c("0 mg C","400 mg C", "800 mg C","1600 mg C"), 
lwd=c(3,3,3,3), lty=c(1,2,3,4))


###################################################
### code chunk number 39: R4Bayes.rnw:795-796
###################################################
#  make p from 0.3 to 0.9 by 100 points 
n.pts =100
p      = seq(0.3,0.9,length=n.pts)
# the prior and the distributions from 4 treatments
pb0  = dbeta(p,1,1)
pb1  = dbeta(p,x[1],n[1]-x[1])
pb2  = dbeta(p,x[2],n[2]-x[2])
pb3  = dbeta(p,x[3],n[3]-x[3])
pb4  = dbeta(p,x[4],n[4]-x[4])

# the maximum to set the yaxis limit
ymax= max(pb1,pb2,pb3,pb4)
# plot the prior and posterior
plot(p,pb0, lwd=1, lty=8, las=1,type="l",xlab="Healing Probability",
 ylab="Density",ylim=c(0,ymax))
lines(p,pb1, lwd=3, lty=1)
lines(p,pb2, lwd=3, lty=2)
lines(p,pb3, lwd=3, lty=3)
lines(p,pb4, lwd=3, lty=4)
legend("topleft", c("0 mg C","400 mg C", "800 mg C","1600 mg C"), 
lwd=c(3,3,3,3), lty=c(1,2,3,4))


