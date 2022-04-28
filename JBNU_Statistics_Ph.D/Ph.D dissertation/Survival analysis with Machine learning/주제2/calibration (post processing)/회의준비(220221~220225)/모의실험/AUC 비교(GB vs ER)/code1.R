library(survival)

uis2 <- read.csv("C:/Users/HSY/Desktop/uis2.csv", header=T, sep=",")

cox0 <- coxph(Surv(time, censor) ~ age+beck+hercoc+ivhx+ndrugtx+race+treat+lot,
              data=uis2)
cox1 <- coxph(Surv(time, censor) ~ age+beck+hercoc+ivhx+ndrugtx+race+treat+lot, 
              data=uis2, iter.max=0, init=c(cox0$coef,rep(0,11-1)))
I.inv <- vcov(cox1)
Sigma.hat <- solve(I.inv[c(6,7,8),c(6,7,8)])

library(survMisc)

gof(cox0, G=11) # GB test


# Kaplan-Meier estimator
library(survival)
uis2 <- read.csv("C:/Users/HSY/Desktop/uis2.csv", header=T, sep=",")
km <- with(uis2, Surv(time, censor))
km
km_fit <- survfit(Surv(time, censor) ~ 1, data=uis2)
summary(km_fit)
km_fit[["cumhaz"]] # cumulative hazard function
km_fit[["surv"]] # survival function
summary(km_fit, times = c(1,5*(1:20)))

# time dependent ROC curve
library(timeROC)
