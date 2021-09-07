### R code from vignette source 'R4Meta.rnw'


###################################################
### code chunk number 2: R4Meta.rnw:70-72
###################################################
library(rmeta)
data(cochrane)


###################################################
### code chunk number 3: R4Meta.rnw:75-79
###################################################
library(xtable)
print(xtable(cochrane, digits=0, align=rep("c", 1+ncol(cochrane)),
caption="Data for Cochrane Collaboration Logo.",label = "meta.data.steriod"),
table.placement = "htbp",caption.placement = "top",include.rownames=FALSE)


###################################################
### code chunk number 4: R4Meta.rnw:298-301
###################################################
library(flexmix)
data("betablocker")
beta=betablocker


###################################################
### code chunk number 5: rmeta
###################################################
library(rmeta)


###################################################
### code chunk number 6: R4Meta.rnw:316-326
###################################################
# Get the data from the ``beta"
n.trt   = betablocker[betablocker$Treatment=="Treated",]$Total
n.ctrl  = betablocker[betablocker$Treatment=="Control",]$Total
ev.trt  = betablocker[betablocker$Treatment=="Treated",]$Deaths
ev.ctrl = betablocker[betablocker$Treatment=="Control",]$Deaths
# call the meta.MH for calculations
betaOR  = meta.MH(n.trt,n.ctrl,ev.trt,ev.ctrl,
	names=paste("Center",1:22,sep=" "))
# print the summary from meta.MH
summary(betaOR)


###################################################
### code chunk number 7: meta.beta.fig4fixed
###################################################
plot(betaOR, ylab="Center")


###################################################
### code chunk number 8: R4Meta.rnw:341-342
###################################################
plot(betaOR, ylab="Center")


###################################################
### code chunk number 9: meta.beta.fig4fixed.funnel
###################################################
funnelplot(betaOR)


###################################################
### code chunk number 10: R4Meta.rnw:358-359
###################################################
funnelplot(betaOR)


###################################################
### code chunk number 11: R4Meta.rnw:370-375
###################################################
# Call the meta.DSL for calculations
betaDSL  = meta.DSL(n.trt,n.ctrl,ev.trt,ev.ctrl,
	names=paste("Center",1:22,sep=" "))
# Print the summary from meta.DSL
summary(betaDSL)


###################################################
### code chunk number 12: meta.beta.fig4random
###################################################
plot(betaDSL, ylab="Center")


###################################################
### code chunk number 13: R4Meta.rnw:389-390
###################################################
plot(betaDSL, ylab="Center")


###################################################
### code chunk number 14: meta.beta.forest
###################################################
# Create the ``text" to include all the outputs
text = cbind(c("","Center",betaOR$names,NA,"Summary"),
             c("Deaths","(Betablockers)",ev.trt,NA,NA),
             c("Deaths","(Placebo)", ev.ctrl, NA,NA),
	     c("","OR",format(exp(betaOR$logOR),digits=2),
		    NA,format(exp(betaOR$logMH),digits=2)))
# Generate the OR and 95\% CI
mean  = c(NA,NA,betaOR$logOR,NA,betaOR$logMH)
sterr = c(NA,NA,betaOR$selogOR,NA,betaOR$selogMH)
l     = mean-1.96*sterr
u     = mean+1.96*sterr
# Call forestplot with a few options
forestplot(text,mean,l,u,zero=0,is.summary=c(TRUE,TRUE,rep(FALSE,22),TRUE),
   clip=c(log(0.1),log(2.5)), xlog=TRUE)


###################################################
### code chunk number 15: R4Meta.rnw:420-421
###################################################
# Create the ``text" to include all the outputs
text = cbind(c("","Center",betaOR$names,NA,"Summary"),
             c("Deaths","(Betablockers)",ev.trt,NA,NA),
             c("Deaths","(Placebo)", ev.ctrl, NA,NA),
	     c("","OR",format(exp(betaOR$logOR),digits=2),
		    NA,format(exp(betaOR$logMH),digits=2)))
# Generate the OR and 95\% CI
mean  = c(NA,NA,betaOR$logOR,NA,betaOR$logMH)
sterr = c(NA,NA,betaOR$selogOR,NA,betaOR$selogMH)
l     = mean-1.96*sterr
u     = mean+1.96*sterr
# Call forestplot with a few options
forestplot(text,mean,l,u,zero=0,is.summary=c(TRUE,TRUE,rep(FALSE,22),TRUE),
   clip=c(log(0.1),log(2.5)), xlog=TRUE)


###################################################
### code chunk number 16: R4Meta.rnw:432-436
###################################################
# Load the data
data(cochrane)
# print it
cochrane


###################################################
### code chunk number 17: R4Meta.rnw:441-446
###################################################
# Fit the fixed-effects model
steroid = meta.MH(n.trt, n.ctrl, ev.trt, ev.ctrl,
                        names=name, data=cochrane)
# Print the model fit
summary(steroid)


###################################################
### code chunk number 18: meta.beta.cochrane
###################################################
# Create the ``tabletext" to include all the outputs
tabletext = cbind(c("","Study",steroid$names,NA,"Summary"),
                 c("Deaths","(Steroid)",cochrane$ev.trt,NA,NA),
                 c("Deaths","(Placebo)",cochrane$ev.ctrl, NA,NA),
		 c("","OR",format(exp(steroid$logOR),digits=2),
			NA,format(exp(steroid$logMH),digits=2)))
# Generate the CI
mean   = c(NA,NA,steroid$logOR,NA,steroid$logMH)
stderr = c(NA,NA,steroid$selogOR,NA,steroid$selogMH)
l      = mean-1.96*stderr
u      = mean+1.96*stderr
# Call forestplot
forestplot(tabletext,mean,l,u,zero=0,
	is.summary=c(TRUE,TRUE,rep(FALSE,8),TRUE),
   	clip=c(log(0.1),log(2.5)), xlog=TRUE)


###################################################
### code chunk number 19: R4Meta.rnw:474-475
###################################################
# Create the ``tabletext" to include all the outputs
tabletext = cbind(c("","Study",steroid$names,NA,"Summary"),
                 c("Deaths","(Steroid)",cochrane$ev.trt,NA,NA),
                 c("Deaths","(Placebo)",cochrane$ev.ctrl, NA,NA),
		 c("","OR",format(exp(steroid$logOR),digits=2),
			NA,format(exp(steroid$logMH),digits=2)))
# Generate the CI
mean   = c(NA,NA,steroid$logOR,NA,steroid$logMH)
stderr = c(NA,NA,steroid$selogOR,NA,steroid$selogMH)
l      = mean-1.96*stderr
u      = mean+1.96*stderr
# Call forestplot
forestplot(tabletext,mean,l,u,zero=0,
	is.summary=c(TRUE,TRUE,rep(FALSE,8),TRUE),
   	clip=c(log(0.1),log(2.5)), xlog=TRUE)


###################################################
### code chunk number 20: R4Meta.rnw:488-489
###################################################
library(meta)


###################################################
### code chunk number 21: R4Meta.rnw:494-495
###################################################
library(help=meta)


###################################################
### code chunk number 22: R4Meta.rnw:503-506
###################################################
angina  = read.csv("angina.csv",header=T)
# Print the data
angina


###################################################
### code chunk number 23: R4Meta.rnw:515-520
###################################################
# Fit fixed-effect model
fixed.angina = metacont(nE, meanE, sqrt(varE),nC,meanC,sqrt(varC),
         data=angina,studlab=Protocol,comb.random=FALSE)
# Print the fitted model
fixed.angina


###################################################
### code chunk number 24: meta.angina.fixed2
###################################################
# Round to 2-digit
fixed.angina$mean.e = round(fixed.angina$mean.e,2)
fixed.angina$sd.e = round(fixed.angina$sd.e,2)
fixed.angina$mean.c = round(fixed.angina$mean.c,2)
fixed.angina$sd.c = round(fixed.angina$sd.c,2)
# Call forest to make plot
forest(fixed.angina)


###################################################
### code chunk number 25: R4Meta.rnw:554-555
###################################################
# Round to 2-digit
fixed.angina$mean.e = round(fixed.angina$mean.e,2)
fixed.angina$sd.e = round(fixed.angina$sd.e,2)
fixed.angina$mean.c = round(fixed.angina$mean.c,2)
fixed.angina$sd.c = round(fixed.angina$sd.c,2)
# Call forest to make plot
forest(fixed.angina)


###################################################
### code chunk number 26: meta.angina.funnel
###################################################
funnel(fixed.angina)


###################################################
### code chunk number 27: R4Meta.rnw:568-569
###################################################
funnel(fixed.angina)


###################################################
### code chunk number 28: R4Meta.rnw:580-581
###################################################
metabias(fixed.angina)


###################################################
### code chunk number 29: R4Meta.rnw:592-597
###################################################
# fit random-effects model
random.angina = metacont(nE, meanE, sqrt(varE),nC,meanC,sqrt(varC),
         data=angina,studlab=Protocol,comb.random=T)
# print the summary fit
random.angina


