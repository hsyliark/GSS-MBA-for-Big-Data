### R code from vignette source 'R4AE.rnw'

###################################################
### code chunk number 1: R4AE.rnw:12-35
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
# lattice setting
library(lattice)
#trellis.par.set(canonical.theme(color=FALSE))
ltheme = canonical.theme(color = FALSE) # in-built B&W theme
ltheme$strip.background$col = "transparent " # change strip bg
lattice.options(default.theme=ltheme) # set as default
# random seed
set.seed(12345678)


###################################################
### code chunk number 2: AE.dat
###################################################
dat = read.csv("datAE.csv",header=T)


###################################################
### code chunk number 3: R4AE.rnw:64-67
###################################################
library(xtable)
print(xtable(dat, caption="AE data for clinical trial.",label = "tab4AE.dat"),
table.placement = "thbp",caption.placement = "top")


###################################################
### code chunk number 4: R4AE.rnw:87-94
###################################################
dat$p2 = dat$fD2/dat$ND2
dat$p1 = dat$fD1/dat$ND1
dat$pC = dat$fC/dat$NC
dat$V2 = dat$p2*(1-dat$p2)/dat$ND2
dat$V1 = dat$p1*(1-dat$p1)/dat$ND1
dat$VC = dat$pC*(1-dat$pC)/dat$NC
len= length(dat[,1])


###################################################
### code chunk number 5: AE.fig4dat
###################################################
plot(p2~Stage, type="n",dat,xlim=c(1,9), 
ylim=c(0,max(dat$p1, dat$p2, dat$pC)), 
xlab="Stage",las=1, ylab="AE Rate")
lines(pC~Stage, dat,lwd=3, lty=8)
lines(p1~Stage, dat,lwd=3, lty=4)
lines(p2~Stage, dat,lwd=3, lty=1)
legend("bottomright", legend = c("C","p1","p2"),
lty = c(8,4,1),lwd=3,title = "Line Types")


###################################################
### code chunk number 6: R4AE.rnw:112-113
###################################################
plot(p2~Stage, type="n",dat,xlim=c(1,9), 
ylim=c(0,max(dat$p1, dat$p2, dat$pC)), 
xlab="Stage",las=1, ylab="AE Rate")
lines(pC~Stage, dat,lwd=3, lty=8)
lines(p1~Stage, dat,lwd=3, lty=4)
lines(p2~Stage, dat,lwd=3, lty=1)
legend("bottomright", legend = c("C","p1","p2"),
lty = c(8,4,1),lwd=3,title = "Line Types")


###################################################
### code chunk number 7: AE.dat (eval = FALSE)
###################################################
## dat = read.csv("datAE.csv", header=T)


###################################################
### code chunk number 8: AE.rate
###################################################
# the AE rate
dat$p2 = dat$fD2/dat$ND2
dat$p1 = dat$fD1/dat$ND1
dat$pC = dat$fC/dat$NC
# The variance
dat$V2 = dat$p2*(1-dat$p2)/dat$ND2
dat$V1 = dat$p1*(1-dat$p1)/dat$ND1
dat$VC = dat$pC*(1-dat$pC)/dat$NC
len    = length(dat[,1])


###################################################
### code chunk number 9: AE.directCI
###################################################
# Function for direct comparison
direct.CI = function(N1,f1,N2,f2,alpha){
	p1 = f1/N1; p2 = f2/N2
	v1 = p1*(1-p1)/N1;v2 = p2*(1-p2)/N2
	z.alpha = qnorm(1-alpha)
	low =(p1-p2)-z.alpha*sqrt(v1+v2)
	up  =(p1-p2)+z.alpha*sqrt(v1+v2)
	data.frame(testp=p1>p2,low=low,diff=p1-p2,
			p1=p1,p2=p2,N1=N1,N2=N2)
}


###################################################
### code chunk number 10: R4AE.rnw:352-360
###################################################
# call ``direct.CI" to compare dose 1 to control
CI1toC = direct.CI(dat$ND1, dat$fD1, dat$NC, dat$fC, 0.025)
# print the calculation
CI1toC
# call ``direct.CI" to compare dose 2 to control
CI2toC = direct.CI(dat$ND2, dat$fD2, dat$NC, dat$fC, 0.025)
# print the calculation
CI2toC


###################################################
### code chunk number 11: AE.fn.indirectCI
###################################################
# Function for indirect comparison: using normal approximation
indirect.CI = function(N,f,alpha){
	p = f/N; v = p*(1-p)/N
	z.alpha = qnorm(1-alpha)
	low     = p-z.alpha*sqrt(v)
	up      = p+z.alpha*sqrt(v)
	data.frame(low=round(low,3),up=round(up,3),N=N,f=f,p=round(p,3))
}


###################################################
### code chunk number 12: AE.indirect.cal
###################################################
CIC = indirect.CI(dat$NC, dat$fC, 0.025)
CIC
CI1 = indirect.CI(dat$ND1, dat$fD1, 0.025)
CI1
CI2 = indirect.CI(dat$ND2, dat$fD2, 0.025)
CI2


###################################################
### code chunk number 13: AE.indirect.cal2
###################################################
# make a dataframe for dose 1 to control
out1toC = data.frame(Stage= dat$Stage, 
	indirect.test = CI1$low > CIC$up,
	low1=CI1$low, upC =CIC$up) 
# print it
print(out1toC)
# make a dataframe for dose 2 to control
out2toC = data.frame(Stage= dat$Stage,
	indirect.test = CI2$low > CIC$up,
	low2=CI2$low, upC =CIC$up) 
# print it
print(out2toC)


###################################################
### code chunk number 14: R4AE.rnw:432-439
###################################################
# Function to calculate the CI for binomial
exact.CI = function(N,f,alpha){
	p   = f/N
	low = qbinom(alpha, N, p)
	up  = qbinom(1-alpha,N,p)
data.frame(N=N,p=p,f=f,low=low,up=up, plow =low/N, pup=up/N)
}


###################################################
### code chunk number 15: R4AE.rnw:443-452
###################################################
# for control
CIC = exact.CI(dat$NC, dat$fC, 0.025)
CIC
# for dose 1
CI1 = exact.CI(dat$ND1, dat$fD1, 0.025)
CI1
# for dose 2
CI2 = exact.CI(dat$ND2, dat$fD2, 0.025)
CI2


###################################################
### code chunk number 16: R4AE.rnw:457-467
###################################################
# dataframe for dose 1 to control
out1toC = data.frame(Stage= dat$Stage, 
	indirect.test = CI1$plow > CIC$pup,
	low1=CI1$plow, upC =CIC$pup) 
print(out1toC)
# dataframe for dose 2 to control
out2toC = data.frame(Stage= dat$Stage, 
	indirect.test = CI2$plow > CIC$pup,
	low2=CI2$plow, upC =CIC$pup) 
print(out2toC)


###################################################
### code chunk number 17: AE.bound
###################################################
# Function to calculate bounds
bound.CI = function(N1,f1,N2,f2){
	p1 = f1/N1
	p2 = f2/N2
	v1 = p1*(1-p1)/N1
	v2 = p2*(1-p2)/N2
data.frame(bound= (p2-p1)/(sqrt(v1)+sqrt(v2)),
	N1=N1,f1=f1,p1=p1,N2=N2,f2=f2,p2=p2)
}


###################################################
### code chunk number 18: R4AE.rnw:503-509
###################################################
# call function ``bound.CI" and make the calculation
d0              = bound.CI(dat$NC,dat$fC,dat$ND2,dat$fD2)
# calculate the alpha from normal approximation
d0$alpha.normal = 1-pnorm(d0$bound)
# print it 
round(d0,4)


###################################################
### code chunk number 19: AE.fignorm
###################################################
plot(dat$Stage, d0$alpha.normal,type="o", xlab="Stage", 
ylab=expression(alpha), las=1, main="")
text(dat$Stage, d0$alpha.normal,round(d0$alpha.normal,3))


###################################################
### code chunk number 20: R4AE.rnw:527-528
###################################################
plot(dat$Stage, d0$alpha.normal,type="o", xlab="Stage", 
ylab=expression(alpha), las=1, main="")
text(dat$Stage, d0$alpha.normal,round(d0$alpha.normal,3))


###################################################
### code chunk number 21: AE.binom
###################################################
# make a funtion for L_D2-UC
fn = function(alpha, stage){
LD2 = exact.CI(dat$ND2,dat$fD2,alpha)$plow[stage]
UC  = exact.CI(dat$NC,dat$fC,alpha)$pup[stage]
LD2-UC
}
# call R function ``uniroot" to solve equation
est.alpha = NULL
for(s in 1:len) 
est.alpha[s]=uniroot(fn, c(0,0.8), stage=s)$root 
est.alpha


###################################################
### code chunk number 22: AE.sample1
###################################################
# set.seed to fix the seed for random number generation
set.seed(123)
# number of simulation
num.sim = 1000
# matrix to hold the output "bound"
bound3 = matrix(0, ncol=len, nrow=num.sim)
for(stage in 1:len){
	ds = d0[stage,] 
     # make the 0 and 1's from the data from control
	x1 = c( rep(1, ds$f1), rep(0, ds$N1-ds$f1)) 
     # from D2
	x2 = c( rep(1, ds$f2), rep(0, ds$N2-ds$f2)) 
     # combine them 
     x  = data.frame(id = 1:(ds$N1+ds$N2),c(x1,x2))

tsim=0
repeat{
	 # sample it and sum the freq
	  f = sample(x$id, ds$N1)
	 f1 = sum(x[f,2])  
	 f2 = sum(x[-f,2])
	 p1 = f1/ds$N1
	 p2 = f2/ds$N2
if(p2 > p1){ 
		tsim=tsim+1
           bound3[tsim,stage] = bound.CI(ds$N1,f1,ds$N2,f2)$bound
		} # end of if
if(tsim ==num.sim) break
	} # end of repeat
} # end of stage


###################################################
### code chunk number 23: R4AE.rnw:625-626
###################################################
dim(bound3)


###################################################
### code chunk number 24: AE.sampling2
###################################################
bound4 = matrix(0, ncol=len, nrow=num.sim)
for(stage in 1:len){
	ds = d0[stage,]
tsim=0
repeat{
	f1 = rbinom(1, ds$N1, (ds$p1+ds$p2)/2) 
	f2 = rbinom(1, ds$N2, (ds$p1+ds$p2)/2) 
	p1 = f1/ds$N1
	p2 = f2/ds$N2
if(p2 > p1){ 
		tsim=tsim+1
            bound4[tsim,stage] = bound.CI(ds$N1,f1,ds$N2,f2)$bound
		} # end of if
if(tsim ==num.sim) break
	} # end of repeat
} # end of stage

dim(bound4)


###################################################
### code chunk number 25: AE.sampling1
###################################################
est.alpha3 = NULL
est.alpha4 = NULL
for(i in 1:len){
	regimen3 = bound3[,i]
	regimen4 = bound4[,i]
	# calculate how many simulations > the obs bound in d0$bound
	est.alpha3[i] = sum(regimen3 > d0$bound[i])/num.sim
	est.alpha4[i] = sum(regimen4 > d0$bound[i])/num.sim
}
est.alpha3
est.alpha4


###################################################
### code chunk number 26: AE.allalpha
###################################################
# make the alpha plot
d0$alpha.binom = est.alpha
d0$alpha.samp1 = est.alpha3
d0$alpha.samp2 = est.alpha4
# keep track of the max values for alpha for y-axis limit
max.alpha = max(d0$alpha.normal, d0$alpha.binom, 
		d0$alpha.samp1, d0$alpha.samp2)
# plot the alpha
plot(dat$Stage, d0$alpha.normal,type="n", ylim=c(0,max.alpha), 
xlab="Stage", ylab=expression(alpha), las=1, main="")
# add lines to the plots
lines(dat$Stage, d0$alpha.normal,lwd=3,lty=1)
lines(dat$Stage, d0$alpha.binom,lwd=3,lty=3)
lines(dat$Stage, d0$alpha.samp1,lwd=3,lty=4)
lines(dat$Stage, d0$alpha.samp2,lwd=3,lty=8)
abline(h=0.025, lwd=4, lty=1)
# add legend
temp <- legend("topright", legend = c(" ", " ","  "," "),
               text.width = strwidth("Normal Approx"),
               lty = c(1,3,4,8),lwd=3,xjust = 1, yjust = 1,
               title = "Line Types")
text(temp$rect$left + temp$rect$w, temp$text$y, c("Normal Approx", 
"Exact Binomial","Resampling 1","Resampling 2"), pos=2)


###################################################
### code chunk number 27: R4AE.rnw:709-710
###################################################
# make the alpha plot
d0$alpha.binom = est.alpha
d0$alpha.samp1 = est.alpha3
d0$alpha.samp2 = est.alpha4
# keep track of the max values for alpha for y-axis limit
max.alpha = max(d0$alpha.normal, d0$alpha.binom, 
		d0$alpha.samp1, d0$alpha.samp2)
# plot the alpha
plot(dat$Stage, d0$alpha.normal,type="n", ylim=c(0,max.alpha), 
xlab="Stage", ylab=expression(alpha), las=1, main="")
# add lines to the plots
lines(dat$Stage, d0$alpha.normal,lwd=3,lty=1)
lines(dat$Stage, d0$alpha.binom,lwd=3,lty=3)
lines(dat$Stage, d0$alpha.samp1,lwd=3,lty=4)
lines(dat$Stage, d0$alpha.samp2,lwd=3,lty=8)
abline(h=0.025, lwd=4, lty=1)
# add legend
temp <- legend("topright", legend = c(" ", " ","  "," "),
               text.width = strwidth("Normal Approx"),
               lty = c(1,3,4,8),lwd=3,xjust = 1, yjust = 1,
               title = "Line Types")
text(temp$rect$left + temp$rect$w, temp$text$y, c("Normal Approx", 
"Exact Binomial","Resampling 1","Resampling 2"), pos=2)


