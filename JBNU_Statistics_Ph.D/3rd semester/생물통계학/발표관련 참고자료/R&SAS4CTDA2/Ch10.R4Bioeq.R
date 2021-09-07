### R code from vignette source 'R4Bioeq.rnw'


###################################################
### code chunk number 2: BE.dat
###################################################
dat = read.csv("ChowLiuTab361.csv", header=T)


###################################################
### code chunk number 3: R4Bioeq.rnw:78-83
###################################################
library(xtable)
print(xtable(dat, digits= c(0,0,0,0,0,3),
caption="AUC data from Chow and Liu (2009) " , 
label = "BE.ChowLiudat"),table.placement = "thbp", 
caption.placement="top",include.rownames=FALSE)


###################################################
### code chunk number 4: BE.dat2
###################################################
dat = read.csv("Cimetidine.csv", header=T)


###################################################
### code chunk number 5: R4Bioeq.rnw:100-105
###################################################
library(xtable)
print(xtable(dat, digits= c(0,0,0,0,0,2,2,1),
caption="AUC, CMAX and TMAX data from Cimetidine" , 
label = "BE.Peacedat"),table.placement = "thbp", 
caption.placement="top",include.rownames=FALSE)


###################################################
### code chunk number 6: BE.dat
###################################################
dat = read.csv("ChowLiuTab361.csv", header=T)


###################################################
### code chunk number 7: BE.meantab
###################################################
# use ``aggregate" to the sample size
tab.n    = aggregate(dat$AUC,list(seq=dat$Sequence,
                              prd=dat$Period),length)
n1       = tab.n[tab.n$seq==1 & tab.n$prd==1,]$x
n2       = tab.n[tab.n$seq==2 & tab.n$prd==1,]$x
n        = n1+n2
# use ``aggregate" to get the mean
tab.mean = aggregate(dat$AUC,list(seq=dat$Sequence,
                              prd=dat$Period),mean)
# use ``aggregate" to get the variance
tab.var  = aggregate(dat$AUC,list(seq=dat$Sequence,
                              prd=dat$Period),var)
# make a dataframe for the summary data
summaryTab = data.frame(Sequence=tab.mean$seq,
Period=tab.mean$prd, numSample = tab.n$x, 
Mean = tab.mean$x, Var=tab.var$x)
# print the summary table
round(summaryTab,2)


###################################################
### code chunk number 8: R4Bioeq.rnw:430-431
###################################################
tapply(dat$AUC, dat[,c("Sequence","Period")], mean)


###################################################
### code chunk number 9: BE.Uik
###################################################
Uik           = aggregate(dat$AUC,
               list(seq = dat$Sequence,sub=dat$Subject), sum) 
colnames(Uik) = c("seq", "sub","Uik")


###################################################
### code chunk number 10: BE.mUk
###################################################
mUk  = aggregate(Uik$Uik, list(seq=Uik$seq), mean) 
colnames(mUk) = c("seq", "mUk")
print(mUk)


###################################################
### code chunk number 11: BE.C
###################################################
hatC = mUk[2,2]-mUk[1,2]
hatC


###################################################
### code chunk number 12: BE.var
###################################################
dU    = merge(Uik, mUk)
sigu2 = sum((dU$Uik-dU$mUk)^2)/(n1+n2-2) 
sigu2


###################################################
### code chunk number 13: R4Bioeq.rnw:501-504
###################################################
se.sigu = sqrt(sigu2*(1/n1+1/n2)) 
TC      = hatC/se.sigu
TC


###################################################
### code chunk number 14: BE.pC
###################################################
pC = 2*(1-pt(abs(TC), n1+n2-2))  
pC


###################################################
### code chunk number 15: BE.dd
###################################################
dik          = aggregate(dat$AUC,
               list(sub=dat$Subject,seq=dat$Sequence),diff)
dik$x        = dik$x/2  
colnames(dik)= c("sub", "seq","dik")
dik


###################################################
### code chunk number 16: R4Bioeq.rnw:537-541
###################################################
mdk           = aggregate(dik$dik, list(seq=dik$seq), mean) 
colnames(mdk) = c("seq", "mdk")
hatF          = mdk[1,2]-mdk[2,2]
hatF


###################################################
### code chunk number 17: BE.vard
###################################################
dF    = merge(dik, mdk)
sigd2 = sum((dF$dik-dF$mdk)^2)/(n1+n2-2) 
sigd2


###################################################
### code chunk number 18: R4Bioeq.rnw:568-571
###################################################
se.sigd = sqrt(sigd2*(1/n1+1/n2)) 
TF      = hatF/se.sigd
TF


###################################################
### code chunk number 19: R4Bioeq.rnw:578-580
###################################################
pF = 2*(1-pt(abs(TF), n1+n2-2))  
pF


###################################################
### code chunk number 20: BE.ANOVA
###################################################
# cat("We first re-format the data into R dataframe","\n")
Data = data.frame(subj = as.factor(dat$Subject),
                   drug = as.factor(dat$Formulation), 
                   seq  = as.factor(dat$Sequence), 
                   prd  = as.factor(dat$Period),
                   AUC  = dat$AUC)
# cat("Then call R function aov for ANOVA Table", "\n")
summary(aov(AUC ~ seq*drug + Error(subj), data = Data))


###################################################
### code chunk number 21: BE.decisionCI
###################################################
# get the mean for AUC by Formulation
mdrug    = tapply(dat$AUC, list(drug=dat$Formulation), mean)
# extract the means
ybarT    = mdrug["T"]
ybarR    = mdrug["R"]
# make the decision CI
dec2.low = theta.L = -0.2*ybarR
dec2.up  = theta.U = 0.25*ybarR
cat("DecisionCI.mean=(",dec2.low,",",dec2.up,")",sep="","\n")


###################################################
### code chunk number 22: BE.CI12
###################################################
# the confidence coefficient: alpha
alphaCI  = .1
# the t-value
qt.alpha = qt(1-alphaCI, n1+n2-2)
qt.alpha
# the lower and upper limits for CI1
low1 = (ybarT-ybarR)-qt.alpha*sqrt(sigd2)*sqrt(1/n1+1/n2)
up1  = (ybarT-ybarR)+qt.alpha*sqrt(sigd2)*sqrt(1/n1+1/n2)
cat("The classical CI1=(", round(low1,3),",", 
round(up1,3),")", sep=" ","\n\n")
# the lower and upper limits for CI2
low2 = (low1/ybarR+1)*100
up2  = (up1/ybarR+1)*100
cat("The Ratio CI2=(", round(low2,3),",", 
round(up2,3),")", sep=" ","\n\n")


###################################################
### code chunk number 23: R4Bioeq.rnw:670-671
###################################################
k12 = 2*(ybarR-ybarT)/sqrt( sigd2*(1/n1+1/n2))


###################################################
### code chunk number 24: R4Bioeq.rnw:676-680
###################################################
k2 = uniroot(function(k2) pt(k12-k2,n1+n2-2)- pt(k2,n1+n2-2)
-(1-alphaCI),lower = -10, upper = 10, tol = 0.0001)$root
k1 =k12-k2
cat("The Westlake k1=",k1," and k2=",k2,sep=" ", "\n\n")


###################################################
### code chunk number 25: R4Bioeq.rnw:684-688
###################################################
low.west = k2*sqrt(sigd2*(1/n1+1/n2))-(ybarR-ybarT)
up.west  = k1*sqrt(sigd2*(1/n1+1/n2))-(ybarR-ybarT)
cat("The Westlake CI for mu_T-mu_A is 
(",low.west,",",up.west,")",sep=" ", "\n\n")


###################################################
### code chunk number 26: R4Bioeq.rnw:696-698
###################################################
TL = (ybarT-ybarR-theta.L)/sqrt(sigd2*(1/n1+1/n2)) 
TU = (ybarT-ybarR-theta.U)/sqrt(sigd2*(1/n1+1/n2))


###################################################
### code chunk number 27: R4Bioeq.rnw:707-710
###################################################
pL     = 1-pt(abs(TL), n1+n2-2)
pU     = pt(TU,n1+n2-2)
p1side = max(pL, pU)


###################################################
### code chunk number 28: BE.bayes
###################################################
tL  = (theta.L -(ybarT-ybarR))/sqrt(sigd2*(1/n1+1/n2)) 
tU  = (theta.U -(ybarT-ybarR))/sqrt(sigd2*(1/n1+1/n2))
pRD = pt(tU, n1+n2-2) - pt(tL, n1+n2-2)
pRD


###################################################
### code chunk number 29: R4Bioeq.rnw:737-744
###################################################
dR   = dat[dat$Formulation=="R",c("Subject","AUC")]
dT   = dat[dat$Formulation=="T",c("Subject","AUC")]
colnames(dR) = c("Subject","AUC4R")
colnames(dT) = c("Subject","AUC4T")
dRT  = merge(dR,dT)
rT2R = dRT$AUC4T/dRT$AUC4R
rT2R


###################################################
### code chunk number 30: R4Bioeq.rnw:748-753
###################################################
k       = 1/sqrt(1-.9)
rbar    = mean(rT2R)
sigrbar = sqrt(var(rT2R)/n)
rbar
sigrbar


###################################################
### code chunk number 31: BE.BTCI
###################################################
low.BT = rbar-k*sigrbar
up.BT  = rbar+k*sigrbar
cat("The Tchebycheff CI for mu_T/mu_A is 
(",low.BT,",",up.BT,")",sep="", "\n\n")


###################################################
### code chunk number 32: BE.boot
###################################################
# B=number of bootstrap
B     = 2000
# boota and bootb to keep track the bootstrap results
boota = bootb = NULL
for(b in 1:B){
# Boottrap the observed individual ratios
boota[b] = mean(sample(rT2R, replace=T))
# boottrap the individuals and calculate the means
tmp = dRT[sample(1:n, replace=T),]
bootb[b] = mean(tmp$AUC4T)/mean(tmp$AUC4R)
}


###################################################
### code chunk number 33: R4Bioeq.rnw:789-794
###################################################

qxa = quantile(boota, c(0.05, 0.95))
qxa
qxb = quantile(bootb, c(0.05, 0.95))
qxb


###################################################
### code chunk number 34: BE.boot1plot
###################################################
hist(boota,nclass=30, freq=F,las=1, 
xlab="Mean of the Ratios", ylab="Density", main="")
box()
den = density(boota)
lines(den, lwd=2)
qya = approx(den$x, den$y, c(qxa,rbar))$y
segments(qxa[1],0,qxa[1],qya[1], lwd=5)
segments(qxa[2],0,qxa[2],qya[2], lwd=5)
segments(rbar,0,rbar,qya[3],lty=4, lwd=5)


###################################################
### code chunk number 35: R4Bioeq.rnw:816-817
###################################################
hist(boota,nclass=30, freq=F,las=1, 
xlab="Mean of the Ratios", ylab="Density", main="")
box()
den = density(boota)
lines(den, lwd=2)
qya = approx(den$x, den$y, c(qxa,rbar))$y
segments(qxa[1],0,qxa[1],qya[1], lwd=5)
segments(qxa[2],0,qxa[2],qya[2], lwd=5)
segments(rbar,0,rbar,qya[3],lty=4, lwd=5)


###################################################
### code chunk number 36: BE.boot2plot
###################################################
hist(bootb,nclass=30, freq=F,las=1, 
xlab="Ratio of the Means", ylab="Density", main="")
box()
den = density(bootb)
lines(den, lwd=2)
rmean = mean(dRT$AUC4T)/mean(dRT$AUC4R)
qyb = approx(den$x, den$y, c(qxb,rmean))$y
segments(qxb[1],0,qxb[1],qyb[1], lwd=5)
segments(qxb[2],0,qxb[2],qyb[2], lwd=5)
segments(rmean,0,rmean,qyb[3],lty=4, lwd=5)


###################################################
### code chunk number 37: R4Bioeq.rnw:839-840
###################################################
hist(bootb,nclass=30, freq=F,las=1, 
xlab="Ratio of the Means", ylab="Density", main="")
box()
den = density(bootb)
lines(den, lwd=2)
rmean = mean(dRT$AUC4T)/mean(dRT$AUC4R)
qyb = approx(den$x, den$y, c(qxb,rmean))$y
segments(qxb[1],0,qxb[1],qyb[1], lwd=5)
segments(qxb[2],0,qxb[2],qyb[2], lwd=5)
segments(rmean,0,rmean,qyb[3],lty=4, lwd=5)


###################################################
### code chunk number 38: R4Bioeq.rnw:849-850
###################################################
shapiro.test(bootb)


###################################################
### code chunk number 39: BE.qqnorm
###################################################
qqnorm(bootb, main="")
qqline(bootb)


###################################################
### code chunk number 40: R4Bioeq.rnw:860-861
###################################################
qqnorm(bootb, main="")
qqline(bootb)


###################################################
### code chunk number 41: BE.pdat
###################################################
datRaw0    = read.csv("CimetidineRaw.csv", header=T)


###################################################
### code chunk number 42: BE.pdat
###################################################
datRaw0    = read.csv("CimetidineRaw.csv", header=T)
print(datRaw0)


###################################################
### code chunk number 43: R4Bioeq.rnw:898-903
###################################################
datRaw           = reshape(datRaw0, direction="long", 
		    varying=-1,idvar = "Subject",sep="")
datRaw$time      = datRaw$time/10
colnames(datRaw) = c("subj","time","conc")
head(datRaw)


###################################################
### code chunk number 44: BE.fig.pdat
###################################################
library(lattice)
print(xyplot(conc~time,group=subj,datRaw,xlab="Time(HR)", 
xlim=c(0,13),auto.key = list(corner=c(1,1),lines = TRUE) ,
ylab="Concentration(mCG/ML)",type=c("p","a")))


###################################################
### code chunk number 45: R4Bioeq.rnw:917-918
###################################################
library(lattice)
print(xyplot(conc~time,group=subj,datRaw,xlab="Time(HR)", 
xlim=c(0,13),auto.key = list(corner=c(1,1),lines = TRUE) ,
ylab="Concentration(mCG/ML)",type=c("p","a")))


###################################################
### code chunk number 46: BE.fig.mpdat
###################################################
# make the mean concentration
dat.mean= aggregate(datRaw$conc, list(time=datRaw$time), mean)
# plot it with a line
plot(conc~time,las=1,type="n",datRaw,xlab="Time",xlim=c(0,13),
ylim=c(0, 4), ylab="Mean Concentration")
lines(x~time,dat.mean, lty=1, lwd=3)


###################################################
### code chunk number 47: R4Bioeq.rnw:939-940
###################################################
# make the mean concentration
dat.mean= aggregate(datRaw$conc, list(time=datRaw$time), mean)
# plot it with a line
plot(conc~time,las=1,type="n",datRaw,xlab="Time",xlim=c(0,13),
ylim=c(0, 4), ylab="Mean Concentration")
lines(x~time,dat.mean, lty=1, lwd=3)


###################################################
### code chunk number 48: BE.fn.beta
###################################################
# function `make.beta' with argument `dt'
make.beta = function(dt){
# terminal elimination beta to find the slope for R2(k+1) <R2(k)
n = length(dt$conc) # get the length of data

# end the loop at tmax
tmax = which.max(dt$conc)

# loop over starting from the last time-conc point to tmax
for(k in n:tmax){ 
dt1 = dt[((k-2):n),] # start with last 3 pts and move on
dt2 = dt[((k-3):n),] # start with last 4 pts and move on
# some date have 0s at the end of t-c curve and make the lm crash
# so make this dataframe at least 3 data points
if( dim(dt1[dt1$conc>0,])[1]>= 3 ){
# fit log(conc) to time and track the r-square
m1    = lm(log(conc)~time, dt1[(dt1$conc>0),]) 
m2    = lm(log(conc)~time, dt2[(dt2$conc>0),]) 
betat = m1$coef[[2]]
#cat("Check=",summary(m1)$r.squared > summary(m2)$r.squared," 
#and Stopped at", k, "with beta=",betat,sep=" ","\n\n")
if(summary(m1)$r.squared > summary(m2)$r.squared) break
	} # end of if-loop
} # end of k-for-loop
#cat("final beta=",betat,"\n\n")
# return
betat
} # end of make-beta function


###################################################
### code chunk number 49: BE.fn.make
###################################################
make   = function(dt){
time   = dt$time; conc = dt$conc 
#  calculate AUC
t.dif  = diff(time) # the t(i)-t(i-1)
c.mean = (conc[-1]+conc[-length(conc)])/2
auc    = sum(t.dif*c.mean)
# Cmax
cmax   = max(conc)
# tmax 
tmax   = dt[which.max(dt$conc),]$time
# terminal elimination beta to find the slope for R2(k+1) <R2(k)
betat  = make.beta(dt)
# terminal halflife
t5     = round(-log(2)/betat*2.303,1)
# AUC infinite
aucinf = auc+ conc[length(conc)]/betat 
# return the results.
c(auc,cmax,tmax, betat, t5, aucinf)
}


###################################################
### code chunk number 50: R4Bioeq.rnw:1004-1015
###################################################
name.subj   = sort(unique(datRaw$subj))
num.subj    = length(name.subj)
endpts      = matrix(0, nrow=num.subj, ncol=7)
colnames(endpts) = c("subj","AUC","CMAX","TMAX",
	"betat","t5","AUCinf") 
for(id in 1:num.subj){
tmp         = datRaw[(datRaw$subj == name.subj[id]),
		c("time","conc")]
endpts[id,] = c(name.subj[id],make(tmp))
}
endpts


###################################################
### code chunk number 51: R4Bioeq.rnw:1027-1028
###################################################
dat = read.csv("Cimetidine.csv",  header=T)


###################################################
### code chunk number 52: BE.pANOVA
###################################################
Data = data.frame(subj = as.factor(dat$Subject),
                   drug = as.factor(dat$Formulation), 
                   seq  = as.factor(dat$Sequence), 
                   prd  = as.factor(dat$Period),
                   AUC  = dat$AUC, lAUC= log(dat$AUC),
		   CMAX = dat$CMAX, lCMAX= log(dat$CMAX))


###################################################
### code chunk number 53: R4Bioeq.rnw:1044-1048
###################################################
nsubj = tapply(dat$Subject, list(dat$Sequence), length)/2
n1    = nsubj[1]
n2    = nsubj[2]
n     = n1+n2


###################################################
### code chunk number 54: R4Bioeq.rnw:1061-1067
###################################################
# the fixed model using lm for "formulation" and "period" 
mdAUC   = lm(AUC ~ seq + subj:seq + prd + drug,data = Data)
print(anova(mdAUC))
# the random effect model using aov for carryover and other effects
mdAUC.R = aov(AUC ~ prd * drug + Error(subj), data = Data)
print(summary(mdAUC.R))


###################################################
### code chunk number 55: R4Bioeq.rnw:1076-1080
###################################################
# the fixed model using lm for "formulation" and "period" 
mdlAUC  = lm(lAUC ~ seq + subj:seq + prd + drug,data = Data)
print(anova(mdlAUC))
# the random effect model using aov for carryover and other effects mdlAUC.R = aov(lAUC ~ prd * drug + Error(subj), data = Data) print(summary(mdlAUC.R))


###################################################
### code chunk number 56: R4Bioeq.rnw:1084-1090
###################################################
# the fixed model using lm for "formulation" and "period" 
mdCMAX  = lm(CMAX ~ seq + subj:seq + prd + drug,data = Data)
print(anova(mdCMAX))
# the random effect model using aov for carryover and other effects
mdCMAX.R = aov(CMAX ~ prd * drug + Error(subj), data = Data)
print(summary(mdCMAX.R))


###################################################
### code chunk number 57: R4Bioeq.rnw:1094-1100
###################################################
# the fixed model using lm for "formulation" and "period" 
mdlCMAX  = lm(lCMAX ~ seq + subj:seq + prd + drug,data = Data)
print(anova(mdlCMAX))
# the random effect model using aov for carryover and other effects
mdlCMAX.R = aov(lCMAX ~ prd * drug + Error(subj), data = Data)
print(summary(mdlCMAX.R))


###################################################
### code chunk number 58: BE.decisionCI
###################################################
mdrug    = tapply(dat$AUC, list(drug=dat$Formulation), mean)
ybarT    = mdrug["T"]
ybarR    = mdrug["R"]
dec2.low = theta.L = -0.2*ybarR
dec2.up  = theta.U = 0.25*ybarR
cat("DecisionCI.mean=(",dec2.low,",",dec2.up,")",sep="","\n")


###################################################
### code chunk number 59: BE.CI12
###################################################
# the confidence coefficient: alpha
alphaCI  = .1
# the t-value
qt.alpha = qt(1-alphaCI, n1+n2-2)
qt.alpha
# the sigma using the ANOVA model instead
sigd2 = anova(mdAUC)[5,3]/2
# the lower and upper limits for CI1
low1 = (ybarT-ybarR)-qt.alpha*sqrt(sigd2)*sqrt(1/n1+1/n2)
up1  = (ybarT-ybarR)+qt.alpha*sqrt(sigd2)*sqrt(1/n1+1/n2)
cat("The classical CI1=(", round(low1,3),",", 
round(up1,3),")", sep=" ","\n\n")
# the lower and upper limits for CI2
low2 = (low1/ybarR+1)*100
up2  = (up1/ybarR+1)*100
cat("The Ratio CI2=(", round(low2,3),",", 
round(up2,3),")", sep=" ","\n\n")


###################################################
### code chunk number 60: R4Bioeq.rnw:1161-1162
###################################################
k12 = 2*(ybarR-ybarT)/sqrt( sigd2*(1/n1+1/n2))


###################################################
### code chunk number 61: R4Bioeq.rnw:1167-1171
###################################################
k2 = uniroot(function(k2) pt(k12-k2,n1+n2-2)- pt(k2,n1+n2-2)
-(1-alphaCI),lower = -10, upper = 10, tol = 0.0001)$root
k1 =k12-k2
cat("The Westlake k1=",k1," and k2=",k2,sep=" ", "\n\n")


###################################################
### code chunk number 62: R4Bioeq.rnw:1175-1179
###################################################
low.west = k2*sqrt(sigd2*(1/n1+1/n2))-(ybarR-ybarT)
up.west  = k1*sqrt(sigd2*(1/n1+1/n2))-(ybarR-ybarT)
cat("The Westlake CI for mu_T-mu_A is 
(",low.west,",",up.west,")",sep=" ", "\n\n")


###################################################
### code chunk number 63: R4Bioeq.rnw:1187-1189
###################################################
TL = (ybarT-ybarR-theta.L)/sqrt(sigd2*(1/n1+1/n2)) 
TU = (ybarT-ybarR-theta.U)/sqrt(sigd2*(1/n1+1/n2))


###################################################
### code chunk number 64: R4Bioeq.rnw:1198-1201
###################################################
pL     = 1-pt(abs(TL), n1+n2-2)
pU     = pt(TU,n1+n2-2)
p1side = max(pL, pU)


###################################################
### code chunk number 65: BE.bayes
###################################################
tL  = (theta.L -(ybarT-ybarR))/sqrt(sigd2*(1/n1+1/n2)) 
tU  = (theta.U -(ybarT-ybarR))/sqrt(sigd2*(1/n1+1/n2))
pRD = pt(tU, n1+n2-2) - pt(tL, n1+n2-2)
pRD


###################################################
### code chunk number 66: R4Bioeq.rnw:1228-1235
###################################################
dR   = dat[dat$Formulation=="R",c("Subject","AUC")]
dT   = dat[dat$Formulation=="T",c("Subject","AUC")]
colnames(dR) = c("Subject","AUC4R")
colnames(dT) = c("Subject","AUC4T")
dRT  = merge(dR,dT)
rT2R = dRT$AUC4T/dRT$AUC4R
rT2R


###################################################
### code chunk number 67: R4Bioeq.rnw:1239-1244
###################################################
k       = 1/sqrt(1-.9)
rbar    = mean(rT2R)
sigrbar = sqrt(var(rT2R)/n)
rbar
sigrbar


###################################################
### code chunk number 68: BE.BTCI
###################################################
low.BT = rbar-k*sigrbar
up.BT  = rbar+k*sigrbar
cat("The Tchebycheff CI for mu_T/mu_A is 
(",low.BT,",",up.BT,")",sep="", "\n\n")


###################################################
### code chunk number 69: BE.boot
###################################################
# B=number of bootstrap
B     = 2000
# boota and bootb to keep track the bootstrap results
boota = bootb = NULL
for(b in 1:B){
# Boottrap the observed individual ratios
boota[b] = mean(sample(rT2R, replace=T))
# boottrap the individuals and calculate the means
tmp = dRT[sample(1:n, replace=T),]
bootb[b] = mean(tmp$AUC4T)/mean(tmp$AUC4R)
}


###################################################
### code chunk number 70: R4Bioeq.rnw:1280-1284
###################################################
qxa = quantile(boota, c(0.05, 0.95))
qxa
qxb = quantile(bootb, c(0.05, 0.95))
qxb


###################################################
### code chunk number 71: BE.boot1plot
###################################################
hist(boota,nclass=30, freq=F,las=1, 
xlab="Mean of the Ratios", ylab="Density", main="")
box()
den = density(boota)
lines(den, lwd=2)
qya = approx(den$x, den$y, c(qxa,rbar))$y
segments(qxa[1],0,qxa[1],qya[1], lwd=5)
segments(qxa[2],0,qxa[2],qya[2], lwd=5)
segments(rbar,0,rbar,qya[3],lty=4, lwd=5)


###################################################
### code chunk number 72: R4Bioeq.rnw:1307-1308
###################################################
hist(boota,nclass=30, freq=F,las=1, 
xlab="Mean of the Ratios", ylab="Density", main="")
box()
den = density(boota)
lines(den, lwd=2)
qya = approx(den$x, den$y, c(qxa,rbar))$y
segments(qxa[1],0,qxa[1],qya[1], lwd=5)
segments(qxa[2],0,qxa[2],qya[2], lwd=5)
segments(rbar,0,rbar,qya[3],lty=4, lwd=5)


###################################################
### code chunk number 73: BE.boot2plot
###################################################
hist(bootb,nclass=30, freq=F,las=1, 
xlab="Ratio of the Means", ylab="Density", main="")
box()
den = density(bootb)
lines(den, lwd=2)
rmean = mean(dRT$AUC4T)/mean(dRT$AUC4R)
qyb = approx(den$x, den$y, c(qxb,rmean))$y
segments(qxb[1],0,qxb[1],qyb[1], lwd=5)
segments(qxb[2],0,qxb[2],qyb[2], lwd=5)
segments(rmean,0,rmean,qyb[3],lty=4, lwd=5)


###################################################
### code chunk number 74: R4Bioeq.rnw:1330-1331
###################################################
hist(bootb,nclass=30, freq=F,las=1, 
xlab="Ratio of the Means", ylab="Density", main="")
box()
den = density(bootb)
lines(den, lwd=2)
rmean = mean(dRT$AUC4T)/mean(dRT$AUC4R)
qyb = approx(den$x, den$y, c(qxb,rmean))$y
segments(qxb[1],0,qxb[1],qyb[1], lwd=5)
segments(qxb[2],0,qxb[2],qyb[2], lwd=5)
segments(rmean,0,rmean,qyb[3],lty=4, lwd=5)


