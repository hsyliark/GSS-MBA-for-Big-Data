library(fda)
library(Hmisc)

str(growth)
str(growth$hgtf)
matplot(growth$hgtf)
set.seed(1)
s <- sample(1:54,3)

h <- Data2fd(growth$age, growth$hgtf,
               fdnames=list("age", "girls", "height"))
plot(h)
#daybasis65 <- create.bspline.basis(c(0, 365), nbasis=35, norder=4)
#daytempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"],
#                          daybasis65, fdnames=list("Day", "Station", "Deg C"))$fd
#plot(daytempfd)

th <- h[-s,]
vh <- h[s,]


idx <- growth$age <= 9
h1 <- Data2fd(growth$age[idx], growth$hgtf[idx,],
             fdnames=list("age", "girls", "height"))
h2 <- Data2fd(growth$age[!idx], growth$hgtf[!idx,],
              fdnames=list("age", "girls", "height"))

th1 <- h1[-s,]
vh1 <- h1[s,]
th2 <- h2[-s,]
vh2 <- h2[s,]

npc1 = 5
h1.pca <- pca.fd(th1, nharm=npc1)   
#plot.pca.fd(h1.pca)  
h1.pca$varprop
s1 <- h1.pca$scores

npc2 = 5
h2.pca <- pca.fd(th2, nharm=npc2)  
#plot.pca.fd(h2.pca)  
cumsum(h2.pca$varprop)
spc2 <- sum(cumsum(h2.pca$varprop) < 0.98) + 1
s2 <- h2.pca$scores[,1:spc2]

#cor(cbind(s1,s2))[(npc1+1):(npc1+spc2),1:npc1]
rcorr(as.matrix(cbind(s1,s2)))$P[(npc1+1):(npc1+spc2),1:npc1]<0.1

d = 5
eta0 <- lm(s2[,1]~s1[,1]-1)$coefficients
eta1 <- lm(s2[,1]~s1[,1:d]-1)$coefficients
eta2 <- lm(s2[,2]~s1[,1:d]-1)$coefficients

mse=mse0=0
mse1=0
mse2=0
for (i in 1:length(s))
{
nd <- vh1[i,]
xi1=c()
for (j in 1:d)
{
xi1[j] <- inprod.bspline(nd-h1.pca$meanfd,h1.pca$harmonics[j,])
}
pd0 <- h2.pca$meanfd + eta0*xi1[1]*h2.pca$harmonics[1,]
pd1 <- h2.pca$meanfd + sum(eta1*xi1)*h2.pca$harmonics[1,]
pd2 <- h2.pca$meanfd + sum(eta1*xi1)*h2.pca$harmonics[1,] + sum(eta2*xi1)*h2.pca$harmonics[2,]
#plot(vh2[i,],col=i,ylim=c(50,180))
#lines(pd1,col=i+1,lty=3)
#lines(pd2,col=i+2,lty=4)
mse0 = mse0 + sum(eval.fd(growth$age[!idx],pd0-vh2[i,])^2)
mse1 = mse1 + sum(eval.fd(growth$age[!idx],pd1-vh2[i,])^2)
mse2 = mse2 + sum(eval.fd(growth$age[!idx],pd2-vh2[i,])^2)
plot(growth$age,growth$hgtf[,s[i]])
lines(pd2,col=2)
}

