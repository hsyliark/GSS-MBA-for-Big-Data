n=1000
x = runif(n)
vc1 = function(x)
{
  1 + 5*x^2 + sin(3*pi*x)
}
vc2 = function(x)
{
  2 + 2*x^3
}
z1 = - runif(n)
z2 = runif(n)
y = 1 + z1*vc1(x) + z2*vc2(x) + rnorm(n)

g=seq(0,1,0.01)
plot(g,vc1(g))
plot(g,vc2(g))

library(mgcv)
fit1 = gam(y ~ s(x,by=z1)+s(x,by=z2))
summary(fit1)
ls(fit1)


### estimated functions at data point2 (1 component)
#aa = predict(fit1,type="terms",se.fit=T)
#aa1 = predict(fit1,type="terms",newdata=data.frame(x=1,z=1))
#aa2 = predict(fit1,type="terms",newdata=data.frame(x=1,z=2))

#ind = order(x)

#plot(fit1)
#lines(g,vc(g),col=4)
#lines(x[ind],aa$fit[ind]/z[ind],col=2)

#hs = aa$fit/z
#lines(x[ind],hs[ind]-2*aa$se.fit[ind]/z[ind],col=2)
#lines(x[ind],hs[ind]+2*aa$se.fit[ind]/z[ind],col=2)



### estimated functions at grid g

aa = predict(fit1,type="terms",se.fit=T, newdata=data.frame(x=g,z1=1,z2=1))

par(mfrow=c(1,2))
plot(fit1)

hs = aa$fit

a1 = range(hs[,1])
l11=a1[1]-diff(a1)/2
l12=a1[2]+diff(a1)/2
plot(g,hs[,1],type="l",ylim=c(l11,l12))
rug(x)
lines(g,vc1(g),col=4)
lines(g,hs[,1]-2*aa$se.fit[,1],col=2,lty="dashed")
lines(g,hs[,1]+2*aa$se.fit[,1],col=2,lty="dashed")

a2 = range(hs[,2])
l21=a2[1]-diff(a2)/2
l22=a2[2]+diff(a2)/2
plot(g,hs[,2],type="l",ylim=c(l21,l22))
rug(x)
lines(g,vc2(g),col=4)
lines(g,hs[,2]-2*aa$se.fit[,2],col=2,lty="dashed")
lines(g,hs[,2]+2*aa$se.fit[,2],col=2,lty="dashed")

