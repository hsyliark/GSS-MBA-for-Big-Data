#### Assignment (E-M Algorithm)




### Q1

set.seed(1)
n <- 100
p <- 0.5
u <- rbinom(n,1,p)
x1 <- rnorm(n,162,6)
x2 <- rnorm(n,173,6)
x <- u*x1+(1-u)*x2
plot(density(x))

mu1 <- mean(x)-sd(x)
mu2 <- mean(x)+sd(x)
s <- var(x)
pi=0.5
newmu1=newmu2=news=newpi=ll=c()
for (i in 1:5000)
{
  hatdelta <- (pi*dnorm(x,mu1,sqrt(s)))/(pi*dnorm(x,mu1,sqrt(s))+(1-pi)*dnorm(x,mu2,sqrt(s)))
  newmu1[i] <- sum(hatdelta*x)/sum(hatdelta)
  newmu2[i]=sum((1-hatdelta)*x)/sum((1-hatdelta))
  news[i]=sum(hatdelta*(x-newmu1[i])^2+(1-hatdelta)*(x-newmu2[i])^2)/n
  newpi[i]=sum(hatdelta)/n
  ll[i]=sum(log(newpi[i]*dnorm(x,newmu1[i],sqrt(news[i]))+(1-newpi[i])*dnorm(x,newmu2[i],sqrt(news[i]))))
  if (max(abs(mu1-newmu1[i]),abs(mu2-newmu2[i]),abs(s-news[i]),abs(pi-newpi[i]))<10E-5) break 
  mu1=newmu1[i]
  mu2=newmu2[i]
  s=news[i]
  pi=newpi[i]
  #print(c(i,mu1,mu2,s,pi))
}
c(mu1,mu2,s,pi)


### Q2

newtheta <- theta <- 0.3
tol <- 10^{-8}
t <- c(125,18,20,34)
for (i in 1:100)
{
  newtheta <- (theta*t[1]+(2+theta)*t[4])/(theta*t[1]+(2+theta)*(t[2]+t[3]+t[4]))
  err <- max(abs(newtheta-theta))
  print(c(i,theta))
  if (err<tol) break
  theta <- newtheta
}
