## 1. Gamma distribution

# Newton's method

theta <- c(1,1)
dl1 <- function(x,alpha,lambda){length(x)*log(lambda)-length(x)*digamma(alpha)+sum(log(x))}
dl2 <- function(x,alpha,lambda){length(x)*alpha/lambda-sum(x)}
ddl11<- function(x,alpha,lambda){-length(x)*trigamma(alpha)}
ddl12<- function(x,alpha,lambda){length(x)/lambda}
ddl21<- function(x,alpha,lambda){length(x)/lambda}
ddl22<- function(x,alpha,lambda){-length(x)*alpha/lambda^2}
tol <- 10^(-8)

set.seed(1)
x <- rgamma(100,shape=3,rate=5)
for (i in 1:100)
{
  ddl <- matrix(c(ddl11(x,theta[1],theta[2]),ddl21(x,theta[1],theta[2]),ddl12(x,theta[1],theta[2]),ddl22(x,theta[1],theta[2])),2,2)
  newtheta <- theta - solve(ddl)%*%c(dl1(x,theta[1],theta[2]),dl2(x,theta[1],theta[2]))
  err <- sum((newtheta-theta)^2)
  print(c(i,theta))
  if (err<tol) break
  theta<-newtheta
}

theta


## 2. Geometric distribution (with re-parametrization)

dl <- function(x,beta){sum(x)/beta-(length(x)+sum(x))/(1+beta)}
ddl <- function(x,beta){-sum(x)/beta^2+(length(x)+sum(x))/(1+beta)^2}
I <- function(x,beta){length(x)/(beta*(1+beta))}
tol <- 10^(-8)

# Newton's method

set.seed(1)
beta <- 1
x <- rgeom(100,1/(1+3))

for (i in 1:100)
{
  newbeta <- beta - dl(x,beta)/ddl(x,beta)
  err <- abs(newbeta-beta)
  print(c(i,beta))
  if (err<tol) break
  beta<-newbeta
}

beta

# Scoring method

set.seed(1)
beta <- 1
x <- rgeom(100,1/(1+3))

for (i in 1:100)
{
  newbeta <- beta + dl(x,beta)/I(x,beta)
  err <- abs(newbeta-beta)
  print(c(i,beta))
  if (err<tol) break
  beta<-newbeta
}

beta

