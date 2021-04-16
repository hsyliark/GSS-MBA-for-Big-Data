### 중간시험 예상문제 답안



## Problem 1
n <- 100000
x <- runif(n)*2-1 ; y <- runif(n)*2-1
d <- sqrt(x^2+y^2)
count <- sum(d<1)
(pi_cal <- (count/n)*4) # answer
pi


## Problem 2

# (a)
n <- 100000
x <- runif(n)
f <- function(x) 4*sqrt(1-x^2)
mean(f(x)) # answer

# (b)
n <- 100000
x <- runif(n)
mean(sin(x)) # answer

# (c)
n <- 100000
x <- runif(n)
f <- function(x) x^3
mean(f(x)) # answer

# (d)
n <- 100000
x <- runif(n)
mean(exp(x)) # answer

# (e)
n <- 100000
x <- rexp(n,0.1)
f <- function(x) 0.1*exp(-x/10)
g <- function(x) (x+1)^(-2)
mean(g(x)/f(x)) # answer


## Problem 3

# (a)
n <- 100
U <- runif(n)
Finv <- function(x) 1-log(x)/log(0.7)
(bernoulli <- as.numeric(Finv(U)>0)) # answer
par(mfrow=c(1,2))
hist(bernoulli,freq=T) 
hist(rbinom(100,1,0.3)) # compare histogram
par(mfrow=c(1,1))

# (b)
n <- 1000
U <- runif(n)
Finv <- function(x) 1-log(x)/log(0.7)
X <- as.numeric(Finv(U)>0)
dim(X)<-c(10,100)
(binomial <- apply(X,2,sum)) # answer
par(mfrow=c(1,2))
hist(binomial,freq=T) 
hist(rbinom(100,10,0.3)) # compare histogram
par(mfrow=c(1,1))

# (c)
n <- 100
U <- runif(n)
Finv <- function(x) -(1/2.5)*log(1-x)
(exponential <- Finv(U)) # answer
par(mfrow=c(1,2))
hist(exponential,freq=T) 
hist(rexp(100,2.5)) # compare histogram
par(mfrow=c(1,1))

# (d)
n <- 1000
U <- runif(n)
Finv <- function(x) -log(1-x)
dim(U)<-c(10,100)
(Gam <- apply(U,2,sum)) # answer
par(mfrow=c(1,2))
hist(Gam,freq=T) 
hist(rgamma(100,shape=10,scale=1)) # compare histogram
par(mfrow=c(1,1))


## Problem 4

# Graph of f(x)
x<-1:10000/10000
f<-function(x){
  -4*(x-1)*(x>=0.5) + 4*x*(x<0.5)
}
plot(x,f(x),type='l')

# Sampling
n <- 10000
a <- 1
b <- 2
U1 <- runif(n)*a
U2 <- runif(n)*b
plot(U1,U2)
lines(x,f(x),col=2,lwd=7,lty=2)
index_ <- U2<f(U1)
plot(U1,U2)
points(U1[index_],U2[index_],col=2)

# answer
n <- 100
a <- 1
b <- 2
U1 <- runif(n)*a
U2 <- runif(n)*b
index_ <- U2<f(U1)
(Y <- U1[index_])
hist(Y)


## Problem 5

# function f and g
x <- 1:9999/10000
f <- function(x) {
  alpha=2.7 ; beta=6.3
  (gamma(alpha+beta)/(gamma(alpha)*gamma(beta)))*(x)^(alpha-1)*(1-x)^(beta-1)
}
g <- function(x) {
  alpha=2 ; beta=6
  (gamma(alpha+beta)/(gamma(alpha)*gamma(beta)))*(x)^(alpha-1)*(1-x)^(beta-1)
}
M <- 1/min(g(x)/f(x))

# Sampling from Gamma(2,1)
n <-10000*2 
Exp <- -log(1-runif(n))
dim(Exp) <-c (2,10000)
Gamma2 <- apply(Exp,2,sum)

# Sampling from Gamma(6,1)
n <- 10000*6
Exp <- -log(1-runif(n))
dim(Exp) <- c(6,10000)
Gamma6 <- apply(Exp,2,sum)

# Making Samples (Beta(2,6))
Beta26 <- Gamma2/(Gamma2+Gamma6)

# Acceptance/Rejection Sampling
Y <- Beta26
U <- runif(10000)
index_ <- c()
for (i in 10000) {
  index_[i] <- U[i] < f(Y[i])/(M*g(Y[i]))
}
X <- Y[index_] # answer
par(mfrow=c(1,2))
hist(X,freq=T) 
hist(rbeta(length(X),shape1=2.7,shape2=6.3)) # compare histogram
par(mfrow=c(1,1))


## Problem 6

# (a)
n <- 10000
X <- rnorm(n)
index_ <- (X > -1.96)&(X < 1.96) 
mean(index_) # answer

# (b)
n <- 10000
X <- rnorm(n)
index_ <- (X > 4.5)
mean(index_)

# (c)
# -> 표준정규분포에서 생성된 난수가 4.5를 넘는 경우는 
# 난수를 많이 뽑지 않는 한 거의 발생하지 않기 때문에 
# 결과를 신뢰하기 어렵다.

# (d)
# importance sampling
n <- 10000
Y <- rexp(n)+4.5
f <- function(x) exp(-(x^2)/2)/sqrt(2*pi)
g <- function(x) exp(-(x-4.5))
mean(f(Y)/g(Y)) # answer


## Problem 7

# (a)
n <- 10000*2
Exp <- rexp(n,rate=1)
dim(Exp) <- c(2,10000)

