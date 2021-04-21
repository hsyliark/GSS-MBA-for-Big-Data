### midterm (202055364 황성윤)

#-----------------------------------------------------------------#

## Problem 1
n <- 100000
x <- runif(n)*2-1
y <- runif(n)*2-1
d <- sqrt(x^2+y^2)
count <- sum(d<1)
(pi_estimate <- (count/n)*4) # answer (계산결과 3.1415 근방의 값이 산출됨을 확인할 수 있음.)
pi # checking real answer

#-----------------------------------------------------------------#

## Problem 2

# (a)
n <- 100000
x <- runif(n)
f <- function(x) 4*sqrt(1-x^2)
mean(f(x)) # answer (약 3.1415)

# (b)
n <- 100000
x <- runif(n)
mean(sin(x)) # answer (약 0.4610)

# (c)
n <- 100000
x <- runif(n)
f <- function(x) x^3
mean(f(x)) # answer (약 0.2501)

# (d)
n <- 100000
x <- runif(n)
mean(exp(x)) # answer (약 1.7170)

# (e)
n <- 100000
x <- rexp(n,0.1)
f <- function(x) 0.1*exp(-x/10)
g <- function(x) (x+1)^(-2)
mean(g(x)/f(x)) # answer (약 0.9998)

#-----------------------------------------------------------------#

## Problem 3

# (a)
n <- 100
U <- runif(n)
(Ber <- as.numeric(U<0.3)) # answer
par(mfrow=c(1,2))
hist(Ber,freq=F) 
hist(rbinom(100,1,0.3),freq=F) # compare histogram
par(mfrow=c(1,1))

# (b)
n <- 100*10
U <- runif(n)
Ber <- as.numeric(U<0.3)
dim(Ber)<-c(10,100)
(Bin <- apply(Ber,2,sum)) # answer
par(mfrow=c(1,2))
hist(Bin,freq=F) 
hist(rbinom(100,10,0.3),freq=F) # compare histogram
par(mfrow=c(1,1))

# (c)
n <- 100
U <- runif(n)
Finv <- function(x) -(1/2.5)*log(1-x)
(Exp <- Finv(U)) # answer
par(mfrow=c(1,2))
hist(Exp,freq=F) 
hist(rexp(100,2.5),freq=F) # compare histogram
par(mfrow=c(1,1))

# (d)
n <- 100*10
U <- runif(n)
Finv <- function(x) -log(1-x)
Exp <- Finv(U)
dim(Exp)<-c(10,100)
(Gam <- apply(Exp,2,sum)) # answer
par(mfrow=c(1,2))
hist(Gam,freq=F) 
hist(rgamma(100,shape=10,scale=1),freq=F) # compare histogram
par(mfrow=c(1,1))

#-----------------------------------------------------------------#

## Problem 4

# Graph of f(x)
x <- 1:10000/10000
f <- function(x) 2-4*abs(x-0.5)
plot(x,f(x),type='l')

# Sampling
n <- 200
a <- 1
b <- 2
U1 <- runif(n)*a
U2 <- runif(n)*b
index_ <- U2<f(U1)
plot(U1,U2)
points(U1[index_],U2[index_],col=2)
lines(x,f(x),col=3,lwd=2,lty=2)
X <- U1[index_] # answer
hist(X,freq=F)

#-----------------------------------------------------------------#

## Problem 5
n <- 100
U <- runif(n)
multi <- 1*(U>=0.1 & U<0.3)+2*(U>=0.3 & U<0.5)+3*(U>=0.5 & U<0.7)+4*(U>=0.7) # answer 
hist(multi,freq=F) 

#-----------------------------------------------------------------#

## Problem 6

# (a)
# function f and g
x <- 1:9999/10000
f <- function(x) {
  alpha=2.7 ; beta=6.3
  gamma(alpha+beta)/gamma(alpha)/gamma(beta) * x^(alpha-1) * (1-x)^(beta-1)
}
g <- function(x) {
  alpha=2 ; beta=6
  gamma(alpha+beta)/gamma(alpha)/gamma(beta) * x^(alpha-1) * (1-x)^(beta-1)
}
(M <- 1/min(g(x)/f(x)))
# Sampling from Gamma(2,1)
n <-200*2 
Exp <- -log(1-runif(n))
dim(Exp) <-c (2,200)
Gamma2 <- apply(Exp,2,sum)
# Sampling from Gamma(6,1)
n <- 200*6
Exp <- -log(1-runif(n))
dim(Exp) <- c(6,200)
Gamma6 <- apply(Exp,2,sum)
# Making Samples (Beta(2,6))
(Beta26 <- Gamma2/(Gamma2+Gamma6)) # answer
par(mfrow=c(1,2))
hist(Beta26,freq=F) 
hist(rbeta(200,shape1=2,shape2=6),freq=F) # compare histogram
par(mfrow=c(1,1))

# (b)
# Acceptance/Rejection Sampling
Y <- Beta26
U <- runif(200)
index_ <- c()
for (i in 1:200) {
  index_[i] <- (U[i] < f(Y[i])/(M*g(Y[i])))
}
(X <- Y[index_]) # answer
par(mfrow=c(1,2))
hist(X,freq=F) 
hist(rbeta(length(X),shape1=2.7,shape2=6.3),freq=F) # compare histogram
par(mfrow=c(1,1))

#-----------------------------------------------------------------#

## Problem 7

# (a)
n <- 10000
X <- rnorm(n)
index_ <- (X > -1.96)&(X < 1.96) 
mean(index_) # answer (약 0.95)

# (b)
n <- 10000
X <- rnorm(n)
index_ <- (X > 4.5)
mean(index_) # answer (약 0)

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
1-pnorm(4.5) # checking real answer

#-----------------------------------------------------------------#

## Problem 8

# (a)
n <- 1000*2
Exp <- rexp(n,rate=1)
dim(Exp) <- c(2,1000)
(Gamma2 <- apply(Exp,2,sum)) # answer
par(mfrow=c(1,2))
hist(Gamma2,freq=F)
hist(rgamma(1000,shape=2,scale=1),freq=F)
par(mfrow=c(1,1))

# (b)
Y <- Gamma2
f <- function(x) x*(x-1)*(x-2)*(x-3)*exp(-x)
g <- function(x) x*exp(-x)
mean(f(Y)/g(Y)) # answer (약 3.3886)
integrate(f,lower=0,upper=Inf) # checking real answer