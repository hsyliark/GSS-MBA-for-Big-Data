## Importance Sampling


#23.8.1

set.seed(1)
n <- 1000
U <- runif(n)
X <- ((((2*cos(pi*U/2))/(3*(1-U^2)))-(2/pi))^2)*(3*(1-U^2)/2)
res <- mean(X)
res


#23.8.3

set.seed(1)
Y <- rbinom(1000, 100, 0.75)
prob <- mean(I(Y>=75)*3^(100-2*Y))
prob
1-pbinom(74, 100, 0.25)


## Gibbs Sampling

set.seed(1)
mu0 <- 1 ; tau0 <- 1 ; n <- 100
sigma20 <- 1/tau0
mu <- mu0 ; tau <- tau0
y <- rnorm(n,mu0,sqrt(sigma20))
N <- 4999
for (i in 2:(N+1)) {
  mu[i] <- rnorm(1, mean=mean(y), sd=sqrt(1/(n*tau[i-1])))
  tau[i] <- rgamma(1, shape=n/2, scale=sum((y-mu[i])^2)/2)
}
sigma2 <- 1/tau
mean(mu[501:5000])
mean(sigma2[501:5000])
par(mfrow=c(1,2))
hist(mu[501:5000])
hist(sigma2[501:5000])

  

