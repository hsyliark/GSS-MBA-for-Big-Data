### Generate random sample from normal distribution

x <- rnorm(100, 0, 2)
mean(x)
var(x)


### Shortest confidence interval


## 1. Bisection method 

lambda <- 0.00102464
a <- 0 ; b <- 99
f <- function(x) {dchisq(x, 99)-lambda}
sign(f(0)*f(99))

tol <- 10^{-8}
for (i in 1:100) {
  
  c <- (a+b)/2
  newa <- c*(f(c)*f(b)<0)+a*(1-(f(c)*f(b)<0))
  newb <- c*(f(c)*f(a)<0)+b*(1-(f(c)*f(a)<0))
  a <- newa ; b <- newb
  err <- abs(a-b)
  print(c(i,a,b))
  
  if (err<tol) break
  
}

(c <- (a+b)/2)


lambda <- 0.00102464
a <- 99 ; b <- 150
f <- function(x) {dchisq(x, 99)-lambda}
sign(f(99)*f(150))

tol <- 10^{-8}
for (i in 1:100) {
  
  d <- (a+b)/2
  newa <- d*(f(d)*f(b)<0)+a*(1-(f(d)*f(b)<0))
  newb <- d*(f(d)*f(a)<0)+b*(1-(f(d)*f(a)<0))
  a <- newa ; b <- newb
  err <- abs(a-b)
  print(c(i,a,b))
  
  if (err<tol) break
  
}

(d <- (a+b)/2)

pchisq(d, 99) - pchisq(c, 99)


## 2. Result

# Shortest confidence interval
S.C.I <- c((99*var(x))/d, (99*var(x))/c) 

# General confidence interval
G.C.I <- c((99*var(x))/qchisq(0.995, 99), (99*var(x))/qchisq(0.005, 99)) 

# Compare size

d - c
S.C.I[2] - S.C.I[1]

qchisq(0.995, 99) - qchisq(0.005, 99)
G.C.I[2] - G.C.I[1]


## 3. Graph

x1 <- rchisq(100, 99)
y <- seq(0, 150, by=0.001)
hist(x1, main="Chi-square distribution with df=99", col="yellow", probability=T)
lines(y, dchisq(y, 99), col="black", lwd=2)
abline(h=lambda, col="green")
abline(v=c, col="red") ; abline(v=d, col="red")
abline(v=qchisq(0.005, 99), col="blue") ; abline(v=qchisq(0.995, 99), col="blue")
legend("topright", c("S.C.I","G.C.I"), col=c("red","blue"), lty=c(1,1))
