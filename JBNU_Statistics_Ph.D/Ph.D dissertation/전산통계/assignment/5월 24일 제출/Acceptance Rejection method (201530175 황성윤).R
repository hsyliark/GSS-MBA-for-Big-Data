## Acceptance-Rejection method

set.seed(1)

i <- 0 ; X=c() ; n <- 0
pdf <- function(x) { 60*x^3*(1-x)^2 }
g <- function(x) { 0.5*4*x^3 + 0.5*3*(1-x)^2 } 

while(n < 500) {
  
  i <- i+1
  
  w <- rbinom(1, 1, 0.5)
  x1 <- rbeta(1, 4, 1) ; x2 <- rbeta(1, 1, 3)
  Y <- w*x1+(1-w)*x2 ; U <- runif(1)
  if (U <= pdf(Y)/(3.2*g(Y)) ) X[i] <- Y
  else X[i] <- -1
  
  n <- sum(X > 0)
}

rX <- X[X>0] ; n/length(X) ; 1/3.2


par(mfrow=c(1,2))
hist(rX) ; hist(rbeta(500,4,3))


par(mfrow=c(1,1))
x <- seq(0, 1, by=0.001)
y <- 3.2*g(x)
plot(y~x, main="Acceptance-Rejection method", type="n", xlim=c(0,1), ylim=c(0,6))
lines(x, pdf(x), lty=1, col="red", lwd=2)
lines(x, 3.2*g(x), lty=2, col="blue", lwd=2)
lines(x, rep(2.0736, length(x)), lty=3, col="black", lwd=2)
legend("topleft", legend=c("pdf f(x)", "envelop c*g(x)", "envelop c*1"), col=c("red", "blue", "black"), lty=c(1,2,3), lwd=c(2,2,2))
