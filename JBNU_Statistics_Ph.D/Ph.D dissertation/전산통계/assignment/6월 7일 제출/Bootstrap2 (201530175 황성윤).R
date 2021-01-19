#### Bootstrap

### Q1

install.packages("mvtnorm")
library(mvtnorm)

n <- 10
B <- 5000 ; M <- 1000
set.seed(1)
Norm.I=Piv.I=Per.I=matrix(0,M,2)

## Bootstrap
for (j in 1:M) {
  Y <- rmvnorm(n,mean=c(0,2),sigma=matrix(c(1,1.5,1.5,4),ncol=2)) 
  x1 <- Y[,1] ; x2 <- Y[,2]
  hatrow <- cor(x1,x2) 
  brow <- c()
  x1b <- matrix(0,B,n) ; x2b <- matrix(0,B,n)
  for (i in 1:B) {
    idx <- sample(1:n,n,replace=TRUE)
    x1b[i,] <- x1[idx] ; x2b[i,] <- x2[idx]
    brow[i] <- cor(x1b[i,],x2b[i,])
  }
  ## C.I.  
  # Normal interval
  Norm.I[j,] <-c(hatrow-qnorm(0.975)*sqrt(var(brow)),hatrow+qnorm(0.975)*sqrt(var(brow)))
  # Pivotal interval
  Piv.I[j,] <- c(2*hatrow-quantile(brow,0.975),2*hatrow-quantile(brow,0.025))
  # Percentile interval
  Per.I[j,] <- c(quantile(brow,0.025),quantile(brow,0.975))
}

## Average length of intervals
c(diff(colMeans(Norm.I)),diff(colMeans(Piv.I)),diff(colMeans(Per.I)))
# coverage proportions
c(mean(apply(Norm.I-0.75,1,prod)<0),mean(apply(Piv.I-0.75,1,prod)<0),mean(apply(Per.I-0.75,1,prod)<0))

### Q2

n <- 10
B <- 5000 ; M <- 1000
set.seed(1)
Norm.I=Piv.I=Per.I=matrix(0,M,2)

## Bootstrap
for (j in 1:M) {
  x1 <- rnorm(n) ; x2 <- x1^2
  hatrow <- cor(x1,x2) 
  brow <- c()
  x1b <- matrix(0,B,n) ; x2b <- matrix(0,B,n)
  for (i in 1:B) {
    idx <- sample(1:n,n,replace=TRUE)
    x1b[i,] <- x1[idx] ; x2b[i,] <- x2[idx]
    brow[i] <- cor(x1b[i,],x2b[i,])
  }
  ## C.I.  
  # Normal interval
  Norm.I[j,] <-c(hatrow-qnorm(0.975)*sqrt(var(brow)),hatrow+qnorm(0.975)*sqrt(var(brow)))
  # Pivotal interval
  Piv.I[j,] <- c(2*hatrow-quantile(brow,0.975),2*hatrow-quantile(brow,0.025))
  # Percentile interval
  Per.I[j,] <- c(quantile(brow,0.025),quantile(brow,0.975))
}

## Average length of intervals
c(diff(colMeans(Norm.I)),diff(colMeans(Piv.I)),diff(colMeans(Per.I)))
# coverage proportions
c(mean(apply(Norm.I-0,1,prod)<0),mean(apply(Piv.I-0,1,prod)<0),mean(apply(Per.I-0,1,prod)<0))
  







