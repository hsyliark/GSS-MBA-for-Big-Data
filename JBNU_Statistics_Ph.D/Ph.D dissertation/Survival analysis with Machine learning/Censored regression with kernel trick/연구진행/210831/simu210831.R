# Loading packages
library(survival)
library(KernSmooth)
library(np)
library(locfit)
library(latex2exp)


# Kaplan-Meier estimator
km.surv <- function(time, cens) {
  
  n <- length(time)
  Y <- matrix(rep(time, times = n), ncol = n)
  IY <- (Y > matrix(rep(time, each = n), ncol = n))
  ny <- apply(IY, 2, sum)
  ny <- ((1 + ny)/(2 + ny))
  NY <- matrix(rep(ny, times = n), ncol = n)
  I <- (matrix(rep(cens, times = n), ncol = n) + IY == 0)
  est <- apply(NY^I, 2, prod)
  
  return(est)
  
}  


# Synthetic response (KSV)
x.i <- seq(0, 1, length.out=100)
e.i <- rnorm(100, mean=0, sd=1)
y.i <- sin(2*pi*x.i)+e.i
pred.value <-  sin(2*pi*seq(0,1,length.out=100))
cen.10 <-rnorm(100, mean=sin(2*pi*x.i)+0.84^2, sd=1) # censoring values
pre.T.10 <- pmin(y.i, cen.10)
delta.10 <- (y.i <= cen.10) * 1
g.10 <- km.surv(pre.T.10, delta.10)
y.i.10 <- ifelse(pre.T.10 <= quantile(pre.T.10, probs=0.98), pre.T.10*delta.10/g.10, 0)


# Gaussian (Radial basis) kernel
my.kernel.matrix <- function(dat.train, dat.test) {
  
  # dat.train : Data frame for training with a response variable is appeared in the first column...
  # dat.test : Data frame for testing with a response variable is appeared in the first column...
  
  data1 <- as.data.frame(dat.train)
  data2 <- as.data.frame(dat.test)
  n <- nrow(data1)
  
  # Training data
  X.train <- as.matrix(data1[,-1])
  y.train <- as.matrix(data1[,1])
  
  # Test data
  X.test <- as.matrix(data2[,-1])
  y.test <- as.matrix(data2[,1])
  
  X1 <- rbind(X.train, X.test)
  
  
  p <- ncol(X1)
  sigma <- 1/p
  n1 <- nrow(X1)
  D <- as.matrix(dist(X1,method="euclidean",p=2)) # Euclidean distance 
  K <- exp( -sigma * D^2 )  # Gaussian kernel
  K.train <- K[1:n,1:n]
  K.test <- K[1:n,(n+1):n1]
  
  return(list(K.train=K.train, K.test=K.test, K=K, y.train=y.train, y.test=y.test))
  
}


# 
