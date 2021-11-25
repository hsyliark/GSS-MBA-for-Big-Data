## Loading packages

library(survival)
library(KernSmooth)
library(np)
library(locfit)
library(latex2exp)

## Loading data
dat1 <- read.csv("",sep=",",header=T)


#### Making some functions

## Kaplan-Meier estimator

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





#------------------------#------------------------#------------------------#------------------------



### 1. Making kernel matrix



# 1. Linear kernel
my.kernel.matrix1 <- function(dat.train, dat.test) {
  
  # dat.train : Data frame for training with a response variable is appeared in the first column...
  # dat.test : Data frame for testing with a response variable is appeared in the first column...
  
  data1 <- as.data.frame(dat.train)
  data2 <- as.data.frame(dat.test)
  n <- nrow(data1)
  
  # Training data
  X.train <- as.matrix(data1[,c(-1,-ncol(data1))])
  y.train.s <- as.matrix(data1[,1]) # synthetic
  y.train <- as.matrix(data1[,ncol(data1)]) # original
  
  # Test data
  X.test <- as.matrix(data2[,c(-1,-ncol(data2))])
  y.test.s <- as.matrix(data2[,1])
  y.test <- as.matrix(data2[,ncol(data2)])
  
  X1 <- rbind(X.train, X.test)
  
  
  p <- ncol(X1)
  sigma <- 1/p
  n1 <- nrow(X1)
  K <- as.matrix(X1%*%t(X1)) # Inner product
  K.train <- K[1:n,1:n]
  K.test <- K[1:n,(n+1):n1]
  
  return(list(K.train=K.train, K.test=K.test, K=K, y.train.s=y.train.s, y.test.s=y.test.s,
              y.train=y.train, y.test=y.test))
  
}

# 2. Gaussian kernel (Radial Basis kernel)
my.kernel.matrix2 <- function(dat.train, dat.test) {
  
  # dat.train : Data frame for training with a response variable is appeared in the first column...
  # dat.test : Data frame for testing with a response variable is appeared in the first column...
  
  data1 <- as.data.frame(dat.train)
  data2 <- as.data.frame(dat.test)
  n <- nrow(data1)
  
  # Training data
  X.train <- as.matrix(data1[,c(-1,-ncol(data1))])
  y.train.s <- as.matrix(data1[,1]) # synthetic
  y.train <- as.matrix(data1[,ncol(data1)]) # original
  
  # Test data
  X.test <- as.matrix(data2[,c(-1,-ncol(data2))])
  y.test.s <- as.matrix(data2[,1])
  y.test <- as.matrix(data2[,ncol(data2)])
  
  X1 <- rbind(X.train, X.test)
  
  
  p <- ncol(X1)
  sigma <- 1/p
  n1 <- nrow(X1)
  D <- as.matrix(dist(X1,method="euclidean",p=2)) # Euclidean distance 
  K <- exp( -sigma * D^2 )  # Gaussian kernel
  K.train <- K[1:n,1:n]
  K.test <- K[1:n,(n+1):n1]
  
  return(list(K.train=K.train, K.test=K.test, K=K, y.train.s=y.train.s, y.test.s=y.test.s,
              y.train=y.train, y.test=y.test))
  
}





### 2. Calculate coefficients using theory



my.kernel.regression <- function(y, K, lambda) {
  
  # y : response variable
  # K : matrix from Linear kernel or Gaussian kernel transformation
  # lambda : penalty parameter
  
  if(lambda <= 0) 
    stop("Lambda is non-positive value. Please insert positive value of lambda.")
  
  y <- as.matrix(y) ; K <-as.matrix(K)
  
  d.hat <- solve(K + lambda*diag(x=1, nrow=length(y), ncol=length(y)))%*%y
  
  y.hat <- K%*%d.hat
  y.hat <- ifelse(y.hat<0, 0, y.hat) # non-zero of survival time
  
  return(list(d.hat=d.hat, y.hat=y.hat))
  
}





### 3. Making function for fitting



fit.kernel <- function(y.train.s, y.train, K.train, lambda) {
  
  # y.train.s : Synthetic response Y* of training data
  # y.train : Original response Y of training data
  # K.train : Kernel matrix from training data 
  # lambda : penalty parameter
  
  if(lambda <= 0) 
    stop("Lambda is non-positive value. Please insert positive value of lambda.")
  
  
  y.train.s <- as.matrix(y.train.s) # Synthetic
  y.train <- as.matrix(y.train) # Original
  K.train <- as.matrix(K.train) 
  
  g <- my.kernel.regression(y.train.s, K.train, lambda)
  
  d.hat <- g$d.hat
  y.pred <- g$y.hat
  
  # rmse <- sqrt(sum((y.train - y.pred)^2)/length(y.train))
  rmse1 <- sqrt(sum((y.train.s - y.pred)^2)/length(y.train.s))
  rmse2 <- sqrt(sum((y.train - y.pred)^2)/length(y.train))
  
  return(list(d.hat=d.hat, y.train.s=y.train.s, y.train=y.train, y.pred=y.pred, 
              rmse1=rmse1, rmse2=rmse2))
  
}






### 4. Making function for predict



pred.kernel <- function(y.test.s, y.test, K.test, d.hat) {
  
  # y.test : Synthetic response Y* of test data
  # y.test.pred : Original response Y of training data
  # K.test : Kernel matrix from test data 
  # d.hat : Estimator of vector d from training data
  
  
  y.test.s <- as.matrix(y.test.s)
  y.test <- as.matrix(y.test)
  K.test <- as.matrix(K.test)
  d.hat <- as.matrix(d.hat)
  
  
  y.hat <- t(K.test)%*%d.hat
  y.hat <- ifelse(y.hat<0, 0, y.hat) # non-zero of survival time
  
  # rmse <- sqrt(sum((y.test - y.hat)^2)/length(y.test))
  rmse1 <- sqrt(sum((y.test.s - y.hat)^2)/length(y.test.s))
  rmse2 <- sqrt(sum((y.test - y.hat)^2)/length(y.test))
  
  return(list(y.test.s=y.test.s, y.test=y.test, rmse1=rmse1, rmse2=rmse2,
              y.hat=y.hat))
  
}






### 5. Making function for K-fold crossvalidation



cv.kernel <- function(y.train.s, y.train, K.train, k, grid.l) {
  
  
  # y.train.s : Synthetic response Y* of training data
  # y.train : Original response Y of training data
  # K.train : Kernel matrix from training data 
  # k : number of criterion for K-fold crossvalidation
  # grid.l : The row of penalty parameter lambda
  
  check <- (grid.l > 0)
  n.check <- length(check)
  
  if(sum(check) != n.check)
    stop("Some of lambda's values are non-positive.
         Please insert positive values of lambda vector...","\n")
  
  
  lambda <- grid.l
  r <- length(lambda)
  
  
  K.sim <- as.matrix(K.train)
  y.sim.s <- as.matrix(y.train.s)
  y.sim <- as.matrix(y.train)
  n <- nrow(K.sim)
  
  cv.index <- sample(1:n,n,replace=F)  
  cv.rmse1 <- NULL  
  cv.rmse2 <- NULL
  
  cat("K-fold crossvalidation is start...","\n")
  
  
  for (j in 1:r) {
    
    rmse1 <- NULL # Root mean squared error with synthetic response
    rmse2 <- NULL # Root mean squared error with original response
    
    
    for (i in 0:(k-1)) {
      
      
      test.index <- cv.index[(1:n)%/%k==i]
      
      K.sim.train <- K.sim[-test.index, -test.index] ; K.sim.test <- K.sim[-test.index, test.index]
      y.sim.train.s <- y.sim.s[-test.index,] ; y.sim.test.s <- y.sim.s[test.index,]
      y.sim.train <- y.sim[-test.index,] ; y.sim.test <- y.sim[test.index,]
      test.size <- length(test.index)
      
      
      a1 <- fit.kernel(y.sim.train.s, y.sim.train, K.sim.train, lambda[j])
      train.d.hat <- a1$d.hat
      
      # a2 <- pred.kernel(y.sim.test, K.sim.test, train.d.hat)
      a2 <- pred.kernel(y.sim.test.s, y.sim.test, K.sim.test, train.d.hat)
      test.y.hat <- a2$y.hat      
      
      
      # rmse <- c(rmse, sqrt(sum((y.sim.test - test.y.hat)^2)/length(y.sim.test)) )
      rmse1 <- c(rmse1, sqrt(sum((y.sim.test.s - test.y.hat)^2)/length(y.sim.test.s)) )
      rmse2 <- c(rmse2, sqrt(sum((y.sim.test - test.y.hat)^2)/length(y.sim.test)) )
      
      
    }
    
    cv.rmse1 <- rbind(cv.rmse1, rmse1)
    cv.rmse2 <- rbind(cv.rmse2, rmse2)
    cat(j, ",")
  }
  
  cat("\n","K-fold crossvalidation complete...")
  
  
  return(list(lambda=grid.l, cv.rmse1=cv.rmse1, cv.rmse2=cv.rmse2))
  
  
}





#------------------------#------------------------#------------------------#------------------------
#------------------------#------------------------#------------------------#------------------------




#### 4 method simulation

# LKR1 : Linear Kernel Regression with Synthetic Response Y*
# LKRS1 : Linear Kernel Regression with Sub-sampling and Synthetic Response Y*
# LKRB1 : Linear Kernel Regression with Bagging and Synthetic Response Y*
# LKRR1 : Linear Kernel Regression with Random Forest and Synthetic Response Y*
  
# LKR2 : Linear Kernel Regression with Generated(Original) Response Y
# LKRS2 : Linear Kernel Regression with Sub-sampling and Generated(Original) Response Y
# LKRB2 : Linear Kernel Regression with Bagging and Generated(Original) Response Y
# LKRR2 : Linear Kernel Regression with Random Forest and Generated(Original) Response Y

# GKR1 : Gaussian Kernel Regression with Synthetic Response Y*
# GKRS1 : Gaussian Kernel Regression with Sub-sampling and Synthetic Response Y*
# GKRB1 : Gaussian Kernel Regression with Bagging and Synthetic Response Y*
# GKRR1 : Gaussian Kernel Regression with Random Forest and Synthetic Response Y*
  
# GKR2 : Gaussian Kernel Regression with Generated(Original) Response Y 
# GKRS2 : Gaussian Kernel Regression with Sub-sampling and Generated(Original) Response Y
# GKRB2 : Gaussian Kernel Regression with Bagging and Generated(Original) Response Y
# GKRR2 : Gaussian Kernel Regression with Random Forest and Generated(Original) Response Y





### Function for 1 independent variable

fit.ftn1 <- function(dat.sim) {
  
  
  
  ### 1. Kernel ridge regression (KR)
  
  LKR1 <- c(rep(0,100))
  LKR2 <- c(rep(0,100))
  
  for (i in 1:100) {
    
    # Making simulation data
    
    set.seed(i)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10)
    
    
    # 5-fold crossvalidation
    
    u <- my.kernel.matrix1(train.sim, test.sim)
    K.train <- u$K.train ; K.test <- u$K.test  
    y.train.s <- u$y.train.s ; y.test.s <- u$y.test.s
    y.train <- u$y.train ; y.test <- u$y.test
    
    grid.l <- 10^seq(-3,2,length=10)
    
    h <- cv.kernel(y.train.s, y.train, K.train, 5, grid.l) 
    
    # Fitting 
    
    mean.rmse1 <- rowMeans(h$cv.rmse1)
    idx1 <- which.min(mean.rmse1)
    best.lam1 <- max(h$lambda[ mean.rmse1 == mean.rmse1[idx1] ])
    
    h1_1 <- fit.kernel(y.train.s, y.train, K.train, best.lam1)
    
    mean.rmse2 <- rowMeans(h$cv.rmse2)
    idx2 <- which.min(mean.rmse2)
    best.lam2 <- max(h$lambda[ mean.rmse2 == mean.rmse2[idx2] ])
    
    h1_2 <- fit.kernel(y.train.s, y.train, K.train, best.lam2)
    
    # Calculate test RMSE
    
    sim.d.hat1 <- h1_1$d.hat
    h2_1 <- pred.kernel(y.test.s, y.test, K.test, sim.d.hat1)
    
    sim.d.hat2 <- h1_2$d.hat
    h2_2 <- pred.kernel(y.test.s, y.test, K.test, sim.d.hat2)
    
    LKR1[i] <- h2_1$rmse1
    LKR2[i] <- h2_2$rmse2
    
  }
  
  LKR1
  LKR2
  boxplot(LKR1)
  boxplot(LKR2)
  dat1_1 <- data.frame(RMSE=LKR1, method=rep("a.LKR1",100), number=as.character(rep(n.train,100)))
  dat1_2 <- data.frame(RMSE=LKR2, method=rep("i.LKR2",100), number=as.character(rep(n.train,100)))
  dat1 <- rbind(dat1_1, dat1_2)
  
  GKR1 <- c(rep(0,100))
  GKR2 <- c(rep(0,100))
  
  for (i in 1:100) {
    
    # Making simulation data
    
    set.seed(i)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10) 
    
    
    # 5-fold crossvalidation
    
    u <- my.kernel.matrix2(train.sim, test.sim)
    K.train <- u$K.train ; K.test <- u$K.test  
    y.train.s <- u$y.train.s ; y.test.s <- u$y.test.s
    y.train <- u$y.train ; y.test <- u$y.test
    
    grid.l <- 10^seq(-3,2,length=10)
    
    h <- cv.kernel(y.train.s, y.train, K.train, 5, grid.l) 
    
    # Fitting 
    
    mean.rmse1 <- rowMeans(h$cv.rmse1)
    idx1 <- which.min(mean.rmse1)
    best.lam1 <- max(h$lambda[ mean.rmse1 == mean.rmse1[idx1] ])
    
    h1_1 <- fit.kernel(y.train.s, y.train, K.train, best.lam1)
    
    mean.rmse2 <- rowMeans(h$cv.rmse2)
    idx2 <- which.min(mean.rmse2)
    best.lam2 <- max(h$lambda[ mean.rmse2 == mean.rmse2[idx2] ])
    
    h1_2 <- fit.kernel(y.train.s, y.train, K.train, best.lam2)
    
    # Calculate test RMSE
    
    sim.d.hat1 <- h1_1$d.hat
    h2_1 <- pred.kernel(y.test.s, y.test, K.test, sim.d.hat1)
    
    sim.d.hat2 <- h1_2$d.hat
    h2_2 <- pred.kernel(y.test.s, y.test, K.test, sim.d.hat2)
    
    GKR1[i] <- h2_1$rmse1
    GKR2[i] <- h2_2$rmse2
    
  }
  
  GKR1
  GKR2
  boxplot(GKR1)
  boxplot(GKR2)
  dat2_1 <- data.frame(RMSE=GKR1, method=rep("e.GKR1",100), number=as.character(rep(n.train,100)))
  dat2_2 <- data.frame(RMSE=GKR2, method=rep("m.GKR2",100), number=as.character(rep(n.train,100)))
  dat2 <- rbind(dat2_1, dat2_2)
  
  
  
  
  ### 2. Kernel ridge regression using Sub-sampling (KRS)
  
  LKRS1 <- c(rep(0,100))
  LKRS2 <- c(rep(0,100))
  
  for (i in 1:100) {
    
    # Making simulation data
    
    set.seed(i)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10)
    
    u <- my.kernel.matrix1(train.sim, test.sim)
    K.train <- u$K.train ; y.train.s <- u$y.train.s 
    K <- u$K ; K.test <- u$K.test ; y.test.s <- u$y.test.s
    y.train <- u$y.train ; y.test <- u$y.test
    
    # Choosing optimal lambda
    
    
    n <- nrow(train.sim) ; n1 <- sum(nrow(train.sim), nrow(test.sim))
    res.rmse1.c <- NULL ; res.lam1.c <- NULL ; res.index1.c <- NULL
    res.rmse2.c <- NULL ; res.lam2.c <- NULL ; res.index2.c <- NULL
    
    for (j in 1:50) {
      
      n <- nrow(train.sim)
      index1.Kc <- sample(1:n, round(0.7*n), replace=F)
      index2.Kc <- index1.Kc
      Kc.train <- K.train[index1.Kc,index1.Kc] ; Kc.test <- K.train[index1.Kc,-index1.Kc]
      yc.train.s <- y.train.s[index1.Kc] ; yc.test.s <- y.train.s[-index1.Kc]
      yc.train <- y.train[index1.Kc] ; yc.test <- y.train[-index1.Kc]
      
      grid.l <- 10^seq(-3,2,length=10)
      
      hc <- cv.kernel(yc.train.s, yc.train, Kc.train, 5, grid.l) # 5-fold crossvalidation 
      
      # synthetic
      mean.rmse1.c <- rowMeans(hc$cv.rmse1)
      idx1.c <- which.min(mean.rmse1.c)
      best.lam1.c <- max(hc$lambda[ mean.rmse1.c == mean.rmse1.c[idx1.c] ])
      
      h1_1.c <- fit.kernel(yc.train.s, yc.train, Kc.train, best.lam1.c) # Fitting
      
      sim.d.hat1.c <- h1_1.c$d.hat
      
      h2_1.c <- pred.kernel(yc.test.s, yc.test, Kc.test, sim.d.hat1.c) # Testing
      rmse1.c <- h2_1.c$rmse1
      
      res.rmse1.c <- c(res.rmse1.c, rmse1.c)
      res.lam1.c <- c(res.lam1.c, best.lam1.c)
      res.index1.c <- cbind(res.index1.c, index1.Kc)
      
      # original
      mean.rmse2.c <- rowMeans(hc$cv.rmse2)
      idx2.c <- which.min(mean.rmse2.c)
      best.lam2.c <- max(hc$lambda[ mean.rmse2.c == mean.rmse2.c[idx2.c] ])
      
      h1_2.c <- fit.kernel(yc.train.s, yc.train, Kc.train, best.lam2.c) # Fitting
      
      sim.d.hat2.c <- h1_2.c$d.hat
      
      h2_2.c <- pred.kernel(yc.test.s, yc.test, Kc.test, sim.d.hat2.c) # Testing
      rmse2.c <- h2_2.c$rmse1
      
      res.rmse2.c <- c(res.rmse2.c, rmse2.c)
      res.lam2.c <- c(res.lam2.c, best.lam2.c)
      res.index2.c <- cbind(res.index2.c, index2.Kc)
      
      
    }
    
    # Fitting  
    
    # synthetic
    res.lambda1 <- res.lam1.c[which.min(res.rmse1.c)]
    res.index1 <- res.index1.c[, which.min(res.rmse1.c)]
    res.K.train1 <- K[res.index1, res.index1]
    res.y.train.s1 <- y.train[res.index1]
    res.y.train1 <- y.train[res.index1]
    K.test1 <- K[res.index1, (n+1):n1]
    
    h1_1 <- fit.kernel(res.y.train.s1, res.y.train1, res.K.train1, res.lambda1)
    
    res.d.hat1 <- h1_1$d.hat
    
    # Calculate test mean square error
    
    h2_1 <- pred.kernel(y.test.s, y.test, K.test1, res.d.hat1)
    LKRS1[i] <- h2_1$rmse1
    
    # original
    res.lambda2 <- res.lam2.c[which.min(res.rmse2.c)]
    res.index2 <- res.index2.c[, which.min(res.rmse2.c)]
    res.K.train2 <- K[res.index2, res.index2]
    res.y.train.s2 <- y.train[res.index2]
    res.y.train2 <- y.train[res.index2]
    K.test2 <- K[res.index2, (n+1):n1]
    
    h1_2 <- fit.kernel(res.y.train.s2, res.y.train2, res.K.train2, res.lambda2)
    
    res.d.hat2 <- h1_2$d.hat
    
    # Calculate test mean square error
    
    h2_2 <- pred.kernel(y.test.s, y.test, K.test2, res.d.hat2)
    LKRS2[i] <- h2_2$rmse2
    
  }
  
  LKRS1
  LKRS2
  boxplot(LKRS1)
  boxplot(LKRS2)
  dat3_1 <- data.frame(RMSE=LKRS1, method=rep("b.LKRS1",100), number=as.character(rep(n.train,100)))
  dat3_2 <- data.frame(RMSE=LKRS2, method=rep("j.LKRS2",100), number=as.character(rep(n.train,100)))
  dat3 <- rbind(dat3_1, dat3_2)
  
  GKRS1 <- c(rep(0,100))
  GKRS2 <- c(rep(0,100))
  
  for (i in 1:100) {
    
    # Making simulation data
    
    set.seed(i)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10)
    
    u <- my.kernel.matrix2(train.sim, test.sim)
    K.train <- u$K.train ; y.train.s <- u$y.train.s 
    K <- u$K ; K.test <- u$K.test ; y.test.s <- u$y.test.s
    y.train <- u$y.train ; y.test <- u$y.test
    
    # Choosing optimal lambda
    
    
    n <- nrow(train.sim) ; n1 <- sum(nrow(train.sim), nrow(test.sim))
    res.rmse1.c <- NULL ; res.lam1.c <- NULL ; res.index1.c <- NULL
    res.rmse2.c <- NULL ; res.lam2.c <- NULL ; res.index2.c <- NULL
    
    for (j in 1:50) {
      
      n <- nrow(train.sim)
      index1.Kc <- sample(1:n, round(0.7*n), replace=F)
      index2.Kc <- index1.Kc
      Kc.train <- K.train[index1.Kc,index1.Kc] ; Kc.test <- K.train[index1.Kc,-index1.Kc]
      yc.train.s <- y.train.s[index1.Kc] ; yc.test.s <- y.train.s[-index1.Kc]
      yc.train <- y.train[index1.Kc] ; yc.test <- y.train[-index1.Kc]
      
      grid.l <- 10^seq(-3,2,length=10)
      
      hc <- cv.kernel(yc.train.s, yc.train, Kc.train, 5, grid.l) # 5-fold crossvalidation 
      
      # synthetic
      mean.rmse1.c <- rowMeans(hc$cv.rmse1)
      idx1.c <- which.min(mean.rmse1.c)
      best.lam1.c <- max(hc$lambda[ mean.rmse1.c == mean.rmse1.c[idx1.c] ])
      
      h1_1.c <- fit.kernel(yc.train.s, yc.train, Kc.train, best.lam1.c) # Fitting
      
      sim.d.hat1.c <- h1_1.c$d.hat
      
      h2_1.c <- pred.kernel(yc.test.s, yc.test, Kc.test, sim.d.hat1.c) # Testing
      rmse1.c <- h2_1.c$rmse1
      
      res.rmse1.c <- c(res.rmse1.c, rmse1.c)
      res.lam1.c <- c(res.lam1.c, best.lam1.c)
      res.index1.c <- cbind(res.index1.c, index1.Kc)
      
      # original
      mean.rmse2.c <- rowMeans(hc$cv.rmse2)
      idx2.c <- which.min(mean.rmse2.c)
      best.lam2.c <- max(hc$lambda[ mean.rmse2.c == mean.rmse2.c[idx2.c] ])
      
      h1_2.c <- fit.kernel(yc.train.s, yc.train, Kc.train, best.lam2.c) # Fitting
      
      sim.d.hat2.c <- h1_2.c$d.hat
      
      h2_2.c <- pred.kernel(yc.test.s, yc.test, Kc.test, sim.d.hat2.c) # Testing
      rmse2.c <- h2_2.c$rmse1
      
      res.rmse2.c <- c(res.rmse2.c, rmse2.c)
      res.lam2.c <- c(res.lam2.c, best.lam2.c)
      res.index2.c <- cbind(res.index2.c, index2.Kc)
      
      
    }
    
    # Fitting  
    
    # synthetic
    res.lambda1 <- res.lam1.c[which.min(res.rmse1.c)]
    res.index1 <- res.index1.c[, which.min(res.rmse1.c)]
    res.K.train1 <- K[res.index1, res.index1]
    res.y.train.s1 <- y.train[res.index1]
    res.y.train1 <- y.train[res.index1]
    K.test1 <- K[res.index1, (n+1):n1]
    
    h1_1 <- fit.kernel(res.y.train.s1, res.y.train1, res.K.train1, res.lambda1)
    
    res.d.hat1 <- h1_1$d.hat
    
    # Calculate test mean square error
    
    h2_1 <- pred.kernel(y.test.s, y.test, K.test1, res.d.hat1)
    GKRS1[i] <- h2_1$rmse1
    
    # original
    res.lambda2 <- res.lam2.c[which.min(res.rmse2.c)]
    res.index2 <- res.index2.c[, which.min(res.rmse2.c)]
    res.K.train2 <- K[res.index2, res.index2]
    res.y.train.s2 <- y.train[res.index2]
    res.y.train2 <- y.train[res.index2]
    K.test2 <- K[res.index2, (n+1):n1]
    
    h1_2 <- fit.kernel(res.y.train.s2, res.y.train2, res.K.train2, res.lambda2)
    
    res.d.hat2 <- h1_2$d.hat
    
    # Calculate test mean square error
    
    h2_2 <- pred.kernel(y.test.s, y.test, K.test2, res.d.hat2)
    GKRS2[i] <- h2_2$rmse2
    
  }
  
  GKRS1
  GKRS2
  boxplot(GKRS1)
  boxplot(GKRS2)
  dat4_1 <- data.frame(RMSE=GKRS1, method=rep("f.GKRS1",100), number=as.character(rep(n.train,100)))
  dat4_2 <- data.frame(RMSE=GKRS2, method=rep("n.GKRS2",100), number=as.character(rep(n.train,100)))
  dat4 <- rbind(dat4_1, dat4_2)
  
  
  
  
  ### 3. Kernel ridge regression using Bagging (KRB)
  
  LKRB1 <- c(rep(0,100))
  LKRB2 <- c(rep(0,100))
  
  for (s in 1:100) {
    
    # Making simulation data
    
    set.seed(s)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    test.sim.y <- test.sim[,1] 
    n.train <- round(nrow(dat.sim)*7/10)
    
    u <- my.kernel.matrix1(train.sim, test.sim)
    K.train <- u$K.train ; y.train.s <- u$y.train.s 
    K <- u$K ; K.test <- u$K.test ; y.test.s <- u$y.test.s
    y.train <- u$y.train ; y.test <- u$y.test
    n <- nrow(train.sim)
    n1 <- sum(nrow(train.sim), nrow(test.sim))
    
    
    # Choosing best lambda
    
    res.rmse1 <- NULL
    res.rmse2 <- NULL
    grid.l <- 10^seq(-3,2,length=10)
    
    for (i in 1:50) {
      
      boot.index <- sample(1:n, n, replace=T)
      boot.K.train <- K.train[boot.index, boot.index]
      boot.K.test <- K.train[boot.index, -boot.index]
      boot.y.train.s <- y.train.s[boot.index]
      boot.y.test.s <- y.train.s[-boot.index]
      boot.y.train <- y.train[boot.index]
      boot.y.test <- y.train[-boot.index]
      
      # synthetic
      boot.rmse1 <- c(rep(0,10))
      # original
      boot.rmse2 <- c(rep(0,10))
      
      for (j in 1:10) {
        
        h1.boot <- fit.kernel(boot.y.train.s, boot.y.train, boot.K.train, grid.l[j])
        boot.d.hat <- h1.boot$d.hat
        h2.boot <- pred.kernel(boot.y.test.s, boot.y.test, boot.K.test, boot.d.hat)
        boot.rmse1[j] <- h2.boot$rmse1
        boot.rmse2[j] <- h2.boot$rmse2
        
      }
      
      res.rmse1 <- rbind(res.rmse1, boot.rmse1)
      res.rmse2 <- rbind(res.rmse2, boot.rmse2)
      
    }
    
    #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
    best.lambda1 <- max(grid.l[colMeans(res.rmse1) == min(colMeans(res.rmse1))])
    best.lambda2 <- max(grid.l[colMeans(res.rmse2) == min(colMeans(res.rmse2))])
    
    # Bootstrapping
    
    boots.y1 <- NULL
    boots.y2 <- NULL
    
    for (r in 1:100) {
      
      boots.index <- sample(1:n, n, replace=T)
      boots.K.train <- K[boots.index, boots.index]
      boots.K.test <- K[boots.index, (n+1):n1]
      boots.y.train.s <- y.train.s[boots.index]
      boots.y.test.s <- y.test.s
      boots.y.train <- y.train[boots.index]
      boots.y.test <- y.test
      
      # synthetic
      h1_1.boots <- fit.kernel(boots.y.train.s, boots.y.train, boots.K.train, best.lambda1)
      boots.d.hat1 <- h1_1.boots$d.hat
      h2_1.boots <- pred.kernel(boots.y.test.s, boots.y.test, boots.K.test, boots.d.hat1)
      boots.y.hat1 <- h2_1.boots$y.hat
      boots.y1 <- cbind(boots.y1, boots.y.hat1)
      
      # original
      h1_2.boots <- fit.kernel(boots.y.train.s, boots.y.train, boots.K.train, best.lambda2)
      boots.d.hat2 <- h1_2.boots$d.hat
      h2_2.boots <- pred.kernel(boots.y.test.s, boots.y.test, boots.K.test, boots.d.hat2)
      boots.y.hat2 <- h2_2.boots$y.hat
      boots.y2 <- cbind(boots.y2, boots.y.hat2)
      
    }
    cat(s, ",")
    
    # Bagging estimate 
    
    y.hat1.bag <- rowMeans(boots.y1)
    y.hat1.bag <- ifelse(y.hat1.bag<0, 0, y.hat1.bag)
    y.hat2.bag <- rowMeans(boots.y2)
    y.hat2.bag <- ifelse(y.hat2.bag<0, 0, y.hat2.bag)
    
    # KRB[s] <- sqrt(sum((y.test - y.hat.bag)^2)/length(y.test))
    LKRB1[s] <- sqrt(sum((y.test.s - y.hat1.bag)^2)/length(y.test.s))
    LKRB2[s] <- sqrt(sum((y.test - y.hat2.bag)^2)/length(y.test))
    
  }
  
  LKRB1
  LKRB2
  boxplot(LKRB1)
  boxplot(LKRB2)
  dat5_1 <- data.frame(RMSE=LKRB1, method=rep("c.LKRB1",100), number=as.character(rep(n.train,100)))
  dat5_2 <- data.frame(RMSE=LKRB2, method=rep("k.LKRB2",100), number=as.character(rep(n.train,100)))
  dat5 <- rbind(dat5_1, dat5_2)
  
  GKRB1 <- c(rep(0,100))
  GKRB2 <- c(rep(0,100))
  
  for (s in 1:100) {
    
    # Making simulation data
    
    set.seed(s)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    test.sim.y <- test.sim[,1] 
    n.train <- round(nrow(dat.sim)*7/10)
    
    u <- my.kernel.matrix2(train.sim, test.sim)
    K.train <- u$K.train ; y.train.s <- u$y.train.s 
    K <- u$K ; K.test <- u$K.test ; y.test.s <- u$y.test.s
    y.train <- u$y.train ; y.test <- u$y.test
    n <- nrow(train.sim)
    n1 <- sum(nrow(train.sim), nrow(test.sim))
    
    
    # Choosing best lambda
    
    res.rmse1 <- NULL
    res.rmse2 <- NULL
    grid.l <- 10^seq(-3,2,length=10)
    
    for (i in 1:50) {
      
      boot.index <- sample(1:n, n, replace=T)
      boot.K.train <- K.train[boot.index, boot.index]
      boot.K.test <- K.train[boot.index, -boot.index]
      boot.y.train.s <- y.train.s[boot.index]
      boot.y.test.s <- y.train.s[-boot.index]
      boot.y.train <- y.train[boot.index]
      boot.y.test <- y.train[-boot.index]
      
      # synthetic
      boot.rmse1 <- c(rep(0,10))
      # original
      boot.rmse2 <- c(rep(0,10))
      
      for (j in 1:10) {
        
        h1.boot <- fit.kernel(boot.y.train.s, boot.y.train, boot.K.train, grid.l[j])
        boot.d.hat <- h1.boot$d.hat
        h2.boot <- pred.kernel(boot.y.test.s, boot.y.test, boot.K.test, boot.d.hat)
        boot.rmse1[j] <- h2.boot$rmse1
        boot.rmse2[j] <- h2.boot$rmse2
        
      }
      
      res.rmse1 <- rbind(res.rmse1, boot.rmse1)
      res.rmse2 <- rbind(res.rmse2, boot.rmse2)
      
    }
    
    #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
    best.lambda1 <- max(grid.l[colMeans(res.rmse1) == min(colMeans(res.rmse1))])
    best.lambda2 <- max(grid.l[colMeans(res.rmse2) == min(colMeans(res.rmse2))])
    
    # Bootstrapping
    
    boots.y1 <- NULL
    boots.y2 <- NULL
    
    for (r in 1:100) {
      
      boots.index <- sample(1:n, n, replace=T)
      boots.K.train <- K[boots.index, boots.index]
      boots.K.test <- K[boots.index, (n+1):n1]
      boots.y.train.s <- y.train.s[boots.index]
      boots.y.test.s <- y.test.s
      boots.y.train <- y.train[boots.index]
      boots.y.test <- y.test
      
      # synthetic
      h1_1.boots <- fit.kernel(boots.y.train.s, boots.y.train, boots.K.train, best.lambda1)
      boots.d.hat1 <- h1_1.boots$d.hat
      h2_1.boots <- pred.kernel(boots.y.test.s, boots.y.test, boots.K.test, boots.d.hat1)
      boots.y.hat1 <- h2_1.boots$y.hat
      boots.y1 <- cbind(boots.y1, boots.y.hat1)
      
      # original
      h1_2.boots <- fit.kernel(boots.y.train.s, boots.y.train, boots.K.train, best.lambda2)
      boots.d.hat2 <- h1_2.boots$d.hat
      h2_2.boots <- pred.kernel(boots.y.test.s, boots.y.test, boots.K.test, boots.d.hat2)
      boots.y.hat2 <- h2_2.boots$y.hat
      boots.y2 <- cbind(boots.y2, boots.y.hat2)
      
    }
    cat(s, ",")
    
    # Bagging estimate 
    
    y.hat1.bag <- rowMeans(boots.y1)
    y.hat1.bag <- ifelse(y.hat1.bag<0, 0, y.hat1.bag)
    y.hat2.bag <- rowMeans(boots.y2)
    y.hat2.bag <- ifelse(y.hat2.bag<0, 0, y.hat2.bag)
    
    # KRB[s] <- sqrt(sum((y.test - y.hat.bag)^2)/length(y.test))
    GKRB1[s] <- sqrt(sum((y.test.s - y.hat1.bag)^2)/length(y.test.s))
    GKRB2[s] <- sqrt(sum((y.test - y.hat2.bag)^2)/length(y.test))
    
  }
  
  GKRB1
  GKRB2
  boxplot(GKRB1)
  boxplot(GKRB2)
  dat6_1 <- data.frame(RMSE=GKRB1, method=rep("g.GKRB1",100), number=as.character(rep(n.train,100)))
  dat6_2 <- data.frame(RMSE=GKRB2, method=rep("o.GKRB2",100), number=as.character(rep(n.train,100)))
  dat6 <- rbind(dat6_1, dat6_2)
  
  
  
  
  
  
  
  ### 4. Kernel ridge regression using Random Forest (KRR)
  
  LKRR1 <- c(rep(0,100))
  LKRR2 <- c(rep(0,100))
  
  for (s in 1:100) {
    
    # Making simulation data
    
    set.seed(s)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10)
    
    train.sim.y.s <- train.sim[,1] ; test.sim.y.s <- test.sim[,1]
    train.sim.y <- train.sim[,ncol(train.sim)]
    test.sim.y <- test.sim[,ncol(test.sim)]
    train.sim.X <- train.sim[,c(-1,-ncol(train.sim))]  
    test.sim.X <- test.sim[,c(-1,-ncol(test.sim))]
    
    n <- nrow(train.sim)
    n1 <- sum(nrow(train.sim), nrow(test.sim))
    # p <- ncol(train.sim) - 1
    p <- ncol(train.sim) - 2
    
    # Choosing best lambda
    
    res.rmse1 <- NULL
    res.rmse2 <- NULL
    grid.l <- 10^seq(-3,2,length=10)
    
    for (i in 1:50) {
      
      rf.index <- sample(1:p, round(sqrt(p)), replace=F)
      boot.index <- sample(1:n, n, replace=T)
      train.sim1 <- cbind(train.sim.y.s[boot.index], train.sim.X[boot.index], 
                          train.sim.y[boot.index])
      test.sim1 <- cbind(train.sim.y.s[-boot.index], train.sim.X[-boot.index],
                         train.sim.y[-boot.index])
      
      u1 <- my.kernel.matrix1(train.sim1, test.sim1)
      boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
      boot.y.train.s <- u1$y.train.s ; boot.y.test.s <- u1$y.test.s
      boot.K <- u1$K
      boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
      
      # synthetic
      boot.rmse1 <- c(rep(0,10))
      # original
      boot.rmse2 <- c(rep(0,10))
      
      for (j in 1:10) {
        
        h1.boot <- fit.kernel(boot.y.train.s, boot.y.train, boot.K.train, grid.l[j])
        boot.d.hat <- h1.boot$d.hat
        h2.boot <- pred.kernel(boot.y.test.s, boot.y.test, boot.K.test, boot.d.hat)
        boot.rmse1[j] <- h2.boot$rmse1
        boot.rmse2[j] <- h2.boot$rmse2
        
      }
      
      res.rmse1 <- rbind(res.rmse1, boot.rmse1)
      res.rmse2 <- rbind(res.rmse2, boot.rmse2)
      
    }
    
    #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
    best.lambda1 <- max(grid.l[colMeans(res.rmse1)==min(colMeans(res.rmse1))])
    best.lambda2 <- max(grid.l[colMeans(res.rmse2)==min(colMeans(res.rmse2))])
    
    # Bootstrapping
    
    boots.y1 <- NULL
    boots.y2 <- NULL
    
    for (r in 1:100) {
      
      rfs.index <- sample(1:p, round(sqrt(p)), replace=F)
      boots.index <- sample(1:n, n, replace=T)
      train.sim2 <- cbind(train.sim.y.s[boots.index], train.sim.X[boots.index],
                          train.sim.y[boots.index])
      test.sim2 <- cbind(test.sim.y.s, test.sim.X, test.sim.y)
      
      u2 <- my.kernel.matrix1(train.sim2, test.sim2)
      boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
      boots.y.train.s <- u2$y.train.s ; boots.y.test.s <- u2$y.test.s
      boots.K <- u2$K
      boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
      
      # synthetic
      h1_1.boots <- fit.kernel(boots.y.train, boots.y.train.pred, boots.K.train, best.lambda1)
      boots.d.hat1 <- h1_1.boots$d.hat
      h2_1.boots <- pred.kernel(boots.y.test, boots.y.test.pred, boots.K.test, boots.d.hat1)
      boots.y.hat1 <- h2_1.boots$y.hat
      boots.y1 <- cbind(boots.y1, boots.y.hat1)
      
      # original
      h1_2.boots <- fit.kernel(boots.y.train, boots.y.train.pred, boots.K.train, best.lambda2)
      boots.d.hat1 <- h1_1.boots$d.hat
      h2_2.boots <- pred.kernel(boots.y.test, boots.y.test.pred, boots.K.test, boots.d.hat2)
      boots.y.hat2 <- h2_2.boots$y.hat
      boots.y2 <- cbind(boots.y2, boots.y.hat2)
      
    }
    cat(s, ",")
    
    # Random Forest estimate 
    
    y.hat1.rf <- rowMeans(boots.y1)
    y.hat1.rf <- ifelse(y.hat1.rf<0, 0, y.hat1.rf)
    y.hat2.rf <- rowMeans(boots.y2)
    y.hat2.rf <- ifelse(y.hat2.rf<0, 0, y.hat2.rf)
    
    # KRR[s] <- sqrt(sum((test.sim.y - y.hat.rf)^2)/length(test.sim.y))
    LKRR1[s] <- sqrt(sum((test.sim.y.s - y.hat1.rf)^2)/length(test.sim.y.s))
    LKRR2[s] <- sqrt(sum((test.sim.y - y.hat2.rf)^2)/length(test.sim.y))
    
  }
  
  LKRR1
  LKRR2
  boxplot(LKRR1)
  boxplot(LKRR2)
  dat7_1 <- data.frame(RMSE=LKRR1, method=rep("d.LKRR1",100), number=as.character(rep(n.train,100)))
  dat7_2 <- data.frame(RMSE=LKRR2, method=rep("l.LKRR2",100), number=as.character(rep(n.train,100)))
  dat7 <- rbind(dat7_1, dat7_2)
  
  GKRR1 <- c(rep(0,100))
  GKRR2 <- c(rep(0,100))
  
  for (s in 1:100) {
    
    # Making simulation data
    
    set.seed(s)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10)
    
    train.sim.y.s <- train.sim[,1] ; test.sim.y.s <- test.sim[,1]
    train.sim.y <- train.sim[,ncol(train.sim)]
    test.sim.y <- test.sim[,ncol(test.sim)]
    train.sim.X <- train.sim[,c(-1,-ncol(train.sim))]  
    test.sim.X <- test.sim[,c(-1,-ncol(test.sim))]
    
    n <- nrow(train.sim)
    n1 <- sum(nrow(train.sim), nrow(test.sim))
    # p <- ncol(train.sim) - 1
    p <- ncol(train.sim) - 2
    
    # Choosing best lambda
    
    res.rmse1 <- NULL
    res.rmse2 <- NULL
    grid.l <- 10^seq(-3,2,length=10)
    
    for (i in 1:50) {
      
      rf.index <- sample(1:p, round(sqrt(p)), replace=F)
      boot.index <- sample(1:n, n, replace=T)
      train.sim1 <- cbind(train.sim.y.s[boot.index], train.sim.X[boot.index], 
                          train.sim.y[boot.index])
      test.sim1 <- cbind(train.sim.y.s[-boot.index], train.sim.X[-boot.index],
                         train.sim.y[-boot.index])
      
      u1 <- my.kernel.matrix2(train.sim1, test.sim1)
      boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
      boot.y.train.s <- u1$y.train.s ; boot.y.test.s <- u1$y.test.s
      boot.K <- u1$K
      boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
      
      # synthetic
      boot.rmse1 <- c(rep(0,10))
      # original
      boot.rmse2 <- c(rep(0,10))
      
      for (j in 1:10) {
        
        h1.boot <- fit.kernel(boot.y.train.s, boot.y.train, boot.K.train, grid.l[j])
        boot.d.hat <- h1.boot$d.hat
        h2.boot <- pred.kernel(boot.y.test.s, boot.y.test, boot.K.test, boot.d.hat)
        boot.rmse1[j] <- h2.boot$rmse1
        boot.rmse2[j] <- h2.boot$rmse2
        
      }
      
      res.rmse1 <- rbind(res.rmse1, boot.rmse1)
      res.rmse2 <- rbind(res.rmse2, boot.rmse2)
      
    }
    
    #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
    best.lambda1 <- max(grid.l[colMeans(res.rmse1)==min(colMeans(res.rmse1))])
    best.lambda2 <- max(grid.l[colMeans(res.rmse2)==min(colMeans(res.rmse2))])
    
    # Bootstrapping
    
    boots.y1 <- NULL
    boots.y2 <- NULL
    
    for (r in 1:100) {
      
      rfs.index <- sample(1:p, round(sqrt(p)), replace=F)
      boots.index <- sample(1:n, n, replace=T)
      train.sim2 <- cbind(train.sim.y.s[boots.index], train.sim.X[boots.index],
                          train.sim.y[boots.index])
      test.sim2 <- cbind(test.sim.y.s, test.sim.X, test.sim.y)
      
      u2 <- my.kernel.matrix2(train.sim2, test.sim2)
      boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
      boots.y.train.s <- u2$y.train.s ; boots.y.test.s <- u2$y.test.s
      boots.K <- u2$K
      boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
      
      # synthetic
      h1_1.boots <- fit.kernel(boots.y.train, boots.y.train.pred, boots.K.train, best.lambda1)
      boots.d.hat1 <- h1_1.boots$d.hat
      h2_1.boots <- pred.kernel(boots.y.test, boots.y.test.pred, boots.K.test, boots.d.hat1)
      boots.y.hat1 <- h2_1.boots$y.hat
      boots.y1 <- cbind(boots.y1, boots.y.hat1)
      
      # original
      h1_2.boots <- fit.kernel(boots.y.train, boots.y.train.pred, boots.K.train, best.lambda2)
      boots.d.hat1 <- h1_1.boots$d.hat
      h2_2.boots <- pred.kernel(boots.y.test, boots.y.test.pred, boots.K.test, boots.d.hat2)
      boots.y.hat2 <- h2_2.boots$y.hat
      boots.y2 <- cbind(boots.y2, boots.y.hat2)
      
    }
    cat(s, ",")
    
    # Random Forest estimate 
    
    y.hat1.rf <- rowMeans(boots.y1)
    y.hat1.rf <- ifelse(y.hat1.rf<0, 0, y.hat1.rf)
    y.hat2.rf <- rowMeans(boots.y2)
    y.hat2.rf <- ifelse(y.hat2.rf<0, 0, y.hat2.rf)
    
    # KRR[s] <- sqrt(sum((test.sim.y - y.hat.rf)^2)/length(test.sim.y))
    GKRR1[s] <- sqrt(sum((test.sim.y.s - y.hat1.rf)^2)/length(test.sim.y.s))
    GKRR2[s] <- sqrt(sum((test.sim.y - y.hat2.rf)^2)/length(test.sim.y))
    
  }
  
  GKRR1
  GKRR2
  boxplot(GKRR1)
  boxplot(GKRR2)
  dat8_1 <- data.frame(RMSE=GKRR1, method=rep("h.GKRR1",100), number=as.character(rep(n.train,100)))
  dat8_2 <- data.frame(RMSE=GKRR2, method=rep("p.GKRR2",100), number=as.character(rep(n.train,100)))
  dat8 <- rbind(dat8_1, dat8_2)
  
  
  
  
  # Final result
  
  dat.res <- rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8)
  
  return(dat.res)
  
}



### Function for many independent variables



fit.ftn <- function(dat.sim) {
  
  
  
  
  ### 1. Kernel ridge regression (KR)
  
  LKR1 <- c(rep(0,100))
  LKR2 <- c(rep(0,100))
  
  for (i in 1:100) {
    
    # Making simulation data
    
    set.seed(i)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10)  
    
    
    # 5-fold crossvalidation
    
    u <- my.kernel.matrix1(train.sim, test.sim)
    K.train <- u$K.train ; K.test <- u$K.test  
    y.train.s <- u$y.train.s ; y.test.s <- u$y.test.s
    y.train <- u$y.train ; y.test <- u$y.test
    
    grid.l <- 10^seq(-3,2,length=10)
    
    h <- cv.kernel(y.train.s, y.train, K.train, 5, grid.l) 
    
    # Fitting 
    
    mean.rmse1 <- rowMeans(h$cv.rmse1)
    idx1 <- which.min(mean.rmse1)
    best.lam1 <- max(h$lambda[ mean.rmse1 == mean.rmse1[idx1] ])
    
    h1_1 <- fit.kernel(y.train.s, y.train, K.train, best.lam1)
    
    mean.rmse2 <- rowMeans(h$cv.rmse2)
    idx2 <- which.min(mean.rmse2)
    best.lam2 <- max(h$lambda[ mean.rmse2 == mean.rmse2[idx2] ])
    
    h1_2 <- fit.kernel(y.train.s, y.train, K.train, best.lam2)
    
    # Calculate test RMSE
    
    sim.d.hat1 <- h1_1$d.hat
    h2_1 <- pred.kernel(y.test.s, y.test, K.test, sim.d.hat1)
    
    sim.d.hat2 <- h1_2$d.hat
    h2_2 <- pred.kernel(y.test.s, y.test, K.test, sim.d.hat2)
    
    LKR1[i] <- h2_1$rmse1
    LKR2[i] <- h2_2$rmse2
    
  }
  
  LKR1
  LKR2
  boxplot(LKR1)
  boxplot(LKR2)
  dat1_1 <- data.frame(RMSE=LKR1, method=rep("a.LKR1",100), number=as.character(rep(n.train,100)))
  dat1_2 <- data.frame(RMSE=LKR2, method=rep("i.LKR2",100), number=as.character(rep(n.train,100)))
  dat1 <- rbind(dat1_1, dat1_2)
  
  GKR1 <- c(rep(0,100))
  GKR2 <- c(rep(0,100))
  
  for (i in 1:100) {
    
    # Making simulation data
    
    set.seed(i)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10) 
    
    
    # 5-fold crossvalidation
    
    u <- my.kernel.matrix2(train.sim, test.sim)
    K.train <- u$K.train ; K.test <- u$K.test  
    y.train.s <- u$y.train.s ; y.test.s <- u$y.test.s
    y.train <- u$y.train ; y.test <- u$y.test
    
    grid.l <- 10^seq(-3,2,length=10)
    
    h <- cv.kernel(y.train.s, y.train, K.train, 5, grid.l) 
    
    # Fitting 
    
    mean.rmse1 <- rowMeans(h$cv.rmse1)
    idx1 <- which.min(mean.rmse1)
    best.lam1 <- max(h$lambda[ mean.rmse1 == mean.rmse1[idx1] ])
    
    h1_1 <- fit.kernel(y.train.s, y.train, K.train, best.lam1)
    
    mean.rmse2 <- rowMeans(h$cv.rmse2)
    idx2 <- which.min(mean.rmse2)
    best.lam2 <- max(h$lambda[ mean.rmse2 == mean.rmse2[idx2] ])
    
    h1_2 <- fit.kernel(y.train.s, y.train, K.train, best.lam2)
    
    # Calculate test RMSE
    
    sim.d.hat1 <- h1_1$d.hat
    h2_1 <- pred.kernel(y.test.s, y.test, K.test, sim.d.hat1)
    
    sim.d.hat2 <- h1_2$d.hat
    h2_2 <- pred.kernel(y.test.s, y.test, K.test, sim.d.hat2)
    
    GKR1[i] <- h2_1$rmse1
    GKR2[i] <- h2_2$rmse2
    
  }
  
  GKR1
  GKR2
  boxplot(GKR1)
  boxplot(GKR2)
  dat2_1 <- data.frame(RMSE=GKR1, method=rep("e.GKR1",100), number=as.character(rep(n.train,100)))
  dat2_2 <- data.frame(RMSE=GKR2, method=rep("m.GKR2",100), number=as.character(rep(n.train,100)))
  dat2 <- rbind(dat2_1, dat2_2)
  
  
  
  
  ### 2. Kernel ridge regression using Sub-sampling (KRS)
  
  LKRS1 <- c(rep(0,100))
  LKRS2 <- c(rep(0,100))
  
  for (i in 1:100) {
    
    # Making simulation data
    
    set.seed(i)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10)
    
    u <- my.kernel.matrix1(train.sim, test.sim)
    K.train <- u$K.train ; y.train.s <- u$y.train.s 
    K <- u$K ; K.test <- u$K.test ; y.test.s <- u$y.test.s
    y.train <- u$y.train ; y.test <- u$y.test
    
    # Choosing optimal lambda
    
    
    n <- nrow(train.sim) ; n1 <- sum(nrow(train.sim), nrow(test.sim))
    res.rmse1.c <- NULL ; res.lam1.c <- NULL ; res.index1.c <- NULL
    res.rmse2.c <- NULL ; res.lam2.c <- NULL ; res.index2.c <- NULL
    
    for (j in 1:50) {
      
      n <- nrow(train.sim)
      index1.Kc <- sample(1:n, round(0.7*n), replace=F)
      index2.Kc <- index1.Kc
      Kc.train <- K.train[index1.Kc,index1.Kc] ; Kc.test <- K.train[index1.Kc,-index1.Kc]
      yc.train.s <- y.train.s[index1.Kc] ; yc.test.s <- y.train.s[-index1.Kc]
      yc.train <- y.train[index1.Kc] ; yc.test <- y.train[-index1.Kc]
      
      grid.l <- 10^seq(-3,2,length=10)
      
      hc <- cv.kernel(yc.train.s, yc.train, Kc.train, 5, grid.l) # 5-fold crossvalidation 
      
      # synthetic
      mean.rmse1.c <- rowMeans(hc$cv.rmse1)
      idx1.c <- which.min(mean.rmse1.c)
      best.lam1.c <- max(hc$lambda[ mean.rmse1.c == mean.rmse1.c[idx1.c] ])
      
      h1_1.c <- fit.kernel(yc.train.s, yc.train, Kc.train, best.lam1.c) # Fitting
      
      sim.d.hat1.c <- h1_1.c$d.hat
      
      h2_1.c <- pred.kernel(yc.test.s, yc.test, Kc.test, sim.d.hat1.c) # Testing
      rmse1.c <- h2_1.c$rmse1
      
      res.rmse1.c <- c(res.rmse1.c, rmse1.c)
      res.lam1.c <- c(res.lam1.c, best.lam1.c)
      res.index1.c <- cbind(res.index1.c, index1.Kc)
      
      # original
      mean.rmse2.c <- rowMeans(hc$cv.rmse2)
      idx2.c <- which.min(mean.rmse2.c)
      best.lam2.c <- max(hc$lambda[ mean.rmse2.c == mean.rmse2.c[idx2.c] ])
      
      h1_2.c <- fit.kernel(yc.train.s, yc.train, Kc.train, best.lam2.c) # Fitting
      
      sim.d.hat2.c <- h1_2.c$d.hat
      
      h2_2.c <- pred.kernel(yc.test.s, yc.test, Kc.test, sim.d.hat2.c) # Testing
      rmse2.c <- h2_2.c$rmse1
      
      res.rmse2.c <- c(res.rmse2.c, rmse2.c)
      res.lam2.c <- c(res.lam2.c, best.lam2.c)
      res.index2.c <- cbind(res.index2.c, index2.Kc)
      
      
    }
    
    # Fitting  
    
    # synthetic
    res.lambda1 <- res.lam1.c[which.min(res.rmse1.c)]
    res.index1 <- res.index1.c[, which.min(res.rmse1.c)]
    res.K.train1 <- K[res.index1, res.index1]
    res.y.train.s1 <- y.train[res.index1]
    res.y.train1 <- y.train[res.index1]
    K.test1 <- K[res.index1, (n+1):n1]
    
    h1_1 <- fit.kernel(res.y.train.s1, res.y.train1, res.K.train1, res.lambda1)
    
    res.d.hat1 <- h1_1$d.hat
    
    # Calculate test mean square error
    
    h2_1 <- pred.kernel(y.test.s, y.test, K.test1, res.d.hat1)
    LKRS1[i] <- h2_1$rmse1
    
    # original
    res.lambda2 <- res.lam2.c[which.min(res.rmse2.c)]
    res.index2 <- res.index2.c[, which.min(res.rmse2.c)]
    res.K.train2 <- K[res.index2, res.index2]
    res.y.train.s2 <- y.train[res.index2]
    res.y.train2 <- y.train[res.index2]
    K.test2 <- K[res.index2, (n+1):n1]
    
    h1_2 <- fit.kernel(res.y.train.s2, res.y.train2, res.K.train2, res.lambda2)
    
    res.d.hat2 <- h1_2$d.hat
    
    # Calculate test mean square error
    
    h2_2 <- pred.kernel(y.test.s, y.test, K.test2, res.d.hat2)
    LKRS2[i] <- h2_2$rmse2
    
  }
  
  LKRS1
  LKRS2
  boxplot(LKRS1)
  boxplot(LKRS2)
  dat3_1 <- data.frame(RMSE=LKRS1, method=rep("b.LKRS1",100), number=as.character(rep(n.train,100)))
  dat3_2 <- data.frame(RMSE=LKRS2, method=rep("j.LKRS2",100), number=as.character(rep(n.train,100)))
  dat3 <- rbind(dat3_1, dat3_2)
  
  GKRS1 <- c(rep(0,100))
  GKRS2 <- c(rep(0,100))
  
  for (i in 1:100) {
    
    # Making simulation data
    
    set.seed(i)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10)
    
    u <- my.kernel.matrix2(train.sim, test.sim)
    K.train <- u$K.train ; y.train.s <- u$y.train.s 
    K <- u$K ; K.test <- u$K.test ; y.test.s <- u$y.test.s
    y.train <- u$y.train ; y.test <- u$y.test
    
    # Choosing optimal lambda
    
    
    n <- nrow(train.sim) ; n1 <- sum(nrow(train.sim), nrow(test.sim))
    res.rmse1.c <- NULL ; res.lam1.c <- NULL ; res.index1.c <- NULL
    res.rmse2.c <- NULL ; res.lam2.c <- NULL ; res.index2.c <- NULL
    
    for (j in 1:50) {
      
      n <- nrow(train.sim)
      index1.Kc <- sample(1:n, round(0.7*n), replace=F)
      index2.Kc <- index1.Kc
      Kc.train <- K.train[index1.Kc,index1.Kc] ; Kc.test <- K.train[index1.Kc,-index1.Kc]
      yc.train.s <- y.train.s[index1.Kc] ; yc.test.s <- y.train.s[-index1.Kc]
      yc.train <- y.train[index1.Kc] ; yc.test <- y.train[-index1.Kc]
      
      grid.l <- 10^seq(-3,2,length=10)
      
      hc <- cv.kernel(yc.train.s, yc.train, Kc.train, 5, grid.l) # 5-fold crossvalidation 
      
      # synthetic
      mean.rmse1.c <- rowMeans(hc$cv.rmse1)
      idx1.c <- which.min(mean.rmse1.c)
      best.lam1.c <- max(hc$lambda[ mean.rmse1.c == mean.rmse1.c[idx1.c] ])
      
      h1_1.c <- fit.kernel(yc.train.s, yc.train, Kc.train, best.lam1.c) # Fitting
      
      sim.d.hat1.c <- h1_1.c$d.hat
      
      h2_1.c <- pred.kernel(yc.test.s, yc.test, Kc.test, sim.d.hat1.c) # Testing
      rmse1.c <- h2_1.c$rmse1
      
      res.rmse1.c <- c(res.rmse1.c, rmse1.c)
      res.lam1.c <- c(res.lam1.c, best.lam1.c)
      res.index1.c <- cbind(res.index1.c, index1.Kc)
      
      # original
      mean.rmse2.c <- rowMeans(hc$cv.rmse2)
      idx2.c <- which.min(mean.rmse2.c)
      best.lam2.c <- max(hc$lambda[ mean.rmse2.c == mean.rmse2.c[idx2.c] ])
      
      h1_2.c <- fit.kernel(yc.train.s, yc.train, Kc.train, best.lam2.c) # Fitting
      
      sim.d.hat2.c <- h1_2.c$d.hat
      
      h2_2.c <- pred.kernel(yc.test.s, yc.test, Kc.test, sim.d.hat2.c) # Testing
      rmse2.c <- h2_2.c$rmse1
      
      res.rmse2.c <- c(res.rmse2.c, rmse2.c)
      res.lam2.c <- c(res.lam2.c, best.lam2.c)
      res.index2.c <- cbind(res.index2.c, index2.Kc)
      
      
    }
    
    # Fitting  
    
    # synthetic
    res.lambda1 <- res.lam1.c[which.min(res.rmse1.c)]
    res.index1 <- res.index1.c[, which.min(res.rmse1.c)]
    res.K.train1 <- K[res.index1, res.index1]
    res.y.train.s1 <- y.train[res.index1]
    res.y.train1 <- y.train[res.index1]
    K.test1 <- K[res.index1, (n+1):n1]
    
    h1_1 <- fit.kernel(res.y.train.s1, res.y.train1, res.K.train1, res.lambda1)
    
    res.d.hat1 <- h1_1$d.hat
    
    # Calculate test mean square error
    
    h2_1 <- pred.kernel(y.test.s, y.test, K.test1, res.d.hat1)
    GKRS1[i] <- h2_1$rmse1
    
    # original
    res.lambda2 <- res.lam2.c[which.min(res.rmse2.c)]
    res.index2 <- res.index2.c[, which.min(res.rmse2.c)]
    res.K.train2 <- K[res.index2, res.index2]
    res.y.train.s2 <- y.train[res.index2]
    res.y.train2 <- y.train[res.index2]
    K.test2 <- K[res.index2, (n+1):n1]
    
    h1_2 <- fit.kernel(res.y.train.s2, res.y.train2, res.K.train2, res.lambda2)
    
    res.d.hat2 <- h1_2$d.hat
    
    # Calculate test mean square error
    
    h2_2 <- pred.kernel(y.test.s, y.test, K.test2, res.d.hat2)
    GKRS2[i] <- h2_2$rmse2
    
  }
  
  GKRS1
  GKRS2
  boxplot(GKRS1)
  boxplot(GKRS2)
  dat4_1 <- data.frame(RMSE=GKRS1, method=rep("f.GKRS1",100), number=as.character(rep(n.train,100)))
  dat4_2 <- data.frame(RMSE=GKRS2, method=rep("n.GKRS2",100), number=as.character(rep(n.train,100)))
  dat4 <- rbind(dat4_1, dat4_2)
  
  
  
  
  
  ### 3. Kernel ridge regression using Bagging (KRB)
  
  LKRB1 <- c(rep(0,100))
  LKRB2 <- c(rep(0,100))
  
  for (s in 1:100) {
    
    # Making simulation data
    
    set.seed(s)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    test.sim.y <- test.sim[,1] 
    n.train <- round(nrow(dat.sim)*7/10)
    
    u <- my.kernel.matrix1(train.sim, test.sim)
    K.train <- u$K.train ; y.train.s <- u$y.train.s 
    K <- u$K ; K.test <- u$K.test ; y.test.s <- u$y.test.s
    y.train <- u$y.train ; y.test <- u$y.test
    n <- nrow(train.sim)
    n1 <- sum(nrow(train.sim), nrow(test.sim))
    
    
    # Choosing best lambda
    
    res.rmse1 <- NULL
    res.rmse2 <- NULL
    grid.l <- 10^seq(-3,2,length=10)
    
    for (i in 1:50) {
      
      boot.index <- sample(1:n, n, replace=T)
      boot.K.train <- K.train[boot.index, boot.index]
      boot.K.test <- K.train[boot.index, -boot.index]
      boot.y.train.s <- y.train.s[boot.index]
      boot.y.test.s <- y.train.s[-boot.index]
      boot.y.train <- y.train[boot.index]
      boot.y.test <- y.train[-boot.index]
      
      # synthetic
      boot.rmse1 <- c(rep(0,10))
      # original
      boot.rmse2 <- c(rep(0,10))
      
      for (j in 1:10) {
        
        h1.boot <- fit.kernel(boot.y.train.s, boot.y.train, boot.K.train, grid.l[j])
        boot.d.hat <- h1.boot$d.hat
        h2.boot <- pred.kernel(boot.y.test.s, boot.y.test, boot.K.test, boot.d.hat)
        boot.rmse1[j] <- h2.boot$rmse1
        boot.rmse2[j] <- h2.boot$rmse2
        
      }
      
      res.rmse1 <- rbind(res.rmse1, boot.rmse1)
      res.rmse2 <- rbind(res.rmse2, boot.rmse2)
      
    }
    
    #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
    best.lambda1 <- max(grid.l[colMeans(res.rmse1) == min(colMeans(res.rmse1))])
    best.lambda2 <- max(grid.l[colMeans(res.rmse2) == min(colMeans(res.rmse2))])
    
    # Bootstrapping
    
    boots.y1 <- NULL
    boots.y2 <- NULL
    
    for (r in 1:100) {
      
      boots.index <- sample(1:n, n, replace=T)
      boots.K.train <- K[boots.index, boots.index]
      boots.K.test <- K[boots.index, (n+1):n1]
      boots.y.train.s <- y.train.s[boots.index]
      boots.y.test.s <- y.test.s
      boots.y.train <- y.train[boots.index]
      boots.y.test <- y.test
      
      # synthetic
      h1_1.boots <- fit.kernel(boots.y.train.s, boots.y.train, boots.K.train, best.lambda1)
      boots.d.hat1 <- h1_1.boots$d.hat
      h2_1.boots <- pred.kernel(boots.y.test.s, boots.y.test, boots.K.test, boots.d.hat1)
      boots.y.hat1 <- h2_1.boots$y.hat
      boots.y1 <- cbind(boots.y1, boots.y.hat1)
      
      # original
      h1_2.boots <- fit.kernel(boots.y.train.s, boots.y.train, boots.K.train, best.lambda2)
      boots.d.hat2 <- h1_2.boots$d.hat
      h2_2.boots <- pred.kernel(boots.y.test.s, boots.y.test, boots.K.test, boots.d.hat2)
      boots.y.hat2 <- h2_2.boots$y.hat
      boots.y2 <- cbind(boots.y2, boots.y.hat2)
      
    }
    cat(s, ",")
    
    # Bagging estimate 
    
    y.hat1.bag <- rowMeans(boots.y1)
    y.hat1.bag <- ifelse(y.hat1.bag<0, 0, y.hat1.bag)
    y.hat2.bag <- rowMeans(boots.y2)
    y.hat2.bag <- ifelse(y.hat2.bag<0, 0, y.hat2.bag)
    
    # KRB[s] <- sqrt(sum((y.test - y.hat.bag)^2)/length(y.test))
    LKRB1[s] <- sqrt(sum((y.test.s - y.hat1.bag)^2)/length(y.test.s))
    LKRB2[s] <- sqrt(sum((y.test - y.hat2.bag)^2)/length(y.test))
    
  }
  
  LKRB1
  LKRB2
  boxplot(LKRB1)
  boxplot(LKRB2)
  dat5_1 <- data.frame(RMSE=LKRB1, method=rep("c.LKRB1",100), number=as.character(rep(n.train,100)))
  dat5_2 <- data.frame(RMSE=LKRB2, method=rep("k.LKRB2",100), number=as.character(rep(n.train,100)))
  dat5 <- rbind(dat5_1, dat5_2)
  
  GKRB1 <- c(rep(0,100))
  GKRB2 <- c(rep(0,100))
  
  for (s in 1:100) {
    
    # Making simulation data
    
    set.seed(s)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    test.sim.y <- test.sim[,1] 
    n.train <- round(nrow(dat.sim)*7/10)
    
    u <- my.kernel.matrix2(train.sim, test.sim)
    K.train <- u$K.train ; y.train.s <- u$y.train.s 
    K <- u$K ; K.test <- u$K.test ; y.test.s <- u$y.test.s
    y.train <- u$y.train ; y.test <- u$y.test
    n <- nrow(train.sim)
    n1 <- sum(nrow(train.sim), nrow(test.sim))
    
    
    # Choosing best lambda
    
    res.rmse1 <- NULL
    res.rmse2 <- NULL
    grid.l <- 10^seq(-3,2,length=10)
    
    for (i in 1:50) {
      
      boot.index <- sample(1:n, n, replace=T)
      boot.K.train <- K.train[boot.index, boot.index]
      boot.K.test <- K.train[boot.index, -boot.index]
      boot.y.train.s <- y.train.s[boot.index]
      boot.y.test.s <- y.train.s[-boot.index]
      boot.y.train <- y.train[boot.index]
      boot.y.test <- y.train[-boot.index]
      
      # synthetic
      boot.rmse1 <- c(rep(0,10))
      # original
      boot.rmse2 <- c(rep(0,10))
      
      for (j in 1:10) {
        
        h1.boot <- fit.kernel(boot.y.train.s, boot.y.train, boot.K.train, grid.l[j])
        boot.d.hat <- h1.boot$d.hat
        h2.boot <- pred.kernel(boot.y.test.s, boot.y.test, boot.K.test, boot.d.hat)
        boot.rmse1[j] <- h2.boot$rmse1
        boot.rmse2[j] <- h2.boot$rmse2
        
      }
      
      res.rmse1 <- rbind(res.rmse1, boot.rmse1)
      res.rmse2 <- rbind(res.rmse2, boot.rmse2)
      
    }
    
    #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
    best.lambda1 <- max(grid.l[colMeans(res.rmse1) == min(colMeans(res.rmse1))])
    best.lambda2 <- max(grid.l[colMeans(res.rmse2) == min(colMeans(res.rmse2))])
    
    # Bootstrapping
    
    boots.y1 <- NULL
    boots.y2 <- NULL
    
    for (r in 1:100) {
      
      boots.index <- sample(1:n, n, replace=T)
      boots.K.train <- K[boots.index, boots.index]
      boots.K.test <- K[boots.index, (n+1):n1]
      boots.y.train.s <- y.train.s[boots.index]
      boots.y.test.s <- y.test.s
      boots.y.train <- y.train[boots.index]
      boots.y.test <- y.test
      
      # synthetic
      h1_1.boots <- fit.kernel(boots.y.train.s, boots.y.train, boots.K.train, best.lambda1)
      boots.d.hat1 <- h1_1.boots$d.hat
      h2_1.boots <- pred.kernel(boots.y.test.s, boots.y.test, boots.K.test, boots.d.hat1)
      boots.y.hat1 <- h2_1.boots$y.hat
      boots.y1 <- cbind(boots.y1, boots.y.hat1)
      
      # original
      h1_2.boots <- fit.kernel(boots.y.train.s, boots.y.train, boots.K.train, best.lambda2)
      boots.d.hat2 <- h1_2.boots$d.hat
      h2_2.boots <- pred.kernel(boots.y.test.s, boots.y.test, boots.K.test, boots.d.hat2)
      boots.y.hat2 <- h2_2.boots$y.hat
      boots.y2 <- cbind(boots.y2, boots.y.hat2)
      
    }
    cat(s, ",")
    
    # Bagging estimate 
    
    y.hat1.bag <- rowMeans(boots.y1)
    y.hat1.bag <- ifelse(y.hat1.bag<0, 0, y.hat1.bag)
    y.hat2.bag <- rowMeans(boots.y2)
    y.hat2.bag <- ifelse(y.hat2.bag<0, 0, y.hat2.bag)
    
    # KRB[s] <- sqrt(sum((y.test - y.hat.bag)^2)/length(y.test))
    GKRB1[s] <- sqrt(sum((y.test.s - y.hat1.bag)^2)/length(y.test.s))
    GKRB2[s] <- sqrt(sum((y.test - y.hat2.bag)^2)/length(y.test))
    
  }
  
  GKRB1
  GKRB2
  boxplot(GKRB1)
  boxplot(GKRB2)
  dat6_1 <- data.frame(RMSE=GKRB1, method=rep("g.GKRB1",100), number=as.character(rep(n.train,100)))
  dat6_2 <- data.frame(RMSE=GKRB2, method=rep("o.GKRB2",100), number=as.character(rep(n.train,100)))
  dat6 <- rbind(dat6_1, dat6_2)
  
  
  
  
  ### 4. Kernel ridge regression using Random Forest (KRR)
  
  LKRR1 <- c(rep(0,100))
  LKRR2 <- c(rep(0,100))
  
  for (s in 1:100) {
    
    # Making simulation data
    
    set.seed(s)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10)
    
    train.sim.y.s <- train.sim[,1] ; test.sim.y.s <- test.sim[,1]
    train.sim.y <- train.sim[,ncol(train.sim)]
    test.sim.y <- test.sim[,ncol(test.sim)]
    train.sim.X <- train.sim[,c(-1,-ncol(train.sim))]  
    test.sim.X <- test.sim[,c(-1,-ncol(test.sim))]
    
    n <- nrow(train.sim)
    n1 <- sum(nrow(train.sim), nrow(test.sim))
    # p <- ncol(train.sim) - 1
    p <- ncol(train.sim) - 2
    
    # Choosing best lambda
    
    res.rmse1 <- NULL
    res.rmse2 <- NULL
    grid.l <- 10^seq(-3,2,length=10)
    
    for (i in 1:50) {
      
      rf.index <- sample(1:p, round(sqrt(p)), replace=F)
      boot.index <- sample(1:n, n, replace=T)
      train.sim1 <- cbind(train.sim.y.s[boot.index], train.sim.X[boot.index, rf.index], 
                          train.sim.y[boot.index])
      test.sim1 <- cbind(train.sim.y.s[-boot.index], train.sim.X[-boot.index, rf.index],
                         train.sim.y[-boot.index])
      
      u1 <- my.kernel.matrix1(train.sim1, test.sim1)
      boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
      boot.y.train.s <- u1$y.train.s ; boot.y.test.s <- u1$y.test.s
      boot.K <- u1$K
      boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
      
      boot.rmse <- c(rep(0,10))
      
      for (j in 1:10) {
        
        h1.boot <- fit.kernel(boot.y.train.s, boot.y.train, boot.K.train, grid.l[j])
        boot.d.hat <- h1.boot$d.hat
        h2.boot <- pred.kernel(boot.y.test.s, boot.y.test, boot.K.test, boot.d.hat)
        boot.rmse1[j] <- h2.boot$rmse1
        boot.rmse2[j] <- h2.boot$rmse2
        
      }
      
      res.rmse1 <- rbind(res.rmse1, boot.rmse1)
      res.rmse2 <- rbind(res.rmse2, boot.rmse2)
      
    }
    
    #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
    best.lambda1 <- max(grid.l[colMeans(res.rmse1)==min(colMeans(res.rmse1))])
    best.lambda2 <- max(grid.l[colMeans(res.rmse2)==min(colMeans(res.rmse2))])
    
    # Bootstrapping
    
    boots.y1 <- NULL
    boots.y2 <- NULL
    
    for (r in 1:100) {
      
      rfs.index <- sample(1:p, round(sqrt(p)), replace=F)
      boots.index <- sample(1:n, n, replace=T)
      train.sim2 <- cbind(train.sim.y.s[boots.index], train.sim.X[boots.index, rfs.index],
                          train.sim.y[boots.index])
      test.sim2 <- cbind(test.sim.y.s, test.sim.X[, rfs.index], test.sim.y)
      
      u2 <- my.kernel.matrix1(train.sim2, test.sim2)
      boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
      boots.y.train.s <- u2$y.train.s ; boots.y.test.s <- u2$y.test.s
      boots.K <- u2$K
      boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
      
      # synthetic
      h1_1.boots <- fit.kernel(boots.y.train.s, boots.y.train, boots.K.train, best.lambda1)
      boots.d.hat1 <- h1_1.boots$d.hat
      h2_1.boots <- pred.kernel(boots.y.test.s, boots.y.test, boots.K.test, boots.d.hat1)
      boots.y.hat1 <- h2_1.boots$y.hat
      boots.y1 <- cbind(boots.y1, boots.y.hat1)
      
      # original
      h1_2.boots <- fit.kernel(boots.y.train.s, boots.y.train, boots.K.train, best.lambda2)
      boots.d.hat2 <- h1_2.boots$d.hat
      h2_2.boots <- pred.kernel(boots.y.test.s, boots.y.test, boots.K.test, boots.d.hat2)
      boots.y.hat2 <- h2_2.boots$y.hat
      boots.y2 <- cbind(boots.y2, boots.y.hat2)
      
    }
    cat(s, ",")
    
    # Random Forest estimate 
    
    y.hat1.rf <- rowMeans(boots.y1)
    y.hat1.rf <- ifelse(y.hat1.rf<0, 0, y.hat1.rf)
    y.hat2.rf <- rowMeans(boots.y2)
    y.hat2.rf <- ifelse(y.hat2.rf<0, 0, y.hat2.rf)
    
    # KRR[s] <- sqrt(sum((test.sim.y - y.hat.rf)^2)/length(test.sim.y))
    LKRR1[s] <- sqrt(sum((test.sim.y.s - y.hat1.rf)^2)/length(test.sim.y.s))
    LKRR2[s] <- sqrt(sum((test.sim.y - y.hat2.rf)^2)/length(test.sim.y))
    
  }
  
  LKRR1
  LKRR2
  boxplot(LKRR1)
  boxplot(LKRR2)
  dat7_1 <- data.frame(RMSE=LKRR1, method=rep("d.LKRR1",100), number=as.character(rep(n.train,100)))
  dat7_2 <- data.frame(RMSE=LKRR2, method=rep("l.LKRR2",100), number=as.character(rep(n.train,100)))
  dat7 <- rbind(dat7_1, dat7_2)
  
  GKRR1 <- c(rep(0,100))
  GKRR2 <- c(rep(0,100))
  
  for (s in 1:100) {
    
    # Making simulation data
    
    set.seed(s)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10)
    
    train.sim.y.s <- train.sim[,1] ; test.sim.y.s <- test.sim[,1]
    train.sim.y <- train.sim[,ncol(train.sim)]
    test.sim.y <- test.sim[,ncol(test.sim)]
    train.sim.X <- train.sim[,c(-1,-ncol(train.sim))]  
    test.sim.X <- test.sim[,c(-1,-ncol(test.sim))]
    
    n <- nrow(train.sim)
    n1 <- sum(nrow(train.sim), nrow(test.sim))
    # p <- ncol(train.sim) - 1
    p <- ncol(train.sim) - 2
    
    # Choosing best lambda
    
    res.rmse1 <- NULL
    res.rmse2 <- NULL
    grid.l <- 10^seq(-3,2,length=10)
    
    for (i in 1:50) {
      
      rf.index <- sample(1:p, round(sqrt(p)), replace=F)
      boot.index <- sample(1:n, n, replace=T)
      train.sim1 <- cbind(train.sim.y.s[boot.index], train.sim.X[boot.index, rf.index], 
                          train.sim.y[boot.index])
      test.sim1 <- cbind(train.sim.y.s[-boot.index], train.sim.X[-boot.index, rf.index],
                         train.sim.y[-boot.index])
      
      u1 <- my.kernel.matrix2(train.sim1, test.sim1)
      boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
      boot.y.train.s <- u1$y.train.s ; boot.y.test.s <- u1$y.test.s
      boot.K <- u1$K
      boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
      
      boot.rmse <- c(rep(0,10))
      
      for (j in 1:10) {
        
        h1.boot <- fit.kernel(boot.y.train.s, boot.y.train, boot.K.train, grid.l[j])
        boot.d.hat <- h1.boot$d.hat
        h2.boot <- pred.kernel(boot.y.test.s, boot.y.test, boot.K.test, boot.d.hat)
        boot.rmse1[j] <- h2.boot$rmse1
        boot.rmse2[j] <- h2.boot$rmse2
        
      }
      
      res.rmse1 <- rbind(res.rmse1, boot.rmse1)
      res.rmse2 <- rbind(res.rmse2, boot.rmse2)
      
    }
    
    #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
    best.lambda1 <- max(grid.l[colMeans(res.rmse1)==min(colMeans(res.rmse1))])
    best.lambda2 <- max(grid.l[colMeans(res.rmse2)==min(colMeans(res.rmse2))])
    
    # Bootstrapping
    
    boots.y1 <- NULL
    boots.y2 <- NULL
    
    for (r in 1:100) {
      
      rfs.index <- sample(1:p, round(sqrt(p)), replace=F)
      boots.index <- sample(1:n, n, replace=T)
      train.sim2 <- cbind(train.sim.y.s[boots.index], train.sim.X[boots.index, rfs.index],
                          train.sim.y[boots.index])
      test.sim2 <- cbind(test.sim.y.s, test.sim.X[, rfs.index], test.sim.y)
      
      u2 <- my.kernel.matrix2(train.sim2, test.sim2)
      boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
      boots.y.train.s <- u2$y.train.s ; boots.y.test.s <- u2$y.test.s
      boots.K <- u2$K
      boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
      
      # synthetic
      h1_1.boots <- fit.kernel(boots.y.train.s, boots.y.train, boots.K.train, best.lambda1)
      boots.d.hat1 <- h1_1.boots$d.hat
      h2_1.boots <- pred.kernel(boots.y.test.s, boots.y.test, boots.K.test, boots.d.hat1)
      boots.y.hat1 <- h2_1.boots$y.hat
      boots.y1 <- cbind(boots.y1, boots.y.hat1)
      
      # original
      h1_2.boots <- fit.kernel(boots.y.train.s, boots.y.train, boots.K.train, best.lambda2)
      boots.d.hat2 <- h1_2.boots$d.hat
      h2_2.boots <- pred.kernel(boots.y.test.s, boots.y.test, boots.K.test, boots.d.hat2)
      boots.y.hat2 <- h2_2.boots$y.hat
      boots.y2 <- cbind(boots.y2, boots.y.hat2)
      
    }
    cat(s, ",")
    
    # Random Forest estimate 
    
    y.hat1.rf <- rowMeans(boots.y1)
    y.hat1.rf <- ifelse(y.hat1.rf<0, 0, y.hat1.rf)
    y.hat2.rf <- rowMeans(boots.y2)
    y.hat2.rf <- ifelse(y.hat2.rf<0, 0, y.hat2.rf)
    
    # KRR[s] <- sqrt(sum((test.sim.y - y.hat.rf)^2)/length(test.sim.y))
    GKRR1[s] <- sqrt(sum((test.sim.y.s - y.hat1.rf)^2)/length(test.sim.y.s))
    GKRR2[s] <- sqrt(sum((test.sim.y - y.hat2.rf)^2)/length(test.sim.y))
    
  }
  
  GKRR1
  GKRR2
  boxplot(GKRR1)
  boxplot(GKRR2)
  dat8_1 <- data.frame(RMSE=GKRR1, method=rep("h.GKRR1",100), number=as.character(rep(n.train,100)))
  dat8_2 <- data.frame(RMSE=GKRR2, method=rep("p.GKRR2",100), number=as.character(rep(n.train,100)))
  dat8 <- rbind(dat8_1, dat8_2)
  
  
  
  
  
  # Final result
  
  dat.res <- rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8)
  
  return(dat.res)
  
}

dat.sim <- dat.sim[dat.sim$method==c('a.LKR1','b.LKRS1','c.LKRB1','d.LKRR1',
                                     'e.GKR1','f.GKRS1','g.GKRB1','h.GKRR1'),]

library(ggplot2)

# Boxplot for one dataset
ggplot(dat.res6_1, aes(x = method, y = RMSE, fill = method)) + geom_boxplot() 
ggplot(dat.res6_2, aes(x = method, y = RMSE, fill = method)) + geom_boxplot()
ggplot(dat.res6_3, aes(x = method, y = RMSE, fill = method)) + geom_boxplot()
# Boxplot for row binded dataset
dat.res6[dat.res6$number==50,]$number <- "1(50)"
dat.res6[dat.res6$number==100,]$number <- "2(100)"
dat.res6[dat.res6$number==200,]$number <- "3(200)"
ggplot(dat.res6, aes(x = number, y = RMSE, fill = number)) + geom_boxplot() +
  facet_wrap(~ method, ncol=16) + theme(axis.text.x=element_text(angle=45, hjust=1))
write.csv(dat.res6, "C:/Users/Administrator/Desktop/211005/p5cen10.csv")


# Print result
mean(dat.res6_1$RMSE[dat.res6_1$method=="a.LKR1"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="a.LKR1"])
mean(dat.res6_1$RMSE[dat.res6_1$method=="b.LKRS1"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="b.LKRS1"])
mean(dat.res6_1$RMSE[dat.res6_1$method=="c.LKRB1"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="c.LKRB1"])
mean(dat.res6_1$RMSE[dat.res6_1$method=="d.LKRR1"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="d.LKRR1"])

mean(dat.res6_1$RMSE[dat.res6_1$method=="e.LKR2"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="e.LKR2"])
mean(dat.res6_1$RMSE[dat.res6_1$method=="f.LKRS2"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="f.LKRS2"])
mean(dat.res6_1$RMSE[dat.res6_1$method=="g.LKRB2"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="g.LKRB2"])
mean(dat.res6_1$RMSE[dat.res6_1$method=="h.LKRR2"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="h.LKRR2"])

mean(dat.res6_1$RMSE[dat.res6_1$method=="i.GKR1"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="i.GKR1"])
mean(dat.res6_1$RMSE[dat.res6_1$method=="j.GKRS1"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="j.GKRS1"])
mean(dat.res6_1$RMSE[dat.res6_1$method=="k.GKRB1"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="k.GKRB1"])
mean(dat.res6_1$RMSE[dat.res6_1$method=="l.GKRR1"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="l.GKRR1"])

mean(dat.res6_1$RMSE[dat.res6_1$method=="m.GKR2"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="m.GKR2"])
mean(dat.res6_1$RMSE[dat.res6_1$method=="n.GKRS2"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="n.GKRS2"])
mean(dat.res6_1$RMSE[dat.res6_1$method=="o.GKRB2"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="o.GKRB2"])
mean(dat.res6_1$RMSE[dat.res6_1$method=="p.GKRR2"])
sd(dat.res6_1$RMSE[dat.res6_1$method=="p.GKRR2"])
#-----------------------------------------------------#

mean(dat.res6_2$RMSE[dat.res6_2$method=="a.LKR1"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="a.LKR1"])
mean(dat.res6_2$RMSE[dat.res6_2$method=="b.LKRS1"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="b.LKRS1"])
mean(dat.res6_2$RMSE[dat.res6_2$method=="c.LKRB1"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="c.LKRB1"])
mean(dat.res6_2$RMSE[dat.res6_2$method=="d.LKRR1"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="d.LKRR1"])

mean(dat.res6_2$RMSE[dat.res6_2$method=="e.LKR2"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="e.LKR2"])
mean(dat.res6_2$RMSE[dat.res6_2$method=="f.LKRS2"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="f.LKRS2"])
mean(dat.res6_2$RMSE[dat.res6_2$method=="g.LKRB2"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="g.LKRB2"])
mean(dat.res6_2$RMSE[dat.res6_2$method=="h.LKRR2"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="h.LKRR2"])

mean(dat.res6_2$RMSE[dat.res6_2$method=="i.GKR1"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="i.GKR1"])
mean(dat.res6_2$RMSE[dat.res6_2$method=="j.GKRS1"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="j.GKRS1"])
mean(dat.res6_2$RMSE[dat.res6_2$method=="k.GKRB1"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="k.GKRB1"])
mean(dat.res6_2$RMSE[dat.res6_2$method=="l.GKRR1"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="l.GKRR1"])

mean(dat.res6_2$RMSE[dat.res6_2$method=="m.GKR2"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="m.GKR2"])
mean(dat.res6_2$RMSE[dat.res6_2$method=="n.GKRS2"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="n.GKRS2"])
mean(dat.res6_2$RMSE[dat.res6_2$method=="o.GKRB2"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="o.GKRB2"])
mean(dat.res6_2$RMSE[dat.res6_2$method=="p.GKRR2"])
sd(dat.res6_2$RMSE[dat.res6_2$method=="p.GKRR2"])
#-----------------------------------------------------#

mean(dat.res6_3$RMSE[dat.res6_3$method=="a.LKR1"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="a.LKR1"])
mean(dat.res6_3$RMSE[dat.res6_3$method=="b.LKRS1"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="b.LKRS1"])
mean(dat.res6_3$RMSE[dat.res6_3$method=="c.LKRB1"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="c.LKRB1"])
mean(dat.res6_3$RMSE[dat.res6_3$method=="d.LKRR1"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="d.LKRR1"])

mean(dat.res6_3$RMSE[dat.res6_3$method=="e.LKR2"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="e.LKR2"])
mean(dat.res6_3$RMSE[dat.res6_3$method=="f.LKRS2"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="f.LKRS2"])
mean(dat.res6_3$RMSE[dat.res6_3$method=="g.LKRB2"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="g.LKRB2"])
mean(dat.res6_3$RMSE[dat.res6_3$method=="h.LKRR2"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="h.LKRR2"])

mean(dat.res6_3$RMSE[dat.res6_3$method=="i.GKR1"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="i.GKR1"])
mean(dat.res6_3$RMSE[dat.res6_3$method=="j.GKRS1"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="j.GKRS1"])
mean(dat.res6_3$RMSE[dat.res6_3$method=="k.GKRB1"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="k.GKRB1"])
mean(dat.res6_3$RMSE[dat.res6_3$method=="l.GKRR1"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="l.GKRR1"])

mean(dat.res6_3$RMSE[dat.res6_3$method=="m.GKR2"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="m.GKR2"])
mean(dat.res6_3$RMSE[dat.res6_3$method=="n.GKRS2"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="n.GKRS2"])
mean(dat.res6_3$RMSE[dat.res6_3$method=="o.GKRB2"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="o.GKRB2"])
mean(dat.res6_3$RMSE[dat.res6_3$method=="p.GKRR2"])
sd(dat.res6_3$RMSE[dat.res6_3$method=="p.GKRR2"])
#-----------------------------------------------------#




#### Summary (p : number of independent variables / n : size of training data)
## Caution -> Please use appropriate 'dat.gen' function


### p=3

## Data generating
# 3 explanatory variables
# censoring --> ( 0% : 20 / 10% : 13  / 30% : 11.76 / 50% : 10.94 )
dat.gen <- function(n, a, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  x1.i <- runif(n,-10,10)
  x2.i <- runif(n,-10,10)
  x3.i <- runif(n,-10,10)
  # e.i <- rnorm(n, mean=0, sd=1)
  y.i <- rnorm(n, mean=10+(x1.i/10)^2+(x2.i/10)^2+(x3.i/10)^2, sd=1)
  y.i <- ifelse(y.i<0, 0, y.i)
  # pred.value <- sin(2*pi*(x1.i*x2.i*x3.i*x4.i*x5.i))+cos(2*pi*(x6.i*x7.i*x8.i*x9.i*x10.i))
  cen <- rnorm(n, mean=a, sd=1) # censoring values
  cen <- ifelse(cen<0, 0, cen)
  pre.T <- pmin(y.i, cen)
  delta <- (y.i <= cen) * 1
  g <- km.surv(pre.T, delta)
  y.i.s <- ifelse(pre.T <= quantile(pre.T, probs=0.98), pre.T*delta/g, 0) # synthetic response (KSV)
  y.i.s <- ifelse(y.i.s<0, 0, y.i.s)
  dat.sim <- data.frame(ys=y.i.s, x1=x1.i, x2=x2.i, x3=x3.i, y=y.i)
  return( dat.sim )
}

## p=3, censoring 0%
dat.res1_1 <- fit.ftn(50, 1000, 20) # n=50
dat.res1_2 <- fit.ftn(100, 1000, 20) # n=100
dat.res1_3 <- fit.ftn(200, 1000, 20) # n=200
dat.res1 <- rbind(dat.res1_1, dat.res1_2, dat.res1_3)

## p=3, censoring 10%
dat.res2_1 <- fit.ftn(50, 1000, 13) # n=50
dat.res2_2 <- fit.ftn(100, 1000, 13) # n=100
dat.res2_3 <- fit.ftn(200, 1000, 13) # n=200
dat.res2 <- rbind(dat.res2_1, dat.res2_2, dat.res2_3)

## p=3, censoring 30%
dat.res3_1 <- fit.ftn(50, 1000, 11.76) # n=50
dat.res3_2 <- fit.ftn(100, 1000, 11.76) # n=100
dat.res3_3 <- fit.ftn(200, 1000, 11.76) # n=200
dat.res3 <- rbind(dat.res3_1, dat.res3_2, dat.res3_3)

## p=3, censoring 50%
dat.res4_1 <- fit.ftn(50, 1000, 10.94) # n=50
dat.res4_2 <- fit.ftn(100, 1000, 10.94) # n=100
dat.res4_3 <- fit.ftn(200, 1000, 10.94) # n=200
dat.res4 <- rbind(dat.res4_1, dat.res4_2, dat.res4_3)


### p=5

## Data generating
# 5 explanatory variables
# censoring --> ( 0% : 20 / 10% : 13.65 / 30% : 12.5 / 50% : 11.61 )
dat.gen <- function(n, a, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  x1.i <- runif(n,-10,10)
  x2.i <- runif(n,-10,10)
  x3.i <- runif(n,-10,10)
  x4.i <- runif(n,-10,10)
  x5.i <- runif(n,-10,10)
  # e.i <- rnorm(n, mean=0, sd=1)
  y.i <- rnorm(n, mean=10+(x1.i/10)^2+(x2.i/10)^2+(x3.i/10)^2+(x4.i/10)^2+(x5.i/10)^2,
               sd=1)
  y.i <- ifelse(y.i<0, 0, y.i)
  # pred.value <- sin(2*pi*(x1.i*x2.i*x3.i*x4.i*x5.i))+cos(2*pi*(x6.i*x7.i*x8.i*x9.i*x10.i))
  cen <- rnorm(n, mean=a, sd=1) # censoring values
  cen <- ifelse(cen<0, 0, cen)
  pre.T <- pmin(y.i, cen)
  delta <- (y.i <= cen) * 1
  g <- km.surv(pre.T, delta)
  y.i.s <- ifelse(pre.T <= quantile(pre.T, probs=0.98), pre.T*delta/g, 0) # synthetic response (KSV)
  y.i.s <- ifelse(y.i.s<0, 0, y.i.s)
  dat.sim <- data.frame(ys=y.i.s, x1=x1.i, x2=x2.i, x3=x3.i, x4=x4.i, x5=x5.i, y=y.i)
  return( dat.sim )
}

## p=5, censoring 0%
dat.res5_1 <- fit.ftn(50, 1000, 20) # n=50
dat.res5_2 <- fit.ftn(100, 1000, 20) # n=100
dat.res5_3 <- fit.ftn(200, 1000, 20) # n=200
dat.res5 <- rbind(dat.res5_1, dat.res5_2, dat.res5_3)

## p=5, censoring 10%
dat.res6_1 <- fit.ftn(50, 1000, 13.65) # n=50
dat.res6_2 <- fit.ftn(100, 1000, 13.65) # n=100
dat.res6_3 <- fit.ftn(200, 1000, 13.65) # n=200
dat.res6 <- rbind(dat.res6_1, dat.res6_2, dat.res6_3)

## p=5, censoring 30%
dat.res7_1 <- fit.ftn(50, 1000, 12.5) # n=50
dat.res7_2 <- fit.ftn(100, 1000, 12.5) # n=100
dat.res7_3 <- fit.ftn(200, 1000, 12.5) # n=200
dat.res7 <- rbind(dat.res7_1, dat.res7_2, dat.res7_3)

## p=5, censoring 50%
dat.res8_1 <- fit.ftn(50, 1000, 11.61) # n=50
dat.res8_2 <- fit.ftn(100, 1000, 11.61) # n=100
dat.res8_3 <- fit.ftn(200, 1000, 11.61) # n=200
dat.res8 <- rbind(dat.res8_1, dat.res8_2, dat.res8_3)


### p=10

## Data generating
# 10 explanatory variables
# censoring --> ( 0% : 20 / 10% : 15.5 / 30% : 14.2 / 50% : 13.31 )
dat.gen <- function(n, a, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  x1.i <- runif(n,-10,10)
  x2.i <- runif(n,-10,10)
  x3.i <- runif(n,-10,10)
  x4.i <- runif(n,-10,10)
  x5.i <- runif(n,-10,10)
  x6.i <- runif(n,-10,10)
  x7.i <- runif(n,-10,10)
  x8.i <- runif(n,-10,10)
  x9.i <- runif(n,-10,10)
  x10.i <- runif(n,-10,10)
  # e.i <- rnorm(n, mean=0, sd=1)
  y.i <- rnorm(n, mean=10+(x1.i/10)^2+(x2.i/10)^2+(x3.i/10)^2+(x4.i/10)^2+(x5.i/10)^2+
                 (x6.i/10)^2+(x7.i/10)^2+(x8.i/10)^2+(x9.i/10)^2+(x10.i/10)^2,
               sd=1)
  y.i <- ifelse(y.i<0, 0, y.i)
  # pred.value <- sin(2*pi*(x1.i*x2.i*x3.i*x4.i*x5.i))+cos(2*pi*(x6.i*x7.i*x8.i*x9.i*x10.i))
  cen <- rnorm(n, mean=a, sd=1) # censoring values
  cen <- ifelse(cen<0, 0, cen)
  pre.T <- pmin(y.i, cen)
  delta <- (y.i <= cen) * 1
  g <- km.surv(pre.T, delta)
  y.i.s <- ifelse(pre.T <= quantile(pre.T, probs=0.98), pre.T*delta/g, 0) # synthetic response (KSV)
  y.i.s <- ifelse(y.i.s<0, 0, y.i.s)
  dat.sim <- data.frame(ys=y.i.s, x1=x1.i, x2=x2.i, x3=x3.i, x4=x4.i, x5=x5.i,
                        x6=x6.i, x7=x7.i, x8=x8.i, x9=x9.i, x10=x10.i, y=y.i)
  return( dat.sim )
}

## p=10, censoring 0%
dat.res9_1 <- fit.ftn(50, 1000, 20) # n=50
dat.res9_2 <- fit.ftn(100, 1000, 20) # n=100
dat.res9_3 <- fit.ftn(200, 1000, 20) # n=200
dat.res9 <- rbind(dat.res9_1, dat.res9_2, dat.res9_3)

## p=10, censoring 10%
dat.res10_1 <- fit.ftn(50, 1000, 15.5) # n=50
dat.res10_2 <- fit.ftn(100, 1000, 15.5) # n=100
dat.res10_3 <- fit.ftn(200, 1000, 15.5) # n=200
dat.res10 <- rbind(dat.res10_1, dat.res10_2, dat.res10_3)

## p=10, censoring 30%
dat.res11_1 <- fit.ftn(50, 1000, 14.2) # n=50
dat.res11_2 <- fit.ftn(100, 1000, 14.2) # n=100
dat.res11_3 <- fit.ftn(200, 1000, 14.2) # n=200
dat.res11 <- rbind(dat.res11_1, dat.res11_2, dat.res11_3)

## p=10, censoring 50%
dat.res12_1 <- fit.ftn(50, 1000, 13.31) # n=50
dat.res12_2 <- fit.ftn(100, 1000, 13.31) # n=100
dat.res12_3 <- fit.ftn(200, 1000, 13.31) # n=200
dat.res12 <- rbind(dat.res12_1, dat.res12_2, dat.res12_3)


### p=20

## Data generating
# 20 explanatory variables
# censoring --> ( 0% : 30 / 10% : 19.1 / 30% : 17.63 / 50% : 16.64 )
dat.gen <- function(n, a, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  x1.i <- runif(n,-10,10)
  x2.i <- runif(n,-10,10)
  x3.i <- runif(n,-10,10)
  x4.i <- runif(n,-10,10)
  x5.i <- runif(n,-10,10)
  x6.i <- runif(n,-10,10)
  x7.i <- runif(n,-10,10)
  x8.i <- runif(n,-10,10)
  x9.i <- runif(n,-10,10)
  x10.i <- runif(n,-10,10)
  x11.i <- runif(n,-10,10)
  x12.i <- runif(n,-10,10)
  x13.i <- runif(n,-10,10)
  x14.i <- runif(n,-10,10)
  x15.i <- runif(n,-10,10)
  x16.i <- runif(n,-10,10)
  x17.i <- runif(n,-10,10)
  x18.i <- runif(n,-10,10)
  x19.i <- runif(n,-10,10)
  x20.i <- runif(n,-10,10)
  # e.i <- rnorm(n, mean=0, sd=1)
  y.i <- rnorm(n, mean=10+(x1.i/10)^2+(x2.i/10)^2+(x3.i/10)^2+(x4.i/10)^2+(x5.i/10)^2+
                 (x6.i/10)^2+(x7.i/10)^2+(x8.i/10)^2+(x9.i/10)^2+(x10.i/10)^2+
                 (x11.i/10)^2+(x12.i/10)^2+(x13.i/10)^2+(x14.i/10)^2+(x15.i/10)^2+
                 (x16.i/10)^2+(x17.i/10)^2+(x18.i/10)^2+(x19.i/10)^2+(x20.i/10)^2,
               sd=1)
  y.i <- ifelse(y.i<0, 0, y.i)
  # pred.value <- sin(2*pi*(x1.i*x2.i*x3.i*x4.i*x5.i))+cos(2*pi*(x6.i*x7.i*x8.i*x9.i*x10.i))
  cen <- rnorm(n, mean=a, sd=1) # censoring values
  cen <- ifelse(cen<0, 0, cen)
  pre.T <- pmin(y.i, cen)
  delta <- (y.i <= cen) * 1
  g <- km.surv(pre.T, delta)
  y.i.s <- ifelse(pre.T <= quantile(pre.T, probs=0.98), pre.T*delta/g, 0) # synthetic response (KSV)
  y.i.s <- ifelse(y.i.s<0, 0, y.i.s)
  dat.sim <- data.frame(ys=y.i.s, x1=x1.i, x2=x2.i, x3=x3.i, x4=x4.i, x5=x5.i,
                        x6=x6.i, x7=x7.i, x8=x8.i, x9=x9.i, x10=x10.i,
                        x11=x11.i, x12=x12.i, x13=x13.i, x14=x14.i, x15=x15.i,
                        x16=x16.i, x17=x17.i, x18=x18.i, x19=x19.i, x20=x20.i, y=y.i)
  return( dat.sim )
}

## p=20, censoring 0%
dat.res13_1 <- fit.ftn(50, 1000, 30) # n=50
dat.res13_2 <- fit.ftn(100, 1000, 30) # n=100
dat.res13_3 <- fit.ftn(200, 1000, 30) # n=200
dat.res13 <- rbind(dat.res13_1, dat.res13_2, dat.res13_3)

## p=20, censoring 10%
dat.res14_1 <- fit.ftn(50, 1000, 19.1) # n=50
dat.res14_2 <- fit.ftn(100, 1000, 19.1) # n=100
dat.res14_3 <- fit.ftn(200, 1000, 19.1) # n=200
dat.res14 <- rbind(dat.res14_1, dat.res14_2, dat.res14_3)

## p=20, censoring 30%
dat.res15_1 <- fit.ftn(50, 1000, 17.63) # n=50
dat.res15_2 <- fit.ftn(100, 1000, 17.63) # n=100
dat.res15_3 <- fit.ftn(200, 1000, 17.63) # n=200
dat.res15 <- rbind(dat.res15_1, dat.res15_2, dat.res15_3)

## p=20, censoring 50%
dat.res16_1 <- fit.ftn(50, 1000, 16.64) # n=50
dat.res16_2 <- fit.ftn(100, 1000, 16.64) # n=100
dat.res16_3 <- fit.ftn(200, 1000, 16.64) # n=200
dat.res16 <- rbind(dat.res16_1, dat.res16_2, dat.res16_3)





