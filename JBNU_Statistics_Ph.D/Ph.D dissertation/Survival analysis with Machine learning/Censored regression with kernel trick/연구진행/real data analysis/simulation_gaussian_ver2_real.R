## Loading packages

library(survival)
library(KernSmooth)
library(np)
library(locfit)
library(latex2exp)






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



# 1. Polynomial kernel
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
  a <- 1/(p^2) ; b <- 3
  n1 <- nrow(X1)
  D <- as.matrix(X1%*%t(X1)) # Inner product
  K <- (a*D + 1)^b # polynomial kernel
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
  
  # y.test.s : Synthetic response Y* of test data
  # y.test : Original response Y of training data
  # K.test : Kernel matrix from test data 
  # d.hat : Estimator of vector d from training data
  
  
  y.test.s <- as.matrix(y.test.s)
  y.test <- as.matrix(y.test)
  K.test <- as.matrix(K.test)
  d.hat <- as.matrix(d.hat)
  
  
  y.hat <- t(K.test)%*%d.hat
  
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




#### 16 method simulation

# PKR1 : Polynomial Kernel Ridge Regression with Synthetic Response Y*
# PKRS1 : Polynomial Kernel Ridge Regression with Sub-sampling and Synthetic Response Y*
# PKRB1 : Polynomial Kernel Ridge Regression with Bagging and Synthetic Response Y*
# PKRR1 : Polynomial Kernel Ridge Regression with Random Forest and Synthetic Response Y*
  
# GKR1 : Gaussian Kernel Ridge Regression with Synthetic Response Y*
# GKRS1 : Gaussian Kernel Ridge Regression with Sub-sampling and Synthetic Response Y*
# GKRB1 : Gaussian Kernel Ridge Regression with Bagging and Synthetic Response Y*
# GKRR1 : Gaussian Kernel Ridge Regression with Random Forest and Synthetic Response Y*

# PKR2 : Polynomial Kernel Ridge Regression with Generated(Original) Response Y
# PKRS2 : Polynomial Kernel Ridge Regression with Sub-sampling and Generated(Original) Response Y
# PKRB2 : Polynomial Kernel Ridge Regression with Bagging and Generated(Original) Response Y
# PKRR2 : Polynomial Kernel Ridge Regression with Random Forest and Generated(Original) Response Y

# GKR2 : Gaussian Kernel Ridge Regression with Generated(Original) Response Y 
# GKRS2 : Gaussian Kernel Ridge Regression with Sub-sampling and Generated(Original) Response Y
# GKRB2 : Gaussian Kernel Ridge Regression with Bagging and Generated(Original) Response Y
# GKRR2 : Gaussian Kernel Ridge Regression with Random Forest and Generated(Original) Response Y





### Function for 1 independent variable

fit.ftn1 <- function(dat.sim) {
  
  
  
  ### 1. Kernel ridge regression (KR)
  
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
  dat2_1 <- data.frame(RMSE=GKR1, method=rep("a.GKR1",100), number=as.character(rep(n.train,100)))
  dat2_2 <- data.frame(RMSE=GKR2, method=rep("e.GKR2",100), number=as.character(rep(n.train,100)))
  dat2 <- rbind(dat2_1, dat2_2)
  
  
  
  
  ### 2. Kernel ridge regression using Sub-sampling (KRS)
  
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
  dat4_1 <- data.frame(RMSE=GKRS1, method=rep("b.GKRS1",100), number=as.character(rep(n.train,100)))
  dat4_2 <- data.frame(RMSE=GKRS2, method=rep("f.GKRS2",100), number=as.character(rep(n.train,100)))
  dat4 <- rbind(dat4_1, dat4_2)
  
  
  
  
  ### 3. Kernel ridge regression using Bagging (KRB)
  
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
    y.hat2.bag <- rowMeans(boots.y2)
    
    # KRB[s] <- sqrt(sum((y.test - y.hat.bag)^2)/length(y.test))
    GKRB1[s] <- sqrt(sum((y.test.s - y.hat1.bag)^2)/length(y.test.s))
    GKRB2[s] <- sqrt(sum((y.test - y.hat2.bag)^2)/length(y.test))
    
  }
  
  GKRB1
  GKRB2
  boxplot(GKRB1)
  boxplot(GKRB2)
  dat6_1 <- data.frame(RMSE=GKRB1, method=rep("c.GKRB1",100), number=as.character(rep(n.train,100)))
  dat6_2 <- data.frame(RMSE=GKRB2, method=rep("g.GKRB2",100), number=as.character(rep(n.train,100)))
  dat6 <- rbind(dat6_1, dat6_2)
  
  
  
  
  
  
  
  ### 4. Kernel ridge regression using Random Forest (KRR)
  
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
      h1_1.boots <- fit.kernel(boots.y.train.s, boots.y.train, boots.K.train, best.lambda1)
      boots.d.hat1 <- h1_1.boots$d.hat
      h2_1.boots <- pred.kernel(boots.y.test.s, boots.y.test, boots.K.test, boots.d.hat1)
      boots.y.hat1 <- h2_1.boots$y.hat
      boots.y1 <- cbind(boots.y1, boots.y.hat1)
      
      # original
      h1_2.boots <- fit.kernel(boots.y.train.s, boots.y.train, boots.K.train, best.lambda2)
      boots.d.hat1 <- h1_1.boots$d.hat
      h2_2.boots <- pred.kernel(boots.y.test.s, boots.y.test, boots.K.test, boots.d.hat2)
      boots.y.hat2 <- h2_2.boots$y.hat
      boots.y2 <- cbind(boots.y2, boots.y.hat2)
      
    }
    cat(s, ",")
    
    # Random Forest estimate 
    
    y.hat1.rf <- rowMeans(boots.y1)
    y.hat2.rf <- rowMeans(boots.y2)
    
    # KRR[s] <- sqrt(sum((test.sim.y - y.hat.rf)^2)/length(test.sim.y))
    GKRR1[s] <- sqrt(sum((test.sim.y.s - y.hat1.rf)^2)/length(test.sim.y.s))
    GKRR2[s] <- sqrt(sum((test.sim.y - y.hat2.rf)^2)/length(test.sim.y))
    
  }
  
  GKRR1
  GKRR2
  boxplot(GKRR1)
  boxplot(GKRR2)
  dat8_1 <- data.frame(RMSE=GKRR1, method=rep("d.GKRR1",100), number=as.character(rep(n.train,100)))
  dat8_2 <- data.frame(RMSE=GKRR2, method=rep("h.GKRR2",100), number=as.character(rep(n.train,100)))
  dat8 <- rbind(dat8_1, dat8_2)
  
  
  
  
  # Final result
  
  dat.res <- rbind(dat2, dat4, dat6, dat8)
  
  return(dat.res)
  
}



### Function for many independent variables



fit.ftn <- function(dat.sim) {
  
  
  
  
  ### 1. Kernel ridge regression (KR)
  
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
  dat2_1 <- data.frame(RMSE=GKR1, method=rep("a.GKR1",100), number=as.character(rep(n.train,100)))
  dat2_2 <- data.frame(RMSE=GKR2, method=rep("e.GKR2",100), number=as.character(rep(n.train,100)))
  dat2 <- rbind(dat2_1, dat2_2)
  
  
  
  
  ### 2. Kernel ridge regression using Sub-sampling (KRS)
  
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
  dat4_1 <- data.frame(RMSE=GKRS1, method=rep("b.GKRS1",100), number=as.character(rep(n.train,100)))
  dat4_2 <- data.frame(RMSE=GKRS2, method=rep("f.GKRS2",100), number=as.character(rep(n.train,100)))
  dat4 <- rbind(dat4_1, dat4_2)
  
  
  
  
  
  ### 3. Kernel ridge regression using Bagging (KRB)
  
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
    y.hat2.bag <- rowMeans(boots.y2)
    
    # KRB[s] <- sqrt(sum((y.test - y.hat.bag)^2)/length(y.test))
    GKRB1[s] <- sqrt(sum((y.test.s - y.hat1.bag)^2)/length(y.test.s))
    GKRB2[s] <- sqrt(sum((y.test - y.hat2.bag)^2)/length(y.test))
    
  }
  
  GKRB1
  GKRB2
  boxplot(GKRB1)
  boxplot(GKRB2)
  dat6_1 <- data.frame(RMSE=GKRB1, method=rep("c.GKRB1",100), number=as.character(rep(n.train,100)))
  dat6_2 <- data.frame(RMSE=GKRB2, method=rep("g.GKRB2",100), number=as.character(rep(n.train,100)))
  dat6 <- rbind(dat6_1, dat6_2)
  
  
  
  
  ### 4. Kernel ridge regression using Random Forest (KRR)
  
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
    y.hat2.rf <- rowMeans(boots.y2)
    
    # KRR[s] <- sqrt(sum((test.sim.y - y.hat.rf)^2)/length(test.sim.y))
    GKRR1[s] <- sqrt(sum((test.sim.y.s - y.hat1.rf)^2)/length(test.sim.y.s))
    GKRR2[s] <- sqrt(sum((test.sim.y - y.hat2.rf)^2)/length(test.sim.y))
    
  }
  
  GKRR1
  GKRR2
  boxplot(GKRR1)
  boxplot(GKRR2)
  dat8_1 <- data.frame(RMSE=GKRR1, method=rep("d.GKRR1",100), number=as.character(rep(n.train,100)))
  dat8_2 <- data.frame(RMSE=GKRR2, method=rep("h.GKRR2",100), number=as.character(rep(n.train,100)))
  dat8 <- rbind(dat8_1, dat8_2)
  
  
  
  
  
  # Final result
  
  dat.res <- rbind(dat2, dat4, dat6, dat8)
  
  return(dat.res)
  
}


library(ggplot2)




uis2 <- read.csv("C:/Users/Hi/Desktop/uis2.csv",
                    sep=",",header=T)
uis2 <-uis2[,-1]
g <- km.surv(uis2$time, uis2$censor)
uis2$y.s <- ifelse(uis2$time <= quantile(uis2$time, probs=0.98), 
                   uis2$time*uis2$censor/g, 0) 
dat.sim <- data.frame(ys1=uis2$y.s, x1=uis2$age, x2=uis2$beck, 
                      x3=uis2$hercoc, x4=uis2$ivhx, x5=uis2$ndrugtx,
                      x6=uis2$race, x7=uis2$treat, x8=uis2$lot,
                      ys2=uis2$y.s) 
dat.res <- fit.ftn(dat.sim)
dat.res1 <- dat.res[dat.res$method=="a.GKR1",]
dat.res2 <- dat.res[dat.res$method=="b.GKRS1",]
dat.res3 <- dat.res[dat.res$method=="c.GKRB1",]
dat.res4 <- dat.res[dat.res$method=="d.GKRR1",]
dat.res_final <- rbind(dat.res1, dat.res2, dat.res3, dat.res4)
write.csv(dat.res_final, "C:/Users/Hi/Desktop/real data analysis/uis2_res.csv")

ggplot(dat.res_final, aes(x = method, y = RMSE, fill = method)) + geom_boxplot()

mean(dat.res_final$RMSE[dat.res_final$method=="a.GKR1"])
sd(dat.res_final$RMSE[dat.res_final$method=="a.GKR1"])
mean(dat.res_final$RMSE[dat.res_final$method=="b.GKRS1"])
sd(dat.res_final$RMSE[dat.res_final$method=="b.GKRS1"])
mean(dat.res_final$RMSE[dat.res_final$method=="c.GKRB1"])
sd(dat.res_final$RMSE[dat.res_final$method=="c.GKRB1"])
mean(dat.res_final$RMSE[dat.res_final$method=="d.GKRR1"])
sd(dat.res_final$RMSE[dat.res_final$method=="d.GKRR1"])


#-----------------------------------------------------#


pbc <- read.csv("C:/Users/Hi/Desktop/pbc.csv",
                 sep=",",header=T)
g <- km.surv(pbc$time, pbc$status)
pbc$y.s <- ifelse(pbc$time <= quantile(pbc$time, probs=0.98), 
                   pbc$time*pbc$status/g, 0) 
dat.sim <- data.frame(ys1=pbc$y.s, x1=pbc$trt, x2=pbc$age, x3=pbc$sex,
                      x4=pbc$ascites, x5=pbc$hepato, x6=pbc$spiders,
                      x7=pbc$edema, x8=pbc$bili, x9=pbc$chol,
                      x10=pbc$albumin, x11=pbc$copper, x12=pbc$alk.phos,
                      x13=pbc$ast, x14=pbc$trig, x15=pbc$platelet,
                      x16=pbc$protime, x17=pbc$stage,
                      ys2=pbc$y.s) 
dat.res <- fit.ftn(dat.sim)
dat.res1 <- dat.res[dat.res$method=="a.GKR1",]
dat.res2 <- dat.res[dat.res$method=="b.GKRS1",]
dat.res3 <- dat.res[dat.res$method=="c.GKRB1",]
dat.res4 <- dat.res[dat.res$method=="d.GKRR1",]
dat.res_final <- rbind(dat.res1, dat.res2, dat.res3, dat.res4)
write.csv(dat.res_final, "C:/Users/Hi/Desktop/real data analysis/pbc_res.csv")

ggplot(dat.res_final, aes(x = method, y = RMSE, fill = method)) + geom_boxplot()

mean(dat.res_final$RMSE[dat.res_final$method=="a.GKR1"])
sd(dat.res_final$RMSE[dat.res_final$method=="a.GKR1"])
mean(dat.res_final$RMSE[dat.res_final$method=="b.GKRS1"])
sd(dat.res_final$RMSE[dat.res_final$method=="b.GKRS1"])
mean(dat.res_final$RMSE[dat.res_final$method=="c.GKRB1"])
sd(dat.res_final$RMSE[dat.res_final$method=="c.GKRB1"])
mean(dat.res_final$RMSE[dat.res_final$method=="d.GKRR1"])
sd(dat.res_final$RMSE[dat.res_final$method=="d.GKRR1"])


#-----------------------------------------------------#


cancer <- read.csv("C:/Users/Hi/Desktop/cancer.csv",
                sep=",",header=T)
g <- km.surv(cancer$time, cancer$status)
cancer$y.s <- ifelse(cancer$time <= quantile(cancer$time, probs=0.98), 
                  cancer$time*cancer$status/g, 0) 
dat.sim <- data.frame(ys1=cancer$y.s, x1=cancer$age, x2=cancer$sex,
                      x3=cancer$ph.ecog, x4=cancer$ph.karno, 
                      x5=cancer$pat.karno, x6=cancer$meal.cal,
                      x7=cancer$wt.loss,
                      ys2=cancer$y.s) 
dat.res <- fit.ftn(dat.sim)
dat.res1 <- dat.res[dat.res$method=="a.GKR1",]
dat.res2 <- dat.res[dat.res$method=="b.GKRS1",]
dat.res3 <- dat.res[dat.res$method=="c.GKRB1",]
dat.res4 <- dat.res[dat.res$method=="d.GKRR1",]
dat.res_final <- rbind(dat.res1, dat.res2, dat.res3, dat.res4)
write.csv(dat.res_final, "C:/Users/Hi/Desktop/real data analysis/cancer_res.csv")

ggplot(dat.res_final, aes(x = method, y = RMSE, fill = method)) + geom_boxplot()

mean(dat.res_final$RMSE[dat.res_final$method=="a.GKR1"])
sd(dat.res_final$RMSE[dat.res_final$method=="a.GKR1"])
mean(dat.res_final$RMSE[dat.res_final$method=="b.GKRS1"])
sd(dat.res_final$RMSE[dat.res_final$method=="b.GKRS1"])
mean(dat.res_final$RMSE[dat.res_final$method=="c.GKRB1"])
sd(dat.res_final$RMSE[dat.res_final$method=="c.GKRB1"])
mean(dat.res_final$RMSE[dat.res_final$method=="d.GKRR1"])
sd(dat.res_final$RMSE[dat.res_final$method=="d.GKRR1"])


#-----------------------------------------------------#









#### Summary (p : number of independent variables / n : size of training data)
## Caution -> Please use appropriate 'dat.gen' function


### p=3

## Data generating
# 3 explanatory variables
# censoring --> ( 0% : 20 / 10% : 3.4  / 30% : 2.15 / 50% : 1.33 )
dat.gen <- function(n, a, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  x1.i <- runif(n,-10,10)
  x2.i <- runif(n,-10,10)
  x3.i <- runif(n,-10,10)
  # e.i <- rnorm(n, mean=0, sd=1)
  y.i <- rnorm(n, mean=1+(x1.i/10)^1+(x2.i/10)^2+(x3.i/10)^3, sd=1)
  # pred.value <- sin(2*pi*(x1.i*x2.i*x3.i*x4.i*x5.i))+cos(2*pi*(x6.i*x7.i*x8.i*x9.i*x10.i))
  cen <- rnorm(n, mean=a, sd=1) # censoring values
  pre.T <- pmin(y.i, cen)
  delta <- (y.i <= cen) * 1
  g <- km.surv(pre.T, delta)
  y.i.s <- ifelse(pre.T <= quantile(pre.T, probs=0.98), pre.T*delta/g, 0) # synthetic response (KSV)
  dat.sim <- data.frame(ys=y.i.s, x1=x1.i, x2=x2.i, x3=x3.i, y=y.i)
  return( dat.sim )
}

## p=3, censoring 0%
dat.res1_1 <- fit.ftn(50, 1000, 20) # n=50
dat.res1_2 <- fit.ftn(100, 1000, 20) # n=100
dat.res1_3 <- fit.ftn(200, 1000, 20) # n=200
dat.res1 <- rbind(dat.res1_1, dat.res1_2, dat.res1_3)

## p=3, censoring 10%
dat.res2_1 <- fit.ftn(50, 1000, 3.4) # n=50
dat.res2_2 <- fit.ftn(100, 1000, 3.4) # n=100
dat.res2_3 <- fit.ftn(200, 1000, 3.4) # n=200
dat.res2 <- rbind(dat.res2_1, dat.res2_2, dat.res2_3)

## p=3, censoring 30%
dat.res3_1 <- fit.ftn(50, 1000, 2.15) # n=50
dat.res3_2 <- fit.ftn(100, 1000, 2.15) # n=100
dat.res3_3 <- fit.ftn(200, 1000, 2.15) # n=200
dat.res3 <- rbind(dat.res3_1, dat.res3_2, dat.res3_3)

## p=3, censoring 50%
dat.res4_1 <- fit.ftn(50, 1000, 1.33) # n=50
dat.res4_2 <- fit.ftn(100, 1000, 1.33) # n=100
dat.res4_3 <- fit.ftn(200, 1000, 1.33) # n=200
dat.res4 <- rbind(dat.res4_1, dat.res4_2, dat.res4_3)


### p=5

## Data generating
# 5 explanatory variables
# censoring --> ( 0% : 20 / 10% : 3.7 / 30% : 2.4 / 50% : 1.5 )
dat.gen <- function(n, a, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  x1.i <- runif(n,-10,10)
  x2.i <- runif(n,-10,10)
  x3.i <- runif(n,-10,10)
  x4.i <- runif(n,-10,10)
  x5.i <- runif(n,-10,10)
  # e.i <- rnorm(n, mean=0, sd=1)
  y.i <- rnorm(n, mean=1+(x1.i/10)^1+(x2.i/10)^2+(x3.i/10)^3+(x4.i/10)^4+(x5.i/10)^5,
               sd=1)
  # pred.value <- sin(2*pi*(x1.i*x2.i*x3.i*x4.i*x5.i))+cos(2*pi*(x6.i*x7.i*x8.i*x9.i*x10.i))
  cen <- rnorm(n, mean=a, sd=1) # censoring values
  pre.T <- pmin(y.i, cen)
  delta <- (y.i <= cen) * 1
  g <- km.surv(pre.T, delta)
  y.i.s <- ifelse(pre.T <= quantile(pre.T, probs=0.98), pre.T*delta/g, 0) # synthetic response (KSV)
  dat.sim <- data.frame(ys=y.i.s, x1=x1.i, x2=x2.i, x3=x3.i, x4=x4.i, x5=x5.i, y=y.i)
  return( dat.sim )
}

## p=5, censoring 0%
dat.res5_1 <- fit.ftn(50, 1000, 20) # n=50
dat.res5_2 <- fit.ftn(100, 1000, 20) # n=100
dat.res5_3 <- fit.ftn(200, 1000, 20) # n=200
dat.res5 <- rbind(dat.res5_1, dat.res5_2, dat.res5_3)

## p=5, censoring 10%
dat.res6_1 <- fit.ftn(50, 1000, 3.7) # n=50
dat.res6_2 <- fit.ftn(100, 1000, 3.7) # n=100
dat.res6_3 <- fit.ftn(200, 1000, 3.7) # n=200
dat.res6 <- rbind(dat.res6_1, dat.res6_2, dat.res6_3)

## p=5, censoring 30%
dat.res7_1 <- fit.ftn(50, 1000, 2.4) # n=50
dat.res7_2 <- fit.ftn(100, 1000, 2.4) # n=100
dat.res7_3 <- fit.ftn(200, 1000, 2.4) # n=200
dat.res7 <- rbind(dat.res7_1, dat.res7_2, dat.res7_3)

## p=5, censoring 50%
dat.res8_1 <- fit.ftn(50, 1000, 1.5) # n=50
dat.res8_2 <- fit.ftn(100, 1000, 1.5) # n=100
dat.res8_3 <- fit.ftn(200, 1000, 1.5) # n=200
dat.res8 <- rbind(dat.res8_1, dat.res8_2, dat.res8_3)


### p=7

## Data generating
# 7 explanatory variables
# censoring --> ( 0% : 30 / 10% : 3.8 / 30% : 2.55 / 50% : 1.65 )
dat.gen <- function(n, a, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  x1.i <- runif(n,-10,10)
  x2.i <- runif(n,-10,10)
  x3.i <- runif(n,-10,10)
  x4.i <- runif(n,-10,10)
  x5.i <- runif(n,-10,10)
  x6.i <- runif(n,-10,10)
  x7.i <- runif(n,-10,10)
  # e.i <- rnorm(n, mean=0, sd=1)
  y.i <- rnorm(n, mean=1+(x1.i/10)^1+(x2.i/10)^2+(x3.i/10)^3+(x4.i/10)^4+(x5.i/10)^5+
                 (x6.i/10)^6+(x7.i/10)^7,
               sd=1)
  # pred.value <- sin(2*pi*(x1.i*x2.i*x3.i*x4.i*x5.i))+cos(2*pi*(x6.i*x7.i*x8.i*x9.i*x10.i))
  cen <- rnorm(n, mean=a, sd=1) # censoring values
  pre.T <- pmin(y.i, cen)
  delta <- (y.i <= cen) * 1
  g <- km.surv(pre.T, delta)
  y.i.s <- ifelse(pre.T <= quantile(pre.T, probs=0.98), pre.T*delta/g, 0) # synthetic response (KSV)
  dat.sim <- data.frame(ys=y.i.s, x1=x1.i, x2=x2.i, x3=x3.i, x4=x4.i, x5=x5.i,
                        x6=x6.i, x7=x7.i, y=y.i)
  return( dat.sim )
}

## p=7, censoring 0%
dat.res9_1 <- fit.ftn(50, 1000, 30) # n=50
dat.res9_2 <- fit.ftn(100, 1000, 30) # n=100
dat.res9_3 <- fit.ftn(200, 1000, 30) # n=200
dat.res9 <- rbind(dat.res9_1, dat.res9_2, dat.res9_3)

## p=7, censoring 10%
dat.res10_1 <- fit.ftn(50, 1000, 3.8) # n=50
dat.res10_2 <- fit.ftn(100, 1000, 3.8) # n=100
dat.res10_3 <- fit.ftn(200, 1000, 3.8) # n=200
dat.res10 <- rbind(dat.res10_1, dat.res10_2, dat.res10_3)

## p=7, censoring 30%
dat.res11_1 <- fit.ftn(50, 1000, 2.55) # n=50
dat.res11_2 <- fit.ftn(100, 1000, 2.55) # n=100
dat.res11_3 <- fit.ftn(200, 1000, 2.55) # n=200
dat.res11 <- rbind(dat.res11_1, dat.res11_2, dat.res11_3)

## p=7, censoring 50%
dat.res12_1 <- fit.ftn(50, 1000, 1.65) # n=50
dat.res12_2 <- fit.ftn(100, 1000, 1.65) # n=100
dat.res12_3 <- fit.ftn(200, 1000, 1.65) # n=200
dat.res12 <- rbind(dat.res12_1, dat.res12_2, dat.res12_3)


### p=9

## Data generating
# 9 explanatory variables
# censoring --> ( 0% : 30 / 10% : 4 / 30% : 2.7 / 50% : 1.75 )
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
  # e.i <- rnorm(n, mean=0, sd=1)
  y.i <- rnorm(n, mean=1+(x1.i/10)^1+(x2.i/10)^2+(x3.i/10)^3+(x4.i/10)^4+(x5.i/10)^5+
                 (x6.i/10)^6+(x7.i/10)^7+(x8.i/10)^8+(x9.i/10)^9,
               sd=1)
  # pred.value <- sin(2*pi*(x1.i*x2.i*x3.i*x4.i*x5.i))+cos(2*pi*(x6.i*x7.i*x8.i*x9.i*x10.i))
  cen <- rnorm(n, mean=a, sd=1) # censoring values
  pre.T <- pmin(y.i, cen)
  delta <- (y.i <= cen) * 1
  g <- km.surv(pre.T, delta)
  y.i.s <- ifelse(pre.T <= quantile(pre.T, probs=0.98), pre.T*delta/g, 0) # synthetic response (KSV)
  dat.sim <- data.frame(ys=y.i.s, x1=x1.i, x2=x2.i, x3=x3.i, x4=x4.i, x5=x5.i,
                        x6=x6.i, x7=x7.i, x8=x8.i, x9=x9.i, y=y.i)
  return( dat.sim )
}

## p=9, censoring 0%
dat.res13_1 <- fit.ftn(50, 1000, 30) # n=50
dat.res13_2 <- fit.ftn(100, 1000, 30) # n=100
dat.res13_3 <- fit.ftn(200, 1000, 30) # n=200
dat.res13 <- rbind(dat.res13_1, dat.res13_2, dat.res13_3)

## p=9, censoring 10%
dat.res14_1 <- fit.ftn(50, 1000, 4) # n=50
dat.res14_2 <- fit.ftn(100, 1000, 4) # n=100
dat.res14_3 <- fit.ftn(200, 1000, 4) # n=200
dat.res14 <- rbind(dat.res14_1, dat.res14_2, dat.res14_3)

## p=9, censoring 30%
dat.res15_1 <- fit.ftn(50, 1000, 2.7) # n=50
dat.res15_2 <- fit.ftn(100, 1000, 2.7) # n=100
dat.res15_3 <- fit.ftn(200, 1000, 2.7) # n=200
dat.res15 <- rbind(dat.res15_1, dat.res15_2, dat.res15_3)

## p=9, censoring 50%
dat.res16_1 <- fit.ftn(50, 1000, 1.75) # n=50
dat.res16_2 <- fit.ftn(100, 1000, 1.75) # n=100
dat.res16_3 <- fit.ftn(200, 1000, 1.75) # n=200
dat.res16 <- rbind(dat.res16_1, dat.res16_2, dat.res16_3)





