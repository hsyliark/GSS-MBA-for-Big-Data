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

# Synthetic response (KSV)
# x1.i <- seq(0, 1, length.out=500)
# x2.i <- seq(0, 1, length.out=500)
# x3.i <- seq(0, 1, length.out=500)
# x4.i <- seq(0, 1, length.out=500)
# x5.i <- seq(0, 1, length.out=500)
x1.i <- runif(100,-1,1)
x2.i <- runif(100,-1,1)
x3.i <- runif(100,-1,1)
x4.i <- runif(100,-1,1)
x5.i <- runif(100,-1,1)
e.i <- rnorm(100, mean=0, sd=1)
y.i <- sin(2*pi*x1.i*x2.i)+cos(2*pi*x3.i*x4.i*x5.i)+e.i
# pred.value <-  sin(2*pi*seq(0,1,length.out=500)^2)+cos(2*pi*seq(0,1,length.out=500)^3)
pred.value <- sin(2*pi*x1.i*x2.i)+cos(2*pi*x3.i*x4.i*x5.i)
cen.10 <-rnorm(100, mean=sin(2*pi*x1.i*x2.i)+cos(2*pi*x3.i*x4.i*x5.i)+0.84^2, sd=1) # censoring values
pre.T.10 <- pmin(y.i, cen.10)
delta.10 <- (y.i <= cen.10) * 1
g.10 <- km.surv(pre.T.10, delta.10)
y.i.10 <- ifelse(pre.T.10 <= quantile(pre.T.10, probs=0.98), pre.T.10*delta.10/g.10, 0) # synthetic response

dat1 <- data.frame(y=y.i, ys=y.i.10, c=cen.10, d=delta.10, yp=pred.value, x1=x1.i, x2=x2.i, x3=x3.i, x4=x4.i, x5=x5.i)
dat.sim <- data.frame(ys=y.i.10, x1=x1.i, x2=x2.i, x3=x3.i, x4=x4.i, x5=x5.i, y=y.i)


#------------------------#------------------------#------------------------#------------------------



### 1. Making kernel matrix




my.kernel.matrix <- function(dat.train, dat.test) {
  
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
  # K : matrix from Gaussian kernel transformation
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
  
  # y.test : Synthetic response Y* of test data
  # y.test.pred : Original response Y of training data
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
  cv.rmse <- NULL   
  
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




#### 5 method simulation


### Function for 1 independent variable

fit.ftn1 <- function(dat.sim) {
  
  
  
  
  ### 1. Kernel ridge regression (KR)
  
  KR1 <- c(rep(0,100))
  KR2 <- c(rep(0,100))
  
  for (i in 1:100) {
    
    # Making simulation data
    
    set.seed(i)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10) 
    
    
    # 5-fold crossvalidation
    
    u <- my.kernel.matrix(train.sim, test.sim)
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
    
    KR1[i] <- h2_1$rmse1
    KR2[i] <- h2_2$rmse2
    
  }
  
  KR1
  KR2
  boxplot(KR1)
  boxplot(KR2)
  dat1_1 <- data.frame(RMSE=KR1, method=rep("1.KR1",100), number=as.character(rep(n.train,100)))
  dat1_2 <- data.frame(RMSE=KR2, method=rep("2.KR2",100), number=as.character(rep(n.train,100)))
  dat1 <- rbind(dat1_1, dat1_2)
  
  
  
  
  ### 2. Kernel ridge regression using Sub-sampling (KRS)
  
  KRS1 <- c(rep(0,100))
  KRS2 <- c(rep(0,100))
  
  for (i in 1:100) {
    
    # Making simulation data
    
    set.seed(i)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10)
    
    u <- my.kernel.matrix(train.sim, test.sim)
    K.train <- u$K.train ; y.train.s <- u$y.train.s 
    K <- u$K ; K.test <- u$K.test ; y.test.s <- u$y.test.s
    y.train <- u$y.train ; y.test <- u$y.test
    
    # Choosing optimal lambda
    
    
    n <- nrow(train.sim) ; n1 <- sum(nrow(train.sim), nrow(test.sim))
    res.rmse1.c <- NULL ; res.lam1.c <- NULL ; res.index1.c <- NULL
    res.rmse2.c <- NULL ; res.lam2.c <- NULL ; res.index2.c <- NULL
    
    for (j in 1:50) {
      
      n <- nrow(train.sim)
      index.Kc <- sample(1:n, round(0.7*n), replace=F)
      Kc.train <- K.train[index.Kc,index.Kc] ; Kc.test <- K.train[index.Kc,-index.Kc]
      yc.train.s <- y.train.s[index.Kc] ; yc.test.s <- y.train.s[-index.Kc]
      yc.train <- y.train[index.Kc] ; yc.test <- y.train[-index.Kc]
      
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
      
      res.index.c <- cbind(res.index.c, index.Kc)
      
    }
    
    # Fitting  
    
    # synthetic
    res.lambda1 <- res.lam1.c[which.min(res.rmse1.c)]
    res.index1 <- res.index.c[, which.min(res.rmse1.c)]
    res.K.train1 <- K[res.index1, res.index1]
    res.y.train.s1 <- y.train[res.index1]
    res.y.train1 <- y.train[res.index1]
    K.test1 <- K[res.index1, (n+1):n1]
    
    h1_1 <- fit.kernel(res.y.train.s1, res.y.train1, res.K.train1, res.lambda1)
    
    res.d.hat1 <- h1_1$d.hat
    
    # Calculate test mean square error
    
    h2_1 <- pred.kernel(y.test.s, y.test, K.test, res.d.hat1)
    KRS1[i] <- h2_1$rmse1
    
    # original
    res.lambda2 <- res.lam2.c[which.min(res.rmse2.c)]
    res.index2 <- res.index.c[, which.min(res.rmse2.c)]
    res.K.train2 <- K[res.index2, res.index2]
    res.y.train.s2 <- y.train[res.index2]
    res.y.train2 <- y.train[res.index2]
    K.test2 <- K[res.index2, (n+1):n1]
    
    h1_2 <- fit.kernel(res.y.train.s2, res.y.train2, res.K.train2, res.lambda2)
    
    res.d.hat2 <- h1_2$d.hat
    
    # Calculate test mean square error
    
    h2_2 <- pred.kernel(y.test.s, y.test, K.test, res.d.hat2)
    KRS2[i] <- h2_2$rmse2
    
  }
  
  KRS1
  KRS2
  boxplot(KRS1)
  boxplot(KRS2)
  dat2_1 <- data.frame(RMSE=KRS1, method=rep("3.KRS1",100), number=as.character(rep(n.train,100)))
  dat2_2 <- data.frame(RMSE=KRS2, method=rep("4.KRS2",100), number=as.character(rep(n.train,100)))
  dat2 <- rbind(dat2_1, dat2_2)
  
  
  
  
  ### 3. Kernel ridge regression using Bagging (KRB)
  
  KRB1 <- c(rep(0,100))
  KRB2 <- c(rep(0,100))
  
  for (s in 1:100) {
    
    # Making simulation data
    
    set.seed(s)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    test.sim.y <- test.sim[,1] 
    n.train <- round(nrow(dat.sim)*7/10)
    
    u <- my.kernel.matrix(train.sim, test.sim)
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
    KRB1[s] <- sqrt(sum((y.test.s - y.hat1.bag)^2)/length(y.test.s))
    KRB2[s] <- sqrt(sum((y.test - y.hat2.bag)^2)/length(y.test))
    
  }
  
  KRB1
  KRB2
  boxplot(KRB1)
  boxplot(KRB2)
  dat3_1 <- data.frame(RMSE=KRB1, method=rep("5.KRB1",100), number=as.character(rep(n.train,100)))
  dat3_2 <- data.frame(RMSE=KRB2, method=rep("6.KRB2",100), number=as.character(rep(n.train,100)))
  dat3 <- rbind(dat3_1, dat3_2)
  
  
  
  
  
  ### 4. Kernel ridge regression using Random Forest (KRR)
  
  KRR1 <- c(rep(0,100))
  KRR2 <- c(rep(0,100))
  
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
      
      u1 <- my.kernel.matrix(train.sim1, test.sim1)
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
      
      u2 <- my.kernel.matrix(train.sim2, test.sim2)
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
    y.hat2.rf <- rowMeans(boots.y2)
    
    # KRR[s] <- sqrt(sum((test.sim.y - y.hat.rf)^2)/length(test.sim.y))
    KRR1[s] <- sqrt(sum((test.sim.y.s - y.hat1.rf)^2)/length(test.sim.y.s))
    KRR2[s] <- sqrt(sum((test.sim.y - y.hat2.rf)^2)/length(test.sim.y))
    
  }
  
  KRR1
  KRR2
  boxplot(KRR1)
  boxplot(KRR2)
  dat4_1 <- data.frame(RMSE=KRR1, method=rep("7.KRR1",100), number=as.character(rep(n.train,100)))
  dat4_2 <- data.frame(RMSE=KRR2, method=rep("8.KRR2",100), number=as.character(rep(n.train,100)))
  dat4 <- rbind(dat4_1, dat4_2)
  
  
  
  # Final result
  
  dat.res <- rbind(dat1, dat2, dat3, dat4)
  
  return(dat.res)
  
}



### Function for many independent variables



fit.ftn <- function(dat.sim) {
  
  
  
  
  ### 1. Kernel ridge regression (KR)
  
  KR1 <- c(rep(0,100))
  KR2 <- c(rep(0,100))
  
  for (i in 1:100) {
    
    # Making simulation data
    
    set.seed(i)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10) 
    
    
    # 5-fold crossvalidation
    
    u <- my.kernel.matrix(train.sim, test.sim)
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
    
    KR1[i] <- h2_1$rmse1
    KR2[i] <- h2_2$rmse2
    
  }
  
  KR1
  KR2
  boxplot(KR1)
  boxplot(KR2)
  dat1_1 <- data.frame(RMSE=KR1, method=rep("1.KR1",100), number=as.character(rep(n.train,100)))
  dat1_2 <- data.frame(RMSE=KR2, method=rep("2.KR2",100), number=as.character(rep(n.train,100)))
  dat1 <- rbind(dat1_1, dat1_2)
  
  
  
  
  ### 2. Kernel ridge regression using Sub-sampling (KRS)
  
  KRS1 <- c(rep(0,100))
  KRS2 <- c(rep(0,100))
  
  for (i in 1:100) {
    
    # Making simulation data
    
    set.seed(i)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    n.train <- round(nrow(dat.sim)*7/10)
    
    u <- my.kernel.matrix(train.sim, test.sim)
    K.train <- u$K.train ; y.train.s <- u$y.train.s 
    K <- u$K ; K.test <- u$K.test ; y.test.s <- u$y.test.s
    y.train <- u$y.train ; y.test <- u$y.test
    
    # Choosing optimal lambda
    
    
    n <- nrow(train.sim) ; n1 <- sum(nrow(train.sim), nrow(test.sim))
    res.rmse1.c <- NULL ; res.lam1.c <- NULL ; res.index1.c <- NULL
    res.rmse2.c <- NULL ; res.lam2.c <- NULL ; res.index2.c <- NULL
    
    for (j in 1:50) {
      
      n <- nrow(train.sim)
      index.Kc <- sample(1:n, round(0.7*n), replace=F)
      Kc.train <- K.train[index.Kc,index.Kc] ; Kc.test <- K.train[index.Kc,-index.Kc]
      yc.train.s <- y.train.s[index.Kc] ; yc.test.s <- y.train.s[-index.Kc]
      yc.train <- y.train[index.Kc] ; yc.test <- y.train[-index.Kc]
      
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
      
      res.index.c <- cbind(res.index.c, index.Kc)
      
    }
    
    # Fitting  
    
    # synthetic
    res.lambda1 <- res.lam1.c[which.min(res.rmse1.c)]
    res.index1 <- res.index.c[, which.min(res.rmse1.c)]
    res.K.train1 <- K[res.index1, res.index1]
    res.y.train.s1 <- y.train[res.index1]
    res.y.train1 <- y.train[res.index1]
    K.test1 <- K[res.index1, (n+1):n1]
    
    h1_1 <- fit.kernel(res.y.train.s1, res.y.train1, res.K.train1, res.lambda1)
    
    res.d.hat1 <- h1_1$d.hat
    
    # Calculate test mean square error
    
    h2_1 <- pred.kernel(y.test.s, y.test, K.test, res.d.hat1)
    KRS1[i] <- h2_1$rmse1
    
    # original
    res.lambda2 <- res.lam2.c[which.min(res.rmse2.c)]
    res.index2 <- res.index.c[, which.min(res.rmse2.c)]
    res.K.train2 <- K[res.index2, res.index2]
    res.y.train.s2 <- y.train[res.index2]
    res.y.train2 <- y.train[res.index2]
    K.test2 <- K[res.index2, (n+1):n1]
    
    h1_2 <- fit.kernel(res.y.train.s2, res.y.train2, res.K.train2, res.lambda2)
    
    res.d.hat2 <- h1_2$d.hat
    
    # Calculate test mean square error
    
    h2_2 <- pred.kernel(y.test.s, y.test, K.test, res.d.hat2)
    KRS2[i] <- h2_2$rmse2
    
  }
  
  KRS1
  KRS2
  boxplot(KRS1)
  boxplot(KRS2)
  dat2_1 <- data.frame(RMSE=KRS1, method=rep("3.KRS1",100), number=as.character(rep(n.train,100)))
  dat2_2 <- data.frame(RMSE=KRS2, method=rep("4.KRS2",100), number=as.character(rep(n.train,100)))
  dat2 <- rbind(dat2_1, dat2_2)
  
  
  
  
  ### 3. Kernel ridge regression using Bagging (KRB)
  
  KRB1 <- c(rep(0,100))
  KRB2 <- c(rep(0,100))
  
  for (s in 1:100) {
    
    # Making simulation data
    
    set.seed(s)
    train.index <- sample(1:nrow(dat.sim), round(nrow(dat.sim)*7/10), replace=F)
    
    train.sim <- dat.sim[train.index,]
    test.sim <- dat.sim[-train.index,]
    test.sim.y <- test.sim[,1] 
    n.train <- round(nrow(dat.sim)*7/10)
    
    u <- my.kernel.matrix(train.sim, test.sim)
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
    KRB1[s] <- sqrt(sum((y.test.s - y.hat1.bag)^2)/length(y.test.s))
    KRB2[s] <- sqrt(sum((y.test - y.hat2.bag)^2)/length(y.test))
    
  }
  
  KRB1
  KRB2
  boxplot(KRB1)
  boxplot(KRB2)
  dat3_1 <- data.frame(RMSE=KRB1, method=rep("5.KRB1",100), number=as.character(rep(n.train,100)))
  dat3_2 <- data.frame(RMSE=KRB2, method=rep("6.KRB2",100), number=as.character(rep(n.train,100)))
  dat3 <- rbind(dat3_1, dat3_2)
  
  
  
  
  
  
  ### 4. Kernel ridge regression using Random Forest (KRR)
  
  KRR1 <- c(rep(0,100))
  KRR2 <- c(rep(0,100))
  
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
      
      u1 <- my.kernel.matrix(train.sim1, test.sim1)
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
      
      u2 <- my.kernel.matrix(train.sim2, test.sim2)
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
    KRR1[s] <- sqrt(sum((test.sim.y.s - y.hat1.rf)^2)/length(test.sim.y.s))
    KRR2[s] <- sqrt(sum((test.sim.y - y.hat2.rf)^2)/length(test.sim.y))
    
  }
  
  KRR1
  KRR2
  boxplot(KRR1)
  boxplot(KRR2)
  dat4_1 <- data.frame(RMSE=KRR1, method=rep("7.KRR1",100), number=as.character(rep(n.train,100)))
  dat4_2 <- data.frame(RMSE=KRR2, method=rep("8.KRR2",100), number=as.character(rep(n.train,100)))
  dat4 <- rbind(dat4_1, dat4_2)
  
  
  
  # Final result
  
  dat.res <- rbind(dat1, dat2, dat3, dat4)
  
  return(dat.res)
  
}


library(ggplot2)

# dat.res1 <- fit.ftn1(dat.sim)
dat.res2 <- fit.ftn(dat.sim)
ggplot(dat.res2, aes(x = method, y = RMSE, fill = method)) + geom_boxplot()  

mean(dat.res2$RMSE[dat.res2$method=="1.KR1"])
sd(dat.res2$RMSE[dat.res2$method=="1.KR1"])
mean(dat.res2$RMSE[dat.res2$method=="2.KR2"])
sd(dat.res2$RMSE[dat.res2$method=="2.KR2"])
mean(dat.res2$RMSE[dat.res2$method=="3.KRS1"])
sd(dat.res2$RMSE[dat.res2$method=="3.KRS1"])
mean(dat.res2$RMSE[dat.res2$method=="4.KRS2"])
sd(dat.res2$RMSE[dat.res2$method=="4.KRS2"])
mean(dat.res2$RMSE[dat.res2$method=="5.KRB1"])
sd(dat.res2$RMSE[dat.res2$method=="5.KRB1"])
mean(dat.res2$RMSE[dat.res2$method=="6.KRB2"])
sd(dat.res2$RMSE[dat.res2$method=="6.KRB2"])
mean(dat.res2$RMSE[dat.res2$method=="7.KRR1"])
sd(dat.res2$RMSE[dat.res2$method=="7.KRR1"])
mean(dat.res2$RMSE[dat.res2$method=="8.KRR2"])
sd(dat.res2$RMSE[dat.res2$method=="2.KRR2"])

