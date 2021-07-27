#### Simulation data

dat.gen <- function(n, beta, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  x1 <- runif(n,-10,10) ; x2 <- runif(n,-10,10) ; x3 <- runif(n,-10,10)
  x4 <- runif(n,-10,10) ; x5 <- runif(n,-10,10) ; x6 <- runif(n,-10,10)
  x7 <- runif(n,-10,10) ; x8 <- runif(n,-10,10) ; x9 <- runif(n,-10,10)
  x10 <- runif(n,-10,10) ; x11 <- runif(n,-10,10) ; x12 <- runif(n,-10,10)
  x13 <- runif(n,-10,10) ; x14 <- runif(n,-10,10) ; x15 <- runif(n,-10,10)
  x16 <- runif(n,-10,10) ; x17 <- runif(n,-10,10) ; x18 <- runif(n,-10,10)
  x19 <- runif(n,-10,10) ; x20 <- runif(n,-10,10)
  
  x <- cbind(1, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
             x11, x12, x13, x14, x15, x16, x17, x18, x19, x20)
  linpred <- x%*%beta
  p <- exp(x%*%beta) / (1+exp(x%*%beta))
  y <- apply(as.matrix(p), c(1,2), rbinom, n=1, size=1)
  dat <- data.frame(y=y, x1=x1, x2=x2, x3=x3, x4=x4, x5=x5,
                    x6=x6, x7=x7, x8=x8, x9=x9, x10=x10,
                    x11=x11, x12=x12, x13=x13, x14=x14, x15=x15,
                    x16=x16, x17=x17, x18=x18, x19=x19, x20=x20, linpred=linpred)
  return( dat )
}




#------------------------#------------------------#------------------------#------------------------





#### Making some functions






### 1. Making kernel matrix


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
  G.kernel <- function(m) exp(-sigma*m^2) # Gaussian kernel
  
  K <- G.kernel(D)
  K.train <- K[1:n,1:n]
  K.test <- K[1:n,(n+1):n1]
  
  result <- list(K.train=K.train, K.test=K.test, K=K,
                 y.train=y.train, y.test=y.test)
  
}








### 2. Calculate coefficients using theory


## Maximum Likelihood Estimation + Newton-Raphson algorithm
## Iteratively Reweighted Least Squares (IRLS)


my.kernel.logistic <- function(y, K, lambda) {
  
  # y : response variable
  # K : matrix from Gaussian kernel transformation
  # lambda : penalty parameter
  
  if(lambda <= 0) 
    stop("Lambda is non-positive value. Please insert positive value of lambda.")
  
  y <- as.matrix(y) ; K <-as.matrix(K)
  
  d.old <- rep(1,ncol(K)) # The initial value of coefficient vector
  d.old <- d.old/sqrt(sum(d.old^2))
  
  for (i in 1:100) {
    
    
    p <- matrix(NA,nrow(K),1)
    p1 <- matrix(NA,nrow(K),1)
    
    V <- K%*%d.old
    
    my.prob1 <- function(m) exp(m)/(1+exp(m))
    my.prob2 <- function(m) (exp(m)/(1+exp(m)))*(1-exp(m)/(1+exp(m)))
    
    p <- my.prob1(V)
    p1 <- my.prob2(V)
    
    W.inverse <- diag(1/as.vector(p1))
    
    d.new <- solve(K + lambda*W.inverse)%*%(K%*%d.old + W.inverse%*%(y-p))
    
    new.logit.hat <- K%*%d.new ; old.logit.hat <- K%*%d.old
    
    log.like <- function(y, logit.hat) {
      sum(y*logit.hat - log(1+exp(logit.hat)))
    } # Log-likelihood function
    diff <- abs(log.like(y, new.logit.hat) - log.like(y, old.logit.hat))  
    
    # diff <- sqrt(sum((d.new - d.old)^2))/sqrt(sum(d.old^2))
    
    cat("( iteration , difference ) = (", i, ",", diff, ")\n")
    
    if (diff < 1E-8) break
    
    d.old <- d.new
    
  }
  
  cat("Algorithm converged...","\n\n")
  
  
  logit.hat <- K%*%d.new
  pi.hat <- exp(logit.hat)/(1+exp(logit.hat))
  y.pred <- rep(1,length(logit.hat))
  y.pred[pi.hat < 0.5] <- 0
  
  res <- list(d.hat=d.new, y.pred=y.pred, pi.hat=pi.hat)
  
}







### 3. Making function for fitting


fit.kernel <- function(y.train, K.train, lambda) {
  
  # y.train : Dependent variable of training data
  # K.train : Kernel matrix from training data 
  # lambda : penalty parameter
  
  if(lambda <= 0) 
    stop("Lambda is non-positive value. Please insert positive value of lambda.")
  
  
  y.train <- as.matrix(y.train)
  K.train <- as.matrix(K.train) 
  
  g <- my.kernel.logistic(y.train, K.train, lambda)
  
  d.hat <- g$d.hat
  y.pred <- g$y.pred
  pi.hat <- g$pi.hat
  
  rate.miss <- mean(y.pred != y.train)
  
  
  result <- list(d.hat=d.hat, y.train=y.train, y.pred=y.pred,
                 misclassification.rate=rate.miss, pi.hat=pi.hat)
  
}







### 4. Making function for predict


pred.kernel <- function(y.test, K.test, d.hat) {
  
  # y.test : Dependent variable of test data
  # K.test : Kernel matrix from test data 
  # d.hat : Estimator of vector d from training data
  
  
  y.test <- as.matrix(y.test) 
  K.test <- as.matrix(K.test)
  d.hat <- as.matrix(d.hat)
  
  
  logit.hat <- t(K.test)%*%d.hat
  pi.hat <- exp(logit.hat)/(1+exp(logit.hat))
  
  predict.new <- rep(1,length(logit.hat))
  predict.new[pi.hat < 0.5] <- 0
  
  rate.miss <- mean(predict.new != y.test)
  
  
  result <- list(y.test=y.test, predict.new=predict.new, logit.hat=logit.hat,
                 misclassification.rate=rate.miss, pi.hat.new=pi.hat)
  
}







### 5. Making function for K-fold crossvalidation


cv.kernel <- function(y.train, K.train, k, grid.l) {
  
  
  # y.train : Dependent variable of training data
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
  y.sim <- as.matrix(y.train)
  n <- nrow(K.sim)
  
  cv.index <- sample(1:n,n,replace=F)  
  cv.logL <- NULL   
  
  cat("K-fold crossvalidation is start...","\n")
  
  
  for (j in 1:r) {
    
    logL <- NULL # minus log-likelihood
    
    
    for (i in 0:(k-1)) {
      
      
      test.index <- cv.index[(1:n)%/%k==i]
      
      K.sim.train <- K.sim[-test.index, -test.index] ; K.sim.test <- K.sim[-test.index, test.index]
      y.sim.train <- y.sim[-test.index,] ; y.sim.test <- y.sim[test.index,]
      test.size <- length(test.index)
      
      
      a1 <- fit.kernel(y.sim.train, K.sim.train, lambda[j])
      train.d.hat <- a1$d.hat
      
      a2 <- pred.kernel(y.sim.test, K.sim.test, train.d.hat)
      test.logit.hat <- a2$logit.hat      
      
      
      logL <- c(logL, -sum(y.sim.test*test.logit.hat - log(1+exp(test.logit.hat))) )
      
      
    }
    
    cv.logL <- rbind(cv.logL, logL)
    cat(j,"\n","\n")
  }
  
  cat("\n","K-fold crossvalidation complete...")
  
  se.logL <- apply(cv.logL, 1, sd)
  mean.logL <- rowMeans(cv.logL)
  idx <- which.min(mean.logL)
  
  plot(log(grid.l, base=10),rowMeans(cv.logL),xlab="log10(lambda)",ylab="Minus log-likelihood"
       ,main="K-fold crossvalidation",type="b")
  abline(h=mean.logL[idx]+se.logL[idx], col="red", lty=2)
  abline(v=log(max(lambda[ mean.logL < mean.logL[idx]+se.logL[idx] ]), base=10), col="blue", lty=2)
  
  result <- list(lambda=grid.l, cv.logL=cv.logL)
  
  
}






#------------------------#------------------------#------------------------#------------------------
#------------------------#------------------------#------------------------#------------------------








#### SIMULATION (n=50)





### 1. Kernel ridge logistic regression classification (KRLC)


KRLC <- c(rep(0,50))

for (i in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(50, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  
  # 5-fold crossvalidation
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train ; K.test <- u$K.test  
  y.train <- u$y.train ; y.test <- u$y.test
  
  grid.l <- 10^seq(-6,1,length=100)
  
  h <- cv.kernel(y.train, K.train, 5, grid.l) 
  
  # Fitting 
  
  se.logL <- apply(h$cv.logL, 1, sd)
  mean.logL <- rowMeans(h$cv.logL)
  idx <- which.min(mean.logL)
  best.lam <- max(h$lambda[ mean.logL < mean.logL[idx]+se.logL[idx] ])
  
  h1 <- fit.kernel(y.train, K.train, best.lam)
  
  # Calculate test misclassification rate
  
  sim.d.hat <- h1$d.hat
  
  h2 <- pred.kernel(y.test, K.test, sim.d.hat)
  KRLC[i] <- h2$misclassification.rate
  
}

KRLC
boxplot(KRLC)

dat1.1 <- data.frame(test.MR=KRLC, method=rep("KRLC",50), number=rep("50",50))








### 2. Kernel ridge logistic regression classification using Sub-sampling (KRLCS)


KRLCS <- c(rep(0,50))

for (i in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(50, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train ; y.train <- u$y.train 
  K <- u$K ; y.test <- u$y.test
  
  # Choosing optimal lambda
  
  
  n <- nrow(train.sim) ; n1 <- sum(nrow(train.sim), nrow(test.sim))
  res.rate.c <- NULL ; res.lam.c <- NULL ; res.index.c <- NULL
  
  for (j in 1:20) {
    
    n <- nrow(train.sim)
    index.Kc <- sample(1:n, round(0.7*n), replace=F)
    Kc.train <- K.train[index.Kc,index.Kc] ; Kc.test <- K.train[index.Kc,-index.Kc]
    yc.train <- y.train[index.Kc] ; yc.test <- y.train[-index.Kc]
    
    grid.l <- 10^seq(-6,1,length=100)
    
    hc <- cv.kernel(yc.train, Kc.train, 5, grid.l) # 5-fold crossvalidation 
    
    se.logL.c <- apply(hc$cv.logL, 1, sd)
    mean.logL.c <- rowMeans(hc$cv.logL)
    idx.c <- which.min(mean.logL.c)
    best.lam.c <- max(hc$lambda[ mean.logL.c < mean.logL.c[idx.c]+se.logL.c[idx.c] ])
    
    h1.c <- fit.kernel(yc.train, Kc.train, best.lam.c) # Fitting
    
    sim.d.hat.c <- h1.c$d.hat
    
    h2.c <- pred.kernel(yc.test, Kc.test, sim.d.hat.c) # Testing
    rate.miss.c <- h2.c$misclassification.rate
    
    res.rate.c <- c(res.rate.c, rate.miss.c)
    res.lam.c <- c(res.lam.c, best.lam.c)
    res.index.c <- cbind(res.index.c, index.Kc)
    
  }
  
  # Fitting  
  
  res.lambda <- res.lam.c[which.min(res.rate.c)]
  res.index <- res.index.c[, which.min(res.rate.c)]
  res.K.train <- K[res.index, res.index]
  res.y.train <- y.train[res.index]
  K.test <- K[res.index, (n+1):n1]
  
  h1 <- fit.kernel(res.y.train, res.K.train, res.lambda)
  
  res.d.hat <- h1$d.hat
  
  # Calculate test misclassification rate
  
  h2 <- pred.kernel(y.test, K.test, res.d.hat)
  KRLCS[i] <- h2$misclassification.rate
  
}

KRLCS
boxplot(KRLCS)

dat1.2 <- data.frame(test.MR=KRLCS, method=rep("KRLCS",50), number=rep("50",50))








### 3. Kernel ridge logistic regression classification using Bagging with logit estimate (KRLCB1)



KRLCB1 <- c(rep(0,50))

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(50, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  test.sim.y <- test.sim[,1]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train 
  y.train <- u$y.train ; y.test <- u$y.test 
  K <- u$K
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    boot.index <- sample(1:n, n, replace=T)
    boot.K.train <- K.train[boot.index, boot.index]
    boot.K.test <- K.train[boot.index, -boot.index]
    boot.y.train <- y.train[boot.index]
    boot.y.test <- y.train[-boot.index]
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.logit <- NULL
  
  for (r in 1:100) {
    
    boots.index <- sample(1:n, n, replace=T)
    boots.K.train <- K[boots.index, boots.index]
    boots.K.test <- K[boots.index, (n+1):n1]
    boots.y.train <- y.train[boots.index]
    boots.y.test <- y.test 
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.logit.hat <- h2.boots$logit.hat
    
    boots.logit <- cbind(boots.logit, boots.logit.hat)
    
  }
  
  # Bagging estimate
  
  logit.hat.bag <- rowMeans(boots.logit)
  pi.hat.bag <- exp(logit.hat.bag)/(1+exp(logit.hat.bag))
  
  y.hat.bag <- c(rep(1,n1-n))
  y.hat.bag[pi.hat.bag < 0.5] <- 0
  y.hat.bag[pi.hat.bag == 0.5] <- rbinom(length(y.hat.bag[pi.hat.bag == 0.5]), 1, 0.5)
  
  KRLCB1[s] <- mean(y.hat.bag != y.test)
  
}

KRLCB1
boxplot(KRLCB1)

dat1.3 <- data.frame(test.MR=KRLCB1, method=rep("KRLCB1",50), number=rep("50",50))









### 4. Kernel ridge logistic regression classification using Bagging with probability estimate (KRLCB2)



KRLCB2 <- c(rep(0,50))

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(50, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  test.sim.y <- test.sim[,1]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train 
  y.train <- u$y.train ; y.test <- u$y.test 
  K <- u$K
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    boot.index <- sample(1:n, n, replace=T)
    boot.K.train <- K.train[boot.index, boot.index]
    boot.K.test <- K.train[boot.index, -boot.index]
    boot.y.train <- y.train[boot.index]
    boot.y.test <- y.train[-boot.index]
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.pi <- NULL
  
  for (r in 1:100) {
    
    boots.index <- sample(1:n, n, replace=T)
    boots.K.train <- K[boots.index, boots.index]
    boots.K.test <- K[boots.index, (n+1):n1]
    boots.y.train <- y.train[boots.index]
    boots.y.test <- y.test 
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.pi.hat <- h2.boots$pi.hat.new
    
    boots.pi <- cbind(boots.pi, boots.pi.hat)
    
  }
  
  # Bagging estimate
  
  pi.hat.bag <- rowMeans(boots.pi)
  
  y.hat.bag <- c(rep(1,n1-n))
  y.hat.bag[pi.hat.bag < 0.5] <- 0
  y.hat.bag[pi.hat.bag == 0.5] <- rbinom(length(y.hat.bag[pi.hat.bag == 0.5]), 1, 0.5)
  
  KRLCB2[s] <- mean(y.hat.bag != y.test)
  
}

KRLCB2
boxplot(KRLCB2)

dat1.4 <- data.frame(test.MR=KRLCB2, method=rep("KRLCB2",50), number=rep("50",50))








### 5. Kernel ridge logistic regression classification using Bagging with majority vote (KRLCB3)



KRLCB3 <- c(rep(0,50))

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(50, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  test.sim.y <- test.sim[,1]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train 
  y.train <- u$y.train ; y.test <- u$y.test 
  K <- u$K
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    boot.index <- sample(1:n, n, replace=T)
    boot.K.train <- K.train[boot.index, boot.index]
    boot.K.test <- K.train[boot.index, -boot.index]
    boot.y.train <- y.train[boot.index]
    boot.y.test <- y.train[-boot.index]
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.predict.new <- NULL
  
  for (r in 1:100) {
    
    boots.index <- sample(1:n, n, replace=T)
    boots.K.train <- K[boots.index, boots.index]
    boots.K.test <- K[boots.index, (n+1):n1]
    boots.y.train <- y.train[boots.index]
    boots.y.test <- y.test 
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.predict <- h2.boots$predict.new
    
    boots.predict.new <- rbind(boots.predict.new, boots.predict)
    
  }
  
  # Bagging estimate
  
  mean.hat.bag <- apply(boots.predict.new, 2, mean)
  y.hat.bag <- c(rep(1,n1-n))
  y.hat.bag[mean.hat.bag < 0.5] <- 0
  y.hat.bag[mean.hat.bag == 0.5] <- rbinom(length(y.hat.bag[mean.hat.bag == 0.5]), 1, 0.5)
  
  KRLCB3[s] <- mean(y.hat.bag != y.test)
  
}

KRLCB3
boxplot(KRLCB3)

dat1.5 <- data.frame(test.MR=KRLCB3, method=rep("KRLCB3",50), number=rep("50",50))








### 6. Kernel ridge logistic regression classification using Random Forest with logit estimate (KRLCR1)



KRLCR1 <- c(rep(0,50))  

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(50, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim.y <- train.sim[,1] ; test.sim.y <- test.sim[,1]
  train.sim.X <- train.sim[,2:21] ; test.sim.X <- test.sim[,2:21]
  
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  p <- ncol(train.sim.X)
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    rf.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boot.index <- sample(1:n, n, replace=T)
    train.sim1 <- cbind(train.sim.y[boot.index], train.sim.X[boot.index, rf.index])
    test.sim1 <- cbind(train.sim.y[-boot.index], train.sim.X[-boot.index, rf.index])
    
    u1 <- my.kernel.matrix(train.sim1, test.sim1)
    boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
    boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.logit <- NULL
  
  for (r in 1:100) {
    
    rfs.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boots.index <- sample(1:n, n, replace=T)
    train.sim2 <- cbind(train.sim.y[boots.index], train.sim.X[boots.index, rfs.index])
    test.sim2 <- cbind(test.sim.y, test.sim.X[, rfs.index])
    
    u2 <- my.kernel.matrix(train.sim2, test.sim2)
    boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
    boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.logit.hat <- h2.boots$logit.hat
    
    boots.logit <- cbind(boots.logit, boots.logit.hat)
    
  }
  
  # Random Forest estimate
  
  logit.hat.rf <- rowMeans(boots.logit)
  pi.hat.rf <- exp(logit.hat.rf)/(1+exp(logit.hat.rf))
  
  y.hat.rf <- c(rep(1,n1-n))
  y.hat.rf[pi.hat.rf < 0.5] <- 0
  y.hat.rf[pi.hat.rf == 0.5] <- rbinom(length(y.hat.rf[pi.hat.rf == 0.5]), 1, 0.5)
  
  KRLCR1[s] <- mean(y.hat.rf != test.sim.y)
  
}

KRLCR1
boxplot(KRLCR1)

dat1.6 <- data.frame(test.MR=KRLCR1, method=rep("KRLCR1",50), number=rep("50",50))







### 7. Kernel ridge logistic regression classification using Random Forest with probability estimate (KRLCR2)



KRLCR2 <- c(rep(0,50))  

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(50, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim.y <- train.sim[,1] ; test.sim.y <- test.sim[,1]
  train.sim.X <- train.sim[,2:21] ; test.sim.X <- test.sim[,2:21]
  
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  p <- ncol(train.sim.X)
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    rf.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boot.index <- sample(1:n, n, replace=T)
    train.sim1 <- cbind(train.sim.y[boot.index], train.sim.X[boot.index, rf.index])
    test.sim1 <- cbind(train.sim.y[-boot.index], train.sim.X[-boot.index, rf.index])
    
    u1 <- my.kernel.matrix(train.sim1, test.sim1)
    boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
    boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.pi <- NULL
  
  for (r in 1:100) {
    
    rfs.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boots.index <- sample(1:n, n, replace=T)
    train.sim2 <- cbind(train.sim.y[boots.index], train.sim.X[boots.index, rfs.index])
    test.sim2 <- cbind(test.sim.y, test.sim.X[, rfs.index])
    
    u2 <- my.kernel.matrix(train.sim2, test.sim2)
    boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
    boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.pi.hat <- h2.boots$pi.hat.new
    
    boots.pi <- cbind(boots.pi, boots.pi.hat)
    
  }
  
  # Random Forest estimate
  
  pi.hat.rf <- rowMeans(boots.pi)
  
  y.hat.rf <- c(rep(1,n1-n))
  y.hat.rf[pi.hat.rf < 0.5] <- 0
  y.hat.rf[pi.hat.rf == 0.5] <- rbinom(length(y.hat.rf[pi.hat.rf == 0.5]), 1, 0.5)
  
  KRLCR2[s] <- mean(y.hat.rf != test.sim.y)
  
}

KRLCR2
boxplot(KRLCR2)

dat1.7 <- data.frame(test.MR=KRLCR2, method=rep("KRLCR2",50), number=rep("50",50))







### 8. Kernel ridge logistic regression classification using Random Forest with majority vote (KRLCR3)



KRLCR3 <- c(rep(0,50))  

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(50, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim.y <- train.sim[,1] ; test.sim.y <- test.sim[,1]
  train.sim.X <- train.sim[,2:21] ; test.sim.X <- test.sim[,2:21]
  
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  p <- ncol(train.sim.X)
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    rf.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boot.index <- sample(1:n, n, replace=T)
    train.sim1 <- cbind(train.sim.y[boot.index], train.sim.X[boot.index, rf.index])
    test.sim1 <- cbind(train.sim.y[-boot.index], train.sim.X[-boot.index, rf.index])
    
    u1 <- my.kernel.matrix(train.sim1, test.sim1)
    boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
    boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.predict.new <- NULL
  
  for (r in 1:100) {
    
    rfs.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boots.index <- sample(1:n, n, replace=T)
    train.sim2 <- cbind(train.sim.y[boots.index], train.sim.X[boots.index, rfs.index])
    test.sim2 <- cbind(test.sim.y, test.sim.X[, rfs.index])
    
    u2 <- my.kernel.matrix(train.sim2, test.sim2)
    boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
    boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.predict <- h2.boots$predict.new
    
    boots.predict.new <- rbind(boots.predict.new, boots.predict)
    
  }
  
  # Random Forest estimate
  
  mean.hat.rf <- apply(boots.predict.new, 2, mean)
  y.hat.rf <- c(rep(1,n1-n))
  y.hat.rf[mean.hat.rf < 0.5] <- 0
  y.hat.rf[mean.hat.rf == 0.5] <- rbinom(length(y.hat.rf[mean.hat.rf == 0.5]), 1, 0.5)
  
  KRLCR3[s] <- mean(y.hat.rf != test.sim.y)
  
}

KRLCR3
boxplot(KRLCR3)

dat1.8 <- data.frame(test.MR=KRLCR3, method=rep("KRLCR3",50), number=rep("50",50))





dat1 <- rbind(dat1.1, dat1.2, dat1.3, dat1.4, dat1.5, dat1.6, dat1.7, dat1.8)













#### SIMULATION (n=100)





### 1. Kernel ridge logistic regression classification (KRLC)


KRLC <- c(rep(0,50))

for (i in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(100, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  
  # 5-fold crossvalidation
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train ; K.test <- u$K.test  
  y.train <- u$y.train ; y.test <- u$y.test
  
  grid.l <- 10^seq(-6,1,length=100)
  
  h <- cv.kernel(y.train, K.train, 5, grid.l) 
  
  # Fitting 
  
  se.logL <- apply(h$cv.logL, 1, sd)
  mean.logL <- rowMeans(h$cv.logL)
  idx <- which.min(mean.logL)
  best.lam <- max(h$lambda[ mean.logL < mean.logL[idx]+se.logL[idx] ])
  
  h1 <- fit.kernel(y.train, K.train, best.lam)
  
  # Calculate test misclassification rate
  
  sim.d.hat <- h1$d.hat
  
  h2 <- pred.kernel(y.test, K.test, sim.d.hat)
  KRLC[i] <- h2$misclassification.rate
  
}

KRLC
boxplot(KRLC)

dat2.1 <- data.frame(test.MR=KRLC, method=rep("KRLC",50), number=rep("100",50))








### 2. Kernel ridge logistic regression classification using Sub-sampling (KRLCS)


KRLCS <- c(rep(0,50))

for (i in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(100, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train ; y.train <- u$y.train 
  K <- u$K ; y.test <- u$y.test
  
  # Choosing optimal lambda
  
  
  n <- nrow(train.sim) ; n1 <- sum(nrow(train.sim), nrow(test.sim))
  res.rate.c <- NULL ; res.lam.c <- NULL ; res.index.c <- NULL
  
  for (j in 1:20) {
    
    n <- nrow(train.sim)
    index.Kc <- sample(1:n, round(0.7*n), replace=F)
    Kc.train <- K.train[index.Kc,index.Kc] ; Kc.test <- K.train[index.Kc,-index.Kc]
    yc.train <- y.train[index.Kc] ; yc.test <- y.train[-index.Kc]
    
    grid.l <- 10^seq(-6,1,length=100)
    
    hc <- cv.kernel(yc.train, Kc.train, 5, grid.l) # 5-fold crossvalidation 
    
    se.logL.c <- apply(hc$cv.logL, 1, sd)
    mean.logL.c <- rowMeans(hc$cv.logL)
    idx.c <- which.min(mean.logL.c)
    best.lam.c <- max(hc$lambda[ mean.logL.c < mean.logL.c[idx.c]+se.logL.c[idx.c] ])
    
    h1.c <- fit.kernel(yc.train, Kc.train, best.lam.c) # Fitting
    
    sim.d.hat.c <- h1.c$d.hat
    
    h2.c <- pred.kernel(yc.test, Kc.test, sim.d.hat.c) # Testing
    rate.miss.c <- h2.c$misclassification.rate
    
    res.rate.c <- c(res.rate.c, rate.miss.c)
    res.lam.c <- c(res.lam.c, best.lam.c)
    res.index.c <- cbind(res.index.c, index.Kc)
    
  }
  
  # Fitting  
  
  res.lambda <- res.lam.c[which.min(res.rate.c)]
  res.index <- res.index.c[, which.min(res.rate.c)]
  res.K.train <- K[res.index, res.index]
  res.y.train <- y.train[res.index]
  K.test <- K[res.index, (n+1):n1]
  
  h1 <- fit.kernel(res.y.train, res.K.train, res.lambda)
  
  res.d.hat <- h1$d.hat
  
  # Calculate test misclassification rate
  
  h2 <- pred.kernel(y.test, K.test, res.d.hat)
  KRLCS[i] <- h2$misclassification.rate
  
}

KRLCS
boxplot(KRLCS)

dat2.2 <- data.frame(test.MR=KRLCS, method=rep("KRLCS",50), number=rep("100",50))








### 3. Kernel ridge logistic regression classification using Bagging with logit estimate (KRLCB1)



KRLCB1 <- c(rep(0,50))

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(100, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  test.sim.y <- test.sim[,1]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train 
  y.train <- u$y.train ; y.test <- u$y.test 
  K <- u$K
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    boot.index <- sample(1:n, n, replace=T)
    boot.K.train <- K.train[boot.index, boot.index]
    boot.K.test <- K.train[boot.index, -boot.index]
    boot.y.train <- y.train[boot.index]
    boot.y.test <- y.train[-boot.index]
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.logit <- NULL
  
  for (r in 1:100) {
    
    boots.index <- sample(1:n, n, replace=T)
    boots.K.train <- K[boots.index, boots.index]
    boots.K.test <- K[boots.index, (n+1):n1]
    boots.y.train <- y.train[boots.index]
    boots.y.test <- y.test 
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.logit.hat <- h2.boots$logit.hat
    
    boots.logit <- cbind(boots.logit, boots.logit.hat)
    
  }
  
  # Bagging estimate
  
  logit.hat.bag <- rowMeans(boots.logit)
  pi.hat.bag <- exp(logit.hat.bag)/(1+exp(logit.hat.bag))
  
  y.hat.bag <- c(rep(1,n1-n))
  y.hat.bag[pi.hat.bag < 0.5] <- 0
  y.hat.bag[pi.hat.bag == 0.5] <- rbinom(length(y.hat.bag[pi.hat.bag == 0.5]), 1, 0.5)
  
  KRLCB1[s] <- mean(y.hat.bag != y.test)
  
}

KRLCB1
boxplot(KRLCB1)

dat2.3 <- data.frame(test.MR=KRLCB1, method=rep("KRLCB1",50), number=rep("100",50))









### 4. Kernel ridge logistic regression classification using Bagging with probability estimate (KRLCB2)



KRLCB2 <- c(rep(0,50))

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(100, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  test.sim.y <- test.sim[,1]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train 
  y.train <- u$y.train ; y.test <- u$y.test 
  K <- u$K
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    boot.index <- sample(1:n, n, replace=T)
    boot.K.train <- K.train[boot.index, boot.index]
    boot.K.test <- K.train[boot.index, -boot.index]
    boot.y.train <- y.train[boot.index]
    boot.y.test <- y.train[-boot.index]
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.pi <- NULL
  
  for (r in 1:100) {
    
    boots.index <- sample(1:n, n, replace=T)
    boots.K.train <- K[boots.index, boots.index]
    boots.K.test <- K[boots.index, (n+1):n1]
    boots.y.train <- y.train[boots.index]
    boots.y.test <- y.test 
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.pi.hat <- h2.boots$pi.hat.new
    
    boots.pi <- cbind(boots.pi, boots.pi.hat)
    
  }
  
  # Bagging estimate
  
  pi.hat.bag <- rowMeans(boots.pi)
  
  y.hat.bag <- c(rep(1,n1-n))
  y.hat.bag[pi.hat.bag < 0.5] <- 0
  y.hat.bag[pi.hat.bag == 0.5] <- rbinom(length(y.hat.bag[pi.hat.bag == 0.5]), 1, 0.5)
  
  KRLCB2[s] <- mean(y.hat.bag != y.test)
  
}

KRLCB2
boxplot(KRLCB2)

dat2.4 <- data.frame(test.MR=KRLCB2, method=rep("KRLCB2",50), number=rep("100",50))








### 5. Kernel ridge logistic regression classification using Bagging with majority vote (KRLCB3)



KRLCB3 <- c(rep(0,50))

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(100, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  test.sim.y <- test.sim[,1]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train 
  y.train <- u$y.train ; y.test <- u$y.test 
  K <- u$K
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    boot.index <- sample(1:n, n, replace=T)
    boot.K.train <- K.train[boot.index, boot.index]
    boot.K.test <- K.train[boot.index, -boot.index]
    boot.y.train <- y.train[boot.index]
    boot.y.test <- y.train[-boot.index]
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.predict.new <- NULL
  
  for (r in 1:100) {
    
    boots.index <- sample(1:n, n, replace=T)
    boots.K.train <- K[boots.index, boots.index]
    boots.K.test <- K[boots.index, (n+1):n1]
    boots.y.train <- y.train[boots.index]
    boots.y.test <- y.test 
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.predict <- h2.boots$predict.new
    
    boots.predict.new <- rbind(boots.predict.new, boots.predict)
    
  }
  
  # Bagging estimate
  
  mean.hat.bag <- apply(boots.predict.new, 2, mean)
  y.hat.bag <- c(rep(1,n1-n))
  y.hat.bag[mean.hat.bag < 0.5] <- 0
  y.hat.bag[mean.hat.bag == 0.5] <- rbinom(length(y.hat.bag[mean.hat.bag == 0.5]), 1, 0.5)
  
  KRLCB3[s] <- mean(y.hat.bag != y.test)
  
}

KRLCB3
boxplot(KRLCB3)

dat2.5 <- data.frame(test.MR=KRLCB3, method=rep("KRLCB3",50), number=rep("100",50))








### 6. Kernel ridge logistic regression classification using Random Forest with logit estimate (KRLCR1)



KRLCR1 <- c(rep(0,50))  

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(100, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim.y <- train.sim[,1] ; test.sim.y <- test.sim[,1]
  train.sim.X <- train.sim[,2:21] ; test.sim.X <- test.sim[,2:21]
  
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  p <- ncol(train.sim.X)
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    rf.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boot.index <- sample(1:n, n, replace=T)
    train.sim1 <- cbind(train.sim.y[boot.index], train.sim.X[boot.index, rf.index])
    test.sim1 <- cbind(train.sim.y[-boot.index], train.sim.X[-boot.index, rf.index])
    
    u1 <- my.kernel.matrix(train.sim1, test.sim1)
    boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
    boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.logit <- NULL
  
  for (r in 1:100) {
    
    rfs.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boots.index <- sample(1:n, n, replace=T)
    train.sim2 <- cbind(train.sim.y[boots.index], train.sim.X[boots.index, rfs.index])
    test.sim2 <- cbind(test.sim.y, test.sim.X[, rfs.index])
    
    u2 <- my.kernel.matrix(train.sim2, test.sim2)
    boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
    boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.logit.hat <- h2.boots$logit.hat
    
    boots.logit <- cbind(boots.logit, boots.logit.hat)
    
  }
  
  # Random Forest estimate
  
  logit.hat.rf <- rowMeans(boots.logit)
  pi.hat.rf <- exp(logit.hat.rf)/(1+exp(logit.hat.rf))
  
  y.hat.rf <- c(rep(1,n1-n))
  y.hat.rf[pi.hat.rf < 0.5] <- 0
  y.hat.rf[pi.hat.rf == 0.5] <- rbinom(length(y.hat.rf[pi.hat.rf == 0.5]), 1, 0.5)
  
  KRLCR1[s] <- mean(y.hat.rf != test.sim.y)
  
}

KRLCR1
boxplot(KRLCR1)

dat2.6 <- data.frame(test.MR=KRLCR1, method=rep("KRLCR1",50), number=rep("100",50))







### 7. Kernel ridge logistic regression classification using Random Forest with probability estimate (KRLCR2)



KRLCR2 <- c(rep(0,50))  

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(100, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim.y <- train.sim[,1] ; test.sim.y <- test.sim[,1]
  train.sim.X <- train.sim[,2:21] ; test.sim.X <- test.sim[,2:21]
  
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  p <- ncol(train.sim.X)
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    rf.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boot.index <- sample(1:n, n, replace=T)
    train.sim1 <- cbind(train.sim.y[boot.index], train.sim.X[boot.index, rf.index])
    test.sim1 <- cbind(train.sim.y[-boot.index], train.sim.X[-boot.index, rf.index])
    
    u1 <- my.kernel.matrix(train.sim1, test.sim1)
    boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
    boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.pi <- NULL
  
  for (r in 1:100) {
    
    rfs.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boots.index <- sample(1:n, n, replace=T)
    train.sim2 <- cbind(train.sim.y[boots.index], train.sim.X[boots.index, rfs.index])
    test.sim2 <- cbind(test.sim.y, test.sim.X[, rfs.index])
    
    u2 <- my.kernel.matrix(train.sim2, test.sim2)
    boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
    boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.pi.hat <- h2.boots$pi.hat.new
    
    boots.pi <- cbind(boots.pi, boots.pi.hat)
    
  }
  
  # Random Forest estimate
  
  pi.hat.rf <- rowMeans(boots.pi)
  
  y.hat.rf <- c(rep(1,n1-n))
  y.hat.rf[pi.hat.rf < 0.5] <- 0
  y.hat.rf[pi.hat.rf == 0.5] <- rbinom(length(y.hat.rf[pi.hat.rf == 0.5]), 1, 0.5)
  
  KRLCR2[s] <- mean(y.hat.rf != test.sim.y)
  
}

KRLCR2
boxplot(KRLCR2)

dat2.7 <- data.frame(test.MR=KRLCR2, method=rep("KRLCR2",50), number=rep("100",50))







### 8. Kernel ridge logistic regression classification using Random Forest with majority vote (KRLCR3)



KRLCR3 <- c(rep(0,50))  

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(100, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim.y <- train.sim[,1] ; test.sim.y <- test.sim[,1]
  train.sim.X <- train.sim[,2:21] ; test.sim.X <- test.sim[,2:21]
  
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  p <- ncol(train.sim.X)
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    rf.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boot.index <- sample(1:n, n, replace=T)
    train.sim1 <- cbind(train.sim.y[boot.index], train.sim.X[boot.index, rf.index])
    test.sim1 <- cbind(train.sim.y[-boot.index], train.sim.X[-boot.index, rf.index])
    
    u1 <- my.kernel.matrix(train.sim1, test.sim1)
    boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
    boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.predict.new <- NULL
  
  for (r in 1:100) {
    
    rfs.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boots.index <- sample(1:n, n, replace=T)
    train.sim2 <- cbind(train.sim.y[boots.index], train.sim.X[boots.index, rfs.index])
    test.sim2 <- cbind(test.sim.y, test.sim.X[, rfs.index])
    
    u2 <- my.kernel.matrix(train.sim2, test.sim2)
    boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
    boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.predict <- h2.boots$predict.new
    
    boots.predict.new <- rbind(boots.predict.new, boots.predict)
    
  }
  
  # Random Forest estimate
  
  mean.hat.rf <- apply(boots.predict.new, 2, mean)
  y.hat.rf <- c(rep(1,n1-n))
  y.hat.rf[mean.hat.rf < 0.5] <- 0
  y.hat.rf[mean.hat.rf == 0.5] <- rbinom(length(y.hat.rf[mean.hat.rf == 0.5]), 1, 0.5)
  
  KRLCR3[s] <- mean(y.hat.rf != test.sim.y)
  
}

KRLCR3
boxplot(KRLCR3)

dat2.8 <- data.frame(test.MR=KRLCR3, method=rep("KRLCR3",50), number=rep("100",50))





dat2 <- rbind(dat2.1, dat2.2, dat2.3, dat2.4, dat2.5, dat2.6, dat2.7, dat2.8)






#### SIMULATION (n=200)





### 1. Kernel ridge logistic regression classification (KRLC)


KRLC <- c(rep(0,50))

for (i in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(200, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  
  # 5-fold crossvalidation
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train ; K.test <- u$K.test  
  y.train <- u$y.train ; y.test <- u$y.test
  
  grid.l <- 10^seq(-6,1,length=100)
  
  h <- cv.kernel(y.train, K.train, 5, grid.l) 
  
  # Fitting 
  
  se.logL <- apply(h$cv.logL, 1, sd)
  mean.logL <- rowMeans(h$cv.logL)
  idx <- which.min(mean.logL)
  best.lam <- max(h$lambda[ mean.logL < mean.logL[idx]+se.logL[idx] ])
  
  h1 <- fit.kernel(y.train, K.train, best.lam)
  
  # Calculate test misclassification rate
  
  sim.d.hat <- h1$d.hat
  
  h2 <- pred.kernel(y.test, K.test, sim.d.hat)
  KRLC[i] <- h2$misclassification.rate
  
}

KRLC
boxplot(KRLC)

dat3.1 <- data.frame(test.MR=KRLC, method=rep("KRLC",50), number=rep("200",50))








### 2. Kernel ridge logistic regression classification using Sub-sampling (KRLCS)


KRLCS <- c(rep(0,50))

for (i in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(200, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train ; y.train <- u$y.train 
  K <- u$K ; y.test <- u$y.test
  
  # Choosing optimal lambda
  
  
  n <- nrow(train.sim) ; n1 <- sum(nrow(train.sim), nrow(test.sim))
  res.rate.c <- NULL ; res.lam.c <- NULL ; res.index.c <- NULL
  
  for (j in 1:20) {
    
    n <- nrow(train.sim)
    index.Kc <- sample(1:n, round(0.7*n), replace=F)
    Kc.train <- K.train[index.Kc,index.Kc] ; Kc.test <- K.train[index.Kc,-index.Kc]
    yc.train <- y.train[index.Kc] ; yc.test <- y.train[-index.Kc]
    
    grid.l <- 10^seq(-6,1,length=100)
    
    hc <- cv.kernel(yc.train, Kc.train, 5, grid.l) # 5-fold crossvalidation 
    
    se.logL.c <- apply(hc$cv.logL, 1, sd)
    mean.logL.c <- rowMeans(hc$cv.logL)
    idx.c <- which.min(mean.logL.c)
    best.lam.c <- max(hc$lambda[ mean.logL.c < mean.logL.c[idx.c]+se.logL.c[idx.c] ])
    
    h1.c <- fit.kernel(yc.train, Kc.train, best.lam.c) # Fitting
    
    sim.d.hat.c <- h1.c$d.hat
    
    h2.c <- pred.kernel(yc.test, Kc.test, sim.d.hat.c) # Testing
    rate.miss.c <- h2.c$misclassification.rate
    
    res.rate.c <- c(res.rate.c, rate.miss.c)
    res.lam.c <- c(res.lam.c, best.lam.c)
    res.index.c <- cbind(res.index.c, index.Kc)
    
  }
  
  # Fitting  
  
  res.lambda <- res.lam.c[which.min(res.rate.c)]
  res.index <- res.index.c[, which.min(res.rate.c)]
  res.K.train <- K[res.index, res.index]
  res.y.train <- y.train[res.index]
  K.test <- K[res.index, (n+1):n1]
  
  h1 <- fit.kernel(res.y.train, res.K.train, res.lambda)
  
  res.d.hat <- h1$d.hat
  
  # Calculate test misclassification rate
  
  h2 <- pred.kernel(y.test, K.test, res.d.hat)
  KRLCS[i] <- h2$misclassification.rate
  
}

KRLCS
boxplot(KRLCS)

dat3.2 <- data.frame(test.MR=KRLCS, method=rep("KRLCS",50), number=rep("200",50))








### 3. Kernel ridge logistic regression classification using Bagging with logit estimate (KRLCB1)



KRLCB1 <- c(rep(0,50))

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(200, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  test.sim.y <- test.sim[,1]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train 
  y.train <- u$y.train ; y.test <- u$y.test 
  K <- u$K
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    boot.index <- sample(1:n, n, replace=T)
    boot.K.train <- K.train[boot.index, boot.index]
    boot.K.test <- K.train[boot.index, -boot.index]
    boot.y.train <- y.train[boot.index]
    boot.y.test <- y.train[-boot.index]
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.logit <- NULL
  
  for (r in 1:100) {
    
    boots.index <- sample(1:n, n, replace=T)
    boots.K.train <- K[boots.index, boots.index]
    boots.K.test <- K[boots.index, (n+1):n1]
    boots.y.train <- y.train[boots.index]
    boots.y.test <- y.test 
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.logit.hat <- h2.boots$logit.hat
    
    boots.logit <- cbind(boots.logit, boots.logit.hat)
    
  }
  
  # Bagging estimate
  
  logit.hat.bag <- rowMeans(boots.logit)
  pi.hat.bag <- exp(logit.hat.bag)/(1+exp(logit.hat.bag))
  
  y.hat.bag <- c(rep(1,n1-n))
  y.hat.bag[pi.hat.bag < 0.5] <- 0
  y.hat.bag[pi.hat.bag == 0.5] <- rbinom(length(y.hat.bag[pi.hat.bag == 0.5]), 1, 0.5)
  
  KRLCB1[s] <- mean(y.hat.bag != y.test)
  
}

KRLCB1
boxplot(KRLCB1)

dat3.3 <- data.frame(test.MR=KRLCB1, method=rep("KRLCB1",50), number=rep("200",50))









### 4. Kernel ridge logistic regression classification using Bagging with probability estimate (KRLCB2)



KRLCB2 <- c(rep(0,50))

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(200, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  test.sim.y <- test.sim[,1]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train 
  y.train <- u$y.train ; y.test <- u$y.test 
  K <- u$K
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    boot.index <- sample(1:n, n, replace=T)
    boot.K.train <- K.train[boot.index, boot.index]
    boot.K.test <- K.train[boot.index, -boot.index]
    boot.y.train <- y.train[boot.index]
    boot.y.test <- y.train[-boot.index]
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.pi <- NULL
  
  for (r in 1:100) {
    
    boots.index <- sample(1:n, n, replace=T)
    boots.K.train <- K[boots.index, boots.index]
    boots.K.test <- K[boots.index, (n+1):n1]
    boots.y.train <- y.train[boots.index]
    boots.y.test <- y.test 
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.pi.hat <- h2.boots$pi.hat.new
    
    boots.pi <- cbind(boots.pi, boots.pi.hat)
    
  }
  
  # Bagging estimate
  
  pi.hat.bag <- rowMeans(boots.pi)
  
  y.hat.bag <- c(rep(1,n1-n))
  y.hat.bag[pi.hat.bag < 0.5] <- 0
  y.hat.bag[pi.hat.bag == 0.5] <- rbinom(length(y.hat.bag[pi.hat.bag == 0.5]), 1, 0.5)
  
  KRLCB2[s] <- mean(y.hat.bag != y.test)
  
}

KRLCB2
boxplot(KRLCB2)

dat3.4 <- data.frame(test.MR=KRLCB2, method=rep("KRLCB2",50), number=rep("200",50))








### 5. Kernel ridge logistic regression classification using Bagging with majority vote (KRLCB3)



KRLCB3 <- c(rep(0,50))

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(200, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  test.sim.y <- test.sim[,1]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train 
  y.train <- u$y.train ; y.test <- u$y.test 
  K <- u$K
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    boot.index <- sample(1:n, n, replace=T)
    boot.K.train <- K.train[boot.index, boot.index]
    boot.K.test <- K.train[boot.index, -boot.index]
    boot.y.train <- y.train[boot.index]
    boot.y.test <- y.train[-boot.index]
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.predict.new <- NULL
  
  for (r in 1:100) {
    
    boots.index <- sample(1:n, n, replace=T)
    boots.K.train <- K[boots.index, boots.index]
    boots.K.test <- K[boots.index, (n+1):n1]
    boots.y.train <- y.train[boots.index]
    boots.y.test <- y.test 
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.predict <- h2.boots$predict.new
    
    boots.predict.new <- rbind(boots.predict.new, boots.predict)
    
  }
  
  # Bagging estimate
  
  mean.hat.bag <- apply(boots.predict.new, 2, mean)
  y.hat.bag <- c(rep(1,n1-n))
  y.hat.bag[mean.hat.bag < 0.5] <- 0
  y.hat.bag[mean.hat.bag == 0.5] <- rbinom(length(y.hat.bag[mean.hat.bag == 0.5]), 1, 0.5)
  
  KRLCB3[s] <- mean(y.hat.bag != y.test)
  
}

KRLCB3
boxplot(KRLCB3)

dat3.5 <- data.frame(test.MR=KRLCB3, method=rep("KRLCB3",50), number=rep("200",50))








### 6. Kernel ridge logistic regression classification using Random Forest with logit estimate (KRLCR1)



KRLCR1 <- c(rep(0,50))  

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(200, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim.y <- train.sim[,1] ; test.sim.y <- test.sim[,1]
  train.sim.X <- train.sim[,2:21] ; test.sim.X <- test.sim[,2:21]
  
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  p <- ncol(train.sim.X)
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    rf.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boot.index <- sample(1:n, n, replace=T)
    train.sim1 <- cbind(train.sim.y[boot.index], train.sim.X[boot.index, rf.index])
    test.sim1 <- cbind(train.sim.y[-boot.index], train.sim.X[-boot.index, rf.index])
    
    u1 <- my.kernel.matrix(train.sim1, test.sim1)
    boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
    boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.logit <- NULL
  
  for (r in 1:100) {
    
    rfs.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boots.index <- sample(1:n, n, replace=T)
    train.sim2 <- cbind(train.sim.y[boots.index], train.sim.X[boots.index, rfs.index])
    test.sim2 <- cbind(test.sim.y, test.sim.X[, rfs.index])
    
    u2 <- my.kernel.matrix(train.sim2, test.sim2)
    boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
    boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.logit.hat <- h2.boots$logit.hat
    
    boots.logit <- cbind(boots.logit, boots.logit.hat)
    
  }
  
  # Random Forest estimate
  
  logit.hat.rf <- rowMeans(boots.logit)
  pi.hat.rf <- exp(logit.hat.rf)/(1+exp(logit.hat.rf))
  
  y.hat.rf <- c(rep(1,n1-n))
  y.hat.rf[pi.hat.rf < 0.5] <- 0
  y.hat.rf[pi.hat.rf == 0.5] <- rbinom(length(y.hat.rf[pi.hat.rf == 0.5]), 1, 0.5)
  
  KRLCR1[s] <- mean(y.hat.rf != test.sim.y)
  
}

KRLCR1
boxplot(KRLCR1)

dat3.6 <- data.frame(test.MR=KRLCR1, method=rep("KRLCR1",50), number=rep("200",50))







### 7. Kernel ridge logistic regression classification using Random Forest with probability estimate (KRLCR2)



KRLCR2 <- c(rep(0,50))  

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(200, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim.y <- train.sim[,1] ; test.sim.y <- test.sim[,1]
  train.sim.X <- train.sim[,2:21] ; test.sim.X <- test.sim[,2:21]
  
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  p <- ncol(train.sim.X)
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    rf.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boot.index <- sample(1:n, n, replace=T)
    train.sim1 <- cbind(train.sim.y[boot.index], train.sim.X[boot.index, rf.index])
    test.sim1 <- cbind(train.sim.y[-boot.index], train.sim.X[-boot.index, rf.index])
    
    u1 <- my.kernel.matrix(train.sim1, test.sim1)
    boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
    boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.pi <- NULL
  
  for (r in 1:100) {
    
    rfs.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boots.index <- sample(1:n, n, replace=T)
    train.sim2 <- cbind(train.sim.y[boots.index], train.sim.X[boots.index, rfs.index])
    test.sim2 <- cbind(test.sim.y, test.sim.X[, rfs.index])
    
    u2 <- my.kernel.matrix(train.sim2, test.sim2)
    boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
    boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.pi.hat <- h2.boots$pi.hat.new
    
    boots.pi <- cbind(boots.pi, boots.pi.hat)
    
  }
  
  # Random Forest estimate
  
  pi.hat.rf <- rowMeans(boots.pi)
  
  y.hat.rf <- c(rep(1,n1-n))
  y.hat.rf[pi.hat.rf < 0.5] <- 0
  y.hat.rf[pi.hat.rf == 0.5] <- rbinom(length(y.hat.rf[pi.hat.rf == 0.5]), 1, 0.5)
  
  KRLCR2[s] <- mean(y.hat.rf != test.sim.y)
  
}

KRLCR2
boxplot(KRLCR2)

dat3.7 <- data.frame(test.MR=KRLCR2, method=rep("KRLCR2",50), number=rep("200",50))







### 8. Kernel ridge logistic regression classification using Random Forest with majority vote (KRLCR3)



KRLCR3 <- c(rep(0,50))  

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(200, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim.y <- train.sim[,1] ; test.sim.y <- test.sim[,1]
  train.sim.X <- train.sim[,2:21] ; test.sim.X <- test.sim[,2:21]
  
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  p <- ncol(train.sim.X)
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    rf.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boot.index <- sample(1:n, n, replace=T)
    train.sim1 <- cbind(train.sim.y[boot.index], train.sim.X[boot.index, rf.index])
    test.sim1 <- cbind(train.sim.y[-boot.index], train.sim.X[-boot.index, rf.index])
    
    u1 <- my.kernel.matrix(train.sim1, test.sim1)
    boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
    boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.predict.new <- NULL
  
  for (r in 1:100) {
    
    rfs.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boots.index <- sample(1:n, n, replace=T)
    train.sim2 <- cbind(train.sim.y[boots.index], train.sim.X[boots.index, rfs.index])
    test.sim2 <- cbind(test.sim.y, test.sim.X[, rfs.index])
    
    u2 <- my.kernel.matrix(train.sim2, test.sim2)
    boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
    boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.predict <- h2.boots$predict.new
    
    boots.predict.new <- rbind(boots.predict.new, boots.predict)
    
  }
  
  # Random Forest estimate
  
  mean.hat.rf <- apply(boots.predict.new, 2, mean)
  y.hat.rf <- c(rep(1,n1-n))
  y.hat.rf[mean.hat.rf < 0.5] <- 0
  y.hat.rf[mean.hat.rf == 0.5] <- rbinom(length(y.hat.rf[mean.hat.rf == 0.5]), 1, 0.5)
  
  KRLCR3[s] <- mean(y.hat.rf != test.sim.y)
  
}

KRLCR3
boxplot(KRLCR3)

dat3.8 <- data.frame(test.MR=KRLCR3, method=rep("KRLCR3",50), number=rep("200",50))





dat3 <- rbind(dat3.1, dat3.2, dat3.3, dat3.4, dat3.5, dat3.6, dat3.7, dat3.8)






#### SIMULATION (n=400)





### 1. Kernel ridge logistic regression classification (KRLC)


KRLC <- c(rep(0,50))

for (i in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(400, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  
  # 5-fold crossvalidation
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train ; K.test <- u$K.test  
  y.train <- u$y.train ; y.test <- u$y.test
  
  grid.l <- 10^seq(-6,1,length=100)
  
  h <- cv.kernel(y.train, K.train, 5, grid.l) 
  
  # Fitting 
  
  se.logL <- apply(h$cv.logL, 1, sd)
  mean.logL <- rowMeans(h$cv.logL)
  idx <- which.min(mean.logL)
  best.lam <- max(h$lambda[ mean.logL < mean.logL[idx]+se.logL[idx] ])
  
  h1 <- fit.kernel(y.train, K.train, best.lam)
  
  # Calculate test misclassification rate
  
  sim.d.hat <- h1$d.hat
  
  h2 <- pred.kernel(y.test, K.test, sim.d.hat)
  KRLC[i] <- h2$misclassification.rate
  
}

KRLC
boxplot(KRLC)

dat4.1 <- data.frame(test.MR=KRLC, method=rep("KRLC",50), number=rep("400",50))








### 2. Kernel ridge logistic regression classification using Sub-sampling (KRLCS)


KRLCS <- c(rep(0,50))

for (i in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(400, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train ; y.train <- u$y.train 
  K <- u$K ; y.test <- u$y.test
  
  # Choosing optimal lambda
  
  
  n <- nrow(train.sim) ; n1 <- sum(nrow(train.sim), nrow(test.sim))
  res.rate.c <- NULL ; res.lam.c <- NULL ; res.index.c <- NULL
  
  for (j in 1:20) {
    
    n <- nrow(train.sim)
    index.Kc <- sample(1:n, round(0.7*n), replace=F)
    Kc.train <- K.train[index.Kc,index.Kc] ; Kc.test <- K.train[index.Kc,-index.Kc]
    yc.train <- y.train[index.Kc] ; yc.test <- y.train[-index.Kc]
    
    grid.l <- 10^seq(-6,1,length=100)
    
    hc <- cv.kernel(yc.train, Kc.train, 5, grid.l) # 5-fold crossvalidation 
    
    se.logL.c <- apply(hc$cv.logL, 1, sd)
    mean.logL.c <- rowMeans(hc$cv.logL)
    idx.c <- which.min(mean.logL.c)
    best.lam.c <- max(hc$lambda[ mean.logL.c < mean.logL.c[idx.c]+se.logL.c[idx.c] ])
    
    h1.c <- fit.kernel(yc.train, Kc.train, best.lam.c) # Fitting
    
    sim.d.hat.c <- h1.c$d.hat
    
    h2.c <- pred.kernel(yc.test, Kc.test, sim.d.hat.c) # Testing
    rate.miss.c <- h2.c$misclassification.rate
    
    res.rate.c <- c(res.rate.c, rate.miss.c)
    res.lam.c <- c(res.lam.c, best.lam.c)
    res.index.c <- cbind(res.index.c, index.Kc)
    
  }
  
  # Fitting  
  
  res.lambda <- res.lam.c[which.min(res.rate.c)]
  res.index <- res.index.c[, which.min(res.rate.c)]
  res.K.train <- K[res.index, res.index]
  res.y.train <- y.train[res.index]
  K.test <- K[res.index, (n+1):n1]
  
  h1 <- fit.kernel(res.y.train, res.K.train, res.lambda)
  
  res.d.hat <- h1$d.hat
  
  # Calculate test misclassification rate
  
  h2 <- pred.kernel(y.test, K.test, res.d.hat)
  KRLCS[i] <- h2$misclassification.rate
  
}

KRLCS
boxplot(KRLCS)

dat4.2 <- data.frame(test.MR=KRLCS, method=rep("KRLCS",50), number=rep("400",50))








### 3. Kernel ridge logistic regression classification using Bagging with logit estimate (KRLCB1)



KRLCB1 <- c(rep(0,50))

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(400, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  test.sim.y <- test.sim[,1]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train 
  y.train <- u$y.train ; y.test <- u$y.test 
  K <- u$K
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    boot.index <- sample(1:n, n, replace=T)
    boot.K.train <- K.train[boot.index, boot.index]
    boot.K.test <- K.train[boot.index, -boot.index]
    boot.y.train <- y.train[boot.index]
    boot.y.test <- y.train[-boot.index]
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.logit <- NULL
  
  for (r in 1:100) {
    
    boots.index <- sample(1:n, n, replace=T)
    boots.K.train <- K[boots.index, boots.index]
    boots.K.test <- K[boots.index, (n+1):n1]
    boots.y.train <- y.train[boots.index]
    boots.y.test <- y.test 
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.logit.hat <- h2.boots$logit.hat
    
    boots.logit <- cbind(boots.logit, boots.logit.hat)
    
  }
  
  # Bagging estimate
  
  logit.hat.bag <- rowMeans(boots.logit)
  pi.hat.bag <- exp(logit.hat.bag)/(1+exp(logit.hat.bag))
  
  y.hat.bag <- c(rep(1,n1-n))
  y.hat.bag[pi.hat.bag < 0.5] <- 0
  y.hat.bag[pi.hat.bag == 0.5] <- rbinom(length(y.hat.bag[pi.hat.bag == 0.5]), 1, 0.5)
  
  KRLCB1[s] <- mean(y.hat.bag != y.test)
  
}

KRLCB1
boxplot(KRLCB1)

dat4.3 <- data.frame(test.MR=KRLCB1, method=rep("KRLCB1",50), number=rep("400",50))









### 4. Kernel ridge logistic regression classification using Bagging with probability estimate (KRLCB2)



KRLCB2 <- c(rep(0,50))

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(400, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  test.sim.y <- test.sim[,1]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train 
  y.train <- u$y.train ; y.test <- u$y.test 
  K <- u$K
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    boot.index <- sample(1:n, n, replace=T)
    boot.K.train <- K.train[boot.index, boot.index]
    boot.K.test <- K.train[boot.index, -boot.index]
    boot.y.train <- y.train[boot.index]
    boot.y.test <- y.train[-boot.index]
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.pi <- NULL
  
  for (r in 1:100) {
    
    boots.index <- sample(1:n, n, replace=T)
    boots.K.train <- K[boots.index, boots.index]
    boots.K.test <- K[boots.index, (n+1):n1]
    boots.y.train <- y.train[boots.index]
    boots.y.test <- y.test 
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.pi.hat <- h2.boots$pi.hat.new
    
    boots.pi <- cbind(boots.pi, boots.pi.hat)
    
  }
  
  # Bagging estimate
  
  pi.hat.bag <- rowMeans(boots.pi)
  
  y.hat.bag <- c(rep(1,n1-n))
  y.hat.bag[pi.hat.bag < 0.5] <- 0
  y.hat.bag[pi.hat.bag == 0.5] <- rbinom(length(y.hat.bag[pi.hat.bag == 0.5]), 1, 0.5)
  
  KRLCB2[s] <- mean(y.hat.bag != y.test)
  
}

KRLCB2
boxplot(KRLCB2)

dat4.4 <- data.frame(test.MR=KRLCB2, method=rep("KRLCB2",50), number=rep("400",50))








### 5. Kernel ridge logistic regression classification using Bagging with majority vote (KRLCB3)



KRLCB3 <- c(rep(0,50))

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(400, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim <- train.sim[,-22] ; test.sim <- test.sim[,-22]
  test.sim.y <- test.sim[,1]
  
  u <- my.kernel.matrix(train.sim, test.sim)
  K.train <- u$K.train 
  y.train <- u$y.train ; y.test <- u$y.test 
  K <- u$K
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    boot.index <- sample(1:n, n, replace=T)
    boot.K.train <- K.train[boot.index, boot.index]
    boot.K.test <- K.train[boot.index, -boot.index]
    boot.y.train <- y.train[boot.index]
    boot.y.test <- y.train[-boot.index]
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.predict.new <- NULL
  
  for (r in 1:100) {
    
    boots.index <- sample(1:n, n, replace=T)
    boots.K.train <- K[boots.index, boots.index]
    boots.K.test <- K[boots.index, (n+1):n1]
    boots.y.train <- y.train[boots.index]
    boots.y.test <- y.test 
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.predict <- h2.boots$predict.new
    
    boots.predict.new <- rbind(boots.predict.new, boots.predict)
    
  }
  
  # Bagging estimate
  
  mean.hat.bag <- apply(boots.predict.new, 2, mean)
  y.hat.bag <- c(rep(1,n1-n))
  y.hat.bag[mean.hat.bag < 0.5] <- 0
  y.hat.bag[mean.hat.bag == 0.5] <- rbinom(length(y.hat.bag[mean.hat.bag == 0.5]), 1, 0.5)
  
  KRLCB3[s] <- mean(y.hat.bag != y.test)
  
}

KRLCB3
boxplot(KRLCB3)

dat4.5 <- data.frame(test.MR=KRLCB3, method=rep("KRLCB3",50), number=rep("400",50))








### 6. Kernel ridge logistic regression classification using Random Forest with logit estimate (KRLCR1)



KRLCR1 <- c(rep(0,50))  

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(400, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim.y <- train.sim[,1] ; test.sim.y <- test.sim[,1]
  train.sim.X <- train.sim[,2:21] ; test.sim.X <- test.sim[,2:21]
  
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  p <- ncol(train.sim.X)
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    rf.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boot.index <- sample(1:n, n, replace=T)
    train.sim1 <- cbind(train.sim.y[boot.index], train.sim.X[boot.index, rf.index])
    test.sim1 <- cbind(train.sim.y[-boot.index], train.sim.X[-boot.index, rf.index])
    
    u1 <- my.kernel.matrix(train.sim1, test.sim1)
    boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
    boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.logit <- NULL
  
  for (r in 1:100) {
    
    rfs.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boots.index <- sample(1:n, n, replace=T)
    train.sim2 <- cbind(train.sim.y[boots.index], train.sim.X[boots.index, rfs.index])
    test.sim2 <- cbind(test.sim.y, test.sim.X[, rfs.index])
    
    u2 <- my.kernel.matrix(train.sim2, test.sim2)
    boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
    boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.logit.hat <- h2.boots$logit.hat
    
    boots.logit <- cbind(boots.logit, boots.logit.hat)
    
  }
  
  # Random Forest estimate
  
  logit.hat.rf <- rowMeans(boots.logit)
  pi.hat.rf <- exp(logit.hat.rf)/(1+exp(logit.hat.rf))
  
  y.hat.rf <- c(rep(1,n1-n))
  y.hat.rf[pi.hat.rf < 0.5] <- 0
  y.hat.rf[pi.hat.rf == 0.5] <- rbinom(length(y.hat.rf[pi.hat.rf == 0.5]), 1, 0.5)
  KRLCR1[s] <- mean(y.hat.rf != test.sim.y)
  
}

KRLCR1
boxplot(KRLCR1)

dat4.6 <- data.frame(test.MR=KRLCR1, method=rep("KRLCR1",50), number=rep("400",50))







### 7. Kernel ridge logistic regression classification using Random Forest with probability estimate (KRLCR2)



KRLCR2 <- c(rep(0,50))  

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(400, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim.y <- train.sim[,1] ; test.sim.y <- test.sim[,1]
  train.sim.X <- train.sim[,2:21] ; test.sim.X <- test.sim[,2:21]
  
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  p <- ncol(train.sim.X)
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    rf.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boot.index <- sample(1:n, n, replace=T)
    train.sim1 <- cbind(train.sim.y[boot.index], train.sim.X[boot.index, rf.index])
    test.sim1 <- cbind(train.sim.y[-boot.index], train.sim.X[-boot.index, rf.index])
    
    u1 <- my.kernel.matrix(train.sim1, test.sim1)
    boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
    boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.pi <- NULL
  
  for (r in 1:100) {
    
    rfs.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boots.index <- sample(1:n, n, replace=T)
    train.sim2 <- cbind(train.sim.y[boots.index], train.sim.X[boots.index, rfs.index])
    test.sim2 <- cbind(test.sim.y, test.sim.X[, rfs.index])
    
    u2 <- my.kernel.matrix(train.sim2, test.sim2)
    boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
    boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.pi.hat <- h2.boots$pi.hat.new
    
    boots.pi <- cbind(boots.pi, boots.pi.hat)
    
  }
  
  # Random Forest estimate
  
  pi.hat.rf <- rowMeans(boots.pi)
  
  y.hat.rf <- c(rep(1,n1-n))
  y.hat.rf[pi.hat.rf < 0.5] <- 0
  y.hat.rf[pi.hat.rf == 0.5] <- rbinom(length(y.hat.rf[pi.hat.rf == 0.5]), 1, 0.5)
  
  KRLCR2[s] <- mean(y.hat.rf != test.sim.y)
  
}

KRLCR2
boxplot(KRLCR2)

dat4.7 <- data.frame(test.MR=KRLCR2, method=rep("KRLCR2",50), number=rep("400",50))







### 8. Kernel ridge logistic regression classification using Random Forest with majority vote (KRLCR3)



KRLCR3 <- c(rep(0,50))  

for (s in 1:50) {
  
  # Making simulation data
  
  beta <- 0.5^(0:20)
  train.sim <- dat.gen(400, beta)
  test.sim <- dat.gen(1000, beta)
  
  train.sim.y <- train.sim[,1] ; test.sim.y <- test.sim[,1]
  train.sim.X <- train.sim[,2:21] ; test.sim.X <- test.sim[,2:21]
  
  n <- nrow(train.sim)
  n1 <- sum(nrow(train.sim), nrow(test.sim))
  p <- ncol(train.sim.X)
  
  # Choosing best lambda
  
  res.rate.miss <- NULL
  grid.l <- 10^seq(-6,1,length=100)
  
  for (i in 1:50) {
    
    rf.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boot.index <- sample(1:n, n, replace=T)
    train.sim1 <- cbind(train.sim.y[boot.index], train.sim.X[boot.index, rf.index])
    test.sim1 <- cbind(train.sim.y[-boot.index], train.sim.X[-boot.index, rf.index])
    
    u1 <- my.kernel.matrix(train.sim1, test.sim1)
    boot.K.train <- u1$K.train ; boot.K.test <- u1$K.test
    boot.y.train <- u1$y.train ; boot.y.test <- u1$y.test
    
    boot.rate.miss <- c(rep(0,100))
    
    for (j in 1:100) {
      
      h1.boot <- fit.kernel(boot.y.train, boot.K.train, grid.l[j])
      boot.d.hat <- h1.boot$d.hat
      h2.boot <- pred.kernel(boot.y.test, boot.K.test, boot.d.hat)
      boot.rate.miss[j] <- h2.boot$misclassification.rate
      
    }
    
    res.rate.miss <- rbind(res.rate.miss, boot.rate.miss)
    
  }
  
  #best.lambda <- grid.l[which.min(colMeans(res.rate.miss))]
  best.lambda <- max(grid.l[colMeans(res.rate.miss)==min(colMeans(res.rate.miss))])
  
  # Bootstraping
  
  boots.predict.new <- NULL
  
  for (r in 1:100) {
    
    rfs.index <- sample(1:p, ceiling(sqrt(p)), replace=F)
    boots.index <- sample(1:n, n, replace=T)
    train.sim2 <- cbind(train.sim.y[boots.index], train.sim.X[boots.index, rfs.index])
    test.sim2 <- cbind(test.sim.y, test.sim.X[, rfs.index])
    
    u2 <- my.kernel.matrix(train.sim2, test.sim2)
    boots.K.train <- u2$K.train ; boots.K.test <- u2$K.test
    boots.y.train <- u2$y.train ; boots.y.test <- u2$y.test
    
    h1.boots <- fit.kernel(boots.y.train, boots.K.train, best.lambda)
    boots.d.hat <- h1.boots$d.hat
    h2.boots <- pred.kernel(boots.y.test, boots.K.test, boots.d.hat)
    boots.predict <- h2.boots$predict.new
    
    boots.predict.new <- rbind(boots.predict.new, boots.predict)
    
  }
  
  # Random Forest estimate
  
  mean.hat.rf <- apply(boots.predict.new, 2, mean)
  y.hat.rf <- c(rep(1,n1-n))
  y.hat.rf[mean.hat.rf < 0.5] <- 0
  y.hat.rf[mean.hat.rf == 0.5] <- rbinom(length(y.hat.rf[mean.hat.rf == 0.5]), 1, 0.5)
  
  KRLCR3[s] <- mean(y.hat.rf != test.sim.y)
  
}

KRLCR3
boxplot(KRLCR3)

dat4.8 <- data.frame(test.MR=KRLCR3, method=rep("KRLCR3",50), number=rep("400",50))





dat4 <- rbind(dat4.1, dat4.2, dat4.3, dat4.4, dat4.5, dat4.6, dat4.7, dat4.8)





### Final result

final.dat <- rbind(dat1, dat2, dat3, dat4)

install.packages("ggplot2")
library(ggplot2)

ggplot(final.dat, aes(x = number, y = test.MR, fill = number)) + geom_boxplot() + facet_wrap(~ method, ncol=8) 