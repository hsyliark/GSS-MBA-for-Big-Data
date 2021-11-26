uis2 <- read.csv("C:/Users/HSY/Desktop/uis2.csv",
                 sep=",",header=T)
uis2 <-uis2[,-1]
g <- km.surv(uis2$time, uis2$censor)
uis2$y.s <- ifelse(uis2$time <= quantile(uis2$time, probs=0.98), 
                   uis2$time*uis2$censor/g, 0) 
dat.sim <- data.frame(ys1=uis2$y.s, x1=uis2$age, x2=uis2$beck, 
                      x3=uis2$hercoc, x4=uis2$ivhx, x5=uis2$ndrugtx,
                      x6=uis2$race, x7=uis2$treat, x8=uis2$lot, censor=uis2$censor,
                      ys2=uis2$y.s) 








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


ggplot(data=data_final, aes(x=ys, y=yhat, colour=method))+
  geom_point(shape=19, size=2)+
  ggtitle("Scatter Plot by different method (all) (Nonkernel, GKR1, GKRS1, GKRB1, GKRR1)")
ggplot(data=data_final_nocensor, aes(x=ys, y=yhat, colour=method))+
  geom_point(shape=19, size=2)+
  ggtitle("Scatter Plot by different method (no censor) (Nonkernel, GKR1, GKRS1, GKRB1, GKRR1)")
ggplot(data=data_final_censor, aes(x=ys, y=yhat, colour=method))+
  geom_point(shape=19, size=2)+
  ggtitle("Scatter Plot by different method (censor) (Nonkernel, GKR1, GKRS1, GKRB1, GKRR1)")

data1 <- data.frame(ys=y_s, yhat=Nonkernel_hat_y_s, method=rep("a.Nonkernel",length(y_s)), censor=test_censor)
data2 <- data.frame(ys=y_s, yhat=GKR1_hat_y_s, method=rep("b.GKR1",length(y_s)), censor=test_censor)
data3 <- data.frame(ys=y_s, yhat=GKRS1_hat_y_s, method=rep("c.GKRS1 (Sub-sampling)",length(y_s)), censor=test_censor)
data4 <- data.frame(ys=y_s, yhat=GKRB1_hat_y_s, method=rep("d.GKRB1 (Bagging)",length(y_s)), censor=test_censor)
data5 <- data.frame(ys=y_s, yhat=GKRR1_hat_y_s, method=rep("e.GKRR1 (Random Forest)",length(y_s)), censor=test_censor)
data_final <- rbind(data1,data2,data3,data4,data5)
data_final

data_final_nocensor <- data_final[data_final$censor==1,]
data_final_censor <- data_final[data_final$censor==0,]




