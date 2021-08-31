library(survival)
library(KernSmooth)
library(np)
library(locfit)
library(latex2exp)

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


set.seed(2021)

n <- 100

M <- 100

new.x <- data.frame(x.i = seq(0, 1, length.out=n))

new.x.10 <- data.frame(x.cen.10 = seq(0, 1, length.out=n))

result.10 <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))

result.10.constant.3 <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))

result.10.cen <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))

result.10.median <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))

result.10.sin <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))


plot.10 <- matrix(0, nrow = length(n), ncol = nrow(new.x))

plot.10.constant.3 <- matrix(0, nrow = length(n), ncol = nrow(new.x))

plot.10.cen <- matrix(0, nrow = length(n), ncol = nrow(new.x))

plot.10.median <- matrix(0, nrow = length(n), ncol = nrow(new.x))

plot.10.sin <- matrix(0, nrow = length(n), ncol = nrow(new.x))



for (i in 1:length(n)){
  
  mx <- sin(2*pi*seq(0, 1, length.out=n))
  
  # 10% Censoring
  mse.10 <- numeric(M)
  var.10 <- numeric(M)
  bias.10 <- numeric(M)
  predict.10 <- matrix(0, nrow = M, ncol = nrow(new.x))
  mse.1 <- matrix(0, nrow = M, ncol = 1)
  
  mse.10.constant.3 <- numeric(M)
  var.10.constant.3 <- numeric(M)
  bias.10.constant.3 <- numeric(M)
  predict.1.cen <- matrix(0, nrow = M, ncol = n[i])
  predict.10.constant.3 <- matrix(0, nrow = M, ncol = nrow(new.x))
  constant.3.predict.10 <- matrix(0, nrow = M, ncol = nrow(new.x))
  mse.1.constant.3 <- matrix(0, nrow = M, ncol = 1)
  
  mse.10.cen <- numeric(M)
  var.10.cen <- numeric(M)
  bias.10.cen <- numeric(M)
  predict.1.cen <- matrix(0, nrow = M, ncol = n[i])
  predict.10.cen <- matrix(0, nrow = M, ncol = nrow(new.x))
  cen.predict.10 <- matrix(0, nrow = M, ncol = nrow(new.x))
  mse.1.cen <- matrix(0, nrow = M, ncol = 1)
  
  
  mse.10.median <- numeric(M)
  var.10.median <- numeric(M)
  bias.10.median <- numeric(M)
  predict.10.median <- matrix(0, nrow = M, ncol = nrow(new.x))
  mse.1.median <- matrix(0, nrow = M, ncol = 1)
  
  mse.10.sin <- numeric(M)
  var.10.sin <- numeric(M)
  bias.10.sin <- numeric(M)
  predict.10.sin <- matrix(0, nrow = M, ncol = nrow(new.x))
  mse.1.sin <- matrix(0, nrow = M, ncol = 1)
  
  
  for (j in 1:M){
    
    x.i <- seq(0, 1, length.out=n)
    
    e.i <- rnorm(n[i], mean = 0, sd = 1)
    
    y.i <- sin(2*pi*x.i)+e.i
    pred.value <-  sin(2*pi*seq(0,1,length.out=n))
    
    rh <- 0.1
    
    
    # 10% Censoring
    
    cen.10 <-rnorm(n,mean = sin(2*pi*x.i)+0.84^2 , sd=1) # censoring values
    pre.T.10 <- pmin(y.i, cen.10)
    delta.10 <- (y.i <= cen.10) * 1
    
    
    g.10<-km.surv(pre.T.10, delta.10)
    
    y.i.10 <- ifelse(pre.T.10<=quantile(pre.T.10,probs=0.98),pre.T.10 * delta.10 / (g.10)  , 0 )
    
    
    model.np.10 <- npreg(y.i.10 ~ x.i, bws = rh , regtype = "ll")
    
    predict.10[j, ] <- predict(model.np.10, newdata = new.x)
    mse.1[j, ] <- mean((predict.10[j, ] - mx )^2)/M   
    
    cubic.10.3 <- lm(y.i.10 ~ x.i + I(x.i^2)+ I(x.i^3))
    
    y.10.3 <- y.i.10 - predict(cubic.10.3)
    
    model.constant.10.3 <- npreg(y.10.3 ~ x.i, bws = rh , regtype = "ll")
    
    predict.10.constant.3[j, ] <- predict(cubic.10.3, newdata = new.x)
    constant.3.predict.10[j, ] <- predict(model.constant.10.3, newdata = new.x)
    
    mse.1.constant.3[j, ] <- mean(((predict.10.constant.3[j, ] + constant.3.predict.10[j, ]) - mx )^2)/M
    
    y.cen.10 <- y.i[delta.10 == 1]
    x.cen.10 <- x.i[delta.10 == 1]
    
    model.cen.1 <- locfit(y.cen.10 ~ x.cen.10, alpha = seq(0.2, 0.8, by = 0.05))
    predict.1.cen[j, ] <- predict(model.cen.1, newdata = x.i)
    
    T.10.cen <- pre.T.10 - predict.1.cen[j, ]
    
    g.10.cen<-km.surv(T.10.cen, delta.10)
    
    y.i.10.cen <- ifelse(T.10.cen<=quantile(T.10.cen,probs=0.98), T.10.cen * delta.10 / g.10.cen  ,0)
    
    
    model.np.1.cen <- npreg(y.i.10.cen ~ x.i, bws = rh , regtype = "ll")
    
    predict.10.cen[j, ] <- predict(model.np.1.cen, newdata = new.x)
    cen.predict.10[j, ] <-predict(model.cen.1, newdata = new.x.10)
    mse.1.cen[j, ] <- mean(((predict.10.cen[j, ] + cen.predict.10[j, ])- mx )^2)/M
    
    
    fit.10 <- read.table(textConnection(capture.output(survfit(Surv(pre.T.10, delta.10) ~ 1))), skip = 2, header = TRUE)
    fit.median.10 <- fit.10$median
    
    T.10.median <- pre.T.10 - fit.median.10
    
    g.10.median<-km.surv(T.10.median, delta.10)
    
    y.i.10.median <- ifelse(T.10.median<=quantile(T.10.median,probs=0.98),T.10.median * delta.10 / g.10.median , 0)
    
    
    model.np.10.median <- npreg(y.i.10.median ~ x.i, bws = rh , regtype = "ll")
    
    predict.10.median[j, ] <- predict(model.np.10.median, newdata = new.x)
    mse.1.median[j, ] <- mean(((predict.10.median[j, ] + fit.median.10)- mx )^2)/M
    
    
    true.value <- sin(2*pi*x.i)
    T.10.sin <- pre.T.10 - true.value
    
    g.10.sin<-km.surv(T.10.sin, delta.10)
    y.i.10.sin <- ifelse(T.10.sin<=quantile(T.10.sin,probs=0.98),T.10.sin * delta.10 / g.10.sin , 0)
    
    
    model.np.10.sin <- npreg(y.i.10.sin ~ x.i, bws = rh , regtype = "ll")
    
    predict.10.sin[j, ] <- predict(model.np.10.sin, newdata = new.x)
    mse.1.sin[j, ] <- mean(((predict.10.sin[j, ] + pred.value)- mx )^2)/M
    
  }
  # 10% Censoring
  mse.10 <- sum(mse.1)
  var.10 <- sum(apply(predict.10, 2, var)/length(mx))
  bias.10 <- mse.10 - var.10
  plot.10[i, ] <- apply(predict.10, 2, mean)
  
  
  result.10[i, ] <- c(bias.10, var.10, mse.10)
  rownames(result.10) <- c("n = 100")
  colnames(result.10) <- c("Bias^2", "Variance", "Mse")
  
  mse.10.constant.3 <- sum(mse.1.constant.3)
  var.10.constant.3 <- sum(apply((predict.10.constant.3 + constant.3.predict.10), 2, var)/length(mx))
  bias.10.constant.3 <- mse.10.constant.3 - var.10.constant.3
  plot.10.constant.3[i, ] <- apply((predict.10.constant.3 + constant.3.predict.10), 2, mean)
  
  result.10.constant.3[i, ] <- c(bias.10.constant.3, var.10.constant.3, mse.10.constant.3)
  rownames(result.10.constant.3) <- c("n = 100")
  colnames(result.10.constant.3) <- c("Bias^2", "Variance", "Mse")
  
  mse.10.cen <- sum(mse.1.cen)
  var.10.cen <- sum(apply((predict.10.cen + cen.predict.10), 2, var)/length(mx))
  bias.10.cen <- mse.10.cen - var.10.cen
  plot.10.cen[i, ] <- apply((predict.10.cen + cen.predict.10), 2, mean)
  
  result.10.cen[i, ] <- c(bias.10.cen, var.10.cen, mse.10.cen)
  rownames(result.10.cen) <- c("n = 100")
  colnames(result.10.cen) <- c("Bias^2", "Variance", "Mse")
  
  mse.10.median <- sum(mse.1.median)
  var.10.median <- sum(apply(predict.10.median, 2, var)/length(mx))
  bias.10.median <- mse.10.median - var.10.median
  plot.10.median[i, ] <- apply((predict.10.median + fit.median.10), 2, mean)
  
  result.10.median[i, ] <- c(bias.10.median, var.10.median, mse.10.median)
  rownames(result.10.median) <- c("n = 100")
  colnames(result.10.median) <- c("Bias^2", "Variance", "Mse")
  
  mse.10.sin <- sum(mse.1.sin)
  var.10.sin <- sum(apply((predict.10.sin), 2, var)/length(mx))
  bias.10.sin <- mse.10.sin - var.10.sin
  plot.10.sin[i, ] <- apply((predict.10.sin), 2, mean) + pred.value
  
  result.10.sin[i, ] <- c(bias.10.sin, var.10.sin, mse.10.sin)
  rownames(result.10.sin) <- c("n = 100")
  colnames(result.10.sin) <- c("Bias^2", "Variance", "Mse")
}

result.10k <- result.10 ; result.10k
result.10.constant.3k <- result.10.constant.3 ;result.10.constant.3k
result.10.cenk <- result.10.cen ; result.10.cenk
result.10.mediank <- result.10.median ; result.10.mediank
result.10.sink <- result.10.sin ; result.10.sink


plot(x.i,mx, ylim=c(-1.5,1.5),type="l", main = "km surv est(N) = 100", xlab = "new.x", ylab = "True value", lwd = 3, col = 1, cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)
abline(v=x.i[delta.i==0],col=7)
lines(x.i,plot.10[1, ], col = 2, lty = 2, lwd = 3)
lines(x.i,plot.10.constant.3[1, ], col = 3, lty = 3, lwd = 3)
lines(x.i,plot.10.cen[1, ], col = 4, lty = 4, lwd = 3)
lines(x.i,plot.10.median[1, ], col = 5, lty = 5, lwd = 3)
lines(x.i,plot.10.sin[1, ], col = 6, lty = 6, lwd = 3)

legend(0.7,1.3, c("True Function", TeX("$GNR_0"), TeX("$GNR_P"),TeX("$GNR_{NP}"), TeX("$GNR_M"), TeX("$GNR_T")), lty = 1:6, 
       col = 1:6,   cex = 1.2, y.intersp = 0.75, lwd = 3, inset = c(-0.15, -0.06))