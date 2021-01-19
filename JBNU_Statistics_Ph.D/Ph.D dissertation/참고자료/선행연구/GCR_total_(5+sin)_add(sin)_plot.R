library(survival)
library(KernSmooth)
library(np)
library(locfit)
rm(list=ls())
gc()
gc(reset = T)

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

set.seed(2018)

n <- c(100, 200, 400)

M <- 400

new.x <- data.frame(x.i = seq(0, 1, by = 0.01))

new.x.no <- data.frame(x.cen.no = seq(0, 1, by = 0.01))

new.x.10 <- data.frame(x.cen.10 = seq(0, 1, by = 0.01))

new.x.30 <- data.frame(x.cen.30 = seq(0, 1, by = 0.01))

new.x.50 <- data.frame(x.cen.50 = seq(0, 1, by = 0.01))

###########################################################################################
result.no.sin <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))
result.10.sin <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))
result.30.sin <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))
result.50.sin <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))

plot.no.sin <- matrix(0, nrow = length(n), ncol = nrow(new.x))
plot.10.sin <- matrix(0, nrow = length(n), ncol = nrow(new.x))
plot.30.sin <- matrix(0, nrow = length(n), ncol = nrow(new.x))
plot.50.sin <- matrix(0, nrow = length(n), ncol = nrow(new.x))


system.time(
  
  for (i in 1:length(n)){
    
    mx <- 5 + sin(2*pi*seq(0, 1, by = 0.01))
    pred.value <- 5 + sin(2*pi*seq(0, 1, by = 0.01))
    
    # 0% Censoring
    
    mse.no.sin <- numeric(M)
    var.no.sin <- numeric(M)
    bias.no.sin <- numeric(M)
    predict.no.sin <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.0.sin <- matrix(0, nrow = 400, ncol = 1)
    
    # 10% Censoring
    
    mse.10.sin <- numeric(M)
    var.10.sin <- numeric(M)
    bias.10.sin <- numeric(M)
    predict.10.sin <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.1.sin <- matrix(0, nrow = 400, ncol = 1)
    
    # 30% Censoring
    
    mse.30.sin <- numeric(M)
    var.30.sin <- numeric(M)
    bias.30.sin <- numeric(M)
    predict.30.sin <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.3.sin <- matrix(0, nrow = 400, ncol = 1)
    
    # 50% Censoring
    mse.50.sin <- numeric(M)
    var.50.sin <- numeric(M)
    bias.50.sin <- numeric(M)
    predict.50.sin <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.5.sin <- matrix(0, nrow = 400, ncol = 1)
    
    
    for (j in 1:M){
      
      x.i <- runif(n[i], min = 0, max = 1)
      
      e.i <- rnorm(n[i], mean = 0, sd = 1)
      
      y.i <- 5 + sin(2*pi*x.i) + e.i
      
      true.value <- 5 + sin(2*pi*x.i)
      
      
      # NO Censoring
      cen.no <- rnorm(n[i], 0, 1) + 100
      pre.T.no <- pmin(y.i, cen.no) 
      delta.no <- (y.i <= cen.no) * 1
      g.no <- km.surv(pre.T.no, delta.no)
      y.i.no <- pre.T.no * delta.no / g.no
      
      
      T.no.sin <- pre.T.no - true.value
      g.no.sin <- km.surv(T.no.sin, delta.no)
      y.i.no.sin <- T.no.sin * delta.no / g.no.sin
      
      bw.0.sin <- dpill(x.i, y.i.no.sin)
      model.np.0.sin <- npreg(y.i.no.sin ~ x.i, bws = bw.0.sin, regtype = "ll")
      
      predict.no.sin[j, ] <- predict(model.np.0.sin, newdata = new.x)
      mse.0.sin[j, ] <- mean(((predict.no.sin[j, ] + pred.value) - mx )^2)/M 
      
      
      
      
      # 10% Censoring
      
      cen.10 <- rnorm(n[i], 5, 1) + 2
      pre.T.10 <- pmin(y.i, cen.10)
      delta.10 <- (y.i <= cen.10) * 1
      g.10 <- km.surv(pre.T.10, delta.10)
      y.i.10 <- pre.T.10 * delta.10 / g.10
      
      
      T.10.sin <- pre.T.10 - true.value
      g.10.sin <- km.surv(T.10.sin, delta.10)
      y.i.10.sin <- T.10.sin * delta.10 / g.10.sin
      
      bw.1.sin <- dpill(x.i, y.i.10.sin)
      model.np.10.sin <- npreg(y.i.10.sin ~ x.i, bws = bw.1.sin, regtype = "ll")
      
      predict.10.sin[j, ] <- predict(model.np.10.sin, newdata = new.x)
      mse.1.sin[j, ] <- mean(((predict.10.sin[j, ] + pred.value)- mx )^2)/M 
      
      # 30% Censoring
      
      cen.30 <- rnorm(n[i], 5, 1) + 0.85
      pre.T.30 <- pmin(y.i, cen.30)
      delta.30 <- (y.i <= cen.30) * 1
      g.30 <- km.surv(pre.T.30, delta.30)
      y.i.30 <- pre.T.30 * delta.30 / g.30
      
      
      T.30.sin <- pre.T.30 - true.value
      g.30.sin <- km.surv(T.30.sin, delta.30)
      y.i.30.sin <- T.30.sin * delta.30 / g.30.sin
      
      bw.3.sin <- dpill(x.i, y.i.30.sin)
      model.np.30.sin <- npreg(y.i.30.sin ~ x.i, bw.3.sin, regtype = "ll")
      
      predict.30.sin[j, ] <- predict(model.np.30.sin, newdata = new.x)
      mse.3.sin[j, ] <- mean(((predict.30.sin[j, ] + pred.value)- mx )^2)/M 
      
      
      # 50% Censoring
      
      cen.50 <- rnorm(n[i], 5, 1) + 0
      pre.T.50 <- pmin(y.i, cen.50)
      delta.50 <- (y.i <= cen.50) * 1
      g.50 <- km.surv(pre.T.50, delta.50)
      y.i.50 <- pre.T.50 * delta.50 / g.50
      
      
      
      T.50.sin <- pre.T.50 - true.value
      g.50.sin <- km.surv(T.50.sin, delta.50)
      y.i.50.sin <- T.50.sin * delta.50 / g.50.sin
      
      bw.5.sin <- dpill(x.i, y.i.50.sin)
      model.np.50.sin <- npreg(y.i.50.sin ~ x.i, bw.5.sin, regtype = "ll")
      
      predict.50.sin[j, ] <- predict(model.np.50.sin, newdata = new.x)
      mse.5.sin[j, ] <- mean(((predict.50.sin[j, ] + pred.value) - mx )^2)/M 
      
      
    }
    
    # No Censoring
    
    
    mse.no.sin <- sum(mse.0.sin)
    var.no.sin <- sum(apply((predict.no.sin), 2, var)/length(mx))
    bias.no.sin <- mse.no.sin - var.no.sin
    plot.no.sin[i, ] <- apply((predict.no.sin), 2, mean) + pred.value
    
    result.no.sin[i, ] <- c(bias.no.sin, var.no.sin, mse.no.sin)
    rownames(result.no.sin) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.no.sin) <- c("Bias^2", "Variance", "Mse")
    
    # 10% Censoring
    
    
    mse.10.sin <- sum(mse.1.sin)
    var.10.sin <- sum(apply((predict.10.sin), 2, var)/length(mx))
    bias.10.sin <- mse.10.sin - var.10.sin
    plot.10.sin[i, ] <- apply((predict.10.sin), 2, mean) + pred.value
    
    result.10.sin[i, ] <- c(bias.10.sin, var.10.sin, mse.10.sin)
    rownames(result.10.sin) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.10.sin) <- c("Bias^2", "Variance", "Mse")
    
    # 30% Censoring
    
    
    
    mse.30.sin <- sum(mse.3.sin)
    var.30.sin <- sum(apply((predict.30.sin ), 2, var)/length(mx))
    bias.30.sin <- mse.30.sin - var.30.sin
    plot.30.sin[i, ] <- apply((predict.30.sin), 2, mean) + pred.value
    
    result.30.sin[i, ] <- c(bias.30.sin, var.30.sin, mse.30.sin)
    rownames(result.30.sin) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.30.sin) <- c("Bias^2", "Variance", "Mse")
    
    # 50% Censoring
    
    
    
    mse.50.sin <- sum(mse.5.sin)
    var.50.sin <- sum(apply((predict.50.sin), 2, var)/length(mx))
    bias.50.sin <- mse.50.sin - var.50.sin
    plot.50.sin[i, ] <- apply((predict.50.sin), 2, mean) +  pred.value
    
    result.50.sin[i, ] <- c(bias.50.sin, var.50.sin, mse.50.sin)
    rownames(result.50.sin) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.50.sin) <- c("Bias^2", "Variance", "Mse")
    
    
    
  }
)
###########################################################################################

# result


result.no.sin





result.10.sin



result.30.sin



result.50.sin

##########################################################################################

