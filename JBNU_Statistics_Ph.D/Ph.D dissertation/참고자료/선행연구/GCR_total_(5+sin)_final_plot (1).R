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

result.no <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))
result.10 <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))
result.30 <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))
result.50 <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))

result.no.constant.3 <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))
result.10.constant.3 <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))
result.30.constant.3 <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))
result.50.constant.3 <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))

result.no.cen <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))
result.10.cen <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))
result.30.cen <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))
result.50.cen <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))


result.no.median <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))
result.10.median <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))
result.30.median <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))
result.50.median <- data.frame(Bias = numeric(length(n)), Variance = numeric(length(n)), Mse = numeric(length(n)))


plot.no <- matrix(0, nrow = length(n), ncol = nrow(new.x))
plot.10 <- matrix(0, nrow = length(n), ncol = nrow(new.x))
plot.30 <- matrix(0, nrow = length(n), ncol = nrow(new.x))
plot.50 <- matrix(0, nrow = length(n), ncol = nrow(new.x))

plot.no.constant.3 <- matrix(0, nrow = length(n), ncol = nrow(new.x))
plot.10.constant.3 <- matrix(0, nrow = length(n), ncol = nrow(new.x))
plot.30.constant.3 <- matrix(0, nrow = length(n), ncol = nrow(new.x))
plot.50.constant.3 <- matrix(0, nrow = length(n), ncol = nrow(new.x))


plot.no.cen <- matrix(0, nrow = length(n), ncol = nrow(new.x))
plot.10.cen <- matrix(0, nrow = length(n), ncol = nrow(new.x))
plot.30.cen <- matrix(0, nrow = length(n), ncol = nrow(new.x))
plot.50.cen <- matrix(0, nrow = length(n), ncol = nrow(new.x))

plot.no.median <- matrix(0, nrow = length(n), ncol = nrow(new.x))
plot.10.median <- matrix(0, nrow = length(n), ncol = nrow(new.x))
plot.30.median <- matrix(0, nrow = length(n), ncol = nrow(new.x))
plot.50.median <- matrix(0, nrow = length(n), ncol = nrow(new.x))


system.time(
  
  for (i in 1:length(n)){
    
    mx <- 5 + sin(2*pi*seq(0, 1, by = 0.01))
    
    # 0% Censoring
    mse.no <- numeric(M)
    var.no <- numeric(M)
    bias.no <- numeric(M)
    predict.no <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.0 <- matrix(0, nrow = 400, ncol = 1)
    
    mse.no.constant.3 <- numeric(M)
    var.no.constant.3 <- numeric(M)
    bias.no.constant.3 <- numeric(M)
    predict.no.constant.3 <- matrix(0, nrow = 400, ncol = nrow(new.x))
    constant.3.predict.no <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.0.constant.3 <- matrix(0, nrow = 400, ncol = 1)
    
    mse.no.cen <- numeric(M)
    var.no.cen <- numeric(M)
    bias.no.cen <- numeric(M)
    predict.0.cen <- matrix(0, nrow = 400, ncol = n[i])
    predict.no.cen <- matrix(0, nrow = 400, ncol = nrow(new.x))
    cen.predict.no <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.0.cen <- matrix(0, nrow = 400, ncol = 1)
    
    mse.no.median <- numeric(M)
    var.no.median <- numeric(M)
    bias.no.median <- numeric(M)
    predict.no.median <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.0.median <- matrix(0, nrow = 400, ncol = 1)
    
    
    # 10% Censoring
    mse.10 <- numeric(M)
    var.10 <- numeric(M)
    bias.10 <- numeric(M)
    predict.10 <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.1 <- matrix(0, nrow = 400, ncol = 1)
    
    mse.10.constant.3 <- numeric(M)
    var.10.constant.3 <- numeric(M)
    bias.10.constant.3 <- numeric(M)
    predict.1.cen <- matrix(0, nrow = 400, ncol = n[i])
    predict.10.constant.3 <- matrix(0, nrow = 400, ncol = nrow(new.x))
    constant.3.predict.10 <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.1.constant.3 <- matrix(0, nrow = 400, ncol = 1)
    
    mse.10.cen <- numeric(M)
    var.10.cen <- numeric(M)
    bias.10.cen <- numeric(M)
    predict.1.cen <- matrix(0, nrow = 400, ncol = n[i])
    predict.10.cen <- matrix(0, nrow = 400, ncol = nrow(new.x))
    cen.predict.10 <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.1.cen <- matrix(0, nrow = 400, ncol = 1)
    
    
    mse.10.median <- numeric(M)
    var.10.median <- numeric(M)
    bias.10.median <- numeric(M)
    predict.10.median <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.1.median <- matrix(0, nrow = 400, ncol = 1)
    
    # 30% Censoring
    mse.30 <- numeric(M)
    var.30 <- numeric(M)
    bias.30 <- numeric(M)
    predict.30 <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.3 <- matrix(0, nrow = 400, ncol = 1)
    
    mse.30.constant.3 <- numeric(M)
    var.30.constant.3 <- numeric(M)
    bias.30.constant.3 <- numeric(M)
    predict.3.cen <- matrix(0, nrow = 400, ncol = n[i])
    predict.30.constant.3 <- matrix(0, nrow = 400, ncol = nrow(new.x))
    constant.3.predict.30 <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.3.constant.3 <- matrix(0, nrow = 400, ncol = 1)
    
    mse.30.cen <- numeric(M)
    var.30.cen <- numeric(M)
    bias.30.cen <- numeric(M)
    predict.3.cen <- matrix(0, nrow = 400, ncol = n[i])
    predict.30.cen <- matrix(0, nrow = 400, ncol = nrow(new.x))
    cen.predict.30 <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.3.cen <- matrix(0, nrow = 400, ncol = 1)
    
    mse.30.median <- numeric(M)
    var.30.median <- numeric(M)
    bias.30.median <- numeric(M)
    predict.30.median <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.3.median <- matrix(0, nrow = 400, ncol = 1)
    
    # 50% Censoring
    mse.50 <- numeric(M)
    var.50 <- numeric(M)
    bias.50 <- numeric(M)
    predict.50 <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.5 <- matrix(0, nrow = 400, ncol = 1)
    
    mse.50.constant.3 <- numeric(M)
    var.50.constant.3 <- numeric(M)
    bias.50.constant.3 <- numeric(M)
    predict.5.cen <- matrix(0, nrow = 400, ncol = n[i])
    predict.50.constant.3 <- matrix(0, nrow = 400, ncol = nrow(new.x))
    constant.3.predict.50 <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.5.constant.3 <- matrix(0, nrow = 400, ncol = 1)
    
    mse.50.cen <- numeric(M)
    var.50.cen <- numeric(M)
    bias.50.cen <- numeric(M)
    predict.5.cen <- matrix(0, nrow = 400, ncol = n[i])
    predict.50.cen <- matrix(0, nrow = 400, ncol = nrow(new.x))
    cen.predict.50 <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.5.cen <- matrix(0, nrow = 400, ncol = 1)
    
    
    mse.50.median <- numeric(M)
    var.50.median <- numeric(M)
    bias.50.median <- numeric(M)
    predict.50.median <- matrix(0, nrow = 400, ncol = nrow(new.x))
    mse.5.median <- matrix(0, nrow = 400, ncol = 1)
    
    
    for (j in 1:M){
      
      x.i <- runif(n[i], min = 0, max = 1)
      
      e.i <- rnorm(n[i], mean = 0, sd = 1)
      
      y.i <- 5 + sin(2*pi*x.i) + e.i
      
      # NO Censoring
      cen.no <- rnorm(n[i], 0, 1) + 100
      pre.T.no <- pmin(y.i, cen.no) 
      delta.no <- (y.i <= cen.no) * 1
      g.no <- km.surv(pre.T.no, delta.no)
      y.i.no <- pre.T.no * delta.no / g.no
      
      bw.0 <- dpill(x.i, y.i.no)
      model.np.0 <- npreg(y.i.no ~ x.i, bws = bw.0, regtype = "ll")
      
      predict.no[j, ] <- predict(model.np.0, newdata = new.x)
      mse.0[j, ] <- mean((predict.no[j, ] - mx)^2)/M   
      
      cubic.no.3 <- lm(y.i.no ~ x.i + I(x.i^2)+ I(x.i^3))
      
      y.no.3 <- y.i.no - predict(cubic.no.3)
      bw.0.constant.3 <- dpill(x.i, y.no.3)
      model.constant.no.3 <- npreg(y.no.3 ~ x.i, bws = bw.0.constant.3, regtype = "ll")
      
      predict.no.constant.3[j, ] <- predict(cubic.no.3, newdata = new.x) 
      constant.3.predict.no[j, ] <- predict(model.constant.no.3, newdata = new.x) 
      
      mse.0.constant.3[j, ] <- mean(((predict.no.constant.3[j, ] + constant.3.predict.no[j, ]) - mx )^2)/M   
      
      y.cen.no <- y.i[delta.no == 1]
      x.cen.no <- x.i[delta.no == 1]
      
      model.cen.0 <- locfit(y.cen.no ~ x.cen.no, alpha = seq(0.2, 0.8, by = 0.05))
      predict.0.cen[j, ] <- predict(model.cen.0, newdata = x.i)
      
      T.no.cen <- pre.T.no - predict.0.cen[j, ]
      g.no.cen <- km.surv(T.no.cen, delta.no)
      y.i.no.cen <- T.no.cen * delta.no / g.no.cen
      
      bw.0.cen <- dpill(x.i, y.i.no.cen)
      model.np.0.cen <- npreg(y.i.no.cen ~ x.i, bws = bw.0.cen, regtype = "ll")
      
      predict.no.cen[j, ] <- predict(model.np.0.cen, newdata = new.x)
      cen.predict.no[j, ] <-predict(model.cen.0, newdata = new.x.no)
      mse.0.cen[j, ] <- mean(((predict.no.cen[j, ] + cen.predict.no[j, ])- mx )^2)/M 
      
      
      fit.no <- read.table(textConnection(capture.output(survfit(Surv(pre.T.no, delta.no) ~ 1))), skip = 2, header = TRUE)
      fit.median.no <- fit.no$median
      
      T.no.median <- pre.T.no - fit.median.no
      g.no.median <- km.surv(T.no.median, delta.no)
      y.i.no.median <- T.no.median * delta.no / g.no.median
      
      bw.0.median <- dpill(x.i, y.i.no.median)
      model.np.0.median <- npreg(y.i.no.median ~ x.i, bws = bw.0.median, regtype = "ll")
      
      predict.no.median[j, ] <- predict(model.np.0.median, newdata = new.x)
      mse.0.median[j, ] <- mean(((predict.no.median[j, ] + fit.median.no)- mx )^2)/M 
      
      # 10% Censoring
      
      cen.10 <- rnorm(n[i], 5, 1) + 2
      pre.T.10 <- pmin(y.i, cen.10)
      delta.10 <- (y.i <= cen.10) * 1
      g.10 <- km.surv(pre.T.10, delta.10)
      y.i.10 <- pre.T.10 * delta.10 / g.10
      
      bw.1 <- dpill(x.i, y.i.10)
      model.np.10 <- npreg(y.i.10 ~ x.i, bws = bw.1, regtype = "ll")
      
      predict.10[j, ] <- predict(model.np.10, newdata = new.x)
      mse.1[j, ] <- mean((predict.10[j, ] - mx )^2)/M   
      
      cubic.10.3 <- lm(y.i.10 ~ x.i + I(x.i^2)+ I(x.i^3))
      
      y.10.3 <- y.i.10 - predict(cubic.10.3)
      bw.10.constant.3 <- dpill(x.i, y.10.3)
      model.constant.10.3 <- npreg(y.10.3 ~ x.i, bws = bw.10.constant.3, regtype = "ll")
      
      predict.10.constant.3[j, ] <- predict(cubic.10.3, newdata = new.x) 
      constant.3.predict.10[j, ] <- predict(model.constant.10.3, newdata = new.x) 
      
      mse.1.constant.3[j, ] <- mean(((predict.10.constant.3[j, ] + constant.3.predict.10[j, ]) - mx )^2)/M   
      
      y.cen.10 <- y.i[delta.10 == 1]
      x.cen.10 <- x.i[delta.10 == 1]
      
      model.cen.1 <- locfit(y.cen.10 ~ x.cen.10, alpha = seq(0.2, 0.8, by = 0.05))
      predict.1.cen[j, ] <- predict(model.cen.1, newdata = x.i)
      
      T.10.cen <- pre.T.10 - predict.1.cen[j, ]
      g.10.cen <- km.surv(T.10.cen, delta.10)
      y.i.10.cen <- T.10.cen * delta.10 / g.10.cen
      
      bw.1.cen <- dpill(x.i, y.i.10.cen)
      model.np.1.cen <- npreg(y.i.10.cen ~ x.i, bws = bw.1.cen, regtype = "ll")
      
      predict.10.cen[j, ] <- predict(model.np.1.cen, newdata = new.x)
      cen.predict.10[j, ] <-predict(model.cen.1, newdata = new.x.10)
      mse.1.cen[j, ] <- mean(((predict.10.cen[j, ] + cen.predict.10[j, ])- mx )^2)/M 
      
      
      fit.10 <- read.table(textConnection(capture.output(survfit(Surv(pre.T.10, delta.10) ~ 1))), skip = 2, header = TRUE)
      fit.median.10 <- fit.10$median
      
      T.10.median <- pre.T.10 - fit.median.10
      g.10.median <- km.surv(T.10.median, delta.10)
      y.i.10.median <- T.10.median * delta.10 / g.10.median
      
      bw.1.median <- dpill(x.i, y.i.10.median)
      model.np.10.median <- npreg(y.i.10.median ~ x.i, bws = bw.1.median, regtype = "ll")
      
      predict.10.median[j, ] <- predict(model.np.10.median, newdata = new.x)
      mse.1.median[j, ] <- mean(((predict.10.median[j, ] + fit.median.10)- mx )^2)/M 
      
      # 30% Censoring
      
      cen.30 <- rnorm(n[i], 5, 1) + 0.85
      pre.T.30 <- pmin(y.i, cen.30)
      delta.30 <- (y.i <= cen.30) * 1
      g.30 <- km.surv(pre.T.30, delta.30)
      y.i.30 <- pre.T.30 * delta.30 / g.30
      
      bw.3 <- dpill(x.i, y.i.30)
      model.np.30 <- npreg(y.i.30 ~ x.i, bws = bw.3, regtype = "ll")
      
      predict.30[j, ] <- predict(model.np.30, newdata = new.x)
      mse.3[j, ] <- mean((predict.30[j, ] - mx )^2)/M   
      
      cubic.30.3 <- lm(y.i.30 ~ x.i + I(x.i^2) + I(x.i^3))
      
      y.30.3 <- y.i.30 - predict(cubic.30.3)
      bw.30.constant.3 <- dpill(x.i, y.30.3)
      model.constant.30.3 <- npreg(y.30.3 ~ x.i, bws = bw.30.constant.3, regtype = "ll")
      
      predict.30.constant.3[j, ] <- predict(cubic.30.3, newdata = new.x) 
      constant.3.predict.30[j, ] <- predict(model.constant.30.3, newdata = new.x) 
      mse.3.constant.3[j, ] <- mean(((predict.30.constant.3[j, ] + constant.3.predict.30[j, ]) - mx )^2)/M   
      
      y.cen.30 <- y.i[delta.30 == 1]
      x.cen.30 <- x.i[delta.30 == 1]
      
      model.cen.3 <- locfit(y.cen.30 ~ x.cen.30, alpha = seq(0.2, 0.8, by = 0.05))
      predict.3.cen[j, ] <- predict(model.cen.3, newdata = x.i)
      
      T.30.cen <- pre.T.30 - predict.3.cen[j, ]
      g.30.cen <- km.surv(T.30.cen, delta.30)
      y.i.30.cen <- T.30.cen * delta.30 / g.30.cen
      
      bw.3.cen <- dpill(x.i, y.i.30.cen)
      model.np.3.cen <- npreg(y.i.30.cen ~ x.i, bws = bw.3.cen, regtype = "ll")
      
      predict.30.cen[j, ] <- predict(model.np.3.cen, newdata = new.x)
      cen.predict.30[j, ] <-predict(model.cen.3, newdata = new.x.30)
      mse.3.cen[j, ] <- mean(((predict.30.cen[j, ] + cen.predict.30[j, ])- mx )^2)/M 
      
      fit.30 <- read.table(textConnection(capture.output(survfit(Surv(pre.T.30, delta.30) ~ 1))), skip = 2, header = TRUE)
      fit.median.30 <- fit.30$median
      
      T.30.median <- pre.T.30 - fit.median.30
      g.30.median <- km.surv(T.30.median, delta.30)
      y.i.30.median <- T.30.median * delta.30 / g.30.median
      
      bw.3.median <- dpill(x.i, y.i.30.median)
      model.np.30.median <- npreg(y.i.30.median ~ x.i, bw.3.median, regtype = "ll")
      
      predict.30.median[j, ] <- predict(model.np.30.median, newdata = new.x)
      mse.3.median[j, ] <- mean(((predict.30.median[j, ] + fit.median.30)- mx )^2)/M 
      
      
      # 50% Censoring
      
      cen.50 <- rnorm(n[i], 5, 1) + 0
      pre.T.50 <- pmin(y.i, cen.50)
      delta.50 <- (y.i <= cen.50) * 1
      g.50 <- km.surv(pre.T.50, delta.50)
      y.i.50 <- pre.T.50 * delta.50 / g.50
      
      bw.5 <- dpill(x.i, y.i.50)
      model.np.50 <- npreg(y.i.50 ~ x.i, bws = bw.5, regtype = "ll")
      
      predict.50[j, ] <- predict(model.np.50, newdata = new.x)
      mse.5[j, ] <- mean((predict.50[j, ] - mx )^2)/M   
      
      cubic.50.3 <- lm(y.i.50 ~ x.i + I(x.i^2)+ I(x.i^3))
      
      y.50.3 <- y.i.50 - predict(cubic.50.3)
      bw.50.constant.3 <- dpill(x.i, y.50.3)
      model.constant.50.3 <- npreg(y.50.3 ~ x.i, bws = bw.50.constant.3, regtype = "ll")
      
      predict.50.constant.3[j, ] <- predict(cubic.50.3, newdata = new.x) 
      constant.3.predict.50[j, ] <- predict(model.constant.50.3, newdata = new.x) 
      mse.5.constant.3[j, ] <- mean(((predict.50.constant.3[j, ] + constant.3.predict.50[j, ]) - mx )^2)/M   
      
      y.cen.50 <- y.i[delta.50 == 1]
      x.cen.50 <- x.i[delta.50 == 1]
      
      model.cen.5 <- locfit(y.cen.50 ~ x.cen.50, alpha = seq(0.2, 0.8, by = 0.05))
      predict.5.cen[j, ] <- predict(model.cen.5, newdata = x.i)
      
      T.50.cen <- pre.T.50 - predict.5.cen[j, ]
      g.50.cen <- km.surv(T.50.cen, delta.50)
      y.i.50.cen <- T.50.cen * delta.50 / g.50.cen
      
      bw.5.cen <- dpill(x.i, y.i.50.cen)
      model.np.5.cen <- npreg(y.i.50.cen ~ x.i, bws = bw.5.cen, regtype = "ll")
      
      predict.50.cen[j, ] <- predict(model.np.5.cen, newdata = new.x)
      cen.predict.50[j, ] <-predict(model.cen.5, newdata = new.x.50)
      mse.5.cen[j, ] <- mean(((predict.50.cen[j, ] + cen.predict.50[j, ])- mx )^2)/M 
      
      
      fit.50 <- read.table(textConnection(capture.output(survfit(Surv(pre.T.50, delta.50) ~ 1))), skip = 2, header = TRUE)
      fit.median.50 <- fit.50$median
      
      T.50.median <- pre.T.50 - fit.median.50
      g.50.median <- km.surv(T.50.median, delta.50)
      y.i.50.median <- T.50.median * delta.50 / g.50.median
      
      bw.5.median <- dpill(x.i, y.i.50.median)
      model.np.50.median <- npreg(y.i.50.median ~ x.i, bws = bw.5.median, regtype = "ll")
      
      predict.50.median[j, ] <- predict(model.np.50.median, newdata = new.x)
      mse.5.median[j, ] <- mean(((predict.50.median[j, ] + fit.median.50)- mx )^2)/M 
      
    }
    
    # No Censoring
    mse.no <- sum(mse.0)
    var.no <- sum(apply(predict.no, 2, var)/length(mx))
    bias.no <- mse.no - var.no
    plot.no[i, ] <- apply(predict.no, 2, mean)
    
    result.no[i, ] <- c(bias.no, var.no, mse.no)
    rownames(result.no) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.no) <- c("Bias^2", "Variance", "Mse")
    
    
    mse.no.constant.3 <- sum(mse.0.constant.3)
    var.no.constant.3 <- sum(apply((constant.3.predict.no + predict.no.constant.3), 2, var)/length(mx))
    bias.no.constant.3 <- mse.no.constant.3 - var.no.constant.3
    plot.no.constant.3[i, ] <- apply((predict.no.constant.3 + constant.3.predict.no), 2, mean)
    
    result.no.constant.3[i, ] <- c(bias.no.constant.3, var.no.constant.3, mse.no.constant.3)
    rownames(result.no.constant.3) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.no.constant.3) <- c("Bias^2", "Variance", "Mse")
    
    mse.no.cen <- sum(mse.0.cen)
    var.no.cen <- sum(apply((cen.predict.no + predict.no.cen), 2, var)/length(mx))
    bias.no.cen <- mse.no.cen - var.no.cen
    plot.no.cen[i, ] <- apply((predict.no.cen + cen.predict.no), 2, mean)
    
    result.no.cen[i, ] <- c(bias.no.cen, var.no.cen, mse.no.cen)
    rownames(result.no.cen) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.no.cen) <- c("Bias^2", "Variance", "Mse")
    
    
    mse.no.median <- sum(mse.0.median)
    var.no.median <- sum(apply((predict.no.median + fit.median.no), 2, var)/length(mx))
    bias.no.median <- mse.no.median - var.no.median
    plot.no.median[i, ] <- apply((predict.no.median + fit.median.no), 2, mean)
    
    result.no.median[i, ] <- c(bias.no.median, var.no.median, mse.no.median)
    rownames(result.no.median) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.no.median) <- c("Bias^2", "Variance", "Mse")
    
    # 10% Censoring
    mse.10 <- sum(mse.1)
    var.10 <- sum(apply(predict.10, 2, var)/length(mx))
    bias.10 <- mse.10 - var.10
    plot.10[i, ] <- apply(predict.10, 2, mean)
    
    
    result.10[i, ] <- c(bias.10, var.10, mse.10)
    rownames(result.10) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.10) <- c("Bias^2", "Variance", "Mse")
    
    mse.10.constant.3 <- sum(mse.1.constant.3)
    var.10.constant.3 <- sum(apply((predict.10.constant.3 + constant.3.predict.10), 2, var)/length(mx))
    bias.10.constant.3 <- mse.10.constant.3 - var.10.constant.3
    plot.10.constant.3[i, ] <- apply((predict.10.constant.3 + constant.3.predict.10), 2, mean)
    
    result.10.constant.3[i, ] <- c(bias.10.constant.3, var.10.constant.3, mse.10.constant.3)
    rownames(result.10.constant.3) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.10.constant.3) <- c("Bias^2", "Variance", "Mse")
    
    mse.10.cen <- sum(mse.1.cen)
    var.10.cen <- sum(apply((predict.10.cen + cen.predict.10), 2, var)/length(mx))
    bias.10.cen <- mse.10.cen - var.10.cen
    plot.10.cen[i, ] <- apply((predict.10.cen + cen.predict.10), 2, mean)
    
    result.10.cen[i, ] <- c(bias.10.cen, var.10.cen, mse.10.cen)
    rownames(result.10.cen) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.10.cen) <- c("Bias^2", "Variance", "Mse")
    
    mse.10.median <- sum(mse.1.median)
    var.10.median <- sum(apply(predict.10.median, 2, var)/length(mx))
    bias.10.median <- mse.10.median - var.10.median
    plot.10.median[i, ] <- apply((predict.10.median + fit.median.10), 2, mean)
    
    result.10.median[i, ] <- c(bias.10.median, var.10.median, mse.10.median)
    rownames(result.10.median) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.10.median) <- c("Bias^2", "Variance", "Mse")
    
    
    # 30% Censoring
    mse.30 <- sum(mse.3)
    var.30 <- sum(apply(predict.30, 2, var)/length(mx))
    bias.30 <- mse.30 - var.30
    plot.30[i, ] <- apply(predict.30, 2, mean)
    
    result.30[i, ] <- c(bias.30, var.30, mse.30)
    rownames(result.30) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.30) <- c("Bias^2", "Variance", "Mse")
    
    mse.30.constant.3 <- sum(mse.3.constant.3)
    var.30.constant.3 <- sum(apply((predict.30.constant.3 + constant.3.predict.30), 2, var)/length(mx))
    bias.30.constant.3 <- mse.30.constant.3 - var.30.constant.3
    plot.30.constant.3[i, ] <- apply((predict.30.constant.3 + constant.3.predict.30), 2, mean)
    
    result.30.constant.3[i, ] <- c(bias.30.constant.3, var.30.constant.3, mse.30.constant.3)
    rownames(result.30.constant.3) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.30.constant.3) <- c("Bias^2", "Variance", "Mse")
    
    mse.30.cen <- sum(mse.3.cen)
    var.30.cen <- sum(apply((predict.30.cen + cen.predict.30), 2, var)/length(mx))
    bias.30.cen <- mse.30.cen - var.30.cen
    plot.30.cen[i, ] <- apply((predict.30.cen + cen.predict.30), 2, mean)
    
    result.30.cen[i, ] <- c(bias.30.cen, var.30.cen, mse.30.cen)
    rownames(result.30.cen) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.30.cen) <- c("Bias^2", "Variance", "Mse")
    
    mse.30.median <- sum(mse.3.median)
    var.30.median <- sum(apply((predict.30.median + fit.median.30), 2, var)/length(mx))
    bias.30.median <- mse.30.median - var.30.median
    plot.30.median[i, ] <- apply((predict.30.median + fit.median.30), 2, mean)
    
    result.30.median[i, ] <- c(bias.30.median, var.30.median, mse.30.median)
    rownames(result.30.median) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.30.median) <- c("Bias^2", "Variance", "Mse")
    
    # 50% Censoring
    mse.50 <- sum(mse.5)
    var.50 <- sum(apply(predict.50, 2, var)/length(mx))
    bias.50 <- mse.50 - var.50
    plot.50[i, ] <- apply(predict.50, 2, mean)
    
    result.50[i, ] <- c(bias.50, var.50, mse.50)
    rownames(result.50) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.50) <- c("Bias^2", "Variance", "Mse")
    
    
    mse.50.constant.3 <- sum(mse.5.constant.3)
    var.50.constant.3 <- sum(apply((predict.50.constant.3 + constant.3.predict.50), 2, var)/length(mx))
    bias.50.constant.3 <- mse.50.constant.3 - var.50.constant.3
    plot.50.constant.3[i, ] <- apply((predict.50.constant.3 + constant.3.predict.50), 2, mean)
    
    result.50.constant.3[i, ] <- c(bias.50.constant.3, var.50.constant.3, mse.50.constant.3)
    rownames(result.50.constant.3) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.50.constant.3) <- c("Bias^2", "Variance", "Mse")
    
    mse.50.cen <- sum(mse.5.cen)
    var.50.cen <- sum(apply((predict.50.cen + cen.predict.50), 2, var)/length(mx))
    bias.50.cen <- mse.50.cen - var.50.cen
    plot.50.cen[i, ] <- apply((predict.50.cen + cen.predict.50), 2, mean)
    
    result.50.cen[i, ] <- c(bias.50.cen, var.50.cen, mse.50.cen)
    rownames(result.50.cen) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.50.cen) <- c("Bias^2", "Variance", "Mse")
    
    mse.50.median <- sum(mse.5.median)
    var.50.median <- sum(apply((predict.50.median + fit.median.50), 2, var)/length(mx))
    bias.50.median <- mse.50.median - var.50.median
    plot.50.median[i, ] <- apply((predict.50.median + fit.median.50), 2, mean)
    
    result.50.median[i, ] <- c(bias.50.median, var.50.median, mse.50.median)
    rownames(result.50.median) <- c("n = 100", "n = 200", "n = 400")
    colnames(result.50.median) <- c("Bias^2", "Variance", "Mse")
    
    
  }
)
###########################################################################################

# result

result.no
result.no.constant.3
result.no.cen
result.no.median



result.10
result.10.constant.3
result.10.cen
result.10.median


result.30
result.30.constant.3
result.30.cen
result.30.median


result.50
result.50.constant.3
result.50.cen
result.50.median

##########################################################################################

# Mean Plot
library(latex2exp)

plot(mx, type="l", main = "표본의 크기(N) = 100", xlab = "new.x", ylab = "True value", lwd = 3, col = 1, cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)
lines(plot.10[1, ], col = 2, lty = 2, lwd = 3)
lines(plot.10.constant.3[1, ], col = 3, lty = 3, lwd = 3)
lines(plot.10.cen[1, ], col = 4, lty = 4, lwd = 3)
lines(plot.10.median[1, ], col = 5, lty = 5, lwd = 3)
lines(plot.10.sin[1, ], col = 6, lty = 6, lwd = 3)
legend("topright", c("True Function", "KSV", TeX("$GNR_P"),TeX("$GNR_{NP}"), TeX("$GNR_M"), TeX("$GNR_T")), lty = 1:6, 
       col = 1:6,   cex = 1.2, y.intersp = 0.25, lwd = 3, inset = c(-0.15, -0.06))

plot(mx, type="l", main = "표본의 크기(N) = 200", xlab = "new.x", ylab = "True value", lwd = 3, col = 1, cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)
lines(plot.10[2, ], col = 2, lty = 2, lwd = 3)
lines(plot.10.constant.3[2, ], col = 3, lty = 3, lwd = 3)
lines(plot.10.cen[2, ], col = 4, lty = 4, lwd = 3)
lines(plot.10.median[2, ], col = 5, lty = 5, lwd = 3)
lines(plot.10.sin[2, ], col = 6, lty = 6, lwd = 3)
legend("topright", c("True Function", "KSV", TeX("$GNR_P"),TeX("$GNR_{NP}"), TeX("$GNR_M"), TeX("$GNR_T")), lty = 1:6, 
       col = 1:6,   cex = 1.2, y.intersp = 0.25, lwd = 3, inset = c(-0.15, -0.06))



plot(mx, type="l", main = "표본의 크기(N) = 400", xlab = "new.x", ylab = "True value", lwd = 3, col = 1, cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)
lines(plot.10[3, ], col = 2, lty = 2, lwd = 3)
lines(plot.10.constant.3[3, ], col = 3, lty = 3, lwd = 3)
lines(plot.10.cen[3, ], col = 4, lty = 4, lwd = 3)
lines(plot.10.median[3, ], col = 5, lty = 5, lwd = 3)
lines(plot.10.sin[3, ], col = 6, lty = 6, lwd = 3)
legend("topright", c("True Function", "KSV", TeX("$GNR_P"),TeX("$GNR_{NP}"), TeX("$GNR_M"), TeX("$GNR_T")), lty = 1:6, 
       col = 1:6,   cex = 1.2, y.intersp = 0.25, lwd = 3, inset = c(-0.15, -0.06))

