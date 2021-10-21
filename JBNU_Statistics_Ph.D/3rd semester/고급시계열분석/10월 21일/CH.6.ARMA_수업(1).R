##############  package
library(data.table)
library(ggplot2)
library(forecast)
library(gridExtra)


#################################################
##### AR(1)
#################################################

###AR(p) 을 생성하는 방법 :

sim_ar <- function(n, mu, phi){
  ### n : sample size
  ### mu : mean
  ### phi : p-dim coefficients
  p <- length(phi)
  z <- rnorm(n+100)  #epsilon ~ WN(0, sigma^2), iid N(0,1)
  delta <- (1-sum(phi))*mu
  
  for (k in (length(phi)+1):(n+100)){
    
    z[k] <- delta + sum(z[(k-1):(k-p)]*phi) + rnorm(1)
  }
  
  return(z[-(1:100)])
}


z <- sim_ar(100, 0, phi=c(0.5))  ##AR(1)

par(mfrow=c(2,1))
# plot.ts(z)
acf(z)
pacf(z)



##AR(1) phi=0.5
z <- arima.sim(n=100, ##order=c(p,d,q)  ARMA : d=0
               list(order=c(1,0,0), ar= 0.5), 
               rand.gen = rnorm)
z <- arima.sim(n=10000, ##order=c(p,d,q)  ARMA : d=0
               list(ar= -0.5) 
               )

par(mfrow=c(2,1))
plot.ts(z)
acf(z, lag.max = 24)
pacf(z, lag.max = 24)

par(mfrow=c(1,1))
plot(z[-1],z[-length(z)], pch=16,
     xlab="Z(t)", ylab='Z(t-1)')



##AR(2)
z <- arima.sim(n=10000, list(order=c(2,0,0), 
                           ar=c(0.5, -0.4)), rand.gen = rnorm)


par(mfrow=c(2,1))
# plot.ts(z)
acf(z, lag.max = 24)
pacf(z, lag.max = 24)




p1 <- ggAcf(z, lag.max=24) + 
  theme_bw() + ylim(-1,1) +
  theme(plot.title = element_blank())

p2 <- ggPacf(z, lag.max=24) + 
  theme_bw() +ylim(-1,1) +
  theme(plot.title = element_blank())

grid.arrange(p1, p2, nrow = 2)


###MA(1) 을 생성하는 방법 :

sim_ma <- function(n, mu, theta){
  ### n : sample size
  ### mu : mean
  ### theta : q-dim coefficients
  
  q <- length(theta)
  
  ep <- rnorm(n+100)
  z <- ep
  
  for (k in (q+1):(n+100)){
    
    z[k] <- mu + ep[k] - sum(ep[(k-1):(k-q)]*theta)
  }
  
  return(z[-(1:100)])
}

z <- sim_ma(100, 0, theta=0.9)


ts.plot(z)
acf(z)
pacf(z)



##MA(1)
z <- arima.sim(n=10000, list(order=c(0,0,1), ma= -0.9), 
               rand.gen = rnorm)
# z <- arima.sim(n=100, list(ma= 0.9), rand.gen = rnorm)

par(mfrow=c(3,1))
ts.plot(z)
acf(z, lag.max = 24)
pacf(z, lag.max = 24)


##MA(2)
z <- arima.sim(n=10000, list(order=c(0,0,2), 
                           ma=c(0.5, -0.2)))

par(mfrow=c(3,1))
plot.ts(z)
acf(z, lag.max=20)
pacf(z, lag.max=20)


##ARMA(1,1)

z <- arima.sim(n=10000, list(order=c(1,0,1), 
                            ar=-0.5, ma=-0.3), 
               rand.gen = rnorm)

par(mfrow=c(2,1))
# plot.ts(z)
acf(z, lag.max=20)
pacf(z, lag.max=20)



