# Assume parametric basis distribution is Exponential distribution... 
lambda1 <- 0.01 # parameter for generate Y  
lambda2 <- 4/lambda1 # parameter for generate C 
beta0 <- 1:5
beta0 <- beta0/sqrt(sum(beta0^2)) # true coefficients
Y <- c() # real survival time
C <- c() # censoring time
t <- c() # observed survival time
censor <- c() # censoring indicator
X <- matrix(, nrow=100, ncol=length(beta0)) # explanatory variables
U1 <- c() # random variable from U(0,1) for generate Y
U2 <- c() # random variable from U(0,1) for generate C 

set.seed(1004)

for (i in 1:100) {
  U1[i] <- runif(1)
  U2[i] <- runif(1)
  X[i,] <- rnorm(5,mean=0,sd=1)
  Y[i] <- -(1/lambda1)*log(U1[i])*exp(-X[i,]%*%beta0)
  C[i] <- (lambda2)*(1-U2[i])
  t[i] <- min(Y[i],C[i])
  censor[i] <- ifelse(Y[i]<=C[i],1,0)
}

mean(1-censor) # censoring rate

dt5 <- data.frame(X=X, Y=Y, C=C, time=t, censor=censor)
head(dt5,10)


library(survival)
So <- with(dt5, Surv(time,censor==1))
dt5 <- dt5[order(So[,1]),]
So <- So[order(So[,1]),]
X <- model.matrix(~X.1+X.2+X.3+X.4+X.5, data=dt5)[,-1]
eta <- predict(coxph(So ~ X))
model1 <- coxph(So ~ X)
beta_hat <- model1[["coefficients"]]
# time.sort <- sort(So[,1]) # sorting survival time
time.sort1 <- seq(from=5, to=250, length.out=20) # sorting survival time
time.sort2 <- c(50, 150, 250) # sorting survival time
lambda <- seq(0,1,0.01) # parameter for k-fold crossvalidation


