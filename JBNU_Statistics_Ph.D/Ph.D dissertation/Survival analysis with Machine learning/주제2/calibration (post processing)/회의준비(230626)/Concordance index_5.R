# reference : https://www.erikdrysdale.com/survConcordance/

## Function for the concordance index (C-Index)

# Function for the i'th concordance
cindex_i <- function(So,eta,i) {
  tt_i <- So[i,1]
  dd_i <- So[i,2]
  idx.k <- which(So[,1] > tt_i)
  conc <- sum(eta[i] > eta[idx.k]  )
  disc <- sum(eta[i] < eta[idx.k]  )
  return(c(conc,disc))
}

# Wrapper for total concordance
cindex <- function(So,eta) {
  conc.disc <- c(0,0)
  for (i in which(So[,2] == 1)) {
    conc.disc <- conc.disc + cindex_i(So,eta,i)
  }
  names(conc.disc) <- c('concordant','discordant')
  return(conc.disc)
}

# Example of the 'veteran' data
library(survival)
So <- with(veteran,Surv(time,status==1))
X <- model.matrix(~factor(trt)+karno+diagtime+age+factor(prior),data=veteran)[,-1]
eta <- predict(coxph(So ~ X))
Som <- as.matrix(So)
Som[,1] <- Som[,1] + (1-Som[,2])*min(Som[,1])/2
survConcordance.fit(So, eta)[1:2]
cindex(Som, eta)
concordance(coxph(So ~ X))


## Equation (2)

sigmoid <- function(x) { 1/(1+exp(-x)) }

l_eta_i <- function(So,eta,i) {
  tt_i <- So[i,1]
  dd_i <- So[i,2]
  idx.k.i <- which(So[,1] > tt_i)
  loss.i <- sum(1 + log( sigmoid(eta[i] - eta[idx.k.i]) )/log(2) )
  return(loss.i)
}

l_eta <- function(So, eta) {
  loss <- 0
  for (i in which(So[,2] == 1)) {
    loss <- loss + l_eta_i(So,eta,i)
  }
  return(-loss / nrow(So)^2)
}


## Equation (3)

sigmoid2 <- function(x) { 1/(1+exp(x)) }

dl_eta_i <- function(eta,So,i) {
  tt_i <- So[i,1]
  dd_i <- So[i,2]
  idx.k <- which(So[,1] > tt_i)
  idx.j <- which(tt_i > So[,1] & So[,2]==1) 
  res.i <- dd_i*sum(sigmoid2(eta[i] - eta[idx.k])) - sum(sigmoid2(eta[idx.j] - eta[i]))
  return(res.i)
}

dl_eta <- function(X,eta,So) {
  grad <- rep(0, ncol(X))
  for (i in seq_along(eta)) {
    grad <- grad + X[i,] * dl_eta_i(eta, So, i)
  }
  grad <- -1 * grad / nrow(X)^2
  return(grad)
}











#--------------------------------------------------------------------------#

## minimize Q(\beta)

library(survival)
So <- with(veteran,Surv(time,status==1))
X <- model.matrix(~factor(trt)+karno+diagtime+age+factor(prior),data=veteran)[,-1]
eta <- predict(coxph(So ~ X))
model1 <- coxph(So ~ X)
beta_hat <- model1[["coefficients"]]
beta <- as.matrix(rep(0.5,ncol(X)))


# Q(\beta)

sigmoid <- function(x) { 1/(1+exp(-x)) }

C_tilde_beta_i <- function(So,beta,X,i) {
  tt_i <- So[i,1]
  dd_i <- So[i,2]
  loss.i <- 0
  for (j in which(So[,1] > tt_i)) {
    loss.i <- loss.i + (1 + log( sigmoid( t(beta)%*%(as.matrix(X[i,]) - as.matrix(X[j,])) ) )/log(2))
  }
  return(loss.i)
}

C_tilde_beta <- function(So,beta,X) {
  loss <- 0
  for (i in which(So[,2] == 1)) {
    loss <- loss + C_tilde_beta_i(So,beta,X,i)
  }
  return(loss / nrow(So)^2)
}

C_tilde_t_beta_i <- function(So,beta,X,t,i) {
  tt_i <- So[i,1]
  dd_i <- So[i,2]
  loss.i <- 0
  for (j in which(So[,1] > t)) {
    loss.i <- loss.i + (1 + log( sigmoid( t(beta)%*%(as.matrix(X[i,]) - as.matrix(X[j,])) ) )/log(2))
  }
  return(loss.i)
} 

C_tilde_t_beta <- function(So,beta,X,t) {
  loss <- 0
  for ( i in which((So[,2] == 1) & (So[,1] <= t)) ) {
    loss <- loss + C_tilde_t_beta_i(So,beta,X,t,i)
  }
  return(loss / nrow(So)^2)
}

lambda <- seq(0,1,0.01)

Q_beta <- function(So,beta,X,lambda) {
  loss_Q <- -C_tilde_beta(So,beta,X) + 
    lambda*(C_tilde_t_beta(So,beta,X,So[,1][2]) - C_tilde_t_beta(So,beta,X,So[,1][1]))^2
  return(loss_Q)
}


# Take the derivative Q(\beta) by \beta

der_C_tilde_beta_i <- function(So,beta,X,i) {
  tt_i <- So[i,1]
  dd_i <- So[i,2]
  loss.i <- 0
  for (j in which(So[,1] > tt_i)) {
    loss.i <- loss.i + 
      (as.matrix(X[j,]) - as.matrix(X[i,]))%*%sigmoid( t(beta)%*%(as.matrix(X[j,]) - as.matrix(X[i,])) )
  }
  return(-loss.i)
}

der_C_tilde_beta <- function(So,beta,X) {
  loss <- 0
  for (i in which(So[,2] == 1)) {
    loss <- loss + der_C_tilde_beta_i(So,beta,X,i)
  }
  return(loss / nrow(So)^2 / log(2))
}

der_C_tilde_t_beta_i <- function(So,beta,X,t,i) {
  tt_i <- So[i,1]
  dd_i <- So[i,2]
  loss.i <- 0
  for (j in which(So[,1] > t)) {
    loss.i <- loss.i + 
      (as.matrix(X[j,]) - as.matrix(X[i,]))%*%sigmoid( t(beta)%*%(as.matrix(X[j,]) - as.matrix(X[i,])) )
  }
  return(-loss.i)
} 

der_C_tilde_t_beta <- function(So,beta,X,t) {
  loss <- 0
  for ( i in which((So[,2] == 1) & (So[,1] <= t)) ) {
    loss <- loss + der_C_tilde_t_beta_i(So,beta,X,t,i)
  }
  return(loss / nrow(So)^2 / log(2))
}

der_Q_beta <- function(So,beta,X,lambda) {
  der_loss_Q <- -der_C_tilde_beta(So,beta,X) + 
    2*lambda*as.numeric(C_tilde_t_beta(So,beta,X,So[,1][2]) - C_tilde_t_beta(So,beta,X,So[,1][1]))*
    (der_C_tilde_t_beta(So,beta,X,So[,1][2]) - der_C_tilde_t_beta(So,beta,X,So[,1][1]))
  return(der_loss_Q)
}



# Gradient descent algorithm

my.gradient.descent <- function(So, X, alpha, lambda) {
  
  # alpha : learning rate
  # lambda : penalty parameter
  
  if(lambda <= 0) 
    stop("Lambda is non-positive value. Please insert positive value of lambda.")
  
  X <-as.matrix(X)
  
  beta.old <- rep(1,ncol(X)) # The initial value of coefficient vector
  beta.old <- beta.old/sqrt(sum(beta.old^2))
  
  for (i in 1:1000) {
    
    
    beta.new <- beta.old - alpha*der_Q_beta(So,beta.old,X,lambda)
    
    diff <- abs(C_tilde_beta(So,beta.new,X) - C_tilde_beta(So,beta.old,X))  
    
    # diff <- sqrt(sum((beta.new - beta.old)^2))/sqrt(sum(beta.old^2))
    
    cat("( iteration , difference ) = (", i, ",", diff, ")\n")
    
    if (diff < 1E-8) break
    
    beta.old <- beta.new
    
  }
  
  cat("Algorithm converged...","\n\n")
  
  return(beta.new)
  
}








#--------------------------------------------------------------------------#

## Generalization

library(survival)
So <- with(veteran, Surv(time,status==1))
X <- model.matrix(~factor(trt)+karno+diagtime+age+factor(prior), data=veteran)[,-1]
eta <- predict(coxph(So ~ X))
model1 <- coxph(So ~ X)
beta_hat <- model1[["coefficients"]]
beta <- as.matrix(rep(0.5, ncol(X)))
time.sort <- sort(So[,1]) # sorting survival time
lambda <- seq(0,1,0.01) # parameter for k-fold crossvalidation




# U(\beta)

sigmoid <- function(x) { 1/(1+exp(-x)) }

C_tilde_beta_i <- function(So, beta, X, i) {
  tt_i <- So[i,1]
  dd_i <- So[i,2]
  loss.i <- 0
  for (j in which(So[,1] > tt_i)) {
    loss.i <- loss.i + (1 + log( sigmoid( t(beta)%*%(as.matrix(X[i,]) - as.matrix(X[j,])) ) )/log(2))
  }
  return(loss.i)
}

C_tilde_beta <- function(So, beta, X) {
  loss <- 0
  for (i in which(So[,2] == 1)) {
    loss <- loss + C_tilde_beta_i(So,beta,X,i)
  }
  return(loss / nrow(So)^2)
}

C_tilde_t_beta_i <- function(So, beta, X, t, i) {
  tt_i <- So[i,1]
  dd_i <- So[i,2]
  loss.i <- 0
  for (j in which(So[,1] > t)) {
    loss.i <- loss.i + (1 + log( sigmoid( t(beta)%*%(as.matrix(X[i,]) - as.matrix(X[j,])) ) )/log(2))
  }
  return(loss.i)
} 

C_tilde_t_beta <- function(So, beta, X, t) {
  loss <- 0
  for ( i in which((So[,2] == 1) & (So[,1] <= t)) ) {
    loss <- loss + C_tilde_t_beta_i(So, beta, X, t, i)
  }
  return(loss / nrow(So)^2)
}

sum_C_tilde_t_beta <- function(So, beta, X, time.sort) {
  loss <- 0
  for ( k in 2:length(time.sort) ) {
    loss <- loss + (C_tilde_t_beta(So, beta, X, time.sort[k]) - C_tilde_t_beta(So, beta, X, time.sort[k-1]))^2
  }
  return(loss)
}

U_beta <- function(So, beta, X, lambda, time.sort) {
  loss_U <- -C_tilde_beta(So, beta, X) + 
    lambda*sum_C_tilde_t_beta(So, beta, X, time.sort)
  return(loss_U)
}




# Take the derivative U(\beta) by \beta

der_C_tilde_beta_i <- function(So, beta, X, i) {
  tt_i <- So[i,1]
  dd_i <- So[i,2]
  loss.i <- 0
  for (j in which(So[,1] > tt_i)) {
    loss.i <- loss.i + 
      (as.matrix(X[j,]) - as.matrix(X[i,]))%*%sigmoid( t(beta)%*%(as.matrix(X[j,]) - as.matrix(X[i,])) )
  }
  return(-loss.i)
}

der_C_tilde_beta <- function(So, beta, X) {
  loss <- 0
  for (i in which(So[,2] == 1)) {
    loss <- loss + der_C_tilde_beta_i(So, beta, X, i)
  }
  return(loss / nrow(So)^2 / log(2))
}

der_C_tilde_t_beta_i <- function(So, beta, X, t, i) {
  tt_i <- So[i,1]
  dd_i <- So[i,2]
  loss.i <- 0
  for (j in which(So[,1] > t)) {
    loss.i <- loss.i + 
      (as.matrix(X[j,]) - as.matrix(X[i,]))%*%sigmoid( t(beta)%*%(as.matrix(X[j,]) - as.matrix(X[i,])) )
  }
  return(-loss.i)
} 

der_C_tilde_t_beta <- function(So, beta, X, t) {
  loss <- 0
  for ( i in which((So[,2] == 1) & (So[,1] <= t)) ) {
    loss <- loss + der_C_tilde_t_beta_i(So, beta, X, t, i)
  }
  return(loss / nrow(So)^2 / log(2))
}

sum_der_C_tilde_t_beta <- function(So, beta, X, time.sort) {
  loss <- 0
  for ( k in 2:length(time.sort) ) {
    loss <- loss + as.numeric(C_tilde_t_beta(So, beta, X, time.sort[k]) - C_tilde_t_beta(So, beta, X, time.sort[k-1]))*
      (der_C_tilde_t_beta(So, beta, X, time.sort[k]) - der_C_tilde_t_beta(So, beta, X, time.sort[k-1]))
  }
  return(loss)
}

der_U_beta <- function(So, beta, X, lambda, time.sort) {
  der_loss_Q <- -der_C_tilde_beta(So, beta, X) + 
    2*lambda*sum_der_C_tilde_t_beta(So, beta, X, time.sort)
  return(der_loss_Q)
}




# Gradient descent algorithm

my.gradient.descent.U <- function(So, X, alpha, lambda, time.sort) {
  
  # alpha : learning rate
  # lambda : penalty parameter
  
  if(lambda < 0) 
    stop("Lambda is negative value. Please insert non-negative value of lambda.")
  
  X <-as.matrix(X)
  
  beta.old <- rep(0.5,ncol(X)) # The initial value of coefficient vector
  beta.old <- beta.old/sqrt(sum(beta.old^2))
  iteration <- c()
  difference <- c()
  
  for (i in 1:100) {
    
    
    beta.new <- beta.old - alpha*der_U_beta(So, beta.old, X, lambda, time.sort)
    
    diff <- abs(C_tilde_beta(So, beta.new, X) - C_tilde_beta(So, beta.old, X))  
    
    iteration <- c(iteration, i)
    difference <- c(difference, diff)
    
    cat("( iteration , difference ) = (", i, ",", diff, ")\n")
    
    if (diff < 1E-5) break
    
    beta.old <- beta.new
    
  }
  
  cat("Algorithm converged...","\n\n")
  
  return(list(beta.new=beta.new, iteration=iteration, difference=difference))
  
}  




# Mini-batch(Stochastic) gradient descent algorithm
# reference : https://jjeongil.tistory.com/577 (stochastic)
# reference : https://data-science-hi.tistory.com/164 (mini-batch)
# reference : https://brunch.co.kr/@linecard/561 (mini-batch)
# reference : https://sonny-daily-story.tistory.com/5 (mini-batch)
# reference : https://www.youtube.com/watch?v=87Q2LIlMWoY (mini-batch)
# reference : https://light-tree.tistory.com/133 (mini-batch)
# reference : https://89douner.tistory.com/43 (mini-batch)

my.mini.gradient.U <- function(So, X, alpha, lambda, k, time.sort) {
  
  # alpha : learning rate
  # lambda : penalty parameter
  # k : number of batch
  
  if(lambda < 0) 
    stop("Lambda is negative value. Please insert non-negative value of lambda.")
  
  X <-as.matrix(X)
  
  beta.old <- rep(0.5,ncol(X)) # The initial value of coefficient vector
  beta.old <- beta.old/sqrt(sum(beta.old^2))
  iteration <- c()
  difference <- c()
  
  for (i in 1:500) {
    
    idx <- sample(1:nrow(X), nrow(X), replace=F)
    grad <- 0
    
    
    for (j in 0:(k-1)) {
      batch.idx <- idx[(1:nrow(X))%%k==j]
      batch.So <- So[batch.idx,]
      batch.X <- X[batch.idx,]
      
      grad <- grad + der_U_beta(batch.So, beta.old, batch.X, lambda, time.sort)
      
      cat("( batch ) = ", j+1, "\n")
    }
    
    mean_grad <- grad/k
      
    beta.new <- beta.old - alpha*mean_grad
      
    diff <- abs(C_tilde_beta(So, beta.new, X) - C_tilde_beta(So, beta.old, X))
    
    iteration <- c(iteration, i)
    difference <- c(difference, diff)
      
    cat("( iteration , difference ) = (", i, ",", diff, ")\n")
      
    if (diff < 1E-5) break
      
    beta.old <- beta.new
    }
    

  cat("Algorithm converged...","\n\n")
  
  return(list(beta.new=beta.new, iteration=iteration, difference=difference))
  
}  

res <- my.mini.gradient.U(So, X, alpha=0.02, lambda=0.3, k=5, time.sort)

library(ggplot2)
dat <- data.frame(iteration=res$iteration, difference=res$difference)
ggplot(data=dat, aes(x=iteration, y=difference, group=1)) +
  geom_line() +
  geom_point() +
  ggtitle('Mini-batch gradient descent algorithm (alpha=0.02, lambda=0.3, k=5)') +
  theme(plot.title = element_text(hjust = 0.5,size=12,face='bold'))

# result of beta estimate 
# Cox PH model : (0.193053118, -0.034084486, 0.001723026, -0.003882848, -0.077640942) 
# minimize U(\beta) (alpha=0.02, lambda=0.3, k=5) : (0.424385535, -0.050732304, -0.012636304, -0.007454847, 0.403322836)






#---------------------------------------------------------------------------------------#

library(survival)
So <- with(veteran, Surv(time,status==1))
X <- model.matrix(~factor(trt)+karno+diagtime+age+factor(prior), data=veteran)[,-1]
eta <- predict(coxph(So ~ X))
model1 <- coxph(So ~ X)
risk.score_cox <- predict(object=coxph(So ~ X), newdata=veteran, type="risk") # risk score
beta_hat <- model1[["coefficients"]]
beta <- as.matrix(rep(0.5, ncol(X)))
time.sort <- sort(So[,1]) # sorting survival time
lambda <- seq(0,1,0.01) # parameter for k-fold crossvalidation


# alpha=0.01, lambda=0.1, k=5, number of sorting time=2 (case01)
time.sort_01 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=2)
res01 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=0.1, k=5, time.sort_01)
risk.score_01 <- exp(X%*%res01$beta.new)
res01_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=0.1, time.sort_01)
risk.score_01_1 <- exp(X%*%res01_1$beta.new)

# alpha=0.01, lambda=0.1, k=5, number of sorting time=10 (case02)
time.sort_02 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=10)
res02 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=0.1, k=5, time.sort_02)
risk.score_02 <- exp(X%*%res02$beta.new)
res02_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=0.1, time.sort_02)
risk.score_02_1 <- exp(X%*%res02_1$beta.new)

# alpha=0.01, lambda=0.1, k=5, number of sorting time=100 (case03)
time.sort_03 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=100)
res03 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=0.1, k=5, time.sort_03)
risk.score_03 <- exp(X%*%res03$beta.new)
res03_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=0.1, time.sort_03)
risk.score_03_1 <- exp(X%*%res03_1$beta.new)


# alpha=0.01, lambda=1, k=5, number of sorting time=2 (case04)
time.sort_04 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=2)
res04 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=1, k=5, time.sort_04)
risk.score_04 <- exp(X%*%res04$beta.new)
res04_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=1, time.sort_04)
risk.score_04_1 <- exp(X%*%res04_1$beta.new)

# alpha=0.01, lambda=1, k=5, number of sorting time=10 (case05)
time.sort_05 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=10)
res05 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=1, k=5, time.sort_05)
risk.score_05 <- exp(X%*%res05$beta.new)
res05_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=1, time.sort_05)
risk.score_05_1 <- exp(X%*%res05_1$beta.new)

# alpha=0.01, lambda=1, k=5, number of sorting time=100 (case06)
time.sort_06 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=100)
res06 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=1, k=5, time.sort_06)
risk.score_06 <- exp(X%*%res06$beta.new)
res06_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=1, time.sort_06)
risk.score_06_1 <- exp(X%*%res06_1$beta.new)


# alpha=0.01, lambda=5, k=5, number of sorting time=2 (case07)
time.sort_07 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=2)
res07 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=5, k=5, time.sort_07)
risk.score_07 <- exp(X%*%res07$beta.new)
res07_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=5, time.sort_07)
risk.score_07_1 <- exp(X%*%res07_1$beta.new)

# alpha=0.01, lambda=5, k=5, number of sorting time=10 (case08)
time.sort_08 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=10)
res08 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=5, k=5, time.sort_08)
risk.score_08 <- exp(X%*%res08$beta.new)
res08_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=5, time.sort_08)
risk.score_08_1 <- exp(X%*%res08_1$beta.new)

# alpha=0.01, lambda=5, k=5, number of sorting time=100 (case09)
time.sort_09 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=100)
res09 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=5, k=5, time.sort_09)
risk.score_09 <- exp(X%*%res09$beta.new)
res09_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=5, time.sort_09)
risk.score_09_1 <- exp(X%*%res09_1$beta.new)


# alpha=0.01, lambda=10, k=5, number of sorting time=2 (case10)
time.sort_10 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=2)
res10 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=10, k=5, time.sort_10)
risk.score_10 <- exp(X%*%res10$beta.new)
res10_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=10, time.sort_10)
risk.score_10_1 <- exp(X%*%res10_1$beta.new)

# alpha=0.01, lambda=10, k=5, number of sorting time=10 (case11)
time.sort_11 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=10)
res11 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=10, k=5, time.sort_11)
risk.score_11 <- exp(X%*%res11$beta.new)
res11_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=10, time.sort_11)
risk.score_11_1 <- exp(X%*%res11_1$beta.new)

# alpha=0.01, lambda=10, k=5, number of sorting time=100 (case12)
time.sort_12 <- seq(from=min(So[,1]), to=max(So[,1]), length.out=100)
res12 <- my.mini.gradient.U(So, X, alpha=0.01, lambda=10, k=5, time.sort_12)
risk.score_12 <- exp(X%*%res12$beta.new)
res12_1 <- my.gradient.descent.U(So, X, alpha=0.01, lambda=10, time.sort_12)
risk.score_12_1 <- exp(X%*%res12_1$beta.new)




## Time-dependent AUC


# Loading packages

library(survival) 
library(timeROC)
library(timereg)
library(ggplot2)
library(reshape)
library(plyr)
library(tidyverse)
library(coxed)
library(gridExtra)
library(mgcv)


# Time dependent AUC (Cumulative/Dynamic time-dependent ROC curve)

# mini-batch gradient descent

ROC.00 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_cox,
                    cause=1,weighting="marginal",times=So[,1])
case00 <- ROC.00[["AUC"]]

ROC.td01 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_01,
                    cause=1,weighting="marginal",times=So[,1])
case01 <- ROC.td01[["AUC"]]

ROC.td02 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_02,
                    cause=1,weighting="marginal",times=So[,1])
case02 <- ROC.td02[["AUC"]]

ROC.td03 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_03,
                    cause=1,weighting="marginal",times=So[,1])
case03 <- ROC.td03[["AUC"]]

ROC.td04 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_04,
                    cause=1,weighting="marginal",times=So[,1])
case04 <- ROC.td04[["AUC"]]

ROC.td05 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_05,
                    cause=1,weighting="marginal",times=So[,1])
case05 <- ROC.td05[["AUC"]]

ROC.td06 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_06,
                    cause=1,weighting="marginal",times=So[,1])
case06 <- ROC.td06[["AUC"]]

ROC.td07 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_07,
                    cause=1,weighting="marginal",times=So[,1])
case07 <- ROC.td07[["AUC"]]

ROC.td08 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_08,
                    cause=1,weighting="marginal",times=So[,1])
case08 <- ROC.td08[["AUC"]]

ROC.td09 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_09,
                    cause=1,weighting="marginal",times=So[,1])
case09 <- ROC.td09[["AUC"]]

ROC.td10 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_10,
                    cause=1,weighting="marginal",times=So[,1])
case10 <- ROC.td10[["AUC"]]

ROC.td11 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_11,
                    cause=1,weighting="marginal",times=So[,1])
case11 <- ROC.td11[["AUC"]]

ROC.td12 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_12,
                    cause=1,weighting="marginal",times=So[,1])
case12 <- ROC.td12[["AUC"]]

dt <- data.frame(time=ROC.00[["times"]], case00=case00,
                 case01=case01, case02=case02, case03=case03, case04=case04,
                 case05=case05, case06=case06, case07=case07, case08=case08,
                 case09=case09, case10=case10, case11=case11, case12=case12)
dt1 <- data.frame(time=rep(ROC.00[["times"]],13),
                  case=c(rep("case00",137), rep("case01",137), rep("case02",137),
                         rep("case03",137), rep("case04",137),
                         rep("case05",137), rep("case06",137),
                         rep("case07",137), rep("case08",137),
                         rep("case09",137), rep("case10",137),
                         rep("case11",137), rep("case12",137)),
                  AUC=c(case00, case01, case02, case03, case04, case05, case06,
                        case07, case08, case09, case10, case11, case12))

dt2 <- na.omit(dt1)

aggregate(dt2$AUC, list(dt2$case), FUN=mean, na.action = na.omit)

ggplot(dt1, aes(x=time, y=AUC, group=case, color=case)) +
  geom_line() +
  ggtitle("Time dependent AUC with risk score (mini-batch gradient descent)") +
  xlab("time") + ylab("AUC") +
  theme(plot.title = element_text(hjust = 0.5))

# gradient descent

ROC.td00_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_cox,
                  cause=1,weighting="marginal",times=So[,1])
case00_1 <- ROC.td00_1[["AUC"]]

ROC.td01_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_01_1,
                    cause=1,weighting="marginal",times=So[,1])
case01_1 <- ROC.td01_1[["AUC"]]

ROC.td02_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_02_1,
                    cause=1,weighting="marginal",times=So[,1])
case02_1 <- ROC.td02_1[["AUC"]]

ROC.td03_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_03_1,
                    cause=1,weighting="marginal",times=So[,1])
case03_1 <- ROC.td03_1[["AUC"]]

ROC.td04_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_04_1,
                    cause=1,weighting="marginal",times=So[,1])
case04_1 <- ROC.td04_1[["AUC"]]

ROC.td05_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_05_1,
                    cause=1,weighting="marginal",times=So[,1])
case05_1 <- ROC.td05_1[["AUC"]]

ROC.td06_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_06_1,
                    cause=1,weighting="marginal",times=So[,1])
case06_1 <- ROC.td06_1[["AUC"]]

ROC.td07_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_07_1,
                    cause=1,weighting="marginal",times=So[,1])
case07_1 <- ROC.td07_1[["AUC"]]

ROC.td08_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_08_1,
                    cause=1,weighting="marginal",times=So[,1])
case08_1 <- ROC.td08_1[["AUC"]]

ROC.td09_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_09_1,
                    cause=1,weighting="marginal",times=So[,1])
case09_1 <- ROC.td09_1[["AUC"]]

ROC.td10_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_10_1,
                    cause=1,weighting="marginal",times=So[,1])
case10_1 <- ROC.td10_1[["AUC"]]

ROC.td11_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_11_1,
                    cause=1,weighting="marginal",times=So[,1])
case11_1 <- ROC.td11_1[["AUC"]]

ROC.td12_1 <- timeROC(T=So[,1],delta=So[,2],marker=risk.score_12_1,
                      cause=1,weighting="marginal",times=So[,1])
case12_1 <- ROC.td12_1[["AUC"]]

dt_1 <- data.frame(time=ROC.td00_1[["times"]], case00_1=case00_1,
                 case01=case01_1, case02=case02_1, case03=case03_1, case04=case04_1,
                 case05=case05_1, case06=case06_1, case07=case07_1, case08=case08_1,
                 case09=case09_1, case10=case10_1, case11=case11_1, case12=case12_1)
dt1_1 <- data.frame(time=rep(ROC.td00_1[["times"]],13),
                  case=c(rep("case00_1",137), rep("case01_1",137), rep("case02_1",137),
                         rep("case03_1",137), rep("case04_1",137),
                         rep("case05_1",137), rep("case06_1",137),
                         rep("case07_1",137), rep("case08_1",137),
                         rep("case09_1",137), rep("case10_1",137),
                         rep("case11_1",137), rep("case12_1",137)),
                  AUC=c(case00_1, case01_1, case02_1, case03_1, 
                        case04_1, case05_1, case06_1, case07_1, 
                        case08_1, case09_1, case10_1, case11_1, case12_1))

dt2_1 <- na.omit(dt1_1)

aggregate(dt2_1$AUC, list(dt2_1$case), FUN=mean, na.action = na.omit)

ggplot(dt1_1, aes(x=time, y=AUC, group=case, color=case)) +
  geom_line() +
  ggtitle("Time dependent AUC with risk score (gradient descent)") +
  xlab("time") + ylab("AUC") +
  theme(plot.title = element_text(hjust = 0.5))


# ROC curve

TP1 <- ROC.td[["TP"]][,colnames(ROC.td[["TP"]])=="t=11"]
FP1 <- ROC.td[["FP"]][,colnames(ROC.td[["FP"]])=="t=11"]
TP2 <- ROC.td[["TP"]][,colnames(ROC.td[["TP"]])=="t=25"]
FP2 <- ROC.td[["FP"]][,colnames(ROC.td[["FP"]])=="t=25"]
TP3 <- ROC.td[["TP"]][,colnames(ROC.td[["TP"]])=="t=36"]
FP3 <- ROC.td[["FP"]][,colnames(ROC.td[["FP"]])=="t=36"]
TP4 <- ROC.td[["TP"]][,colnames(ROC.td[["TP"]])=="t=411"]
FP4 <- ROC.td[["FP"]][,colnames(ROC.td[["FP"]])=="t=411"]
TP5 <- ROC.td[["TP"]][,colnames(ROC.td[["TP"]])=="t=587"]
FP5 <- ROC.td[["FP"]][,colnames(ROC.td[["FP"]])=="t=587"]
TP <- c(TP1,TP2,TP3,TP4,TP5)
FP <- c(FP1,FP2,FP3,FP4,FP5)

time.survival <- c(rep("t=11",138),rep("t=25",552),rep("t=36",138),
                   rep("t=411",138),rep("t=587",138))
dt.ROC <- data.frame(times=time.survival,TP=TP,FP=FP)

ggplot(data=dt.ROC)+
  geom_line(mapping=aes(x=FP,
                        y=TP,
                        group=times,
                        color=times),
            size=1) +
  ggtitle("Time dependent ROC curve (veteran data)") +
  xlab("False Positive Rate") + ylab("True Positive Rate") +
  theme(plot.title = element_text(hjust = 0.5))


# Make AUC data

# uncensored
uncensored.density <- density(So[,1][So[,2]==1])
data5_1 <- data.frame(time=uncensored.density[["x"]], 
                      density=uncensored.density[["y"]])
data5_1$density.adjusted <- (max(na.omit(ROC.td[["AUC"]]))-min(na.omit(ROC.td[["AUC"]])))*
  (uncensored.density[["y"]]-min(uncensored.density[["y"]]))/
  (max(uncensored.density[["y"]])-min(uncensored.density[["y"]]))+min(na.omit(ROC.td[["AUC"]]))

# censored
censored.density <- density(So[,1][So[,2]==0])
data5_2 <- data.frame(time=censored.density[["x"]], 
                      density=censored.density[["y"]])
data5_2$density.adjusted <- (max(na.omit(ROC.td[["AUC"]]))-min(na.omit(ROC.td[["AUC"]])))*
  (censored.density[["y"]]-min(censored.density[["y"]]))/
  (max(censored.density[["y"]])-min(censored.density[["y"]]))+min(na.omit(ROC.td[["AUC"]]))

res.AUC <- data.frame(time=ROC.td[["times"]],
                      AUC=ROC.td[["AUC"]])
res.AUC <- melt(data = res.AUC, 
                id.vars = "time", 
                measure.vars = c("AUC"))
res.AUC$status <- sort(So)[,2]


# Time dependent AUC plot with ggplot2

res.AUC$censor <- ifelse(So[,2]==1,
                      1,min(na.omit(ROC.td[["AUC"]])))
ggplot() +
  geom_line(data=res.AUC, aes(x=time, y=value,
                           group=variable, color=variable)) +
  geom_point(data=res.AUC, aes(x=time, y=censor, color=as.factor(So[,2]))) +
  geom_line(data=data5_1, aes(x=time, y=density.adjusted), 
            size=1, colour='black') +
  geom_line(data=data5_2, aes(x=time, y=density.adjusted), 
            size=1, colour='orange') +
  scale_x_continuous(limits = c(0, max(So[,1]))) +
  geom_vline(xintercept = min(So[,1][So[,2]==1]),
             linetype="dotted", color = "blue", size=1) +
  ggtitle("Time dependent AUC with risk score (alpha=0.01, lambda=0.1, k=5, time.sort : 2)") +
  xlab("time") + ylab("AUC") +
  theme(plot.title = element_text(hjust = 0.5))
