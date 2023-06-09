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
  for ( k in 2:nrow(X) ) {
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
  for ( k in 2:nrow(X) ) {
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
  
  if(lambda <= 0) 
    stop("Lambda is non-positive value. Please insert positive value of lambda.")
  
  X <-as.matrix(X)
  
  beta.old <- rep(0.5,ncol(X)) # The initial value of coefficient vector
  beta.old <- beta.old/sqrt(sum(beta.old^2))
  
  for (i in 1:100) {
    
    
    beta.new <- beta.old - alpha*der_U_beta(So, beta.old, X, lambda, time.sort)
    
    diff <- abs(C_tilde_beta(So, beta.new, X) - C_tilde_beta(So, beta.old, X))  
    
    # diff <- sqrt(sum((beta.new - beta.old)^2))/sqrt(sum(beta.old^2))
    
    cat("( iteration , difference ) = (", i, ",", diff, ")\n")
    
    if (diff < 1E-6) break
    
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
  
  if(lambda <= 0) 
    stop("Lambda is non-positive value. Please insert positive value of lambda.")
  
  X <-as.matrix(X)
  
  beta.old <- rep(0.5,ncol(X)) # The initial value of coefficient vector
  beta.old <- beta.old/sqrt(sum(beta.old^2))
  iteration <- c()
  difference <- c()
  
  for (i in 1:100) {
    
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
      
    if (diff < 1E-6) break
      
    beta.old <- beta.new
    }
    

  cat("Algorithm converged...","\n\n")
  
  return(list(beta.new=beta.new, iteration=iteration, difference=difference))
  
}  

res <- my.mini.gradient.U(So, X, alpha=0.02, lambda=0.3, k=5, time.sort)

library(ggplot2)
dat <- data.frame(iteration=res$iteration, difference=res$difference)
ggplot(data=dat, aes(x=iteration, y=difference, group=1)) +
  geom_line()+
  geom_point()+
  ggtitle('Mini-batch gradient descent algorithm (alpha=0.02, lambda=0.3, k=5)')+
  theme(plot.title = element_text(hjust = 0.5,size=12,face='bold'))

# result of beta estimate 
# Cox PH model : (0.193053118, -0.034084486, 0.001723026, -0.003882848, -0.077640942) 
# minimize U(\beta) (alpha=0.02, lambda=0.3, k=5) : (0.424385535, -0.050732304, -0.012636304, -0.007454847, 0.403322836)

