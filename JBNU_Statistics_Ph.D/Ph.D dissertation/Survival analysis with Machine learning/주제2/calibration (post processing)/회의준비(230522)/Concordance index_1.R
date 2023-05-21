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

  





