
dat.gen <- function(n, p, beta, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  x <- matrix( runif( n*p , -10 , 10 ), nrow=n , ncol=p )
  linpred <- x%*%beta
  prob <- exp(x%*%beta) / (1+exp(x%*%beta))
  y <- apply(as.matrix(prob), c(1,2), rbinom, n=1, size=1)
  dat <- cbind( y , x , linpred )
  colnames( dat ) <- c( "y" , paste0("x", 1:p) , "linpred")
  return( as.data.frame(dat) )
}





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
  K <- exp( -sigma * D^2 )  # Gaussian kernel
  K.train <- K[1:n,1:n]
  K.test <- K[1:n,(n+1):n1]
  
  return( list(K.train=K.train, K.test=K.test, K=K,
               y.train=y.train, y.test=y.test) )
  
}




