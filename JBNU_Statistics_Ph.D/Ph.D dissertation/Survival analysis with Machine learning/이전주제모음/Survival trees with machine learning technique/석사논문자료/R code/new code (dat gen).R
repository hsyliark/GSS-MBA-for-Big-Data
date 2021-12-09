
dat.gen <- function(n, p, beta, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  x <- matrix( runif( n*p , -2 , 2 ), nrow=n , ncol=p )
  linpred <- colSums(apply(x , 1 , function(tmp) tmp^(1:p)) * beta)
  prob <- exp(linpred) / (1+exp(linpred))
  y <- apply(as.matrix(prob), c(1,2), rbinom, n=1, size=1)
  dat <- cbind( y , x )
  colnames( dat ) <- c( "y" , paste0("x", 1:p) )
  return( as.data.frame(dat) )
}


