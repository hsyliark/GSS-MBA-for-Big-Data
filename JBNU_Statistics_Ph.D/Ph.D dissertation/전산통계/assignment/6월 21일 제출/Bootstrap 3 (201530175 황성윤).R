### Bootstrap regression




## 등분산모형


B <- 1000 ; M <- 1000
set.seed(1)
n <- 10
BT <- matrix(0,M,2)

x <- matrix(runif(n*M,1,2),M,n)
c22 <- 1/rowSums((x-rowMeans(x))^2)
ep <- matrix(rnorm(n*M,mean=0,sd=1),M,n)
y <- 2*x+ep
sxx <- rowSums(x^2)-n*rowMeans(x)^2
sxy <- rowSums(x*y)-n*rowMeans(x)*rowMeans(y)
hb <- cbind(rowMeans(y)-(sxy/sxx)*rowMeans(x),sxy/sxx)
e <- y-hb[,1]-x*hb[,2]
hsigma <- sqrt(rowSums(e^2)/(n-2))
Bm=Bmp=matrix(0,M,B)

for (i in 1:B) {
  hy <- hb[,1]+hb[,2]*x+t(apply(e,1,sample,replace=T))
  sxyb <- rowSums(x*hy)-n*rowMeans(x)*rowMeans(hy)
  hbb <- cbind(rowMeans(hy)-(sxyb/sxx)*rowMeans(x),sxyb/sxx)
  he <- hy-hbb[,1]-hbb[,2]*x
  Bm[,i] <- (hbb[,2] - hb[,2])/(sqrt(c22*rowSums(he^2)/(n-2))) 
}

# Bootstrap t C.I. (residual)
BT <- cbind(hb[,2]-hsigma*sqrt(c22)*apply(Bm,1,quantile,probs=0.975),
            hb[,2]-hsigma*sqrt(c22)*apply(Bm,1,quantile,probs=0.025)) 


for (i in 1:B) {
  hi <- sample(1:n,replace=T)
  hx <- x[,hi]
  hy <- y[,hi]
  sxxb <- rowSums(hx*hx)-n*rowMeans(hx)*rowMeans(hx)
  sxyb <- rowSums(hx*hy)-n*rowMeans(hx)*rowMeans(hy)
  hbb <- cbind(rowMeans(hy)-(sxyb/sxxb)*rowMeans(hx),sxyb/sxxb)
  he <- hy-hbb[,1]-hbb[,2]*hx
  c22b <- 1/rowSums((hx-rowMeans(hx))^2)
  Bmp[,i] <- (hbb[,2] - hb[,2])/( sqrt(c22b*rowSums(he^2)/(n-2))) 
}

# Bootstrap t C.I. (pair)
BTp <- cbind(hb[,2]-hsigma*sqrt(c22)*apply(Bmp,1,quantile,probs=0.975),
             hb[,2]-hsigma*sqrt(c22)*apply(Bmp,1,quantile,probs=0.025)) 

# Asymp. Normal C.I.
AN <- cbind(hb[,2]-qnorm(0.975)*sqrt(hsigma^2*c22),hb[,2]+qnorm(0.975)*sqrt(hsigma^2*c22)) 

# average lengths of intervals
c(diff(colMeans(BT)),diff(colMeans(BTp)),diff(colMeans(AN)))
# coverage proportions
c(mean(apply(BT-2,1,prod)<0),mean(apply(BTp-2,1,prod)<0),mean(apply(AN-2,1,prod)<0))




## 이분산모형


B <- 1000 ; M <- 1000
set.seed(1)
n <- 10
BT <- matrix(0,M,2)

x <- matrix(runif(n*M,1,2),M,n)
c22 <- 1/rowSums((x-rowMeans(x))^2)
ep <- matrix(rnorm(n*M,mean=0,sd=sqrt((x^2)*(exp(x)-exp(1)+1))),M,n)
y <- 2*x+ep
sxx <- rowSums(x^2)-n*rowMeans(x)^2
sxy <- rowSums(x*y)-n*rowMeans(x)*rowMeans(y)
hb <- cbind(rowMeans(y)-(sxy/sxx)*rowMeans(x),sxy/sxx)
e <- y-hb[,1]-x*hb[,2]
hsigma <- sqrt(rowSums(e^2)/(n-2))
Bm=Bmp=matrix(0,M,B)

for (i in 1:B) {
  hy <- hb[,1]+hb[,2]*x+t(apply(e,1,sample,replace=T))
  sxyb <- rowSums(x*hy)-n*rowMeans(x)*rowMeans(hy)
  hbb <- cbind(rowMeans(hy)-(sxyb/sxx)*rowMeans(x),sxyb/sxx)
  he <- hy-hbb[,1]-hbb[,2]*x
  Bm[,i] <- (hbb[,2] - hb[,2])/(sqrt(c22*rowSums(he^2)/(n-2))) 
}

# Bootstrap t C.I. (residual)
BT <- cbind(hb[,2]-hsigma*sqrt(c22)*apply(Bm,1,quantile,probs=0.975),
            hb[,2]-hsigma*sqrt(c22)*apply(Bm,1,quantile,probs=0.025)) 


for (i in 1:B) {
  hi <- sample(1:n,replace=T)
  hx <- x[,hi]
  hy <- y[,hi]
  sxxb <- rowSums(hx*hx)-n*rowMeans(hx)*rowMeans(hx)
  sxyb <- rowSums(hx*hy)-n*rowMeans(hx)*rowMeans(hy)
  hbb <- cbind(rowMeans(hy)-(sxyb/sxxb)*rowMeans(hx),sxyb/sxxb)
  he <- hy-hbb[,1]-hbb[,2]*hx
  c22b <- 1/rowSums((hx-rowMeans(hx))^2)
  Bmp[,i] <- (hbb[,2] - hb[,2])/( sqrt(c22b*rowSums(he^2)/(n-2))) 
}

# Bootstrap t C.I. (pair)
BTp <- cbind(hb[,2]-hsigma*sqrt(c22)*apply(Bmp,1,quantile,probs=0.975),
             hb[,2]-hsigma*sqrt(c22)*apply(Bmp,1,quantile,probs=0.025)) 

# Asymp. Normal C.I.
AN <- cbind(hb[,2]-qnorm(0.975)*sqrt(hsigma^2*c22),hb[,2]+qnorm(0.975)*sqrt(hsigma^2*c22)) 

# average lengths of intervals
c(diff(colMeans(BT)),diff(colMeans(BTp)),diff(colMeans(AN)))
# coverage proportions
c(mean(apply(BT-2,1,prod)<0),mean(apply(BTp-2,1,prod)<0),mean(apply(AN-2,1,prod)<0))
