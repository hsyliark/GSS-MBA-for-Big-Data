x <- seq(-5,5,0.01)
f <- function(x) {ifelse(x>0, 1, 0)}
sig <- function(x) {1/(1+exp(-x))}
g <- function(x) {1+(log(sig(x))/log(2))}
plot(x,f(x),xlim=c(-5,5),ylim=c(-8,3),col="red",lwd=2,lty=1,xlab="x",ylab="f(x) or g(x)",
     main="Indicator function VS Log-sigmoid lower bound")
lines(x,g(x),col="blue",lwd=2,lty=2)
legend("bottomright",legend=c("Indicator function","Log-sigmoid lower bound"),
       fill=c("red","blue"),lwd=c(2,2),lty=c(1,2),cex=1.2)
