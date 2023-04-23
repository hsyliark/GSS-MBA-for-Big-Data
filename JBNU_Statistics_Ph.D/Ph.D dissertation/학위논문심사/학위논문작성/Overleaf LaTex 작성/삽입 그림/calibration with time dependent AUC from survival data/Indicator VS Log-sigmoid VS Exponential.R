x <- seq(-5,5,0.01)
f <- function(x) {ifelse(x>0, 1, 0)}
sig <- function(x) {1/(1+exp(-x))}
g <- function(x) {1+(log(sig(x))/log(2))}
h <- function(x) {1-exp(-x)}
plot(x,f(x),xlim=c(-5,5),ylim=c(-8,3),col="black",lwd=2,lty=1,xlab="x",ylab="f(x) or g(x) or h(x)",
     main="Indicator function VS Log-sigmoid function VS Exponential function")
lines(x,g(x),col="red",lwd=2,lty=2)
lines(x,h(x),col="blue",lwd=2,lty=3)
legend("bottomright",legend=c("Indicator function","Log-sigmoid function","Exponential function"),
       fill=c("black","red","blue"),lwd=c(2,2,2),lty=c(1,2,3),cex=1.2)
