theta <- seq(0,6,0.001)
K_theta <- function(theta) {
  pnorm(sqrt(10)*(6-theta)/2-qnorm(0.95))
}
plot(theta,K_theta(theta),type="l")