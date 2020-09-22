install.packages("tidyverse")
install.packages("gridExtra")
install.packages("MASS")
install.packages("ggplot2")
library(tidyverse)
library(gridExtra)
library(MASS)
library(ggplot2)

## Problem1

n <- 100
alpha <- 2
lambda <- 3

data <- rweibull(n, shape=alpha, scale=1/lambda)
exp_H <- tibble(t= sort(data), St=1-pweibull(t, shape=alpha, scale=1/lambda))

ggplot(exp_H, aes(x=log(t)))+
  geom_line(aes(y=log(-log(St))))+
  xlab("log t")+ylab("log[-log S(t)]")+
  labs(title="Graph of (log t,log[-log S(t)])")

(wei_e <- fitdistr(data, 'weibull'))

## Problem2

alpha <- c(3,7,11)
lambda <- 1
t <- t <- seq(0, 5, length=200) 
gamma_d <- tibble(t, f1=dgamma(t, shape=alpha[1], scale=1/lambda), s1=1-pgamma(t, shape=alpha[1], scale=1/lambda), h1=f1/s1,
                f2=dgamma(t, shape=alpha[2], scale=1/lambda), s2=1-pgamma(t, shape=alpha[2], scale=1/lambda), h2=f2/s2,
                f3=dgamma(t, shape=alpha[3], scale=1/lambda), s3=1-pgamma(t, shape=alpha[3], scale=1/lambda), h3=f3/s3)
g1 <- ggplot(gamma_d, aes(x=t))+
  geom_line(aes(y=f1, color="alpha=3"))+
  geom_line(aes(y=f2, color="alpha=7"))+
  geom_line(aes(y=f3, color="alpha=11"))+
  xlab("t")+ylab("f(t)")+
  labs(title="Gamma density", color=" ")
g2 <-  ggplot(gamma_d, aes(x=t))+
  geom_line(aes(y=s1, color="alpha=3"))+
  geom_line(aes(y=s2, color="alpha=7"))+
  geom_line(aes(y=s3, color="alpha=11"))+
  xlab("t")+ylab("S(t)")+
  labs(title="Gamma survival function", color=" ")
g3 <- ggplot(gamma_d, aes(x=t))+
  geom_line(aes(y=h1, color="alpha=3"))+
  geom_line(aes(y=h2, color="alpha=7"))+
  geom_line(aes(y=h3, color="alpha=11"))+
  xlab("t")+ylab("h(t)")+
  labs(title="Gamma hazard function", color=" ")
grid.arrange(g1,g2,g3)


