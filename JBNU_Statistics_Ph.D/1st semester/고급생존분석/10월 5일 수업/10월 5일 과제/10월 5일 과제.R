library(ggplot2)
library(gridExtra)
t <- seq(0,10,0.001)
N1 <- function(t) {
  res <- ifelse(t>=2, 1, 0)
  return(res)}
N2 <- function(t) {
  res <- 0
  return(res)}
N3 <- function(t) {
  res <- ifelse(t>=6, 1, 0)
  return(res)}
Y1 <- function(t) {
  res <- ifelse(t<=2, 1, 0)
  return(res)}
Y2 <- function(t) {
  res <- ifelse(t<=8, 1, 0)
  return(res)}
Y3 <- function(t) {
  res <- ifelse(t<=6, 1, 0)
  return(res)}
n1 <- c()
n2 <- c()
n3 <- c()
y1 <- c()
y2 <- c()
y3 <- c()
for (i in 1:length(t)) {
  n1 <- c(n1,N1(t[i]))
  n2 <- c(n2,N2(t[i]))
  n3 <- c(n3,N3(t[i]))
  y1 <- c(y1,Y1(t[i]))
  y2 <- c(y2,Y2(t[i]))
  y3 <- c(y3,Y3(t[i]))
}
my_survival_data1 <- data.frame(t=t, n1=n1, n2=n2, n3=n3, y1=y1, y2=y2, y3=y3)
g1 <- ggplot(my_survival_data1, aes(x=t))+
  geom_line(aes(y=n1))+
  xlab("t")+ylab("N1(t)")+
  labs(title="Counting process N1(t)")
g2 <- ggplot(my_survival_data1, aes(x=t))+
  geom_line(aes(y=n2))+
  xlab("t")+ylab("N2(t)")+
  labs(title="Counting process N2(t)")
g3 <- ggplot(my_survival_data1, aes(x=t))+
  geom_line(aes(y=n3))+
  xlab("t")+ylab("N3(t)")+
  labs(title="Counting process N3(t)")
g4 <- ggplot(my_survival_data1, aes(x=t))+
  geom_line(aes(y=y1))+
  xlab("t")+ylab("Y1(t)")+
  labs(title="Counting process Y1(t)")
g5 <- ggplot(my_survival_data1, aes(x=t))+
  geom_line(aes(y=y2))+
  xlab("t")+ylab("Y2(t)")+
  labs(title="Counting process Y2(t)")
g6 <- ggplot(my_survival_data1, aes(x=t))+
  geom_line(aes(y=y3))+
  xlab("t")+ylab("Y3(t)")+
  labs(title="Counting process Y3(t)")
grid.arrange(g1,g4,g2,g5,g3,g6, nrow=3, ncol=2)

t <- seq(0,10,0.001)
Nt <- function(t) {
  res <- ifelse(t>=2, 1, 0)+0+ifelse(t>=6, 1, 0)
  return(res)}
Yt <- function(t) {
  res <- ifelse(t<=2, 1, 0)+ifelse(t<=8, 1, 0)+ifelse(t<=6, 1, 0)
  return(res)}
nt <- c()
yt <- c()
for (i in 1:length(t)) {
  nt <- c(nt,Nt(t[i]))
  yt <- c(yt,Yt(t[i]))
}
my_survival_data2 <- data.frame(t=t, nt=nt, yt=yt)
g11 <- ggplot(my_survival_data2, aes(x=t))+
  geom_line(aes(y=nt))+
  xlab("t")+ylab("N(t)")+
  labs(title="Counting process N(t)")
g12 <- ggplot(my_survival_data2, aes(x=t))+
  geom_line(aes(y=yt))+
  xlab("t")+ylab("Y(t)")+
  labs(title="Counting process Y(t)")
grid.arrange(g11,g12, nrow=2, ncol=1)
