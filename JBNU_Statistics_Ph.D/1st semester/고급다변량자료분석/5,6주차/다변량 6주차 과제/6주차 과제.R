# problem 5.4
# (a) 
S <- matrix(c(2.879,10.010,-1.810,10.010,199.788,-5.640,-1.810,-5.640,3.628),
            3,3,byrow=T)
eigen(S)
qf(0.9,3,17)
# (b)
sweet <- read.csv("C:/Users/stat/Desktop/???޴ٺ????ڷ??м?/?ٺ??? 6???? ??��/sweet data_5_4.csv",
                  sep=",", header=T)
qqnorm(sweet$x1) ; qqline(sweet$x1)
qqnorm(sweet$x2) ; qqline(sweet$x2)
qqnorm(sweet$x3) ; qqline(sweet$x3)
plot(sweet)


# problem 5.20
# (a)
library(ggplot2)
bird <- read.csv("C:/Users/stat/Desktop/???޴ٺ????ڷ??м?/?ٺ??? 6???? ??��/bird data_5_20.csv",
                 sep=",", header=T)
ggplot(bird, aes(x1, x2)) +
  geom_point() +
  stat_ellipse()
# (b) 
mean(bird$x1)
mean(bird$x2)
var(bird)
c(193.6222-sqrt((2*44/43)*qf(0.95,2,43))*sqrt(120.6949/45),193.6222+sqrt((2*44/43)*qf(0.95,2,43))*sqrt(120.6949/45))
c(279.7778-sqrt((2*44/43)*qf(0.95,2,43))*sqrt(208.5404/45),279.7778+sqrt((2*44/43)*qf(0.95,2,43))*sqrt(208.5404/45))
c(193.6222-qt(0.9875,44)*sqrt(120.6949/45),193.6222+qt(0.9875,44)*sqrt(120.6949/45))
c(279.7778-qt(0.9875,44)*sqrt(208.5404/45),279.7778+qt(0.9875,44)*sqrt(208.5404/45))
# (c)
qqnorm(bird$x1) ; qqline(bird$x1)
qqnorm(bird$x2) ; qqline(bird$x2)
