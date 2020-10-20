## problem 6.1
data1 <- read.csv("C:/Users/HSY/Desktop/problem 6_1.csv",sep=",",header=T)
x <- data.matrix(data1)
x
(n <- nrow(x))
alpha <- 0.05
p <- 2
((n-1)*p/(n*(n-p))*qf(alpha,p,n-p,lower.tail=F))
(d.mean <- rbind(c(mean(x[,1]-x[,3])),c(mean(x[,2]-x[,4]))))
(diff <- cbind(x[,1]-x[,3],x[,2]-x[,4]))
(Sd <- cov(diff))
solve(Sd)
t(d.mean)%*%solve(Sd)%*%d.mean

## problem 6.17
parity <- read.csv("C:/Users/HSY/Desktop/problem 6_17.csv",sep=",",header=T)
parity
summary(parity)
fit <- manova(cbind(parity$different,parity$same)~parity$treat)
summary.aov(fit)           # univariate ANOVA tables
summary(fit, test="Wilks") # ANOVA table of Wilk's lambda
summary(fit, test="Pillai") # ANOVA table of Pillai's Trace
summary(fit, test="Hotelling-Lawley") # ANOVA table of Hotelling-Lawley Trace
summary(fit, test="Roy") # ANOVA table of Roy's Greatest Root

## problem 6.20
library(ggplot2)
hook <- read.csv("C:/Users/HSY/Desktop/problem 6_20.csv",sep=",",header=T)
summary(hook)
ggplot(hook, aes(x1,x2)) + geom_point(aes(colour = factor(gender)), size = 2)
fit1 <- manova(cbind(hook$x1,hook$x2)~hook$gender)
summary.aov(fit1)           # univariate ANOVA tables
summary(fit1, test="Wilks") # ANOVA table of Wilk's lambda
summary(fit1, test="Pillai") # ANOVA table of Pillai's Trace
summary(fit1, test="Hotelling-Lawley") # ANOVA table of Hotelling-Lawley Trace
summary(fit1, test="Roy") # ANOVA table of Roy's Greatest Root 

## problem 6.25
crude <- read.csv("C:/Users/HSY/Desktop/problem 6_25.csv",sep=",",header=T)
fit2 <- manova(cbind(crude$x1,crude$x2,crude$x3,crude$x4,crude$x5)~crude$zone)
summary.aov(fit2)           # univariate ANOVA tables
summary(fit2, test="Wilks") # ANOVA table of Wilk's lambda
summary(fit2, test="Pillai") # ANOVA table of Pillai's Trace
summary(fit2, test="Hotelling-Lawley") # ANOVA table of Hotelling-Lawley Trace
summary(fit2, test="Roy") # ANOVA table of Roy's Greatest Root 