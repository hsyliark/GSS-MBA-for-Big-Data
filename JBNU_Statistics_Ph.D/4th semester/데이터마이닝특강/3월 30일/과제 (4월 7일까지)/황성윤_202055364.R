### problem 2

auto <- read.csv("C:/Users/stat/Desktop/데이터마이닝특강/과제 (4월 7일까지)/Auto.csv",
                 sep=",",header=T)
summary(auto)
auto$horsepower <- as.numeric(auto$horsepower)
auto <- na.omit(auto)
summary(auto)

auto.lm <- lm(mpg ~ horsepower, data=auto)
summary(auto.lm)

new.auto <- data.frame(horsepower=98)
predict(auto.lm, new=new.auto)
predict(auto.lm, new=new.auto, interval = "confidence")
anova(auto.lm)

plot(auto$horsepower, auto$mpg, main="Plot of mpg vs horsepower",
     xlab="horsepower", ylab="mpg", col="blue", type="p")
abline(auto.lm, col="red", lwd=2)



### problem 3

set.seed(1)
x1 <- runif(100)
x2 <- 0.5*x1+rnorm(100)/10
y <- 2+2*x1+0.3*x2+rnorm(100)

dt <- data.frame(x1=x1, x2=x2, y=y)
plot(dt$x1, dt$x2, xlab="x1", ylab="x2", main="Plot of x1 vs x2")

dt.lm <- lm(y ~ x1+x2, data=dt)
summary(dt.lm)

dt.lm1 <- lm(y ~ x1, data=dt)
summary(dt.lm1)

dt.lm2 <- lm(y ~ x2, data=dt)
summary(dt.lm2)

x1 <- c(x1, 0.1)
x2 <- c(x2, 0.8)
y <- c(y, 6)
dtw <- data.frame(x1=x1, x2=x2, y=y)  

dtw.lm <- lm(y ~ x1+x2, data=dtw)
summary(dtw.lm)

dtw.lm1 <- lm(y ~ x1, data=dtw)
summary(dtw.lm1)

dtw.lm2 <- lm(y ~ x2, data=dtw)
summary(dtw.lm2)

par(mfrow=c(2,2))
plot(dt.lm)
plot(dtw.lm)
par(mfrow=c(1,1))

