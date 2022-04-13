
################## 1 ##################
set.seed(123)
x <- rnorm(1000)

y1 <- 2*x + rnorm(1000) 
y2 <- 0.5*x + 0.1*x^2 + rnorm(1000)
y3 <- 0.5*x + 0.1*x^2 + 0.05*x^3 + 0.1*x^4 + rnorm(1000,0)


## (a) real : y1
model1 <- lm(y1 ~ x)
summary(model1)

model2 <- lm(y1 ~ x + I(x^2) + I(x^3))
summary(model2)

paste0("model1 : SSE = ", round(sum((y1 - model1$fitted.values)^2),3))
paste0("model2 : SSE = ", round(sum((y1 - model2$fitted.values)^2),3))

## (b) real : y2
model1 <- lm(y2 ~ x)
summary(model1)

model2 <- lm(y2 ~ x + I(x^2) + I(x^3))
summary(model2)

paste0("model1 : SSE = ", round(sum((y2 - model1$fitted.values)^2),3))
paste0("model2 : SSE = ", round(sum((y2 - model2$fitted.values)^2),3))

## (b) real : y3
model1 <- lm(y3 ~ x)
summary(model1)

model2 <- lm(y3 ~ x + I(x^2) + I(x^3))
summary(model2)

paste0("model1 : SSE = ", round(sum((y3 - model1$fitted.values)^2),3))
paste0("model2 : SSE = ", round(sum((y3 - model2$fitted.values)^2),3))


################## 2 ##################
setwd("C:\\R-Project\\DAT\\INtroduction SL")

dt <- read.csv("Auto.csv", na.strings = "?", stringsAsFactors = T)
dim(dt)

dt <- na.omit(dt)
dim(dt)
###

dt_sub <- dt[,c('mpg', 'horsepower')]

dt_fit <- dt_sub[dt_sub$horsepower != '?', ]
dt_fit$horsepower <- as.numeric(dt_fit$horsepower)

##(a)
model <- lm(mpg ~ horsepower, data = dt_fit)
summary(model)

sqrt(summary(model)$r.squared)  #r^2 = R^2

predict(model, 
        newdata = data.frame(horsepower = 98), 
        interval = "confidence", level = 0.95)

##(b)
par(mfrow=c(1,1))
plot(dt_fit$horsepower, dt_fit$mpg, 
     pch=16, col='steelblue', xlab='Horsepower', ylab='MPG')
abline(model, lty=2, col='red', lwd=2)

################## 3 ##################


##(a)
set.seed(1)
x1 <- runif(100)
x2 <- 0.5*x1 + rnorm(100)/10
y <- 2 + 2*x1 + 0.3*x2 + rnorm(100)

##(b)
plot(x1, x2, pch=16, col='darkorange')
cor.test(x1, x2, alternative = 'two.sided')

##(c)
model12 <- lm(y ~ x1 + x2)
summary(model12)

##(d)
model1 <- lm(y ~ x1)
summary(model1)

##(e)
model2 <- lm(y ~ x2)
summary(model2)

##(g)
x11 <- c(x1, 0.1)
x22 <- c(x2, 0.8)
y1 <- c(y, 6)

##(b) new
plot(x1, x2, pch=16, col='darkorange', ylim=c(0,1))
points(0.1,0.8, pch=17, cex=2, col='red')

##(c) new
model12_new <- lm(y1 ~ x11 + x22)
summary(model12_new)
summary(model12)
influence.measures(model12_new)

##(d) new
model1_new <- lm(y1 ~ x11)
summary(model1_new)
summary(model1)
influence.measures(model1_new)

##(e) new
model2_new <- lm(y1 ~ x22)
summary(model2_new)
summary(model2)
influence.measures(model2_new)

##
par(mfrow=c(1,2))

plot(x1, y, pch=16, col='darkorange')
points(0.1,6, pch=17, cex=2, col='red')
abline(model1, lty=2, col='blue')
abline(model1_new, lty=2, col='red')


plot(x2, y, pch=16, col='darkorange', xlim=c(0,0.8))
points(0.8,6, pch=17, cex=2, col='red')
abline(model2, lty=2, col='blue')
abline(model2_new, lty=2, col='red')

##residual
par(mfrow=c(1,3))
plot(model12_new$residuals)
abline(h=0, lty=2)
points(101, model12_new$residuals[101], pch=16, col='red')
plot(model1_new$residuals)
abline(h=0, lty=2)
points(101, model1_new$residuals[101], pch=16, col='red')
plot(model2_new$residuals)
abline(h=0, lty=2)
points(101, model2_new$residuals[101], pch=16, col='red')


