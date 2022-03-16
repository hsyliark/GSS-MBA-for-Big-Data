##############  경로지정
setwd("C:\\R-Project\\DAT\\INtroduction SL")



####################################################
####################################################

adv.data <- read.csv('Advertising.csv', stringsAsFactors = F)
head(adv.data)
adv.data$X <- NULL
head(adv.data)

##### 산점도 (y vs. x1, y vs. x2, y vs. x3)
par(mfrow=c(1,3))
plot(adv.data$TV, adv.data$sales, pch=16, col='steelblue', 
     lwd=2, cex.lab=1.5,
     xlab="TV", ylab="Sales")
abline(lm(sales ~ TV, data= adv.data), col='red', lwd=2)
plot(adv.data$radio, adv.data$sales, pch=16, col='steelblue', 
     lwd=2, cex.lab=1.5,
     xlab="Radio", ylab="Sales")
abline(lm(sales ~ radio, data= adv.data), col='red', lwd=2)
plot(adv.data$newspaper, adv.data$sales, pch=16, col='steelblue', 
     lwd=2, cex.lab=1.5,
     xlab="Newspaper", ylab="Sales")
abline(lm(sales ~ newspaper, data= adv.data), col='red', lwd=2)

##### 산점도 행렬

pairs(adv.data, pch=16, col='steelblue')


panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(adv.data, pch=16, col='steelblue', upper.panel = panel.cor)

#### 단순회귀모형

reg_model <- lm(sales ~ TV, data= adv.data)
reg_model
summary(reg_model)


coef(reg_model) ## 회귀계수
confint(reg_model)  ##95% 신뢰구간 

#y의 추정값
head(reg_model$fitted.values, 10)
head(predict(reg_model), 10)
head(fitted(reg_model), 10)
head(fitted.values(reg_model), 10)
y_hat <- fitted(reg_model)

#잔차
adv.data$sales - reg_model$fitted.values
reg_model$residuals
resid(reg_model)
residual <- resid(reg_model)
s_residual <- rstandard(reg_model)  ## 내적으로 스튜던트화된 잔차

par(mfrow=c(1,2))
plot(y_hat, residual, pch=16, col='steelblue')
abline(h=0, lty=2)

plot(y_hat, s_residual, pch=16, col='steelblue')
abline(h=0, lty=2)


#### 중회귀모형
reg_model <- lm(sales ~ ., data= adv.data)
reg_model <- lm(sales ~ TV + radio + newspaper, data= adv.data)

a <- summary(reg_model)
a$r.squared
a$adj.r.squared
summary(reg_model)$r.squared
summary(reg_model)$adj.r.squared


coef(reg_model) ## 회귀계수
confint(reg_model)  ##95% 신뢰구간 

#y의 추정값
reg_model$fitted.values  
predict(reg_model)
fitted(reg_model)
fitted.values(reg_model)
y_hat <- fitted(reg_model)

#잔차
adv.data$sales - reg_model$fitted.values
resid(reg_model)
residual <- resid(reg_model)
s_residual <- rstandard(reg_model)  ## 내적으로 스튜던트화된 잔차

par(mfrow=c(1,2))
plot(y_hat, residual, pch=16, col='steelblue')
abline(h=0, lty=2)

plot(y_hat, s_residual, pch=16, col='steelblue')
abline(h=0, lty=2)
out <- which.max(abs(s_residual))
text (y_hat[out], s_residual[out], out,adj = c(0,0))


###### 변수선택
###Forward - AIC
model_start = lm(sales ~ 1, data = adv.data)
model_forward = step(
  model_start, 
  scope = sales ~ TV + radio + newspaper, 
  # trace = F,
  direction = "forward")
summary(model_forward)

model_forward = step(
  model_start, 
  scope = reg_model$call, 
  # trace = F,
  direction = "forward")
summary(model_forward)


###Backward - AIC
model_back = step(reg_model, direction = "backward")
summary(model_back)


###Step - AIC
model_step = step(
  model_start, 
  scope = sales ~ TV + radio + newspaper, 
  direction = "both")
summary(model_step)




##########################################

###
head(Carseats)
###

lm.fit <- lm(Sales ~ ., 
             data = Carseats)
summary(lm.fit)

contrasts(Carseats$ShelveLoc)


# 교호작용
lm.fit_inter <- lm(Sales ~ . + Income:Advertising + Price:Age, 
             data = Carseats)
summary(lm.fit_inter)
###



















