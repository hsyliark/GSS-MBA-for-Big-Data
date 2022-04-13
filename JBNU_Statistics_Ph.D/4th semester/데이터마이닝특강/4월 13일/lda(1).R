##########################################
########## Classification : LDA ########## 
##########################################
library(MASS) ## for lda
library(ROCR) ## for ROC

setwd("C:\\R-Project\\DAT\\INtroduction SL")

##########################################
#############  iris
iris
pairs(iris[,-5], col=c('steelblue', 'darkorange', 'darkgreen')[iris$Species])


######## LDA
lda.fit <- lda(Species ~ ., data = iris)

lda.fit

lda.fit$scaling
lda.fit$means


#### predict 
ls(predict(lda.fit))

predict(lda.fit)$class
round(predict(lda.fit)$posterior,4)
head(predict(lda.fit)$x)

plot(lda.fit, abbrev=-1)
points(predict(lda.fit)$x, 
       pch=c(15,16,17)[predict(lda.fit)$class], 
       col=c('steelblue', 'darkorange', 'darkgreen')[iris$Species])
points(predict(lda.fit)$x[iris$Species != predict(lda.fit)$class,], pch=1, col='red', cex=3)
legend(x="top", legend = levels(iris$Species), 
       col=c('steelblue', 'darkorange', 'darkgreen'), 
       pch=c(15,16,17),
       bty = "n", ncol=3)

plot(lda.fit, dimen=1)

#### misclassification rate
table(iris$Species, predict(lda.fit)$class)
prop.table(table(iris$Species, predict(lda.fit)$class))

mean(iris$Species == predict(lda.fit)$class) ## 정분류율
mean(iris$Species != predict(lda.fit)$class) ## 오분류율


######## QDA
qda.fit <- qda(Species ~ ., data = iris)

qda.fit

#### predict 
ls(predict(qda.fit))

predict(qda.fit)$class
round(predict(qda.fit)$posterior,4)


#### misclassification rate
table(iris$Species, predict(qda.fit)$class)

mean(iris$Species == predict(qda.fit)$class) ## 정분류율
mean(iris$Species != predict(qda.fit)$class) ## 오분류율


##########################################
######## LDA vs. Logistic

sub_iris <- iris[iris$Species != 'versicolor',-(1:2)]
sub_iris <- droplevels(sub_iris)

plot(sub_iris$Petal.Length, sub_iris$Petal.Width, 
     col=c('darkorange', 'steelblue')[sub_iris$Species], pch=16,
     xlab = 'Petal.Length',
     ylab = 'Petal.Width')

glm.fits <- glm(
  Species ~ Petal.Length + Petal.Width ,
  data = sub_iris, 
  family = binomial
)
summary(glm.fits)

contrasts(sub_iris$Species)
round(predict(glm.fits, type='response'),4)

lda.fit <- lda(Species ~ ., data = sub_iris)
plot(lda.fit)
round(predict(lda.fit)$posterior,4)



##########
sub_iris <- iris[iris$Species != 'setosa',-(1:2)]
sub_iris <- droplevels(sub_iris)

plot(sub_iris$Petal.Length, sub_iris$Petal.Width, 
     col=c('darkorange', 'steelblue')[sub_iris$Species], pch=16,
     xlab = 'Petal.Length',
     ylab = 'Petal.Width')

glm.fits <- glm(
  Species ~ Petal.Length + Petal.Width ,
  data = sub_iris, 
  family = binomial
)
summary(glm.fits)

contrasts(sub_iris$Species)
round(predict(glm.fits, type='response'),4)
predicted_class <- ifelse(predict(glm.fits, type='response') > 0.5, 'virginica', 'versicolor')
table(sub_iris$Species, predicted_class)
mean(sub_iris$Species !=  predicted_class) ## 오분류율


lda.fit <- lda(Species ~ ., data = sub_iris)
lda.fit
plot(lda.fit)
round(predict(lda.fit)$posterior,4)

qda.fit <- qda(Species ~ ., data = sub_iris)
qda.fit
round(predict(qda.fit)$posterior,4)

table(sub_iris$Species, predict(lda.fit)$class)
table(sub_iris$Species, predict(qda.fit)$class)

mean(sub_iris$Species != predict(lda.fit)$class) ## 오분류율
mean(sub_iris$Species != predict(qda.fit)$class) ## 오분류율


##########################################
dt <- read.csv("Default.csv", stringsAsFactors = T)

lda.fit <- lda(default ~ balance, data = dt)
lda.pred <- predict(lda.fit)

lda.pred$class

addmargins(table( dt$default, lda.pred$class))


##### ROC
head(lda.pred$posterior)

pred_ROCR <- prediction(lda.pred$posterior[,2], dt$default)
roc_ROCR <- performance(pred_ROCR, measure = "tpr", x.measure = "fpr")
plot(roc_ROCR, main = "ROC curve", col='steelblue')
abline(a = 0, b = 1, lty=2, col='grey')


### logistic
glm.fits <- glm(
  default ~ balance,
  data = dt, 
  family = binomial
)

predict(glm.fits, type='response')

pred_ROCR_logsitc <- prediction(predict(glm.fits, type='response'), dt$default)
roc_ROCR_logsitc <- performance(pred_ROCR_logsitc, measure = "tpr", x.measure = "fpr")
plot(roc_ROCR_logsitc, main = "ROC curve ", col='steelblue')
abline(a = 0, b = 1, lty=2, col='grey')


par(mfrow=c(1,2))
plot(roc_ROCR, main = "ROC curve : LDA", col='steelblue')
abline(a = 0, b = 1, lty=2, col='grey')

plot(roc_ROCR_logsitc, main = "ROC curve : Logistic", col='steelblue')
abline(a = 0, b = 1, lty=2, col='grey')

##### AUC
auc_ROCR <- performance(pred_ROCR, measure = "auc")
auc_ROCR@y.values[[1]]

auc_ROCR_logsitc <- performance(pred_ROCR_logsitc, measure = "auc")
auc_ROCR_logsitc@y.values[[1]]


