##############  package
library(randomForest) ##random Forest
library(ipred) ##bagging
library(gbm) ## Boosting
library(xgboost)  ## xgboost  
library(adabag)  ##Adaboosting : boosting
library(rpart)
library(rpart.plot)
library(ggplot2)
library(data.table)


##############  경로지정
setwd("C:\\R-Project\\DAT\\INtroduction SL")


#################################################
##### Boston
#################################################

################ Tree ####################
set.seed(1)

Boston <- read.csv("Boston.csv", 
               stringsAsFactors = T)
dim(dt)
head(dt)


train <- sample(1:nrow(Boston), nrow(Boston) / 2)
Boston_train <- Boston[train,-1]
Boston_test <- Boston[-train,-1]

cart_boston<- rpart(medv ~ ., data=Boston_train)

cart_boston

rpart.plot(cart_boston, 
           main = "Regression using CART")

summary(cart_boston)
cart_boston$variable.importance  ##변수중요도


cart_boston$cptable  ##complexity parameter 
printcp(cart_boston)
plotcp(cart_boston)


#가지치기
prune_cart_boston <- rpart(medv ~ ., data=Boston_train,
                           control = rpart.control(cp = 0.04))


rpart.plot(prune_cart_boston, main = "Regression using CART")



## 예측
yhat <- predict(prune_cart_boston, 
                newdata = Boston_test)
yhat
## type : vector, tree, class, where

## Test MSE
mean((yhat-Boston_test$medv)^2)  ###MSE

## 가지치기 전 tree 에서의 MSE
mean((predict(cart_boston, 
              newdata=Boston_test)-Boston_test$medv)^2)  


################ Bagging ####################
fit.bagg<- ipredbagg(Boston_train$medv, 
                     Boston_train[,-13],
                     nbagg=1000, 
                     coob=T)  
fit.bagg

pred<-predict(fit.bagg, 
              newdata = Boston_test)
mean((pred-Boston_test$medv)^2)  ###MSE


################ Random Forest ####################
rf.boston <- randomForest(medv ~ ., 
                          data = Boston_train,
                          ntree=200,
                          mtry = 4, 
                          importance = TRUE)

rf.boston

yhat.bag <- predict(rf.boston, 
                    newdata = Boston_test)

plot(yhat.bag, Boston_test$medv, pch=16)
abline(0, 1, col='steelblue')

mean((yhat.bag - Boston_test$medv)^2) ##Test MSE

plot(rf.boston$mse, type='l', lwd=2)
plot(rf.boston$rsq, type='l', lwd=2)

rf.boston$importance
rf.boston$importanceSD

importance(rf.boston)
varImpPlot(rf.boston, pch=16)

rf.boston$oob.times

par(mar = c(5.1, 6, 4.1, 2.1))
#c(5.1, 4.1, 4.1, 2.1)## default
par(mfrow=c(1,2))
barplot(sort(importance(rf.boston)[,"%IncMSE"], decreasing = T), 
        horiz = T,las=1,
        col='skyblue',
        xlab = "%IncMSE")

barplot(sort(importance(rf.boston)[,"IncNodePurity"], decreasing = T), 
        horiz = T,las=1,
        col='skyblue',
        xlab = "IncNodePurity")
par(mfrow=c(1,1))



### Bagging
bag.boston <- randomForest(medv ~ ., data = Boston_train, 
                           mtry = 13)
yhat.bag <- predict(bag.boston, newdata = Boston_test)
mean((yhat.bag - Boston_test$medv)^2) ## Test error


### number of mtry
rf.boston <- randomForest(medv ~ ., 
                          data = Boston_train,
                          mtry = 4, 
                          xtest=Boston_test[,-13], 
                          ytest=Boston_test$medv)

mtry_rf <- function(m){
  return(randomForest(medv ~ ., 
                      data = Boston_train,
                      mtry = m)$mse)
}
mtry_rf(1)

tmp_dt <- data.table( num_tree = 1:500,
                      rf_1 = mtry_rf(1),
                      rf_4 = rf.boston$mse,
                      rf_8 = mtry_rf(8),
                      rf_13 = bag.boston$mse,
                      rf_test = rf.boston$test$mse)

melt.tmp <- melt(tmp_dt, id=1)

ggplot(melt.tmp, aes(num_tree, value, col=variable)) + 
  geom_line(lwd=1) +
  labs(y='MSE', col="") + theme_bw()


################ Boosting ####################

boosting_Boston<- gbm(medv~.,
                      data=Boston_train,
                      distribution="gaussian",
                      interaction.depth = 1,
                      n.trees=500)
boosting_Boston
summary(boosting_Boston)

plot(boosting_Boston, i = "rm")
plot(boosting_Boston, i = "lstat")

yhat.boost <- predict(boosting_Boston,
                      newdata = Boston_test, n.trees = 500)

plot(Boston_test$medv, yhat.boost, pch=16)
abline(a=0, b=1, col='blue')

mean((yhat.boost - Boston_test$medv)^2)  ##Test MSE

plot(boosting_Boston$train.error, type='l')
plot(boosting_Boston$oobag.improve,type='l')
abline(h=0)

##CV

boosting_Boston<- gbm(medv~.,data=Boston_train,
                      distribution="gaussian",
                      interaction.depth = 4,
                      cv.folds = 10,
                      n.trees=500)

gbm.perf(boosting_Boston, plot.it = TRUE, 
         oobag.curve = F, overlay = TRUE, method='cv')
gbm.perf(boosting_Boston, plot.it = TRUE, 
         oobag.curve = F, overlay = TRUE, method='OOB')

boosting_Boston$cv.fitted

plot(Boston_train$medv, boosting_Boston$cv.fitted, pch=16)
abline(a=0, b=1, col='blue')



### number of interaction.depth

Boost_ID <- function(d){
  
  tmp_boost_model <- gbm(medv~.,data=Boston_train,
                         distribution="gaussian",
                         interaction.depth = d,
                         n.trees=100)
  
  return(tmp_boost_model$oobag.improve)
  
}

tmp_dt <- data.table( num_tree = 1:100,
                      Boost_1 = Boost_ID(1),
                      Boost_2 = Boost_ID(2),
                      Boost_4 = Boost_ID(4),
                      Boost_12 = Boost_ID(12)
                      )

melt.tmp <- melt(tmp_dt, id=1)

ggplot(melt.tmp, aes(num_tree, value, col=variable)) + 
  geom_line(lwd=2) +
  labs(y='oobag.improve', col="depth")+ theme_bw()



#################################################
##### Heart
#################################################

Heart <- read.csv("Heart.csv", stringsAsFactors = T)
Heart <- na.omit(Heart)
head(Heart)

set.seed (123)
train = sample (1: nrow(Heart), nrow(Heart)*0.7)

Heart_train <- Heart[train,-1]
Heart_test <- Heart[-train,-1]

################ Tree ####################
cart_Heart<- rpart(AHD ~ ., data=Heart_train)
cart_Heart

rpart.plot(cart_Heart, 
           main = "Calssification using CART")

summary(cart_Heart)
cart_Heart$variable.importance  ##변수중요도


cart_Heart$cptable  ##complexity parameter 
printcp(cart_Heart)
plotcp(cart_Heart)


#가지치기
prune_cart_Heart <- rpart(AHD ~ ., data=Heart_train,
                           control = rpart.control(cp = 0.04))


rpart.plot(prune_cart_Heart, main = "Calssification using CART")



## 예측
yhat <- predict(prune_cart_Heart, newdata = Heart_test, type='class')
yhat
## type : vector, tree, class, where

## 오분류율
mis_rate_tree <- mean(yhat!=Heart_test$AHD)  
mis_rate_tree

mean(predict(cart_Heart, 
             newdata=Heart_test, type='class')!=Heart_test$AHD)


################ Bagging ####################
fit.bagg<- ipredbagg(Heart_train$AHD, 
                     Heart_train[, -grep('AHD', names(Heart_train))],
                     nbagg=100, coob=T)  
fit.bagg

pred<-predict(fit.bagg, 
              newdata = Heart_test)
table(Heart_test$AHD, pred)
mis_rate_bagg <- mean(pred!=Heart_test$AHD)    ###test 오분류율 
mis_rate_bagg

################ Random Forest ####################
rf.Heart<- randomForest(AHD ~ ., 
                        data = Heart_train,
                        mtry = 1,
                        importance = TRUE,
                        xtest=Heart_test[,-grep('AHD', names(Heart_train))], 
                        ytest=Heart_test$AHD
                        )
rf.Heart

yhat.rf <- rf.Heart$test$predicted

ls(rf.Heart$test)
rf.Heart$test$votes

table(Heart_test$AHD, yhat.rf)
mis_rate_rf_4 <- mean(yhat.rf!=Heart_test$AHD)    ###오분류율 
mis_rate_rf_4

plot(rf.Heart$err.rate[,1], type='l', 
     ylab='error rate',
     xlab='num_trees')
lines(rf.Heart$test$err.rate[,1], col='red', lty=2)

rf.Heart$importance
rf.Heart$importanceSD

importance(rf.Heart)
varImpPlot(rf.Heart, pch=16)


par(mar = c(5.1, 6, 4.1, 2.1))
#c(5.1, 4.1, 4.1, 2.1)## default
par(mfrow=c(1,2))
barplot(sort(importance(rf.Heart)[,"MeanDecreaseAccuracy"], decreasing = T), 
        horiz = T,las=1,
        col='skyblue',
        xlab = "MeanDecreaseAccuracy")

barplot(sort(importance(rf.Heart)[,"MeanDecreaseGini"], decreasing = T), 
        horiz = T,las=1,
        col='skyblue',
        xlab = "MeanDecreaseGini")
par(mar = c(5.1, 4.1 , 4.1, 2.1))
par(mfrow=c(1,1))

### Bagging
bag.Heart <- randomForest(AHD ~ ., data = Heart_train, 
                          mtry = 13)
yhat.bag <- predict(bag.Heart, newdata = Heart_test)
table(Heart_test$AHD, yhat.bag)
mean(yhat.bag!=Heart_test$AHD)    ###test 오분류율 


#### number of mtry

mtry_rf_heart <- function(m){
  return(randomForest(AHD ~ ., 
                      data = Heart_train,
                      mtry = m)$err.rate[,1])
}

tmp_dt <- data.table( num_tree = 1:500,
                      rf_1 = mtry_rf_heart(1),
                      rf_4 = mtry_rf_heart(4),
                      rf_8 = mtry_rf_heart(8),
                      rf_13 = mtry_rf_heart(13),
                      rf_test = rf.Heart$test$err.rate[,1])

melt.tmp <- melt(tmp_dt, id=1)

ggplot(melt.tmp, aes(num_tree, value, col=variable)) + geom_line(lwd=1) +
  labs(y='error rate', col="") + theme_bw()

################ Boosting ####################
boosting_Heart<- gbm(AHD ~.,data=Heart_train,
                     distribution = 'bernoulli',
                     interaction.depth = 4,
                     n.trees=500)

Heart_train_2 <- Heart_train

Heart_train_2$AHD1 <- as.numeric(Heart_train_2$AHD)-1
Heart_train_2 <- Heart_train_2[,-which(names(Heart_train_2)=='AHD')]

boosting_Heart<- gbm(AHD1 ~.,data=Heart_train_2,
                     distribution = 'bernoulli',
                     interaction.depth = 1,
                     n.trees=500)
boosting_Heart
summary(boosting_Heart)

plot(boosting_Heart, i = "Oldpeak")
plot(boosting_Heart, i = "ChestPain")

yhat.boost_prob <- predict(boosting_Heart,
                           newdata = Heart_test, 
                           n.trees = 500, 
                           type='response')

yhat.boost <- ifelse(yhat.boost_prob > 0.5, 'Yes', 'No')

table(Heart_test$AHD, yhat.boost)
mis_rate_boost <- mean(yhat.boost!=Heart_test$AHD)    ###test 오분류율 
mis_rate_boost

plot(boosting_Heart$train.error, type='l')
plot(boosting_Heart$oobag.improve, type='l')
abline(h=0)

boosting_Heart<- gbm(AHD1 ~.,data=Heart_train_2,
                     distribution = 'bernoulli',
                     interaction.depth = 4,
                     cv.folds = 10,
                     n.trees=500)

gbm.perf(boosting_Heart, plot.it = TRUE, 
         overlay = TRUE, method='cv')
gbm.perf(boosting_Heart, plot.it = TRUE, 
         oobag.curve = T, overlay = TRUE, method='OOB')


### number of interaction.depth

Boost_ID <- function(d){
  
  tmp_boost_model <- gbm(AHD1~.,data=Heart_train_2,
                         distribution="bernoulli",
                         interaction.depth = d,
                         n.trees=50)
  
  return(tmp_boost_model$oobag.improve)
  
}

tmp_dt <- data.table( num_tree = 1:50,
                      Boost_1 = Boost_ID(1),
                      Boost_2 = Boost_ID(2),
                      Boost_4 = Boost_ID(4),
                      Boost_12 = Boost_ID(12)
)

melt.tmp <- melt(tmp_dt, id=1)

ggplot(melt.tmp, aes(num_tree, value, col=variable)) + 
  geom_line(lwd=1) +
  labs(y='oobag.improve', col="depth")+ theme_bw()



####### AdaBoost
fit.boo <- boosting(AHD ~ ., data=Heart_train, 
                    boos=T,
                    mfinal=500)
summary(fit.boo)

ls(fit.boo)

fit.boo$importance

plot(fit.boo$trees[[1]],margin=0.3)
text(fit.boo$trees[[1]], cex=0.8)

#### Predict
pred  <- predict(fit.boo, 
                 newdata = Heart_test)
table(pred$class,Heart_test$AHD)

## 오분류율
mean(pred$class!= Heart_test$AHD)
