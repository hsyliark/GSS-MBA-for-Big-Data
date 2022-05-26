######################################
########## HW03  #####################
######################################
library(rpart)  ##Recursive Partitioning and Regression Trees
library(rpart.plot)
library(ipred) ##bagging
library(randomForest) ##random Forest
library(data.table)
library(ISLR2)
library(ggplot2)

setwd('C:\\R-Project\\DAT\\Introduction SL\\')


######################################
#### 3 
######################################

dt_carset <- read.csv('Carseats.csv')
head(dt_carset)

#### (a)
set.seed(2022)
train_id <- sample(1:nrow(dt_carset), nrow(dt_carset)*0.6)

train_dt <- dt_carset[train_id,]
test_dt <- dt_carset[-train_id,]

#### (b)
cart_car <- rpart(Sales  ~ ., data=train_dt)
summary(cart_car)

rpart.plot(cart_car, 
           main = "Classification using CART")

hat_y <- predict(cart_car, newdata = test_dt)
test_MSE <- mean((hat_y-test_dt$Sales)^2)
test_MSE

#### (c)
cart_car$cptable
plotcp(cart_car)

#가지치기
prune_cart_car <- rpart(Sales ~ ., data=train_dt,
                        control = rpart.control(cp = 0.05))
rpart.plot(prune_cart_car, 
           main = "Classification using CART")

test_MSE <- mean((predict(prune_cart_car, 
                          newdata=test_dt)-test_dt$Sales)^2)
test_MSE

#### (d)
fit.bagg_car<- ipredbagg(train_dt$Sales, 
                     train_dt[,-1],
                     nbagg=1000, 
                     coob=T)  
fit.bagg_car

hat_y_bag<-predict(fit.bagg_car, 
                   newdata = test_dt)
mean((hat_y_bag-test_dt$Sales)^2)  ###MSE

#### (e)
fit.rf_car <- randomForest(Sales ~ ., 
                          data = train_dt,
                          ntree= 200,
                          mtry = 3, ## m ~ p/3
                          # xtest=test_dt[,-1], 
                          # ytest=test_dt$Sales,
                          importance = TRUE)

hat_y_rf <-predict(fit.rf_car, 
                   newdata = test_dt)

mean((hat_y_rf-test_dt$Sales)^2) 

fit.rf_car$importance
varImpPlot(fit.rf_car, pch=16)

mtry_rf <- function(m){
  set.seed(m)
  return(randomForest(Sales ~ ., 
                      data = train_dt,
                      mtry = m,
                      xtest=test_dt[,-1], 
                      ytest=test_dt$Sales))
}

tmp_dt <- data.table( num_tree = 1:500,
                      rf_1 = mtry_rf(1)$mse,
                      rf_3 = mtry_rf(3)$mse,
                      rf_4 = mtry_rf(4)$mse,
                      rf_7 = mtry_rf(7)$mse,
                      rf_10 =  mtry_rf(10)$mse)

melt.tmp <- melt(tmp_dt, id=1)

ggplot(melt.tmp, aes(num_tree, value, col=variable)) + 
  geom_line(lwd=1) +
  labs(y='MSE', col="") + theme_bw()


###########
tmp_dt <- data.table( num_tree = 1:500,
                      rf_1 = mtry_rf(1)$test$mse,
                      rf_3 = mtry_rf(3)$test$mse,
                      rf_4 = mtry_rf(4)$test$mse,
                      rf_7 = mtry_rf(7)$test$mse,
                      rf_10 =  mtry_rf(10)$test$mse)

melt.tmp <- melt(tmp_dt, id=1)

ggplot(melt.tmp, aes(num_tree, value, col=variable)) + 
  geom_line(lwd=1) +
  labs(y='MSE', col="") + theme_bw()

# test_mse
mtry_rf(1)$test$mse[500]
mtry_rf(3)$test$mse[500]
mtry_rf(4)$test$mse[500]
mtry_rf(7)$test$mse[500]
mtry_rf(10)$test$mse[500]



######################################
#### 4
######################################


#### (a)
dt_OJ <- OJ
dim(OJ)

set.seed(1234)
train_id <- sample(1:nrow(dt_OJ), 800)

train_dt <- dt_OJ[train_id,]
test_dt <- dt_OJ[-train_id,]

#### (b)
cart_OJ<- rpart(Purchase  ~ ., data=train_dt)
cart_OJ

mean(predict(cart_OJ, type='class') != train_dt$Purchase)  ##training error

#### (c)
rpart.plot(cart_OJ, 
           main = "Classification using CART")

#### (d)
hat_y_OJ <- predict(cart_OJ, newdata = test_dt, type='class')
table(hat_y_OJ, test_dt$Purchase) ## confusion matrix
mean(hat_y_OJ != test_dt$Purchase) ##test error

#### (e)
cart_OJ$cptable  ##complexity parameter 
printcp(cart_OJ)
plotcp(cart_OJ)


#가지치기
prune_cart_OJ <- rpart(Purchase ~ ., data=train_dt,
                       control = rpart.control(cp = 0.02))


rpart.plot(prune_cart_OJ, main = "Classification using CART")


#### (f)
mean(predict(cart_OJ, type='class') != train_dt$Purchase)
mean(predict(prune_cart_OJ, type='class') != train_dt$Purchase)


#### (g)
mean(hat_y_OJ != test_dt$Purchase) ##test error before pruning 
mean(predict(prune_cart_OJ,newdata = test_dt, type='class') != test_dt$Purchase)  ##test error after pruning
