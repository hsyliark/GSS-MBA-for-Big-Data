## Problem 3

library(rpart)  ##Recursive Partitioning and Regression Trees
library(rpart.plot)
library(tree)

# (a)
carseats <- read.csv("C:/Users/HSY/Desktop/Carseats.csv",sep=",",header=T)
set.seed(55364)
train_id <- sample(1:nrow(carseats), nrow(carseats)*0.6)
train_dt <- carseats[train_id,]
test_dt <- carseats[-train_id,]

# (b)
tree_carseats <- rpart(Sales ~ ., data=train_dt)
tree_carseats

summary(tree_carseats)

rpart.plot(tree_carseats)
plot(tree_carseats)
text(tree_carseats, all=T)

yhat <- predict(tree_carseats, newdata=test_dt)
yhat
mean((yhat-test_dt$Sales)^2)

# (c)
tree_carseats$cptable  
printcp(tree_carseats)
plotcp(tree_carseats)

prune_tree_carseats <- rpart(Sales ~., data=train_dt,
                        control = rpart.control(cp = 0.023))
prune_tree_carseats
summary(prune_tree_carseats)

rpart.plot(prune_tree_carseats, main = "Regression using CART / Carseats data")
plot(prune_tree_carseats)
text(prune_tree_carseats, all=T)

yhat <- predict(prune_tree_carseats, newdata=test_dt)
yhat
mean((yhat-test_dt$Sales)^2)  

# (d)
library(randomForest) ##random Forest
library(ipred) ##bagging
library(gbm) ## Boosting
library(xgboost)  ## xgboost  
library(adabag)  ##Adaboosting : boosting
library(rpart)
library(rpart.plot)
library(ggplot2)
library(data.table)

fit.bagg <- ipredbagg(train_dt$Sales, 
                     train_dt[,-1],
                     nbagg=1000, 
                     coob=T)  
fit.bagg

pred<-predict(fit.bagg, 
              newdata = test_dt)
mean((pred-test_dt$Sales)^2)  

# (e)
rf.carseats <- randomForest(Sales ~ ., 
                          data = train_dt,
                          ntree=1000,
                          mtry = 3, 
                          importance = TRUE)

rf.carseats
yhat.bag <- predict(rf.carseats, 
                    newdata = test_dt)
mean((yhat.bag - test_dt$Sales)^2) ##Test MSE

importance(rf.carseats)
varImpPlot(rf.carseats, pch=16)

mtry_rf <- function(m){
  return(randomForest(Sales ~ ., 
                      data = train_dt,
                      ntree = 1000,
                      mtry = m,
                      xtest=test_dt[,-1], 
                      ytest=test_dt$Sales)$mse)
}
mtry_rf(1)

tmp_dt <- data.table( num_tree = 1:1000,
                      rf_1 = mtry_rf(1),
                      rf_2 = mtry_rf(2),
                      rf_3 = mtry_rf(3),
                      rf_4 = mtry_rf(4),
                      rf_5 = mtry_rf(5),
                      rf_6 = mtry_rf(6),
                      rf_7 = mtry_rf(7),
                      rf_8 = mtry_rf(8),
                      rf_9 = mtry_rf(9),
                      rf_10 = mtry_rf(10) # Bagging
                      )

melt.tmp <- melt(tmp_dt, id=1)

ggplot(melt.tmp, aes(num_tree, value, col=variable)) + 
  geom_line(lwd=1) +
  labs(y='MSE', col="") + theme_bw()



## Problem 4

# (a)
oj <- read.csv("C:/Users/HSY/Desktop/OJ.csv",sep=",",header=T)
set.seed(55364)
train_id <- sample(1:nrow(oj), 800)
train_dt <- oj[train_id,]
test_dt <- oj[-train_id,]

# (b)
tree_oj <- rpart(Purchase ~ ., data=train_dt)
tree_oj

summary(tree_oj)

yhat <- predict(tree_oj, newdata = train_dt, 
                type='class')
yhat
table(yhat, train_dt$Purchase)
mean(yhat!=train_dt$Purchase)

# (c)
rpart.plot(tree_oj, 
           main = "Classification using CART")
par(mar = c(5.1, 6.0, 4.1, 2.1))
plot(tree_oj)
text(tree_oj, all=T)

# (d)
yhat <- predict(tree_oj, newdata = test_dt, 
                type='class')
yhat
table(yhat, test_dt$Purchase)
mean(yhat!=test_dt$Purchase)

# (e)
tree_oj$cptable  
printcp(tree_oj)
plotcp(tree_oj)

prune_tree_oj <- rpart(Purchase ~ ., data=train_dt,
                       control = rpart.control(cp = 0.02))

prune_tree_oj

summary(prune_tree_oj)

plot(prune_tree_oj)
text(prune_tree_oj, all=T)

# (f)
tree_oj
yhat <- predict(tree_oj, newdata = train_dt, 
                type='class')
table(yhat, train_dt$Purchase)
mean(yhat!=train_dt$Purchase)
prune_tree_oj
yhat <- predict(prune_tree_oj, newdata = train_dt, 
                type='class')
table(yhat, train_dt$Purchase)
mean(yhat!=train_dt$Purchase)

# (g)
tree_oj
yhat <- predict(tree_oj, newdata = test_dt, 
                type='class')
table(yhat, test_dt$Purchase)
mean(yhat!=test_dt$Purchase)
prune_tree_oj
yhat <- predict(prune_tree_oj, newdata = test_dt, 
                type='class')
table(yhat, test_dt$Purchase)
mean(yhat!=test_dt$Purchase)



