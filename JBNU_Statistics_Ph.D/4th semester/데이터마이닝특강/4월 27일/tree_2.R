##############  package
library(rpart)  ##Recursive Partitioning and Regression Trees
library(rpart.plot)
library(tree)
#library(ISLR2) #carseats

##############  경로지정
setwd("C:\\R-Project\\DAT\\INtroduction SL")


#################################################
##### IRIS
#################################################
set.seed(1000)
train_id <- sample(1:nrow(iris), nrow(iris)*0.7)

train_dt <- iris[train_id,]
test_dt <- iris[-train_id,]



################# TREE ###########################
tree_iris <- tree(Species ~ ., data=train_dt)
tree_iris

summary(tree_iris)

plot(tree_iris)
text(tree_iris)
text(tree_iris, all=T)


## 가지치기
sub_tree <- snip.tree(tree_iris, nodes=c(6,7))
plot(sub_tree)
text(sub_tree)


tree_p2  <- prune.misclass(tree_iris)
plot(tree_p2)

cv.iris <- cv.tree(tree_iris,FUN = prune.misclass, K=10)
cv.iris
plot(cv.iris)

par(mfrow = c(1, 2))
plot(cv.iris$size, cv.iris$dev, type = "b")
plot(cv.iris$k, cv.iris$dev, type = "b")

fin.tr_iris  <- prune.tree(tree_iris, best=3)
plot(fin.tr_iris)
text(fin.tr_iris)


## 예측
yhat <- predict(fin.tr_iris, newdata = test_dt, 
                type='class')
yhat

## type : vector, tree, class, where

## 정분류율

table(yhat, test_dt$Species)
mean(yhat==test_dt$Species)

mean( predict(tree_iris, newdata = test_dt, 
              type='class')==test_dt$Species)  ##가지치기전 오분류율 



## 재귀

plot(train_dt$Petal.Length, 
     train_dt$Petal.Width,type = "n")
points(train_dt$Petal.Length, 
       train_dt$Petal.Width, pch=16, 
       col = c(2:(length(levels(train_dt$Species)) + 1))[train_dt$Species])
partition.tree(fin.tr_iris, add=TRUE, cex=1.5)

#################################################


##################### Rpart #####################
cart_iris <- rpart(Species ~ ., data=train_dt, 
                   control = rpart.control(cp = 0.00001,
                                           minsplit = 1))

cart_iris

rpart.plot(cart_iris, 
           main = "Classification using CART")

summary(cart_iris)
cart_iris$variable.importance  ##변수중요도


cart_iris$cptable  ##complexity parameter 
printcp(cart_iris)
plotcp(cart_iris)


#가지치기
prune_cart_iris <- rpart(Species ~ ., data=train_dt,
                  control = rpart.control(cp = 0.2))
             

rpart.plot(prune_cart_iris, main = "Classification using CART")



## 예측
yhat <- predict(prune_cart_iris, newdata = test_dt, 
                type='class')
yhat
## type : vector, tree, class, where

## 정분류율

table(yhat, test_dt$Species)
mean(yhat==test_dt$Species)
mean(yhat!=test_dt$Species)






#################################################
##### airquality
#################################################

head(airquality)
dt_air <- subset(airquality, !is.na(Ozone))
dim(dt_air)


#데이터 분할
set.seed(10)
sample.num            <- sample(1:nrow(dt_air), 0.7*nrow(dt_air))
train                 <- dt_air[sample.num,]
test                  <- dt_air[-sample.num,]

cart_air <- rpart(Ozone ~., data=train,
                  control = rpart.control(cp = 0.01,
                                          minsplit = 5))
cart_air
rpart.plot(cart_air, main = "Regression using CART")

printcp(cart_air)
plotcp(cart_air)

prune_cart_air <- rpart(Ozone ~., data=train,
                        control = rpart.control(cp = 0.05,
                                                minsplit = 5))
prune_cart_air
summary(prune_cart_air)

rpart.plot(prune_cart_air, main = "Regression using CART")

prune_cart_air$variable.importance

yhat <- predict(prune_cart_air, newdata=test)
yhat
mean((yhat-test$Ozone)^2)  ###MSE

mean((predict(cart_air, newdata=test)-test$Ozone)^2)  
###가지치기 전 MSE




