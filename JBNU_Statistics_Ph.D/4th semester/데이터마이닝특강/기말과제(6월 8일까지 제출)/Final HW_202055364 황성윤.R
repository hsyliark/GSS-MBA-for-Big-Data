#### Problem 1

### Loading packages
library(arules)
library(arulesViz)

## (a)
data(Income)

## (b)
?Income
inspect(Income)
summary(Income)

## (c)
Income.40000 <- Income[Income %in% "income=$40,000+"]
itemFrequencyPlot(Income.40000,
                  support=0.4,
                  main = "Histogram of Income data with 'income=$40,000+' ")

## (d)
rules <- apriori(Income,
                 parameter = list(supp=0.1, conf=0.8))
inspect(rules)
rules.sub <- subset(rules, rhs %in% "income=$40,000+" & lift > 1.0)
inspect(sort(rules.sub, by="lift"))


#### Problem 2

### Loading packages
library(rpart)  # Recursive Partitioning and Regression Trees
library(rpart.plot)
library(tree)
library(randomForest) # random Forest
library(ipred) # bagging
library(gbm) # Boosting
library(xgboost)  # xgboost  
library(adabag)  # Adaboosting : boosting
library(ggplot2)
library(data.table)

## (a)
DT <- read.csv("C:/Users/HSY/Desktop/sample_DT.csv",sep=",",header=T)
set.seed(1234)
train_id <- sample(1:nrow(DT), nrow(DT)*0.7)
train_dt <- DT[train_id,]
test_dt <- DT[-train_id,]

## (b)
tree_dt <- rpart(DEFECT_TYPE ~ ., data=train_dt)
tree_dt
summary(tree_dt)

rpart.plot(tree_dt)
plot(tree_dt)
text(tree_dt, all=T)

tree_dt$cptable  
printcp(tree_dt)
plotcp(tree_dt)

prune_tree_dt <- rpart(DEFECT_TYPE ~., data=train_dt,
                             control = rpart.control(cp = 0.0125))
prune_tree_dt
summary(prune_tree_dt)

rpart.plot(prune_tree_dt, main = "Classification using CART (After pruning)")
plot(prune_tree_dt)
text(prune_tree_dt, all=T)

yhat <- predict(prune_tree_dt, newdata=test_dt, type="class")
table(yhat, test_dt$DEFECT_TYPE)
mean(yhat!=test_dt$DEFECT_TYPE)

## (c)
DT1 <- DT
DT1$DEFECT_TYPE <- as.factor(DT1$DEFECT_TYPE)
set.seed(1234)
train_id <- sample(1:nrow(DT1), nrow(DT1)*0.7)
train_dt <- DT1[train_id,]
test_dt <- DT1[-train_id,]

fit.bagg <- randomForest(DEFECT_TYPE ~ ., 
                         data = train_dt,
                         ntree=1000,
                         mtry = 20, 
                         importance = TRUE,
                         proximity = TRUE)
fit.bagg$importance

yhat <- predict(fit.bagg, newdata=test_dt, type="class")
table(yhat, test_dt$DEFECT_TYPE)
mean(yhat!=test_dt$DEFECT_TYPE)

## (d)
DT2 <- DT1
DT2$DEFECT_TYPE <- ifelse(DT2$DEFECT_TYPE=="G",1,0)
set.seed(1234)
train_id <- sample(1:nrow(DT2), nrow(DT2)*0.7)
train_dt <- DT2[train_id,]
test_dt <- DT2[-train_id,]

boosting_dt <- gbm(DEFECT_TYPE ~.,data=train_dt,
                     distribution = 'bernoulli',
                     interaction.depth = 2,
                     n.trees=1000)

boosting_dt
summary(boosting_dt)

yhat.boost_prob <- predict(boosting_dt,
                           newdata = test_dt, 
                           n.trees = 1000, 
                           type='response')

yhat <- ifelse(yhat.boost_prob > 0.5, 'G', 'NG')

train_dt <- DT1[train_id,]
test_dt <- DT1[-train_id,]
table(yhat, test_dt$DEFECT_TYPE)
mean(yhat!=test_dt$DEFECT_TYPE)

## (e)
DT1 <- DT
DT1$DEFECT_TYPE <- as.factor(DT1$DEFECT_TYPE)
set.seed(1234)
train_id <- sample(1:nrow(DT1), nrow(DT1)*0.7)
train_dt <- DT1[train_id,]
test_dt <- DT1[-train_id,]

bag.Heart <- randomForest(DEFECT_TYPE ~ ., data = train_dt, 
                          mtry = 4, ntree=1000)
yhat.bag <- predict(bag.Heart, newdata = test_dt)
table(yhat, test_dt$DEFECT_TYPE)
mean(yhat!=test_dt$DEFECT_TYPE)

mtry_rf_dt <- function(m){
  return(randomForest(DEFECT_TYPE ~ ., 
                      data = train_dt,
                      mtry = m)$err.rate[,1])
}

tmp_dt <- data.table( num_tree = 1:500,
                      rf_1 = mtry_rf_dt(1),
                      rf_4 = mtry_rf_dt(4),
                      rf_8 = mtry_rf_dt(7),
                      rf_13 = mtry_rf_dt(15),
                      rf_20 = mtry_rf_dt(20))

melt.tmp <- melt(tmp_dt, id=1)

ggplot(melt.tmp, aes(num_tree, value, col=variable)) + geom_line(lwd=1) +
  labs(y='error rate', col="") + theme_bw()


#### Problem 3

wcd <- read.csv("C:/Users/HSY/Desktop/Wholesale_customers_data.csv",
                sep=",",header=T)
wcd <- wcd[,-(1:2)]

### k-means
data_k <- kmeans(wcd ,centers=3, nstart=5) 
attributes(data_k)

data_k$cluster
data_k$centers
data_k$withinss
data_k$tot.withinss
data_k$size
data_k$iter

par(mfrow=c(2,3))
data_k <- kmeans(wcd ,centers=3, nstart=5) 
plot(wcd[,c(1,2)], pch = 16, col =  data_k$cluster, 
     main = round(data_k$tot.withinss,2),
     xaxt='n', yaxt='n',
     xlab='Fresh', ylab='Milk')
par(mfrow=c(1,1))

plot(wcd$Fresh,wcd$Milk)
wcd_c <- wcd[(wcd$Fresh<=60000)&(wcd$Milk<=30000),]
plot(wcd_c$Fresh,wcd_c$Milk)

par(mfrow=c(2,3))
data_k <- kmeans(wcd_c ,centers=3, nstart=5) 
plot(wcd_c[,c(1,2)], pch = 16, col =  data_k$cluster, 
     main = round(data_k$tot.withinss,2),
     xaxt='n', yaxt='n',
     xlab='Fresh', ylab='Milk')
par(mfrow=c(1,1))

### k-medoids
library(cluster)  # pam
library(factoextra)  # fviz_cluster
library(MASS)
library(cowplot)
library(gridExtra)

fit <- pam(wcd,k=3)
summary(fit)

par(mfrow=c(1,2))
fit1 <- pam(wcd,k=3)
fit2 <- pam(wcd,k=3)
plot(wcd[,c(1,2)], pch = 16, col =  fit1$clustering, 
     xaxt='n', yaxt='n',
     xlab='Fresh', ylab='Milk',
     main="With influential point")
plot(wcd_c[,c(1,2)], pch = 16, col =  fit2$clustering, 
     xaxt='n', yaxt='n',
     xlab='Fresh', ylab='Milk',
     main="Without influential point")
par(mfrow=c(1,1))
