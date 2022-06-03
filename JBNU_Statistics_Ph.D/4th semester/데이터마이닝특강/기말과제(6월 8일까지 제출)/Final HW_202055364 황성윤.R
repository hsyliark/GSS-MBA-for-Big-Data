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
                      rf_8 = mtry_rf_dt(8),
                      rf_13 = mtry_rf_dt(13),
                      rf_20 = mtry_rf_dt(20))

melt.tmp <- melt(tmp_dt, id=1)

ggplot(melt.tmp, aes(num_tree, value, col=variable)) + geom_line(lwd=1) +
  labs(y='error rate', col="") + theme_bw()


#### Problem 3

wcd <- read.csv("C:/Users/HSY/Desktop/Wholesale_customers_data.csv",
                sep=",",header=T)
wcd <- wcd[,-(1:2)]
wcd_scale <- scale(wcd)

### k-means
tot_withinss <- c()
for (i in 1:20) {
  set.seed(5790) # for reproducibility 
  kmeans_cluster <- kmeans(wcd_scale, centers = i, iter.max = 1000)
  tot_withinss[i] <- kmeans_cluster$tot.withinss
}

library(factoextra)
library(cluster)
avg_silhouette <- c(0)
for (i in 2:20) {
  set.seed(5790) # for reproducibility 
  km.res <- kmeans(wcd_scale, centers = i, iter.max = 1000)
  sil <- silhouette(km.res$cluster, dist(wcd_scale))[,3] 
  avg_silhouette[i] <- mean(sil)
}

par(mfrow=c(1,2))
plot(c(1:20), tot_withinss, type="b", 
     main="Optimal number of clusters", 
     xlab="Number of clusters", 
     ylab="Total within-cluster sum of squares")
plot(c(1:20), avg_silhouette, type="b", 
     main="Optimal number of clusters", 
     xlab="Number of clusters", 
     ylab="Average Silhouette width")
par(mfrow=c(1,1))
 
data_k <- kmeans(wcd_scale ,centers=5, nstart=5) 
attributes(data_k)

data_k$cluster
data_k$centers
data_k$withinss
data_k$tot.withinss
data_k$size
data_k$iter

plot(wcd$Fresh,wcd$Milk)
wcd_c <- wcd[(wcd$Fresh<=60000)&(wcd$Milk<=30000),]
plot(wcd_c$Fresh,wcd_c$Milk)
wcd_c_scale <- scale(wcd_c)

par(mfrow=c(1,2))
data_k1 <- kmeans(wcd_scale ,centers=5, nstart=5) 
plot(wcd[,c(1,2)], pch = 16, col =  data_k1$cluster, 
     main = paste0("k-means with influential point, \n tot.withinss=",round(data_k1$tot.withinss,2)),
     xaxt='n', yaxt='n',
     xlab='Fresh', ylab='Milk')
data_k2 <- kmeans(wcd_c_scale ,centers=5, nstart=5) 
plot(wcd_c[,c(1,2)], pch = 16, col =  data_k2$cluster, 
     main = paste0("k-means without influential point, \n tot.withinss=",round(data_k2$tot.withinss,2)),
     xaxt='n', yaxt='n',
     xlab='Fresh', ylab='Milk')
par(mfrow=c(1,1))

### k-medoids
library(cluster)  # pam
library(factoextra)  # fviz_cluster
library(MASS)
library(cowplot)
library(gridExtra)

wss <- fviz_nbclust(wcd_scale, pam, method = "wss")
p1 <- wss+theme(axis.text = element_text(size = 8, color = "red"), 
                title = element_text(size = 8, color = "blue"))
wss$data

sil <- fviz_nbclust(wcd_scale, pam, method = "silhouette")
p2 <- sil+theme(axis.text = element_text(size = 8, color = "red"), 
                title = element_text(size = 8, color = "blue"))
sil$data

grid.arrange(p1, p2, nrow = 1)

fit <- pam(wcd_scale,k=5)
summary(fit)

par(mfrow=c(1,2))
fit1 <- pam(wcd_scale,k=2)
fit2 <- pam(wcd_c_scale,k=2)
plot(wcd[,c(1,2)], pch = 16, col =  fit1$clustering, 
     xaxt='n', yaxt='n',
     xlab='Fresh', ylab='Milk',
     main="PAM with influential point")
plot(wcd_c[,c(1,2)], pch = 16, col =  fit2$clustering, 
     xaxt='n', yaxt='n',
     xlab='Fresh', ylab='Milk',
     main="PAM without influential point")
par(mfrow=c(1,1))

### SOM
library(kohonen)  # som
library(gclus)

set.seed(7)
wcd.som <- som(wcd_scale, 
                grid = somgrid(5, 4, topo = "hexagonal"))
summary(wcd.som)
attributes(wcd.som)

wcd.som$distances
wcd.som$unit.classif

par(mfrow=c(1,2))
plot(wcd.som, main="Wine data")
plot(wcd.som, type="mapping",  col = wcd.som$unit.classif, 
     pch = wcd.som$unit.classif, main="mapping plot")
plot(wcd.som, type="counts", main="wine data: counts")
plot(wcd.som, type="quality", main="wine data: mapping quality")
par(mfrow=c(1,1))

### DBSCAN
library(dbscan)  #dbscan, knndist
library(cluster)
library(ggplot2)

eps <- 0.2
res1 <- dbscan(wcd_scale, eps = eps , minPts = 5)
p1 <- ggplot(wcd, aes(Fresh,Milk, col=as.factor(res1$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('DBscan : epsilon = ', eps))+labs(col = "cluster")+
  theme_bw()

eps <- 0.3
res2 <- dbscan(wcd_scale, eps = eps , minPts = 5)
p2 <- ggplot(wcd, aes(Fresh,Milk, col=as.factor(res2$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('DBscan : epsilon = ', eps))+labs(col = "cluster")+
  theme_bw()

eps <- 0.7
res3 <- dbscan(wcd_scale, eps = eps , minPts = 5)
p3 <- ggplot(wcd, aes(Fresh,Milk, col=as.factor(res3$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('DBscan : epsilon = ', eps))+labs(col = "cluster")+
  theme_bw()

grid.arrange(p1, p2, p3, nrow = 1)

m <- 3
res1 <- dbscan(wcd_scale, eps = 0.3 , minPts = m)
p1 <- ggplot(wcd, aes(Fresh,Milk, col=as.factor(res1$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('BDscan : minPts = ', m))+
  labs(col = "cluster")+
  theme_bw()

m <- 5
res2 <- dbscan(wcd_scale, eps = 0.3 , minPts = m)
p2 <- ggplot(wcd, aes(Fresh,Milk, col=as.factor(res2$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('BDscan : minPts = ', m))+
  labs(col = "cluster")+
  theme_bw()

m <- 10
res3 <- dbscan(wcd_scale, eps = 0.3 , minPts = m)
p3 <- ggplot(wcd, aes(Fresh,Milk, col=as.factor(res3$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('BDscan : minPts = ', m))+
  labs(col = "cluster")+
  theme_bw()

grid.arrange(p1, p2, p3, nrow = 1)

eps <- 0.2
res1 <- dbscan(wcd_c_scale, eps = eps , minPts = 5)
p1 <- ggplot(wcd_c, aes(Fresh,Milk, col=as.factor(res1$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('DBscan : epsilon = ', eps))+labs(col = "cluster")+
  theme_bw()

eps <- 0.3
res2 <- dbscan(wcd_c_scale, eps = eps , minPts = 5)
p2 <- ggplot(wcd_c, aes(Fresh,Milk, col=as.factor(res2$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('DBscan : epsilon = ', eps))+labs(col = "cluster")+
  theme_bw()

eps <- 0.7
res3 <- dbscan(wcd_c_scale, eps = eps , minPts = 5)
p3 <- ggplot(wcd_c, aes(Fresh,Milk, col=as.factor(res3$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('DBscan : epsilon = ', eps))+labs(col = "cluster")+
  theme_bw()

grid.arrange(p1, p2, p3, nrow = 1)

m <- 3
res1 <- dbscan(wcd_c_scale, eps = 0.3 , minPts = m)
p1 <- ggplot(wcd_c, aes(Fresh,Milk, col=as.factor(res1$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('BDscan : minPts = ', m))+
  labs(col = "cluster")+
  theme_bw()

m <- 5
res2 <- dbscan(wcd_c_scale, eps = 0.3 , minPts = m)
p2 <- ggplot(wcd_c, aes(Fresh,Milk, col=as.factor(res2$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('BDscan : minPts = ', m))+
  labs(col = "cluster")+
  theme_bw()

m <- 10
res3 <- dbscan(wcd_c_scale, eps = 0.3 , minPts = m)
p3 <- ggplot(wcd_c, aes(Fresh,Milk, col=as.factor(res3$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('BDscan : minPts = ', m))+
  labs(col = "cluster")+
  theme_bw()

grid.arrange(p1, p2, p3, nrow = 1)
