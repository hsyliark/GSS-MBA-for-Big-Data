####################################################
##### 군집분석
##### K-means
####################################################

crime <- read.csv("C:\\R-Project\\DAT\\data mining\\Crime.csv")
head(crime)

plot(crime$murder, crime$rape, pch=16)

data_k <- kmeans(crime[,c(3,4)] ,centers = 3) 
attributes(data_k)

data_k$cluster
data_k$centers
data_k$withinss
data_k$tot.withinss
data_k$size
data_k$iter

par(mfrow=c(2,3))
data_k <-  kmeans(crime[,c(3,4)] ,
                  centers = 3)
plot(crime[,c(3,4)], pch = 16, col =  data_k$cluster, 
     main = round(data_k$tot.withinss,2),
     xaxt='n', yaxt='n',
     xlab='', ylab='')

par(mfrow=c(2,3))
set.seed(4)
data_k <-  kmeans(crime[,c(3,4)],
                  centers = 3, nstart = 20)
plot(crime[,c(3,4)], pch = 16, col =  data_k$cluster, 
     main = round(data_k$tot.withinss,2),
     xaxt='n', yaxt='n',
     xlab='', ylab='')

par(mfrow=c(1,1))

crime$cluster <- data_k$cluster
plot(crime[,c(3,4)], 
     pch = crime$cluster-1, 
     col = crime$cluster, 
     main = "K-means clustering")
points(data_k$centers, cex=2, col=c(1,2,3), pch=c(15,16,17))
text(crime[,c(3,4)], labels = crime$city, adj = -0.2, cex = 0.8)



