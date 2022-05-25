####################################################
##### 군집분석
##### DBScan
####################################################
library(dbscan)  #dbscan, knndist
library(cluster)
library(ggplot2)

set.seed(1004)
x <- c(rnorm(50,1,0.05), rnorm(300,2,0.4), rnorm(100,1,0.2))
y <- c(rnorm(50,0.5,0.05), rnorm(300,2,0.4), rnorm(100,-1,0.2))

dt <- data.frame(x=x, y=y)

ggplot(dt, aes(x,y))+
  geom_point(col='steelblue')+
  xlab("")+ylab("")+
  theme_bw()

plot(dt, pch=16, col='steelblue')

kNNdistplot(dt, k = 5)
kNNdist(dt, k = 5)
head(kNNdist(dt, k = 5),15)
abline(h = seq(0.2,0.4,0.05), col = "red", lty=2)
abline(h=0.1)
kNNdistplot(dt, k = 10)

eps <- 0.3
res <- dbscan(dt, eps = eps , minPts = 5)
str(res)
res
ggplot(dt, aes(x,y, col=as.factor(res$cluster)))+
    geom_point()+
    xlab("")+ylab("")+ggtitle(paste0('BDscan : epsilon = ', eps))+labs(col = "cluster")+
    theme_bw()


#################################################################### 
eps <- 0.2
res1 <- dbscan(dt, eps = eps , minPts = 5)
p1 <- ggplot(dt, aes(x,y, col=as.factor(res1$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('BDscan : epsilon = ', eps))+labs(col = "cluster")+
  theme_bw()

eps <- 0.3
res2 <- dbscan(dt, eps = eps , minPts = 5)
p2 <- ggplot(dt, aes(x,y, col=as.factor(res2$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('BDscan : epsilon = ', eps))+labs(col = "cluster")+
  theme_bw()

eps <- 0.7
res3 <- dbscan(dt, eps = eps , minPts = 5)
p3 <- ggplot(dt, aes(x,y, col=as.factor(res3$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('BDscan : epsilon = ', eps))+labs(col = "cluster")+
  theme_bw()

grid.arrange(p1, p2, p3, nrow = 1)

#################################################################### 

#################################################################### 
m <- 3
res1 <- dbscan(dt, eps = 0.3 , minPts = m)
p1 <- ggplot(dt, aes(x,y, col=as.factor(res1$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('BDscan : minPts = ', m))+
  labs(col = "cluster")+
  theme_bw()

m <- 5
res2 <- dbscan(dt, eps = 0.3 , minPts = m)
p2 <- ggplot(dt, aes(x,y, col=as.factor(res2$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('BDscan : minPts = ', m))+
  labs(col = "cluster")+
  theme_bw()

m <- 10
res3 <- dbscan(dt, eps = 0.4 , minPts = m)
p3 <- ggplot(dt, aes(x,y, col=as.factor(res3$cluster)))+
  geom_point()+
  xlab("")+ylab("")+ggtitle(paste0('BDscan : minPts = ', m))+
  labs(col = "cluster")+
  theme_bw()

grid.arrange(p1, p2, p3, nrow = 1)

#################################################################### 

res_k <- kmeans(dt, center=3, nstart = 20)
ggplot(dt, aes(x,y, col=as.factor(res_k$cluster)))+geom_point()+
  xlab("")+ylab("")+ggtitle('Kmeans')+
  labs(col = "cluster")+
  theme_bw()




#######################################################
#######################################################
data(ruspini)
head(ruspini)

ggplot(ruspini, aes(x,y))+geom_point(col='steelblue')+
  xlab("")+ylab("")+
  theme_bw()

kNNdistplot(ruspini, k = 5)
abline(h = 15:20, col = "red", lty=2)

######## DBscan
res <- dbscan(ruspini, eps = 15, minPts = 5)
res

ggplot(ruspini, aes(x,y, col=as.factor(res$cluster)))+geom_point()+
  xlab("")+ylab("")+coord_fixed()+labs(col = "cluster")+
  theme_bw()

######## K-means
res_k <- kmeans(ruspini,centers =4) 

ruspini$cluster_k <- as.factor(res_k$cluster)

ggplot(ruspini, aes(x,y, col=cluster_k))+geom_point()+
  xlab("")+ylab("")+coord_fixed()+labs(col = "cluster")+
  theme_bw()



#######################################################
#######################################################

set.seed(-1)
get.sample <- function(n=1000, p=0.7){
  x1 <- rnorm(n)
  y1 <- rnorm(n)
  r2 <- 7 + rnorm(n)
  t2 <- runif(n,0,2*pi)
  x2 <- r2*cos(t2)
  y2 <- r2*sin(t2)
  r <- runif(n)>p
  x <- ifelse(r,x1,x2)
  y <- ifelse(r, y1, y2)
  d <- data.frame(x=x, y=y)
  d
}

dt <- get.sample()

ggplot(dt, aes(x,y))+geom_point(col='steelblue')+
  xlab("")+ylab("")+
  theme_bw()

######### DBScan vs. kmeans
kNNdistplot(dt, k = 5)
abline(h = seq(1,1.5,0.1), col = "red", lty=2)

res_db <- dbscan(dt, eps = 1, minPts = 5)
res_kmeans <- kmeans(dt,centers =2) 

p1 <- ggplot(dt, aes(x,y, col=as.factor(res_db$cluster)))+geom_point()+
  xlab("")+ylab("")+ggtitle('DB scan')+
  labs(col = "cluster")+
  theme_bw()

p2 <- ggplot(dt, aes(x,y, col=as.factor(res_kmeans$cluster)))+geom_point()+
  xlab("")+ylab("")+ggtitle('K menas')+
  labs(col = "cluster")+
  theme_bw()

grid.arrange(p1, p2, nrow = 1)

#######################################################
#######################################################

set.seed(12)

x1 <-  runif(300,-1,1)
y1 <- 2*x1^2-1.5+ rnorm(300,0,0.2)


x2 <- runif(300,-2,0)
y2 <- -2*(x2+1)^2 + 1.5 + rnorm(300,0,0.2)


dt1 <- data.frame(x = c(x1,x2),
                  y = c(y1,y2))



ggplot(dt1, aes(x,y))+geom_point(col='steelblue')+
  xlab("")+ylab("")+
  theme_bw()


######### DBScan
kNNdistplot(dt1, k = 5)
abline(h = 0.2, col = "red", lty=2)

res_db <- dbscan(dt1, eps = 0.2, minPts = 5)
res_kmeans <- kmeans(dt1,centers =2) 

p1 <- ggplot(dt1, aes(x,y, col=as.factor(res_db$cluster)))+geom_point()+
  xlab("")+ylab("")+ggtitle('DB scan')+
  labs(col = "cluster")+
  theme_bw()

p2 <- ggplot(dt1, aes(x,y, col=as.factor(res_kmeans$cluster)))+geom_point()+
  xlab("")+ylab("")+ggtitle('K menas')+
  labs(col = "cluster")+
  theme_bw()

grid.arrange(p1, p2, nrow = 1)



