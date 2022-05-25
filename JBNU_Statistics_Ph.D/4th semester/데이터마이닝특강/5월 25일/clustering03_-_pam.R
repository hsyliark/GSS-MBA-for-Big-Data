####################################################
##### 군집분석
##### K-Medois
####################################################

library(cluster)  ##pam
library(factoextra)  ## fviz_cluster
library(MASS)
library(cowplot)
library(gridExtra)


################################################## k-medoids with pam algorithm
################################################## k-medoids with pam algorithm 

# pam algorithm 
data("USArrests")      # Loading the data set
sam = scale(USArrests) # Scaling the data

fit = pam(sam,k=4)
summary(fit)


fit$id.med # Medoid samples 
row.names(sam)[fit$id.med]
fit$clusinfo # cluster summary
fit$clustering # cluster index 
row.names(sam)[fit$clustering==1] # The first cluster entries 

fviz_cluster(fit,sam) + theme_bw()


################################################## k-means and k-medoids with outliers
################################################## k-means and k-medoids with outliers 

sam <- sam[,1:2]

tail(sam,n=3) 
new.sam = rbind(sam,c(4,4)) # Adding new sample 
row.names(new.sam)[51] = "New City"
tail(new.sam,n=3)

# construct clusers 
kmeans.fit = kmeans(sam,centers=2)
new.kmeans.fit = kmeans(new.sam,centers=2)
pam.fit = pam(sam,2)
new.pam.fit = pam(new.sam,2)

# plots 
list.plot = list() 
list.plot[[1]] = fviz_cluster(kmeans.fit,
                              sam, 
                              main='K-means w/o outlier')+theme_bw()
list.plot[[3]] = fviz_cluster(new.kmeans.fit,
                              new.sam, 
                              main='K-means w/ outlier') +theme_bw()
list.plot[[2]] = fviz_cluster(pam.fit,sam, main='K-medoids w/o outlier')+theme_bw()
list.plot[[4]] = fviz_cluster(new.pam.fit,new.sam, main='K-medoids w/ outlier')+theme_bw()
plot_grid(plotlist = list.plot, ncol = 2)




################################################## 
################################################## 


# validation of clustering: fviz_nbclust
wss = fviz_nbclust(sam, pam, method = "wss")
p1 <- wss+theme(axis.text = element_text(size = 8, color = "red"), 
          title = element_text(size = 8, color = "blue"))
wss$data

sil = fviz_nbclust(sam, pam, method = "silhouette")
p2 <- sil+theme(axis.text = element_text(size = 8, color = "red"), 
          title = element_text(size = 8, color = "blue"))
sil$data


grid.arrange(p1, p2, nrow = 1)


# validation of clustering: fviz_nbclust
wss = fviz_nbclust(sam, kmeans, method = "wss")
p1 <- wss+theme(axis.text = element_text(size = 8, color = "red"), 
                title = element_text(size = 8, color = "blue"))
wss$data

sil = fviz_nbclust(sam, kmeans, method = "silhouette")
p2 <- sil+theme(axis.text = element_text(size = 8, color = "red"), 
                title = element_text(size = 8, color = "blue"))
sil$data


grid.arrange(p1, p2, nrow = 1)

