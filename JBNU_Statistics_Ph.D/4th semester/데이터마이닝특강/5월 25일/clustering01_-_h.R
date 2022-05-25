####################################################
##### 군집분석
##### 계층적 군집분석
####################################################

crime <- read.csv("C:\\R-Project\\DAT\\data mining\\Crime.csv")
head(crime)

x <- crime[,3:4]

plot(x, pch=16)
text(x, labels = crime$city, adj = -0.1, cex = 0.8)

D1 <- dist(x)  # 유클리드거리
round(D1,2)

D2 <- dist(x, method = "manhattan") # 맨하탄거리
round(D2,2)

hc1 <- hclust(dist(x)^2, method = "single")   #최단 연결법
plot(hc1, 
     labels = crime$city, hang = 0.1,
     main = "dendrogram : 최단 연결법")
plot(hc1, 
     hang = 0.1,
     main = "dendrogram : 최단 연결법")

hc1$merge
hc1$height
hc1$dist.method


par(mfrow=c(2,3))
plot(hc1, labels = crime$city, hang = 0.1,
     main = "dendrogram : 최단 연결법", xlab='', sub='')


hc2 <- hclust(dist(x)^2, method = "complete")   #최장 연결법
plot(hc2, labels = crime$city,main = "dendrogram : 최장 연결법", xlab='', sub='')

hc3 <- hclust(dist(x)^2, method = "centroid")   #중심법
plot(hc3, labels = crime$city,main = "dendrogram : 중심법", xlab='', sub='')

hc4 <- hclust(dist(x)^2, method = "ward.D")   #Ward
plot(hc4, labels = crime$city,main = "dendrogram : Ward", xlab='', sub='')

hc5 <- hclust(dist(x)^2, method = "average")   #평균법
plot(hc5, labels = crime$city,main = "dendrogram : 평균법", xlab='', sub='')

hc6 <- hclust(dist(x)^2, method = "median")   #중앙값법
plot(hc6, labels = crime$city,main = "dendrogram : 중앙값법", xlab='', sub='')



############################################################
################## 그룹 별 산점도로 결과 확인
############################################################

c1.num <- 2  #군집 수 설정

############################################################
hc1.result <- cutree(hc1, k = c1.num) #최단 연결법 결과

plot(x, pch = c(15,16,17)[hc1.result], 
     col = hc1.result, main = "single")
text(x, labels = crime$city, adj = -0.1, cex = 0.8)
############################################################


############################################################
hc2.result <- cutree(hc2, k = c1.num) #최장 연결법 결과

# 그룹 별 산점도로 결과 확인
plot(x, pch =  c(15,16,17)[hc2.result], col = hc2.result, main = "complete")
text(x, labels = crime$city, adj = -0.1, cex = 1)
############################################################


############################################################
hc3.result <- cutree(hc3, k = c1.num) 

# 그룹 별 산점도로 결과 확인
plot(x, pch =  c(15,16,17)[hc3.result], col = hc3.result, main = "중심법")
text(x, labels = crime$city, adj = -0.1, cex = 1)
############################################################

############################################################
hc4.result <- cutree(hc4, k = c1.num) 

# 그룹 별 산점도로 결과 확인
plot(x, pch =  c(15,16,17)[hc4.result], col = hc4.result, main = "Ward")
text(x, labels = crime$city, adj = -0.1, cex = 1)
############################################################

############################################################
hc5.result <- cutree(hc5, k = c1.num)

# 그룹 별 산점도로 결과 확인
plot(x, pch =  c(15,16,17)[hc5.result], col = hc5.result, main = "평균법")
text(x, labels = crime$city, adj = -0.1, cex = 1)
############################################################

############################################################
hc6.result <- cutree(hc6, k = c1.num) 

# 그룹 별 산점도로 결과 확인
plot(x, pch =  c(15,16,17)[hc6.result], col = hc6.result, main = "중앙값법")
text(x, labels = crime$city, adj = -0.1, cex = 1)
############################################################



############################################################
################## 표준화
############################################################

sc_x <- scale(x)  ## 표준화 

hc2_sc <- hclust(dist(sc_x)^2, method = "complete")   #최장 연결법

par(mfrow=c(1,2))
plot(hc2, labels = crime$city, hang = 0.1,
     main = "dendrogram : 최장 연결법")

plot(hc2_sc, labels = crime$city, hang = 0.1,
     main = "dendrogram : 최장 연결법, 표준화")

c1.num <- 3
hc2.result <- cutree(hc2, k = c1.num) #최장 연결법 결과
hc2_sc.result <- cutree(hc2_sc, k = c1.num) #최장 연결법 결과

plot(x, pch =  c(15,16,17)[hc2.result], col = hc2.result, main = "최장연결법")
plot(x, pch =  c(15,16,17)[hc2_sc.result], col = hc2_sc.result, main = "최장연결법/표준화")



