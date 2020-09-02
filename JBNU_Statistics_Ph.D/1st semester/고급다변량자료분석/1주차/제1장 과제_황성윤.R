### Chapter 1 Homework

### 1.6
## Reading data
air <- read.csv("C:/Users/stat/Desktop/다변량 1주차 과제 (진행중)/data/air_pollution.csv", sep=",", header=T)
## (a)
install.packages("ggplot2")
library(ggplot2)
install.packages("ggExtra")
library(ggExtra)
p <- ggplot(air, aes(x=O3, y=HC)) +
  geom_point() +
  theme(legend.position="none")
# with marginal histogram
p1 <- ggMarginal(p, type="histogram")
# marginal density
p2 <- ggMarginal(p, type="density")
# marginal boxplot
p3 <- ggMarginal(p, type="boxplot")
## (b)
sapply(air,mean)
var(air)
cor(air)

### 1.15
## Reading data
therapy <- read.csv("C:/Users/stat/Desktop/다변량 1주차 과제 (진행중)/data/radiotherapy.csv", sep=",", header=T)

## (a)
p <- ggplot(therapy, aes(x=x2, y=x3)) +
  geom_point() +
  theme(legend.position="none")
p1 <- ggMarginal(p, type="histogram")
## (b)
sapply(therapy,mean)
var(therapy)
cor(therapy)

### 1.18
## Reading data
track <- read.csv("C:/Users/stat/Desktop/다변량 1주차 과제 (진행중)/data/track.csv", sep=",", header=T)
rownames(track) <- track[,1]
track <- track[,-1]
track[,(4:7)] <- 60*track[,(4:7)]
track[,1] <- 100/track[,1]
track[,2] <- 200/track[,2]
track[,3] <- 400/track[,3]
track[,4] <- 800/track[,4]
track[,5] <- 1500/track[,5]
track[,6] <- 3000/track[,6]
track[,7] <- 42195/track[,7]
sapply(track,mean)
var(track)
cor(track)

### 1.23
## Reading data
cereal <- read.csv("C:/Users/stat/Desktop/다변량 1주차 과제 (진행중)/data/cereal.csv", sep=",", header=T)
rownames(cereal) <- cereal[,1]
cereal <- cereal[,-1]
## Stars
stars(cereal)
## Chernoff faces
install.packages("aplpack")
library(aplpack)
faces(cereal)
