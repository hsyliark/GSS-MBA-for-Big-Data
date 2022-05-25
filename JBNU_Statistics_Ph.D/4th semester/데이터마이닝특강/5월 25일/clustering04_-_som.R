####################################################
##### 군집분석
##### SOM
####################################################

library(kohonen)  #som
library(gclus)  #wine data


data(wine)
head(wine)

wine.sc <- scale(wine)

set.seed(7)
wine.som <- som(wine.sc, 
                grid = somgrid(5, 4, topo = "hexagonal"))
summary(wine.som)
attributes(wine.som)
nunits(wine.som)

wine.som$distances
wine.som$unit.classif

plot(wine.som, main="Wine data")

par(mfrow=c(1,3))
plot(wine.som, type="counts", main="wine data: counts")
plot(wine.som, type="quality", main="wine data: mapping quality")
plot(wine.som, type="mapping",  col = wine.som$unit.classif, 
     pch = wine.som$unit.classif, main="mapping plot")
graphics.off()







