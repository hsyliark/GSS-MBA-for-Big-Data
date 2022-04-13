####################################################
##### Naive Bayes
##### 분류모형 
####################################################
library(mlbench) 
library(e1071)  



###################################################
###################################################
model_iris <- naiveBayes(Species ~ ., data = iris)
model_iris

predict(model_iris, iris)
predict(model_iris, iris, type='raw')
head(round(predict(model_iris, iris, type='raw'),4))

table(iris$Species, predict(model_iris, iris))
mean(predict(model_iris, iris)==iris$Species)
mean(predict(model_iris, iris)!=iris$Species)


###################################################
###################################################
data(HouseVotes84, package = "mlbench")
head(HouseVotes84)

summary(HouseVotes84)

train_id <- sample(1:nrow(HouseVotes84), nrow(HouseVotes84)*0.7)

train_dt <- HouseVotes84[train_id,]
test_dt <- HouseVotes84[-train_id,]

model <- naiveBayes(Class  ~ ., data = train_dt)
model

head(round(predict(model, test_dt, type='raw'),4))
head(predict(model, test_dt))

pred <- predict(model, test_dt)

table(pred, test_dt$Class)
print(paste0("정분류율 = ", 
             round(mean(pred==test_dt$Class)*100,2), "%"))








