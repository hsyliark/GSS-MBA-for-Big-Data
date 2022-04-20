# Problem 2

## (a)
auto <- read.csv("C:/Users/HSY/Desktop/Auto.csv",
                 sep=",",header=T,na.strings="?",stringsAsFactors=T)
auto <- na.omit(auto)
auto$mpg01 <- ifelse(auto$mpg > median(auto$mpg), 1, 0)

## (b)
pairs(mpg01~cylinders+displacement+horsepower+weight+acceleration+year+origin,
      main = "Auto Data", 
      pch = 21, bg = c("red","blue")[as.factor(auto$mpg01)], data=auto)
round(cor(auto[,c(-1,-9)]),4)

## (c)
auto$mpg01 <- as.factor(auto$mpg01)
auto$origin <- as.factor(auto$origin)
set.seed(55364)
id_train <- sample(x=1:nrow(auto),size=round(0.6*nrow(auto),0),replace=F)
auto_train <- auto[id_train,]
auto_test <- auto[-id_train,]

## (d)
library(MASS)
lda.fit <- lda(mpg01~displacement+horsepower+weight, data=auto_train)
lda.fit
table(auto_test$mpg01, predict(lda.fit,auto_test)$class)
mean(auto_test$mpg01 != predict(lda.fit,auto_test)$class)

## (e)
qda.fit <- qda(mpg01~displacement+horsepower+weight, data=auto_train)
qda.fit
table(auto_test$mpg01, predict(qda.fit,auto_test)$class)
mean(auto_test$mpg01 != predict(qda.fit,auto_test)$class)

## (f)
glm.fits <- glm(mpg01~displacement+horsepower+weight, data=auto_train, 
                family = binomial)
summary(glm.fits)
fitted_class <- as.factor(ifelse(predict(glm.fits, auto_test,
                                         type='response')>0.5,1,0))
table(auto_test$mpg01, fitted_class)
mean(auto_test$mpg01 != fitted_class)

## (g)
library(class)
auto_train_x <- auto_train[,c(3,4,5)]
auto_train_y <- auto_train[,10]
auto_test_x <- auto_test[,c(3,4,5)]
auto_test_y <- auto_test[,10] 

pred_test <- knn(auto_train_x, auto_test_x, 
                 cl=auto_train_y, k=1, prob = T)  
table(auto_test_y, pred_test) 
mean(pred_test!=auto_test_y)  

pred_test <- knn(auto_train_x, auto_test_x, 
                 cl=auto_train_y, k=2, prob = T)  
table(auto_test_y, pred_test)  
mean(pred_test!=auto_test_y)  

pred_test <- knn(auto_train_x, auto_test_x, 
                 cl=auto_train_y, k=3, prob = T)   
table(auto_test_y, pred_test)  
mean(pred_test!=auto_test_y)  

k.grid <- 1:21
acc_df <- data.frame(k = k.grid, acc = NA)
for(i in 1:length(k.grid)){
  pred_test <- knn(auto_train_x, auto_test_x, 
                   cl=auto_train_y, k=k.grid[i])  
  acc_df[i,2] <- mean(pred_test!=auto_test_y)  
}
acc_df$k1 <- 1/acc_df$k
plot(acc_df$k1, acc_df$acc, type='b',
     xlab="1/k", ylab="test misclassification rate", main="KNN of Auto data")
plot(acc_df$k, acc_df$acc, type='b',
     xlab="k", ylab="test misclassification rate", main="KNN of Auto data")



# Problem 3

boston <- read.csv("C:/Users/HSY/Desktop/Boston.csv",sep=",",header=T)
boston$crim01 <- as.factor(ifelse(boston$crim > median(boston$crim), 1, 0))
boston <- boston[,-c(1,2)]
set.seed(55364)
id_train <- sample(x=1:nrow(boston),size=round(0.6*nrow(boston),0),replace=F)
boston_train <- boston[id_train,]
boston_test <- boston[-id_train,]

## Logistic regression
glm.fits <- glm(crim01~., data=boston_train, family = binomial)
summary(glm.fits)
fitted_class <- as.factor(ifelse(predict(glm.fits, boston_test,
                                         type='response')>0.5,1,0))
table(boston_test$crim01, fitted_class)
mean(boston_test$crim01 != fitted_class)

## LDA
lda.fit <- lda(crim01~., data=boston_train)
lda.fit
table(boston_test$crim01, predict(lda.fit,boston_test)$class)
mean(boston_test$crim01 != predict(lda.fit,boston_test)$class)

## KNN
boston_train_x <- boston_train[,1:12]
boston_train_y <- boston_train[,13]
boston_test_x <- boston_test[,1:12]
boston_test_y <- boston_test[,13]

k.grid <- 1:21
acc_df <- data.frame(k = k.grid, acc = NA)
for(i in 1:length(k.grid)){
  pred_test <- knn(boston_train_x, boston_test_x, 
                   cl=boston_train_y, k=k.grid[i])  
  acc_df[i,2] <- mean(pred_test!=boston_test_y)  
}
acc_df$k1 <- 1/acc_df$k
plot(acc_df$k1, acc_df$acc, type='b',
     xlab="1/k", ylab="test misclassification rate", main="KNN of Boston data")
plot(acc_df$k, acc_df$acc, type='b',
     xlab="k", ylab="test misclassification rate", main="KNN of Boston data")

pred_test <- knn(boston_train_x, boston_test_x, 
                 cl=boston_train_y, k=3, prob = T)   
table(boston_test_y, pred_test)  
mean(pred_test!=boston_test_y) 
