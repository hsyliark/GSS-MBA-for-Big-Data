###########################################  
################## HW02  ################## 
###########################################  
setwd("C:\\R-Project\\DAT\\INtroduction SL")

library(MASS)  #lda, qda
library(class)  #KNN


#################  2  ############################
dt <- read.csv("Auto.csv", 
               na.strings = "?", 
               stringsAsFactors = T)
dim(dt)

dt <- na.omit(dt)
dim(dt)

### (a)
dt$mpg01 <- as.factor(ifelse(dt$mpg > median(dt$mpg), 1, 0))
dt <- dt[,-which(names(dt) %in% c('mpg', 'name'))]


### (b)
boxplot(cylinders ~ mpg01, data=dt)
boxplot(displacement  ~ mpg01, data=dt)
boxplot(horsepower ~ mpg01, data=dt)
boxplot(weight ~ mpg01, data=dt)
boxplot(acceleration   ~ mpg01, data=dt)
boxplot(year  ~ mpg01, data=dt)
table(dt$origin, dt$mpg01)

library(ggplot2)
library(data.table)

dt <- as.data.table(dt)
dim(dt)

melt_dt <- melt(dt, id=8)
dim(melt_dt)

ggplot(melt_dt, aes(mpg01, value)) + 
  geom_boxplot() +
  facet_wrap(variable~., scales = 'free_y') + labs(y='')+
  theme_bw()

boxplot(value~ mpg01 + variable, melt_dt)

barplot(table( dt$origin, dt$mpg01))

t.test(acceleration~mpg01, dt)
t.test(year~mpg01, dt)


###(c)
set.seed(1)
train_id <- sample(1:nrow(dt), nrow(dt)*0.6)

names(dt)

dt <- as.data.frame(dt)
train_dt <- dt[train_id,]
test_dt <- dt[-train_id,]

###(d)

lda.fit <- lda(mpg01 ~ ., data = train_dt)

lda_fitted_test <- predict(lda.fit, newdata = test_dt)$class
mean(lda_fitted_test != test_dt$mpg01)


###(e)
qda.fit <- qda(mpg01 ~ ., data = train_dt)

qda_fitted_test <- predict(qda.fit, newdata = test_dt)$class
mean(qda_fitted_test != test_dt$mpg01)

###(f)

glm.fit <- glm(mpg01 ~ ., data = train_dt, 
                    family = binomial)

logistic_fitted_test_prob <- predict(glm.fit, newdata = test_dt,
                                type='response')
logistic_fitted_test <- ifelse(logistic_fitted_test_prob>0.5, 1,0)
mean(logistic_fitted_test != test_dt$mpg01)

###(g)
acc <- data.frame(k = seq(1,21,2),
                  miss_class_rate = 0)

for (k in 1:nrow(acc)){
  
  acc[k,2] <- mean( test_dt$mpg01 != knn(train = train_dt[,-grep("mpg01", names(train_dt))], 
                                         test = test_dt[,-grep("mpg01", names(train_dt))], 
                                         cl=train_dt$mpg01, 
                                         k=acc[k,1]) )
}
plot(acc, type='b')

knn_fitted_test <- knn(train = train_dt[,-grep("mpg01", names(train_dt))], 
                       test = test_dt[,-grep("mpg01", names(train_dt))], 
                       cl=train_dt$mpg01, k=acc[which.min(acc$miss_class_rate),1])  ## test data 예측 

mean(knn_fitted_test != test_dt$mpg01)



#################  3  ############################
dt <- read.csv("Boston.csv", 
               stringsAsFactors = T)
dim(dt)
head(dt)

### 변수생성
dt$crim01 <- as.factor(ifelse(dt$crim > median(dt$crim), 1, 0))
head(dt)


### 데이터 분할 
set.seed(1004)
train_id <- sample(1:nrow(dt), nrow(dt)*0.6)

names(dt)

train_dt <- dt[train_id,-c(1,2)]
test_dt <- dt[-train_id,-c(1,2)]


### 시각화 
table(train_dt$chas, train_dt$crim01) #
barplot(table(train_dt$chas, train_dt$crim01)) #

melt_dt <- melt(as.data.table(train_dt)[, -which(names(train_dt)=='chas'), with=F], id='crim01')
ggplot(melt_dt, aes(crim01, value)) + 
  geom_boxplot() +
  facet_wrap(variable~., scales = 'free_y') + labs(y='')+
  theme_bw()

ggplot(melt_dt, aes(value)) + 
  geom_histogram() +
  facet_wrap(variable~crim01, scales = 'free') + labs(y='')+
  theme_bw()

train_dt <- as.data.frame(train_dt)
t.test(rm ~ crim01, train_dt)


### LDA
lda.fit <- lda(crim01 ~ ., data = train_dt)

lda_fitted_test <- predict(lda.fit, newdata = test_dt)$class
mean(lda_fitted_test != test_dt$crim01)


### QDA
qda.fit <- qda(crim01 ~ ., data = train_dt)

qda_fitted_test <- predict(qda.fit, newdata = test_dt)$class
mean(qda_fitted_test != test_dt$crim01)

### Logistic
glm.fit <- glm(crim01 ~ ., data = train_dt, 
               family = binomial)
summary(glm.fit)


logistic_fitted_test_prob <- predict(glm.fit, newdata = test_dt,
                                     type='response')
logistic_fitted_test <- ifelse(logistic_fitted_test_prob>0.5, 1,0)
mean(logistic_fitted_test != test_dt$crim01)

### KNN
acc <- data.frame(k = seq(1,21,2),
                  miss_class_rate = 0)

for (k in 1:nrow(acc)){
  
  acc[k,2] <- mean( test_dt$crim01 != knn(train = train_dt[,-grep("crim01", names(train_dt))], 
                                         test = test_dt[,-grep("crim01", names(train_dt))], 
                                         cl=train_dt$crim01, k=acc[k,1]) )
}
plot(acc, type='b')

acc[which.min(acc$miss_class_rate),1]

knn_fitted_test <- knn(train = train_dt[,-grep("crim01", names(train_dt))], 
                       test = test_dt[,-grep("crim01", names(train_dt))], 
                       cl=train_dt$crim01,
                       k=acc[which.min(acc$miss_class_rate),1])  ## test data 예측 
mean(knn_fitted_test != test_dt$crim01)

