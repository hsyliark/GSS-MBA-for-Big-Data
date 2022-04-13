####################################################
##### Knn
##### 분류모형 
####################################################
library(class) ##knn


###################################################
###################################################
###################################################
###################################################
set.seed(1234)
train_id <- sample(1:nrow(iris), nrow(iris)*0.7)

train_dt_x <- iris[train_id,1:2]
train_dt_y <- iris[train_id,5]
test_dt_x <- iris[-train_id,1:2]
test_dt_y <- iris[-train_id,5]

pred_test <- knn(train_dt_x, test_dt_x, 
                 cl=train_dt_y, k=1, prob = T)  ## test data 예측 

table(test_dt_y, pred_test)  ## confusion matrix
mean(pred_test!=test_dt_y)  # 오분류율

pred_test <- knn(train_dt_x, test_dt_x, 
                 cl=train_dt_y, k=2, prob = T)  ## test data 예측 

table(test_dt_y, pred_test)  ## confusion matrix
mean(pred_test!=test_dt_y)  # 오분류율

pred_test <- knn(train_dt_x, test_dt_x, 
                 cl=train_dt_y, k=5, prob = T)  ## test data 예측 

table(test_dt_y, pred_test)  ## confusion matrix
mean(pred_test!=test_dt_y)  # 오분류율


pred_test <- knn(train_dt_x, test_dt_x, 
                 cl=train_dt_y, k=10, prob = T)  ## test data 예측 

table(test_dt_y, pred_test)  ## confusion matrix
mean(pred_test!=test_dt_y)  # 오분류율


#### 최적의 k찾기

k.grid <- seq(1,21,2)
acc_df <- data.frame(k = k.grid, acc = NA)

for(i in 1:length(k.grid)){
  pred_test <- knn(train_dt_x, test_dt_x, 
                   cl=train_dt_y, k=k.grid[i])  ## test data 예측 
  acc_df[i,2] <- mean(pred_test!=test_dt_y)  # 오분류율
}

acc_df$k1 <- 1/acc_df$k
plot(acc_df$k1, acc_df$acc, type='b')

plot(acc_df$k, acc_df$acc, type='b')















