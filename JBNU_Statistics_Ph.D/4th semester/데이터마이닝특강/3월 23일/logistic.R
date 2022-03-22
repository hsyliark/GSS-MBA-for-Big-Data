##############  경로지정
setwd("C:\\R-Project\\DAT\\INtroduction SL")

######################################################
######## 데이터 불러오기 
######################################################
dt <- read.csv("Default.csv", stringsAsFactors = T)
head(dt)
dim(dt)


#### 산점도와 상자그림 
def.par <- par(no.readonly = TRUE) # save default, for resetting...

layout(matrix(c(1,1,2,3),1,4, byrow = TRUE))
plot(dt$balance, dt$income,type='n',
     xlab='Balance',
     ylab = 'Income', cex.lab=1.5)
points(dt[dt$default=='No',]$balance, dt[dt$default=='No',]$income,col='steelblue')
points(dt[dt$default=='Yes',]$balance, dt[dt$default=='Yes',]$income,col='darkorange', pch=3)


boxplot(balance~default, data=dt,
        xlab = 'Default', ylab='Balance', cex.lab=1.5, col=c('steelblue', 'darkorange'), cex.axis=1.3)
boxplot(income~default, data=dt,
        xlab = 'Default', ylab='Income', cex.lab=1.5, col=c('steelblue', 'darkorange'), cex.axis=1.3)

par(def.par)
pairs(dt)

#### Logistic regression
glm.fits <- glm(
  default ~ balance,
  data = dt, 
  family = binomial
)

summary(glm.fits)

coef(glm.fits)
summary(glm.fits)$coef
summary(glm.fits)$coef[, 4]

contrasts(dt$default)
contrasts(dt$student)

glm.fits$fitted.values
head(glm.fits$fitted.values)
head(predict(glm.fits))
predict(glm.fits, type='response')
head(predict(glm.fits, type='response'))


plot(dt$balance[order(dt$balance)], 
     glm.fits$fitted.values[order(dt$balance)],  
     type='l', col='darkorange', lwd=2,ylim=c(-0.1,1.1),
     xlab='Balance', ylab='Probability of Default')
abline(h=c(0,1), lty=2)


### classification : cut value 
dt$fitted_class <- ifelse(predict(glm.fits, type='response')>0.5,'Yes', 'No')
head(dt)

table(dt$default, dt$fitted_class)
addmargins(table(dt$default, dt$fitted_class))
(9625+100)/10000 * 100
mean(dt$fitted_class==dt$default)


dt$fitted_class <- ifelse(predict(glm.fits, type='response')>0.4,'Yes', 'No')
table(dt$default, dt$fitted_class)
addmargins(table(dt$default, dt$fitted_class))

mean(dt$fitted_class==dt$default)


##################################

glm.fits2 <- glm(
  default ~ student,
  data = dt, family = binomial
)

summary(glm.fits2)

glm.fits3 <- glm(
  default ~ balance + income + student,
  data = dt, family = binomial
)

summary(glm.fits3)




### plot
dt2 <- dt[order(dt$balance),]

glm.fits4 <- glm(
  default ~ balance + income + student,
  data = dt2, family = binomial
)

dt2$fitted <- glm.fits4$fitted.values

par(mfrow=c(1,2))
plot(dt2$balance, dt2$fitted,  
     type='n',ylim=c(-0.1,1.1),
     xlab='Balance', ylab='Probability of Default', cex.axis=0.7)
lines(dt2[dt2$student=='Yes',]$balance, dt2[dt2$student=='Yes',]$fitted, 
      col='steelblue')
lines(dt2[dt2$student=='No',]$balance, dt2[dt2$student=='No',]$fitted, col='darkorange')
abline(h=0)
abline(h=mean(dt2[dt2$student=='Yes',]$fitted), lty=2, col='steelblue')
abline(h=mean(dt2[dt2$student=='No',]$fitted), lty=2, col='darkorange')

boxplot(balance~student, data=dt,
        xlab = 'Student Status', ylab='Credit Card Balnce', 
        col=c('darkorange', 'steelblue'), cex.axis=0.7)
par(mfrow=c(1,1))


###################################



################################################################
# 남미의 심장병 발생 확률이 높은 지역의 남성에 대한 retrospective 조사를 한 데이터
################################################################
SAheart <- read.csv("SAheartdata.csv", stringsAsFactors = T)
head(SAheart)
SAheart$chd <- as.factor(SAheart$chd)
SAheart$famhist <- as.factor(SAheart$famhist)

heart <- SAheart[c(2,5,8:10)]

head(heart)

# tobacco: 누적 흡연량
# famhist: 유전도
# alcohol:알코올
# age:나이
# chd: 심장병 유무

par(mfrow=c(1,3))
# 심장병 유무 v.s 흡연량
boxplot(tobacco~chd, data=heart,
        xlab = 'chd', ylab='tobacco', cex.lab=1.5, col=c('steelblue', 'darkorange'), cex.axis=1.3)
# 심장병 유무 v.s 유전
boxplot(alcohol~chd, data=heart,
        xlab = 'chd', ylab='alcohol', cex.lab=1.5, col=c('steelblue', 'darkorange'), cex.axis=1.3)
# 심장병 유무 v.s 나이
boxplot(age~chd, data=heart,
        xlab = 'chd', ylab='age', cex.lab=1.5, col=c('steelblue', 'darkorange'), cex.axis=1.3)
par(mfrow=c(1,1))

table(heart$famhist, heart$chd)

barplot(table(heart$chd, heart$famhist), 
        col=c('steelblue', 'darkorange'), 
        beside=T,
        xlab='chd', ylab='famhist')

barplot(table(heart$famhist, heart$chd), 
        col=c('steelblue', 'darkorange'), 
        beside=T,
        xlab='famhist', ylab='chd')


#============================================================
# 모형 적합
#============================================================

heartfit <- glm(chd ~ ., data = heart, family = binomial)
round(coef(summary(heartfit)), 3)
# tobacco, famhist1, age는 p-value<0.05이므로 유의수준 5%에서 유의
# alcohol변수는 p-value가 0.984로 0.5보다 크므로 유의수준 5%에서 유의하지 않음

# 변수선택
new_heartfit <- step(heartfit)
# alcohol변수가 제거됨

# 회귀계수 추정
round(coef(summary(new_heartfit)), 3)
# tobacco가 한단위 증가하면 심장병에 걸릴 오즈가 exp(0.08) 증가
# famhist1가 한단위 증가하면 심장병에 걸릴 오즈가 exp(0.97) 증가
# age가 한단위 증가하면 심장병에 걸릴 오즈가 exp(0.048) 증가

# cut-off가 0.7기준일 경우
y.pred1 = ifelse(heartfit$fitted.values > .7, 1, 0)
table(heart$chd, y.pred1)
# 정분류율
round(mean(heart$chd == y.pred1),2)
