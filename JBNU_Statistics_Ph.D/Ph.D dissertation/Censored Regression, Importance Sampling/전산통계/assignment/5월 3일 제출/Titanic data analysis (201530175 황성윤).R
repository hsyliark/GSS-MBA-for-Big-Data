### Statistical Computing third assignment



### 1. Binary regression

titanic <- as.data.frame(Titanic)
titanic.no <- titanic[(1:16),-4]
colnames(titanic.no) <- c("Class","Sex","Age","Death") 
Survive <- titanic[(17:32),5]
titanic.binary <- cbind(titanic.no,Survive)
titanic.binary
attach(titanic.binary)

## logistic link

logit1 <- glm(cbind(Survive,Death)~Class+Sex+Age, family=binomial(link="logit"))
summary(logit1)

## probit link

probit1 <- glm(cbind(Survive,Death)~Class+Sex+Age, family=binomial(link="probit"))
summary(probit1)

## complementary log-log link

cloglog1 <- glm(cbind(Survive,Death)~Class+Sex+Age, family=binomial(link="cloglog"))
summary(cloglog1)

detach(titanic.binary)



### 2. poisson regression

titanic.poisson <- titanic[-c(1,2,4,5,6,8,20,24),]
attach(titanic.poisson)

## Question 1

loglink1 <- glm(Freq~Class*Sex*Age*Survived, family=poisson(link="log"))
summary(loglink1)

## Question 2,3

install.packages('epiDisplay')
library(epiDisplay)

loglink.1 <- glm(Freq~Class+Sex+Age+Survived, family=poisson(link="log"))
loglink.s <- glm(Freq~Class*Sex*Age*Survived, family=poisson(link="log"))
lrtest(loglink.1,loglink.s)

loglink.2 <- glm(Freq~Class+Sex+Age+Survived+Class:Sex+
                   Class:Age+Class:Survived+Sex:Age+Sex:Survived
                 +Age:Survived, family=poisson(link="log"))
lrtest(loglink.2,loglink.s)

loglink.3 <- glm(Freq~Class+Sex+Age+Survived+Class:Sex+
                   Class:Age+Class:Survived+Sex:Age+Sex:Survived
                 +Age:Survived+Class:Sex:Age
                 +Class:Sex:Survived, family=poisson(link="log"))
lrtest(loglink.3,loglink.s)

detach(titanic.poisson)
