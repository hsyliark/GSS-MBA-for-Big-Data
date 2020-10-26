library(survival)
library(survminer) 
library(gridExtra)
library(tidyverse)
library(fitdistrplus)
pro5 <- tibble( time = c(3,4,5,6,6,8,11,14,15,16),
                status= c(1,1,0,1,0,0,1,1,1,0))
pro5_1 <- pro5 %>% 
  mutate(left=ifelse(status==1, time, NA),
         right=time)%>%
  dplyr::select(left,right)
pro5_1 <- as.data.frame(pro5_1)
(fit <- fitdistcens(pro5_1, "weibull"))
c(fit$estimate[1]-1.96*fit$sd[1],fit$estimate[1]+1.96*fit$sd[1])
c(fit$estimate[2]-1.96*fit$sd[2],fit$estimate[2]+1.96*fit$sd[2])

fit1 <- survfit(Surv(time,status)~1,type="kaplan", data=pro5)
summary(fit1)
fit2 <- survfit(Surv(time,status)~1,type="fleming", data=pro5)
summary(fit2)
ggsurvplot(fit1, conf.int=TRUE,
           xlab = 'Time', 
           ylab = 'Survival Probability',
           title='Kaplan-Meier Survival Curve')
ggsurvplot(fit1, conf.int=TRUE,
           xlab = 'Time', 
           ylab = 'Survival Probability',
           title='Nelson-Aalen Survival Curve')






