install.packages("survival")
library(survival)
data(retinopathy)
head(retinopathy, 20)
install.packages("fitdistrplus")
install.packages("tidyverse")
library(tidyverse)
library(fitdistrplus)
install.packages("dplyr")
library(dplyr)

fcdata <- retinopathy %>% 
  mutate(left=ifelse(status==1, futime, NA),
         right=futime)%>%
  dplyr::select(left,right)
head(fcdata)


fcdata <- as.data.frame(fcdata)
(fit <- fitdistcens(fcdata, "weibull"))
(alphahat <- fit$estimate[1])
(lambdahat <- 1/fit$estimate[2])
fit$sd
fit$estimate[1]-1.96*fit$sd[1] 

fit$estimate[1]+1.96*fit$sd[1]
