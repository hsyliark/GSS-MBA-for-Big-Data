#============================================================
# title: "연관규칙분석 - Apriori"
#============================================================

#============================================================
# 패키지 설치
#============================================================

# install.packages("arules")
# install.packages("arulesViz")
library(arules)
library(arulesViz)

#============================================================
# 예제1 - 셔츠, 넥타이 구매
#============================================================
#Input data
Shopping_List <- list(c("넥타이","셔츠","양말"),
                      c("양말","벨트","장갑","셔츠"),
                      c("지갑","넥타이","셔츠"),
                      c("양말","벨트","장갑","바지"))

#A collection of itemsets represented as a binary incidence matrix. 
Shopping_itemMatrix <- as(Shopping_List,"itemMatrix")
Shopping_df<- as(Shopping_itemMatrix, "matrix")
rownames(Shopping_df) <- paste0("고객", 1:4)
Shopping_df
apply(Shopping_df,2,as.numeric)

#Convert to transaction data
Shopping_trans <- as(Shopping_df, "transactions")
inspect(Shopping_trans)

# 지지도가 0.5 이상인 item들의 빈도 막대그래프
itemFrequencyPlot(Shopping_trans, support = 0.5, cex.names=0.8)

# 지지도가 0.5 이상인 연관규칙
rules <- apriori(Shopping_trans,
                 parameter = list(minlen=2, supp=0.5, conf=0.01))
inspect(rules)

# lhs가 셔츠인 연관규칙 
rules.target <- subset(rules, lhs %oin% "셔츠")
inspect(sort(rules.target, by="lift"))

# lhs가 셔츠이고, rhs가 넥타이인 연관규칙
rules.target <- subset(rules, lhs %oin% "셔츠" & rhs %oin% "넥타이")
inspect(rules.target)


#============================================================
# 예제1-1 - 주류 구매
#============================================================
#Input data
Shopping_List <- list(c("소주","콜라","맥주"),
                      c("소주","콜라","와인"),
                      c("소주","주스"),
                      c("콜라","맥주"),
                      c('소주','콜라','맥주','와인'),
                      c('주스'))

#A collection of itemsets represented as a binary incidence matrix. 
Shopping_itemMatrix <- as(Shopping_List,"itemMatrix")
Shopping_df<- as(Shopping_itemMatrix, "matrix")
rownames(Shopping_df) <- paste0("고객", 1:6)
Shopping_df
apply(Shopping_df,2,as.numeric)

#Convert to transaction data
Shopping_trans <- as(Shopping_df, "transactions")
inspect(Shopping_trans)

# 지지도가 0.5 이상인 item들의 빈도 막대그래프
itemFrequencyPlot(Shopping_trans, 
                  support = 0.5, 
                  cex.names=0.8)

# 지지도가 0.5 이상인 연관규칙
rules <- apriori(Shopping_trans,
                 parameter = list(minlen=2, supp=0.5, conf=0.01))
inspect(rules)

# lhs가 맥주인 연관규칙 
rules.target <- subset(rules, lhs %in% "맥주")
inspect(sort(rules.target, by="confidence"))

# lhs가 맥주이고, rhs가 콜라인 연관규칙
rules.target <- subset(rules, lhs %in% "맥주" & rhs %in% "콜라")
inspect(rules.target)


#============================================================
# 예제2 : Adult data
#============================================================
# 미국 Census Bureau의 Census Income데이터 베이스에 추출한 설문조사 자료
# 관측치의 개수 : 48843개
# 나이, 직업군, 교육정도 등의 주로 범주형인 15개의 변수 포함

data("AdultUCI")
dim(AdultUCI)
AdultUCI[1:2,]

#변수제거 
AdultUCI[["fnlwgt"]] <- NULL
AdultUCI[["education-num"]] <- NULL

#연속형 변수를 범주형 변수로 변환 
AdultUCI[["age"]] <- ordered(cut(AdultUCI[[ "age"]], 
                                 c(15,25,45,65,100)),
                             labels = c("Young", "Middle-aged",
                                        "Senior", "Old"))
AdultUCI[["hours-per-week"]] <- ordered(cut(AdultUCI[[ "hours-per-week"]],
                                        c(0,25,40,60,168)),
                                        labels = c("Part-time", "Full-time", 
                                                   "Over-time", "Workaholic"))
AdultUCI[["capital-gain"]] <- ordered(cut(AdultUCI[[ "capital-gain"]],
                                      c(-Inf,0,median(AdultUCI[[ "capital-gain"]][AdultUCI[[ "capital-gain"]]>0]),Inf)),
                                      labels = c("None", "Low", "High"))
AdultUCI[["capital-loss"]] <- ordered(cut(AdultUCI[[ "capital-loss"]],
                                      c(-Inf,0, median(AdultUCI[[ "capital-loss"]][AdultUCI[[ "capital-loss"]]>0]),Inf)),
                                      labels = c("none", "low", "high"))

head(AdultUCI)

Adult <- as(AdultUCI, "transactions")
inspect(Adult[1:2,])



#============================================================
# 예제3 : Adult data
#============================================================

# (Q) 어떤 규칙을 갖는 사람들이 남성일까?
data(Adult)
summary(Adult)

# 지지도(support)가 0.4이상인 item들의 빈도 막대그래프
itemFrequencyPlot(Adult, 
                  support = 0.4, 
                  main = "Item Frequency Plot above support 0.4")

# 지지도 기준 상위 10개 item들의 빈도 막대그래프
itemFrequencyPlot(Adult, 
                  topN = 20, 
                  main = "Histogram of support top 10 items")

## the following example compares the item frequencies
## of Male in the data set
Adult.male <- Adult[Adult %in% "sex=Male"]

## simple plot
itemFrequencyPlot(Adult.male,
                  support = 0.2)   ### lhs='Sex=Male' 인 경우의 신뢰도 

## plot lift ratio (frequency in x / frequency in population)
## for items with a support of 20% in the population
itemFrequencyPlot(Adult.male, 
                  population = Adult, support = 0.2, 
                  lift = TRUE, horiz = TRUE)

#=============================================================
#  연관규칙분석
#============================================================

# 지지도가 0.4 이상인 연관규칙
rules <- apriori(Adult, parameter = list(support = 0.4))
summary(rules)
inspect(rules)

# 지지도가 0.4이상이면서 향상도가 1.3 이상인 것
rules.sub <- subset(rules, 
                    subset = rhs %pin% "sex" & lift > 1.3)
inspect(sort(rules.sub, by="lift"))
# 위 결과로는 “시민권자와 결혼하였으면 남성이다.” , “남편이면 남성이다”, “시민권자와 결혼했고 남편이면 남성이다.”라는 것들을 찾음


rules.sub <- subset(rules, 
                    subset = lift > 1.5)
inspect(sort(rules.sub, by="lift"))


# 모든 규칙의 산점도
# x축: 지지도, y축: 향상도, 점의 색: 신뢰도
plot(rules, 
     measure = c("support", "lift"), 
     shading = "confidence")




