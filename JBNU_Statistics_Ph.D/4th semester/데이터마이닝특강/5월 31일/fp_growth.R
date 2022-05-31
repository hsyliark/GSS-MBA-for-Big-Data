#============================================================
# title: "연관규칙분석 - fp_growth"
#============================================================

#============================================================
# 패키지 설치
#============================================================
#install.packages("rCBA")
library(rCBA)
library(data.table)

#============================================================
# 예제1
#============================================================

Shopping_List <- list(c('A','B','E'),
                      c('B','D'),
                      c('B','C'),
                      c('A','B','D'),
                      c('A','C'),
                      c('B','C'),
                      c('A','C'),
                      c('A','B','C','E'),
                      c('A','B','C'))

#A collection of itemsets represented as a binary incidence matrix. 
Shopping_itemMatrix <- as(Shopping_List,"itemMatrix")
Shopping_df<- as(Shopping_itemMatrix, "matrix")
rownames(Shopping_df) <- paste0("T", seq(10,90,10))
Shopping_df
Shopping_df <- apply(Shopping_df,2,as.numeric)

Shopping_df1 <- as.data.table(Shopping_df)
Shopping_df1[, A := ifelse(A==1, "1", NA)]
Shopping_df1[, B := ifelse(B==1, "1", NA)]
Shopping_df1[, C := ifelse(C==1, "1", NA)]
Shopping_df1[, D := ifelse(D==1, "1", NA)]
Shopping_df1[, E := ifelse(E==1, "1", NA)]
Shopping_df1

#Convert to transaction data
Shopping_trans <- as(Shopping_df1, "transactions")
inspect(Shopping_trans)



# FP-Growth Algorithm
rules = fpgrowth(Shopping_trans, 
                 support = 2/9, 
                 confidence = 0.03, 
                 maxLength=5,
                 consequent='B')
inspect(rules)

itemFrequencyPlot(Shopping_trans, support=2/9)

rules = fpgrowth(Shopping_trans, 
                 support = 2/9, 
                 confidence = 0.03, 
                 maxLength=5,
                 consequent = 'C')
inspect(rules)



