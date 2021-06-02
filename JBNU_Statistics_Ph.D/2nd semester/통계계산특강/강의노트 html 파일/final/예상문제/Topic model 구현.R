## Topic model 구현

install.packages("tidyverse")
library(tidyverse)

doclen <- 5 
word <- c('손흥민','골','골','박지성','패스',
          '골','확률','패스','손흥민','골',
          '골','골','확률','패스','골',
          '골','박지성','통계','확률','골',
          '확률','통계','박지성','통계','골',
          '확률','확률','골','통계','AI',
          '확률','확률','통계','통계','AI',
          '확률','빅데이터','데이터과학','빅데이터','AI',
          '확률','데이터과학','AI','데이터과학','데이터과학')
doc <- c(rep(1,doclen),rep(2,doclen),rep(3,doclen),
         rep(4,doclen),rep(5,doclen),rep(6,doclen),
         rep(7,doclen),rep(8,doclen),rep(9,doclen))
topic <- c(2,1,1,2,2,
           1,2,1,2,2,
           1,2,1,2,2,
           1,2,1,2,2,
           1,2,1,2,2,
           1,2,2,2,2,
           1,2,1,1,2,
           1,2,1,1,2,
           2,1,1,2,2)
topicdata <- tibble(word=word,doc=doc,topic=topic)
topicdata

topic_doc  <-  function(topicdata,k,d,alpha=0.1){
  #k: 토픽의 인덱스, d:doc의 인덱스
  count_ = dim(filter(topicdata,doc==d,topic==k))[1] + alpha 
  total_ = dim(filter(topicdata,doc==d))[1] + alpha * K # K는 전체 토픽의 수 
  count_ / total_ 
  } 

word_topic  <-  function(topicdata,w,k,d,alpha=0.1){ 
  #k: 토픽의 인덱스, w: word의 인덱스 
  word_ = filter(topicdata,doc==d)$word[w]
  count_ = dim(filter(topicdata, topic==k,word==word_))[1] + alpha
  total_ = dim(filter(topicdata, topic==k))[1] + alpha * W  # W=length(unique(data$word))
  count_ / total_ 
  }

topicprob <- function(topicdata,d,w){
  rtn <- c()
  for(k in 1:K){
    rtn[k] <- topic_doc(topicdata,k,d)* word_topic(topicdata,w,k,d)
    } 
  rtn/sum(rtn)
  }

D <- max(topicdata$doc) # 다큐먼트의 수 
doclen <- dim(filter(topicdata,doc==D))[1] # 다큐먼트의 길이
K <- 2 # 토픽의 수 
W <- length(unique(topicdata$word))
for(t in 1:50){
  for(d in 1:D){
    for(w in 1:doclen){
      i <- (d-1)*doclen+w 
      topicdata[i,]$topic <- which(topicprob(topicdata[-i,],d,w)==max(topicprob(topicdata[-i,],d,w)))
    }
  }
}
data

  