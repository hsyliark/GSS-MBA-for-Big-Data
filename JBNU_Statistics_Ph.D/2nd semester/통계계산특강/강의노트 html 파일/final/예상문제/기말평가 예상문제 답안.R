install.packages("tidyverse")
library(tidyverse)

topicdata <- read.csv("C:/Users/HSY/Desktop/data.csv",sep=",",header=T)
# topicdata <- read.csv("C:/Users/HSY/Desktop/utf/data.csv",sep=",",header=T,fileEncoding="utf-8") 

topic_doc  <-  function(topicdata,k,d,alpha=0.1) {
  #k: index of topic, d: index of doc
  count_ = dim(filter(topicdata,doc==d,topic==k))[1] + alpha 
  total_ = dim(filter(topicdata,doc==d))[1] + alpha * K # K : number of all topics 
  count_ / total_ 
} 

word_topic  <-  function(topicdata,w,k,d,alpha=0.1) { 
  #k: index of topic, w: index of word 
  word_ = filter(topicdata,doc==d)$word[w]
  count_ = dim(filter(topicdata,topic==k,word==word_))[1] + alpha
  total_ = dim(filter(topicdata,topic==k))[1] + alpha * W  # W=length(unique(data$word))
  count_ / total_  
}

topicprob <- function(topicdata,d,w) {
  rtn <- c()
  for(k in 1:K) {
    rtn[k] <- topic_doc(topicdata,k,d)* word_topic(topicdata,w,k,d)
  } 
  rtn/sum(rtn)
}

D <- max(topicdata$doc) # number of all documents 
K <- 3 # number of all topics 
W <- length(unique(topicdata$word))
doclen <- c() # length of documents
for(m in 1:D){
  doclen[m] <- dim(filter(topicdata,doc==m))[1]
}
doclen <- c(0,doclen)

# method 1
for(t in 1:400) {
  for(d in 1:D) {
    len <- dim(filter(topicdata,doc==d))[1] 
    for(w in 1:len) {
      i <- cumsum(doclen)[d]+w 
      prob_res <- cumsum(c(0,topicprob(topicdata[-i,],d,w)))
      topic_res <- c()
      U <- runif(1)
      for(j in 1:K) {
        topic_res[j] <- j*(U>=prob_res[j] & U<prob_res[j+1])
      }
      topicdata[i,]$topic <- sum(topic_res)
    }
  }
}
topicdata

# method 2
for(t in 1:100) {
  for(d in 1:D) {
    len <- dim(filter(topicdata,doc==d))[1] 
    for(w in 1:len) {
      i <- cumsum(doclen)[d]+w 
      prob_res <- topicprob(topicdata[-i,],d,w)
      topic_res <- c()
      for(j in 1:K) {
        topic_res[j] <- j*(prob_res[j]==max(prob_res))
      }
      topicdata[i,]$topic <- sum(topic_res)
    }
  }
}
topicdata

# method 3
for(t in 1:500) {
  i <- 0
  for(d in 1:D) {
    len <- dim(filter(topicdata,doc==d))[1] 
    for(w in 1:len) {
      i <- i+1
      prob_res <- cumsum(c(0,topicprob(topicdata[-i,],d,w)))
      topic_res <- c()
      U <- runif(1)
      for(j in 1:K) {
        topic_res[j] <- j*(U>=prob_res[j] & U<prob_res[j+1])
      }
      topicdata[i,]$topic <- sum(topic_res)
    }
  }
}
topicdata

write.csv(topicdata,file="C:/Users/HSY/Desktop/202055364.csv")