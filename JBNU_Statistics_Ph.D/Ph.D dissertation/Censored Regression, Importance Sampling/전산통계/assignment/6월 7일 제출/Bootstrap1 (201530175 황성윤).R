## Bootstrap

LSAT<-c(576,635,558,578,666,580,555,661,651,605,653,575,545,572,594)
GPA<-c(3.39,3.30,2.81,3.03,3.44,3.07,3.00,3.43,3.36,3.13,3.12,2.74,2.76,2.88,2.96)

# Q1

cor(LSAT, GPA)

# Q2

B <- 1000
BL <- matrix(0, B, length(LSAT))
BG <- matrix(0, B, length(GPA))
Bcor <- c()
set.seed(1)
for (i in 1:B) {
  idx <- sample(length(LSAT), replace=T)
  BL[i,] <- LSAT[idx]
  BG[i,] <- GPA[idx]
  Bcor[i] <- cor(BL[i,], BG[i,])
}
print(c(mean(Bcor), var(Bcor)))

# Q3

hist(Bcor, prob=T) ; lines(density(Bcor))

# Q4

TX <- cor(LSAT, GPA) ; Bv <- var(Bcor)
c(TX-qnorm(0.975)*sqrt(Bv), TX+qnorm(0.975)*sqrt(Bv))

bq <- quantile(Bcor, probs=c(0.025, 0.975))
c(2*TX-bq[2], 2*TX-bq[1])

c(bq[1], bq[2])