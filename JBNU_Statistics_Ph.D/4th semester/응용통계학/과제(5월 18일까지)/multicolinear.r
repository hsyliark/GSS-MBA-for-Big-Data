Nrep = 200

n = 1000

hatb1 <- hatb2 <- c()
te1 <- te2 <- c()

for (k in 1:Nrep)
{
x1 = runif(n,-3,3)
x2 = x1 + rnorm(n,0,0.01)
#X = model.matrix(~ x1 + x2)
#solve(t(X)%*%X)
y = 2 + x1 + 2*x2 + rnorm(n)

ind = sample(1:n,500)
tx1 = x1[ind]
tx2 = x2[ind]
test_x1 = x1[-ind]
test_x2 = x2[-ind]
ty  = y[ind]
test_y = y[-ind]

fit1 = lm(ty~tx1+tx2)
fit2 = lm(ty~tx1)

hatb1[k] = fit1$coefficients[2]
hatb2[k] = fit2$coefficients[2]

te1[k] = mean((test_y - predict(fit1,newdata=data.frame(test_x1,test_x2)))^2)
te2[k] = mean((test_y - predict(fit2,newdata=data.frame(test_x1)))^2)
#summary(lm(y~x1+x2))
print(k)
}

c(mean(te1),mean(te2))
c(mean((hatb1-1)^2),mean((hatb2-2)^2))

