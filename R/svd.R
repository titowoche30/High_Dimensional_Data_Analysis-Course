library(rafalib)
library(MASS)
n <- 100
set.seed(42)
y <- t(mvrnorm(n,c(0,0), matrix(c(1,0.95,0.95,1),2,2)))
LIM <- c(-3.5,3.5)
mypar()
plot(y[1,],y[2,],xlim=LIM,ylim=LIM)

s <- svd(y)   #s eh uma lista com 3 elementos, sendo eles d,u e v
PC1 <- s$d[1]*s$v[,1] #PC1 eh o d[1] * v[:,1]
PC2 <- s$d[2]*s$v[,2] #PC2 eh o d[2] * v[:,2]
plot(PC1,PC2,xlim=LIM,ylim=LIM)

library(tissuesGeneExpression)
data(tissuesGeneExpression)
dim(e)
ind <- sample(nrow(e),500)
X <- t(apply(e[ind,], 1, scale))
s <- svd(X)

U <- s$u
V <- s$v
D <- diag(s$d)

X_hat <- U %*% D %*% t(V)
resid <- X - X_hat
max(abs(resid))

plot(s$d)
k <- sum(s$d<=1e-11)

new_k <- ncol(U) - k
X_hat2 <- U[,1:new_k] %*% D[1:new_k,1:new_k] %*% t(V[,1:new_k])
resid <- X - X_hat2
max(abs(resid))

plot(s$d^2 / crossprod(s$d,s$d))

#vamo mandar metade das colunas embora(95)
k <- 95
new_k <- ncol(U) - k
X_hat2 <- U[,1:new_k] %*% D[1:new_k,1:new_k] %*% t(V[,1:new_k])
resid <- X - X_hat2
max(abs(resid))
boxplot(resid,ylim=LIM)

#Porcentagem de variância que mantemos
1-var(as.vector(resid)/var(as.vector(X)))
var(as.vector(X_hat2))
sum(s$d[1:new_k]^2)/crossprod(s$d,s$d)
