library(rafalib)
library(tissuesGeneExpression)
data(tissuesGeneExpression)

#Getting all categories
group <- factor(tissue)

#Multi-Dimensional Scaling Plot by hand - PCoA
mat <- e
A = mat-rowMeans(mat)
s <- svd(A)
U <- as.matrix(s$u)
V <- as.matrix(s$v)
d <- s$d
#Can take any pair of columns, but 1,2 explains most of the variance
PC <- diag(d[1:2]) %*% t(V[,1:2])
#--or--
#PC1<-s$d[1]*s$v[,1]
#PC2<-s$d[2]*s$v[,2]
#------

mypar(1,1)
plot(t(PC),pch=21,bg=as.numeric(group),xlab="First dimension",ylab="Second dimension")
legend("bottomright",levels(group),col=seq(along=levels(group)),pch=19,cex=.6)
  
  #By R Function
dists <- dist(t(mat))
mds <- cmdscale(dists,k=2)
mypar2(1,1)
plot(mds,bg=as.numeric(group),pch=21,xlab="First dimension",ylab="Second dimension")
legend("bottomleft",levels(group),col=seq(along=levels(group)),pch=15,cex=.7)

#Variance explained by each column
plot(d^2/crossprod(d))

#Difference between PCoA and PCA
#PCA is used for quantitative variables, so the axes in graphic have a quantitative weight. 
#And the position of the samples are in relation with those weight. 
#On the other hand, PCoA is used when characters or variables are qualitative or discrete.
#In that case the bi- or tri- plot shows relations among different samples on a smaller 
#number of variables. If the distance matrix of the input of PCoA is calculated by the
#euclidean distance, both methods are equivallent