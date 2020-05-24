library(rafalib)
library(tissuesGeneExpression)
data(tissuesGeneExpression)

#Getting all categories
group <- factor(tissue)


mat <- e
pca <- prcomp(t(mat),center=TRUE)   
s = svd(mat - rowMeans(mat))

pc1_s = s$d[1]*s$v[,1]
pc2_s = s$d[2]*s$v[,2]
mypar(1,1)
plot(pc1_s,pc2_s,pch=21,bg=as.numeric(group),xlab="First dimension",ylab="Second dimension")
legend("bottomright",levels(group),col=seq(along=levels(group)),pch=19,cex=.6)

pc1=pca$x[,1] 
pc2=pca$x[,2] 
mypar(1,1)
plot(pc1,pc2,bg=as.numeric(group),pch=21,xlab="First dimension",ylab="Second dimension")
legend("bottomleft",levels(group),col=seq(along=levels(group)),pch=15,cex=.7)
