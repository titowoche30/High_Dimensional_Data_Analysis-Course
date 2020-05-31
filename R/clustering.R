library(rafalib)
library(tissuesGeneExpression)
data(tissuesGeneExpression)

# Hierarchichal Clustering - Dendogram

d <- dist(t(e))     # distance between genes(dist calculates distance between rows)
hc <- hclust(d)
class(hc)
plot(hc,cex=.5,label=tissue)    #plot the tissue that each gene belongs

mypar2(1,1)
#fumeric turns into factors and then into numbers
myplclust(hc,cex=.5,label=tissue,lab.col=as.fumeric(tissue))
abline(h=120)

clusters = cutree(hc,h=120)
table(true=tissue,cluster=clusters)


#Exercises
set.seed(1)
m = 10000
n = 24
x = matrix(rnorm(m*n),m,n)
colnames(x)=1:n

d_ <- dist(t(x))    
hc_ <- hclust(d_)
class(hc_)
plot(hc_,cex=.5)

#k-means

km <- kmeans( t(e),centers = 7,algorith="Lloyd")    #distance between genes
table(tissue,clusters=km$cluster)

d <- dist(t(e))
mds <- cmdscale(d)
plot(mds[,1],mds[,2],col=km$cluster)

#Heatmaps

#library("BiocManager")
#BiocManager::install("genefilter")
library(genefilter)

#top 40 vars
rv <- rowVars(e)
idx <- order(-rv)[1:40]

library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap(e[idx,],col=hmcol)









