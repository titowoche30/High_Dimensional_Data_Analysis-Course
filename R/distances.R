'''
library(devtools)
install_github("genomicsclass/tissuesGeneExpression")
'''

# library(tissuesGeneExpression)
# data(tissuesGeneExpression)
# dim(e) ##e contains the expression data

tissue
table(tissue)

x <- e[,1]
y <- e[,2]
z <- e[,87]

sqrt(sum((x-y)^2))
sqrt(sum((x-z)^2))

sqrt(crossprod(x-y))
sqrt(crossprod(x-z))

d <- dist(t(e))                 #dist faz a distância entre cada linha, como quero de cada coluna, meto o transposto
d <- as.matrix(d)
d[1,2]
d[1,8]
d[1,87]
image(d)

#-----exercícios-----
sum(tissue == "hippocampus")

d[3,45]

a<- e["210486_at",]
b<- e["200805_at",]

sqrt(crossprod(a-b))


