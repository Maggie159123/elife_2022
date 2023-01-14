#! /usr/bin/Rscript
args <- commandArgs(T)

#library(gplots)

data1 <- t(as.matrix(read.table(file=args[1], header=T, sep='\t', row.names=1)))
data2 <- t(as.matrix(read.table(file=args[2], header=T, sep='\t', row.names=1)))
dim(data1)
dim(data2)
cor_matrix <- cor(data1,data2,method = "pearson")
write.table(cor_matrix,file=args[3],quote=F,sep='\t')
#pheatmap(cor_matrix)
dev.off()
