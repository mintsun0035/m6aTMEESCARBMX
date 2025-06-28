

#???ð?
rm(list = ls())

setwd("C:\\Users\\sun\\Desktop\\m6a-esca\\17.m6aExp")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)           #???ð?
# expFile="merge_GSE13898_GSE53622_GSE53624.txt"      #?????????ļ?
# expFile="merge_GSE13898_GSE53622_GSE53624_raw.txt"      #?????????ļ?
expFile="merge_GSE13898.txt"      #?????????ļ?

geneFile="gene.txt"      #?????б??ļ?
# setwd("D:\\biowolf\\m6aTME\\17.m6aExp")    #???ù???Ŀ¼
# setwd("D:\\biowolf\\m6aTME\\17.m6aExp")    #???ù???Ŀ¼
setwd("C:\\Users\\sun\\Desktop\\m6a-esca\\17.m6aExp")    #???ù???Ŀ¼


# 00.00 -------------------------------------------------------------------


#??ȡ?????ļ??????????ݽ??д???--------------------------------------
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#??ȡm6a?????ı???��
gene=read.table(geneFile, header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

#????????
out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file="m6aGeneExp.txt",sep="\t",quote=F,col.names=F)


