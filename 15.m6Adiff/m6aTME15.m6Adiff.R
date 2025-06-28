######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#???ð?
library(limma)
library(reshape2)
library(ggpubr)
expFile="TCGA.TPM.txt"      #????????????
geneFile="gene.txt"         #m6A?????б?
setwd("C:\\Users\\sun\\Desktop\\m6a-esca\\15.m6Adiff")       #???ù???Ŀ¼

#??ȡ?????ļ?
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#??ȡm6A????????��
gene=read.table(geneFile, header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), row.names(data))
data=data[sameGene,]

#????????????Ŀ
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       #????????Ʒ??Ŀ
treatNum=length(group[group==0])     #????????Ʒ??Ŀ
sampleType=c(rep(1,conNum), rep(2,treatNum))

#??????ת????ggplot2?????ļ?
exp=log2(data+1)
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=sampleType)
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#????????ͼ
p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
	     ylab="Gene expression",
	     xlab="",
	     legend.title="Type",
	     palette = c("blue", "red"),
	     width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#????????ͼ
pdf(file="boxplot.pdf", width=7, height=5)
print(p1)
dev.off()


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056
