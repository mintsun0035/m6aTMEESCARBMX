

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")

# install.packages('Rcpp')

rm(list = ls())
gc()
setwd("C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\28.PCA")
# setwd("C:\\Users\\123456\\Desktop\\ESCA_20220223_2026\\26.GSVA")      #设置工作目录



library(Rcpp)


#引用包
library(limma)
library(ggplot2)
expFile="m6aGeneExp.txt"         #表达输入文件
clusterFile="m6aCluster.txt"     #m6A分型文件
# setwd("D:\\biowolf\\m6aTME\\28.PCA")      #设置工作目录
# setwd("L:/124m6aTME分型_直肠癌m6aTME_雷佩捷202105/28.PCA")



#读取输入文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=t(data)

# sampleGene
dim(data)
# [1] 279  17


#PCA分析01-------------------
data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)
write.table(pcaPredict, file="newTab_scale.xls", quote=F, sep="\t")

#读取分型文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
m6Acluster=as.vector(cluster[,1])

#设置颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
m6aCluCol=bioCol[1:length(levels(factor(m6Acluster)))]

#可视化
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], m6Acluster=m6Acluster)
PCA.mean=aggregate(PCA[,1:2], list(m6Acluster=PCA$m6Acluster), mean)

pdf(file="PCA_scale.pdf", height=5, width=6.5)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = m6Acluster)) +
	scale_colour_manual(name="m6Acluster", values =m6aCluCol)+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$m6Acluster, cex=7)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()





#PCA分析02---------------
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
write.table(pcaPredict, file="newTab_02noScale.xls", quote=F, sep="\t")

#读取分型文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
m6Acluster=as.vector(cluster[,1])

#设置颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
m6aCluCol=bioCol[1:length(levels(factor(m6Acluster)))]

#可视化
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], m6Acluster=m6Acluster)
PCA.mean=aggregate(PCA[,1:2], list(m6Acluster=PCA$m6Acluster), mean)

pdf(file="PCA_02noScale.pdf", height=5, width=6.5)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = m6Acluster)) +
  scale_colour_manual(name="m6Acluster", values =m6aCluCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$m6Acluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()









# ================ -------------------------------------------------------

options(stringsAsFactors = F)
library(Rcpp)
library(stringr)
library(ggsci)
library(ggplot2)
library(Rtsne)
library(umap)
library(cowplot)
library(Rcpp)




# read in expression value
# data = read.table("datExp.normalized.txt",header = T,sep="\t",row.names = 1)
head(data)
# data=t(data)
head(data)

# read in phenotype
# pheno0 = read.table("target.txt",header = T,sep="\t")
# View(head(pheno0))
# pheno0[,2] = str_remove_all(pheno0[,2]," +$")

# pheno = pheno0[match(colnames(data),pheno0[,1]),]
# View(head(pheno))


## PCA
data=as.data.frame(data)
pca <- princomp(data)
# 您可以使用prcomp代替princomp     
# pca <- prcomp(data)

#scatter plot
# Group = pheno[,2]
# Group =m6Acluster
Group = cluster[,1]

# Group =factor(m6Acluster)

# Group =levels(factor(m6Acluster))
# [1] "A" "B" "C"
# Group =as.numeric(factor(m6Acluster))


m6aCluCol=bioCol[1:length(levels(factor(m6Acluster)))]
Group =m6aCluCol

p1 = ggplot(as.data.frame(pca$loadings[,1:2]),aes(Comp.1,Comp.2)) +
  geom_point(size=2, stat="identity") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA,size=0.5),
        panel.grid.minor = element_line(colour = "#cccccc", size=0.1),
        panel.grid.major = element_line(colour = "#cccccc", size=0.1)) +
  labs(title="PCA plot") +
  xlab("PC1") + ylab("PC2") +
  scale_color_lancet()
p1




## tSNE00-----------------
# 需要数据是样本基因矩阵---
# 需要数据是样本基因矩阵---
# tSNE对行来操作的
dim(data)
# [1] 122  20
# [1] 279  17

# Group = pheno[,2]
set.seed(1)

library(Rcpp)
# data=t(data)



# tsne_out <- Rtsne(t(data),theta = 0.1,perplexity = 3) # Run TSNE

# data是样本基因矩阵
tsne_out <- Rtsne(data,theta = 0.1,perplexity = 3) # Run TSNE




#scatter plot
plot.data_tsne <- as.data.frame(tsne_out$Y)
colnames(plot.data_tsne) = c("tSNE1","tSNE2")


# Group =as.numeric(factor(m6Acluster))
Group =factor(m6Acluster)
# Group =factor(m6Acluster,levels = c("#0066FF","#FF9900"))



p2 = ggplot(plot.data_tsne,aes(tSNE1,tSNE2, shape=Group)) + 
  geom_point(aes(colour=Group),size=2, stat="identity") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA,size=0.5), 
        panel.grid.minor = element_line(colour = "#cccccc", size=0.1), 
        panel.grid.major = element_line(colour = "#cccccc", size=0.1)) +
  labs(title="tSNE plot") +
  scale_color_lancet()
p2




## UMAP-------------
# Group = pheno[,2]
# umap.out = umap::umap(t(data))
umap.out = umap::umap(data)
#scatter plot
plot.data_umap <- as.data.frame(umap.out$layout)
colnames(plot.data_umap) = c("UMAP1","UMAP2")
p3 = ggplot(plot.data_umap,aes(UMAP1,UMAP2, shape=Group)) + 
  geom_point(aes(colour=Group),size=2, stat="identity") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA,size=0.5), 
        panel.grid.minor = element_line(colour = "#cccccc", size=0.1), 
        panel.grid.major = element_line(colour = "#cccccc", size=0.1)) +
  labs(title="UMAP plot") +
  scale_color_lancet()
p3




## output plot together
p.out = plot_grid(p1,p2,p3,labels = c("A","B","C"),align = "h",nrow = 1)
p.out
ggsave2("Figure1.pca_tsne_umap.pdf",p.out,device = "pdf",
        height = 4,width = 15,units = "in")
ggsave2("Figure1.pca_tsne_umap.tiff",p.out,device = "tiff",
        height = 4,width = 15,units = "in",dpi = 300)


save.image()


