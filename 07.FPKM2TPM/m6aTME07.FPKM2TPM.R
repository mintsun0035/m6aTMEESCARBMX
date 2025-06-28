
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
library(limma)             #引用包
inputFile="symbol.txt"     #输入文件
setwd("C:\\Users\\sun\\Desktop\\m6a-esca\\07.FPKM2TPM")     #工作目录
#??ȡ?????ļ?,?????????ļ?????
outTab=data.frame()
rt=read.table(inputFile, header=T, sep="\t", check.names=F)#读取输入文件
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)#R语⾔中创建矩阵的基本语法是：
#matrix(data, nrow, ncol, byrow, dimnames)
#以下是所使的参数的说明：
#data - 是这成为矩阵的数据元素输向量。
#nrow - 是要创建的行数。
#ncol - 要被创建的列的数⽬。
#byrow - 是一个合乎逻辑。如果为True，那么输⼊向量元素在安排的⾏。
#dimname - 是分配给行和列名称。
data=avereps(data)#有的基因出现过多行,把出现多行的gene取平均值
#FPKM转化为TPM
fpkmToTpm=function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpm=apply(data, 2, fpkmToTpm)#每一列数据都进行上述操作exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
#tpm=apply(data, 1, fpkmToTpm)每一行数据进行处理
#输出文件
tpmOut=rbind(ID=colnames(tpm), tpm)#进行行的叠加，列数必须相同
write.table(tpmOut, file="TCGA.TPM.txt", sep="\t", col.names=F, quote=F)


