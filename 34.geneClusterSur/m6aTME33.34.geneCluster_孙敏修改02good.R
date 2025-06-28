

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")

rm(list = ls())
gc()

setwd("C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\34.geneClusterSur")


#引用包
library(limma)
library(ConsensusClusterPlus)
expFile="uniSigGeneExp.txt"        #表达输入文件
 # expFile="uniSigGeneExp_p0.01.txt"        #表达输入文件


# workDir="D:\\资源\\课题\\00AAA直肠癌m6aTME_雷佩捷202105\\33.geneCluster"      #设置工作目录
# workDir="L:/124m6aTME分型_直肠癌m6aTME_雷佩捷202105/33.geneCluster_02"      #设置工作目录

workDir="C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\34.geneClusterSur"      #设置工作目录


# setwd("L:/124m6aTME分型_直肠癌m6aTME_雷佩捷202105/33.geneCluster_02")






# 00-------------- --------------------------------------------------------






setwd(workDir)       #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

dim(data)
# [1]  68 279


#聚类
maxK=9
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="km",
              distance="euclidean",
              seed=123456,
              plot="png")







# #####后面写成了for循环00 ---------------------------------------------------------
#输出分型结果
clusterNum=3        #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("geneCluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$geneCluster))
cluster$geneCluster=letter[match(cluster$geneCluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="geneCluster.txt", sep="\t", quote=F, col.names=F)
# #####后面写成了for循环01 ---------------------------------------------------------



# 00.01-------------------------------------------

# rm(list = ls())

# 01.00-------------------------------------------

#输出分型结果
#install.packages("survival")
#install.packages("survminer")


#引用包
library(survival)
library(survminer)
library(dplyr)


################################################################################
########################   try                               ###################
################################################################################

cliFile="time.txt"               #生存数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/30

summary(cli$futime)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   13.32   21.27   38.18   47.38  152.80 


#需要修改01#########################

cutData = seq(0,12,1)
pValueOutput = data.frame()

for (j in cutData) {
  cli = cli %>% dplyr::filter(futime >= j)                            
  #需要修改02#########################
  for (i in 2:3){
    # i =3
    clusterNum=i       #分几类，根据判断标准判断
    cluster=results[[clusterNum]][["consensusClass"]]
    cluster=as.data.frame(cluster)
    colnames(cluster)=c("geneCluster")
    letter=c("A","B","C","D","E","F","G")
    uniqClu=levels(factor(cluster$geneCluster))
    
    cluster$geneCluster=letter[match(cluster$geneCluster, uniqClu)]
    
    clusterOut=rbind(ID=colnames(cluster), cluster)
    write.table(clusterOut, file=paste0("geneCluster_", i,"_",j,".txt"), sep="\t", quote=F, col.names=F)
    
    # clusterFile=paste0("geneCluster_", i,".txt")    #m6A分型文件
    
    # setwd("D:\\biowolf\\m6aTME\\24.geneClusterSur")      #设置工作目录
    
    #读取输入文件
    # clusterOut=clusterOut 
    # cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
    rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
    
    
    
    #数据合并
    sameSample=intersect(row.names(cluster), row.names(cli))
    rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])
    
    #生存差异统计
    length=length(levels(factor(rt$geneCluster)))
    diff=survdiff(Surv(futime, fustat) ~ geneCluster, data = rt)
    pvalue=1-pchisq(diff$chisq, df=length-1)
    pvalueDF = sprintf("%.03f",pvalue)
    if(pvalue<0.001){
      pValue="p<0.001"     #pValue 大写
    }else{
      pValue=paste0("p=",sprintf("%.03f",pvalue))
    }
    fit <- survfit(Surv(futime, fustat) ~ geneCluster, data = rt)
    # print(surv_median(fit))
    # print(pValue)
    
    #绘制生存曲线
    bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
    bioCol=bioCol[1:length]
    print(paste(i,j,sep = "_"))
    surPlot=ggsurvplot(fit, 
                       data=rt,
                       conf.int=F,
                       pval=pValue,
                       pval.size=6,
                       legend.title="geneCluster",
                       legend.labs=levels(factor(rt[,"geneCluster"])),
                       legend = c(0.8, 0.8),
                       font.legend=10,
                       xlab="Time(Months)",
                       # break.time.by = 1,       #?????
                       palette = bioCol,
                       # surv.median.line = "hv",
                       risk.table=T,
                       cumevents=F,
                       risk.table.height=.25)
    
    pdf(file=paste0("survival_",i,"_",j,"_",pvalueDF,".pdf"),onefile = FALSE,width=7,height=5.5)
    print(surPlot)
    dev.off()
    pValueOutput=rbind(pValueOutput,
                       cbind(clusterpart=i,
                             cutoffData=j,
                             pvalue = pvalueDF))
    
  }
}
#输出单因素的结果
# outTab=cbind(outTab, km)
write.table(pValueOutput,file="m6ATME_geneCluster_kmPvalue.txt",sep="\t",row.names=F,col.names = T,quote=F)


# 01.01-------------------------------------------


# 02.00-------------------------------------------

setwd("C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\34.geneClusterSur")

rm(list = ls())
gc()

setwd("C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\34.geneClusterSur")


#引用包
library(limma)
library(ConsensusClusterPlus)
expFile="uniSigGeneExp.txt"        #表达输入文件
# expFile="uniSigGeneExp_p0.01.txt"        #表达输入文件

dir.create("C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\34.geneClusterSur")

workDir="C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\34.geneClusterSur"      #设置工作目录
#读取输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

dim(data)
# [1]  68 279


#聚类
maxK=9
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="png")



cliFile="time.txt"               #生存数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/30

summary(cli$futime)







#需要修改01#########################

setwd("C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\34.geneClusterSur")

cutData = seq(0,12,1)
pValueOutput = data.frame()

for (j in cutData) {
  cli = cli %>% dplyr::filter(futime >= j)                            
  #需要修改02#########################
  for (i in 2:3){
    # i =3
    clusterNum=i       #分几类，根据判断标准判断
    cluster=results[[clusterNum]][["consensusClass"]]
    cluster=as.data.frame(cluster)
    colnames(cluster)=c("geneCluster")
    letter=c("A","B","C","D","E","F","G")
    uniqClu=levels(factor(cluster$geneCluster))
    
    cluster$geneCluster=letter[match(cluster$geneCluster, uniqClu)]
    
    clusterOut=rbind(ID=colnames(cluster), cluster)
    write.table(clusterOut, file=paste0("geneCluster_", i,"_",j,".txt"), sep="\t", quote=F, col.names=F)
    
    # clusterFile=paste0("geneCluster_", i,".txt")    #m6A分型文件
    
    # setwd("D:\\biowolf\\m6aTME\\24.geneClusterSur")      #设置工作目录
    
    #读取输入文件
    # clusterOut=clusterOut 
    # cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
    rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
    
    
    
    #数据合并
    sameSample=intersect(row.names(cluster), row.names(cli))
    rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])
    
    #生存差异统计
    length=length(levels(factor(rt$geneCluster)))
    diff=survdiff(Surv(futime, fustat) ~ geneCluster, data = rt)
    pvalue=1-pchisq(diff$chisq, df=length-1)
    pvalueDF = sprintf("%.03f",pvalue)
    if(pvalue<0.001){
      pValue="p<0.001"     #pValue 大写
    }else{
      pValue=paste0("p=",sprintf("%.03f",pvalue))
    }
    fit <- survfit(Surv(futime, fustat) ~ geneCluster, data = rt)
    # print(surv_median(fit))
    # print(pValue)
    
    #绘制生存曲线
    bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
    bioCol=bioCol[1:length]
    print(paste(i,j,sep = "_"))
    surPlot=ggsurvplot(fit, 
                       data=rt,
                       conf.int=F,
                       pval=pValue,
                       pval.size=6,
                       legend.title="geneCluster",
                       legend.labs=levels(factor(rt[,"geneCluster"])),
                       legend = c(0.8, 0.8),
                       font.legend=10,
                       xlab="Time(Months)",
                       # break.time.by = 1,       #?????
                       palette = bioCol,
                       # surv.median.line = "hv",
                       risk.table=T,
                       cumevents=F,
                       risk.table.height=.25)
    
    pdf(file=paste0("survival_",i,"_",j,"_",pvalueDF,".pdf"),onefile = FALSE,width=7,height=5.5)
    print(surPlot)
    dev.off()
    pValueOutput=rbind(pValueOutput,
                       cbind(clusterpart=i,
                             cutoffData=j,
                             pvalue = pvalueDF))
    
  }
}
#输出单因素的结果
# outTab=cbind(outTab, km)
write.table(pValueOutput,file="m6ATME_geneCluster_kmPvalue.txt",sep="\t",row.names=F,col.names = T,quote=F)


# 02.01 -------------------------------------------------------------------




