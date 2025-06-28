
# 00 ----------------------------------------------------------------------
#多组学数据分析
#用途：

#方法：

# 说明
###1.原step37用的是PCA的积分，主成分1和主成分2
###2.对m6Ascore增加for循环 95line
# 原step38 中可以选中位cutoff，也可以选择最优的cutoff


#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

# BiocManager::install(version = '3.14')

#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")

#install.packages("corrplot")
#install.packages("ggpubr")
#BiocManager::install("maftools")

# BiocManager::install("xfun")


rm(list = ls())   #清空内存  rm=remove
options(stringsAsFactors = F)


#引用包
library(plyr)
library(ggplot2)
library(ggpubr)
library(pacman)
library(survival)
library(survminer)
library(ggalluvial)
library(ggplot2)
library(dplyr)
library(corrplot)
library(maftools)

library(Rcpp)

library(Rcpp)
library(stringr)
library(ggsci)
library(ggplot2)
library(Rtsne)
library(umap)
library(cowplot)
library(Rcpp)
if(!requireNamespace("tidyr",quietly = TRUE)) install.packages("tidyr",update = F,ask = F)
if(!requireNamespace("tidyverse",quietly = TRUE)) install.packages("tidyverse",update = F,ask = F)

if(!requireNamespace("magrittr",quietly = TRUE)) install.packages("magrittr",update = F,ask = F)
if(!requireNamespace("data.table",quietly = TRUE)) install.packages("data.table",update = F,ask = F)

# 安装加载需要的R包
# install.packages("pacman", repos = 'https://mirror.lzu.edu.cn/CRAN/')
library(pacman)
#按F1和F2可以查看p_load这个函数的属性和原代码
# p_load  是函数
# p_load(TCGAbiolinks, tidyverse, magrittr, data.table, biomaRt)
#ctrl+L 清屏

projectPath = "C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\01.data"    #设置主工作目录


# dataPath = paste(projectPath, "01.data", sep = "/")
# if(!dir.exists(dataPath)) dir.create(dataPath)  # 新建一个Data文件夹，用于存放从TCGA上下载的数据
# setwd(dataPath)

# resultPath = paste(projectPath, "01.result", sep = "/")
# if(!dir.exists(resultPath)) dir.create(resultPath)  # 新建一个Data文件夹，用于存放从TCGA上下载的数据


dataPath = paste(projectPath, "02.data", sep = "/")
if(!dir.exists(dataPath)) dir.create(dataPath)  # 新建一个Data文件夹，用于存放从TCGA上下载的数据
setwd(dataPath)

# resultPath = paste(projectPath, "02.result", sep = "/")
# if(!dir.exists(resultPath)) dir.create(resultPath)  # 新建一个Data文件夹，用于存放从TCGA上下载的数据


# resultPath = paste(projectPath, "03.result", sep = "/")
# if(!dir.exists(resultPath)) dir.create(resultPath)  # 新建一个Data文件夹，用于存放从TCGA上下载的数据

resultPath = paste(projectPath, "03.result_umap", sep = "/")
if(!dir.exists(resultPath)) dir.create(resultPath)  # 新建一个Data文件夹，用于存放从TCGA上下载的数据




# 01 ----------------------------------------------------------------------



# 原step37 -----------------------------------------------------------------

expFile="uniSigGeneExp.txt"      #表达输入文件
# 基因样本文件
# setwd("C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\01.data")     #设置工作目录

#读取输入文件
uniSigGeneExp = read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
uniSigGeneExp_t = t(uniSigGeneExp)
# 样本基因数据框

# data = uniSigGeneExp_t
# dim(data)
# [1] 122  82

# 02 ----------------------------------------------------------------------

#PCA分析-----------------------
#PCAnoScale分析-----------------------
#tsne分析----------------
#umap分析

# methods = c("PCA","PCAnoScale","tsne","umap")
# methods = c("PCAnoScale","tsne","umap")
# methods = c("tsne","umap")
methods = c("umap")

for(method in methods){
  
  # method = c("PCA")
  # method = c("tsne")
  method = c("umap")
  
  if (method == c("PCA")){
    pca=prcomp(uniSigGeneExp_t, scale=TRUE)
    value=predict(pca)
    m6Ascore=value[,1]+value[,2]
    m6Ascore=as.data.frame(m6Ascore)
    scoreOut=rbind(id=colnames(m6Ascore), m6Ascore)
    write.table(scoreOut, file=paste0("m6Ascore_",method,".txt"), sep="\t", quote=F, col.names=F)
  }
  
  if (method == c("PCAnoScale")){
    pca=prcomp(uniSigGeneExp_t)
    value=predict(pca)
    m6Ascore=value[,1]+value[,2]
    m6Ascore=as.data.frame(m6Ascore)
    scoreOut=rbind(id=colnames(m6Ascore), m6Ascore)
    write.table(scoreOut, file=paste0("m6Ascore_",method,".txt"), sep="\t", quote=F, col.names=F)
  }

  if (method == c("tsne")){
    set.seed(1)
    # data=t(data)
    # tsne_out <- Rtsne(t(data),theta = 0.1,perplexity = 3) # Run TSNE
    tsne_out <- Rtsne(uniSigGeneExp_t,theta = 0.1,perplexity = 3) # Run TSNE
    #scatter plot
    plot.data_tsne <- as.data.frame(tsne_out$Y)
    colnames(plot.data_tsne) = c("tSNE1","tSNE2")
    m6Ascore=plot.data_tsne[,1]+plot.data_tsne[,2]
    m6Ascore=as.data.frame(m6Ascore)
    rownames(m6Ascore) = rownames(uniSigGeneExp_t)
  }  
  
  if (method == c("umap")){
    umap.out = umap::umap(uniSigGeneExp_t)
    #scatter plot
    plot.data_umap <- as.data.frame(umap.out$layout)
    colnames(plot.data_umap) = c("UMAP1","UMAP2") 
    m6Ascore=plot.data_umap[,1]+plot.data_umap[,2]
    m6Ascore=as.data.frame(m6Ascore)
    rownames(m6Ascore) = rownames(uniSigGeneExp_t)
  }  
  
  # pca=prcomp(uniSigGeneExp_t, scale=TRUE)
  # value=predict(pca)
  # m6Ascore=value[,1]+value[,2]
  # m6Ascore=as.data.frame(m6Ascore)
  # scoreOut=rbind(id=colnames(m6Ascore), m6Ascore)
  # write.table(scoreOut, file="m6Ascore_PCA.txt", sep="\t", quote=F, col.names=F)
  
  # 输出
  # outfile = paste0("LIHC_Portal_RNA_", str_match(data_type, "- ([^\\s]*)$")[,2], ".txt")
  # fwrite(scoreOut, outfile, row.names = F, sep = "\t", quote = F)
  
  
  # 原step38 -----------------------------------------------------------------
  
  
  #install.packages("survival")
  #install.packages("survminer")
  
  
  #引用包
  library(survival)
  library(survminer)
  # scoreFile="m6Ascore.txt"     #m6A打分文件
  cliFile="time.txt"           #生存数据文件
  # setwd("D:\\biowolf\\m6aTME\\38.scoreSur")      #设置工作目录
  
  
  #读取输入文件
  
  
  
  # score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
  
  # m6Ascore
  
  
  sampleType=gsub("(.*?)\\_.*", "\\1", row.names(m6Ascore))
  m6AscoreType=cbind(m6Ascore, sampleType)
  rownames(m6AscoreType)=gsub("(.*?)\\_(.*?)", "\\2", rownames(m6AscoreType))
  
  
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
  colnames(cli)=c("futime", "fustat")
  cli$futime=cli$futime/365
  
  #数据合并
  sameSample=intersect(row.names(m6AscoreType), row.names(cli))
  DataTimeSurM6aScore=cbind(cli[sameSample,], m6AscoreType[sameSample,c("m6Ascore")])
  names(DataTimeSurM6aScore)[3] = c("m6Ascore")
  
  # 第1种可能中位cutoff -------------------------------------------------------------------  
  
  
  
  Type=ifelse(DataTimeSurM6aScore[,"m6Ascore"]<=median(DataTimeSurM6aScore[,"m6Ascore"]), "Low", "High")
  DataTimeSurM6aScore$group=Type
  outTab=rbind(id=colnames(DataTimeSurM6aScore), DataTimeSurM6aScore)
  write.table(outTab, file=paste0("m6Ascore.group_",method,"_median.txt"), sep="\t", quote=F, col.names=F)
  
  #计算高低风险组生存差异
  DataTimeSurM6aScore$group=factor(DataTimeSurM6aScore$group, levels=c("Low", "High"))
  diff=survdiff(Surv(futime, fustat) ~ group, data = DataTimeSurM6aScore)
  length=length(levels(factor(DataTimeSurM6aScore[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ group, data = DataTimeSurM6aScore)
  #print(surv_median(fit))
  
  #绘制生存曲线
  bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit, 
                     data=DataTimeSurM6aScore,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="m6Ascore",
                     legend.labs=levels(factor(DataTimeSurM6aScore[,"group"])),
                     legend = c(0.8, 0.8),
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette = bioCol,
                     surv.median.line = "hv",
                     risk.table=T,
                     cumevents=F,
                     risk.table.height=.25)
  
  
  
  
  #保存图片
  pdf(file=paste0(resultPath,"/","m6AscoreSurvival_",method,"_中位cutoff.pdf"), onefile = FALSE, width=7, height=5.5)
  print(surPlot)
  dev.off()
  
  
  
  # # 第2种可能最优cutoff -------------------------------------------------------------------
  #获取最优cutoff
  res.cut=surv_cutpoint(DataTimeSurM6aScore, time="futime", event="fustat", variables=c("m6Ascore"))
  cutoff=as.numeric(res.cut$cutpoint[1])
  print(cutoff)

  Type=ifelse(DataTimeSurM6aScore[,"m6Ascore"]<=cutoff, "Low", "High")
  DataTimeSurM6aScore$group=Type
  outTab=rbind(id=colnames(DataTimeSurM6aScore), DataTimeSurM6aScore)
  write.table(outTab, file=paste0("m6Ascore.group_",method,".txt"), sep="\t", quote=F, col.names=F)

  #计算高低风险组生存差异
  DataTimeSurM6aScore$group=factor(DataTimeSurM6aScore$group, levels=c("Low", "High"))
  diff=survdiff(Surv(futime, fustat) ~ group, data = DataTimeSurM6aScore)
  length=length(levels(factor(DataTimeSurM6aScore[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ group, data = DataTimeSurM6aScore)
  #print(surv_median(fit))

  #绘制生存曲线
  bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit,
                     data=DataTimeSurM6aScore,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="m6Ascore",
                     legend.labs=levels(factor(DataTimeSurM6aScore[,"group"])),
                     legend = c(0.8, 0.8),
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette = bioCol,
                     surv.median.line = "hv",
                     risk.table=T,
                     cumevents=F,
                     risk.table.height=.25)




  #保存图片
  pdf(file=paste0(resultPath,"/","m6AscoreSurvival_",method,"_最优cutoff.pdf"), onefile = FALSE, width=7, height=5.5)
  print(surPlot)
  dev.off()
  

  
  
  
  
  # 原step39.00-----------------------------------------------------------------
  
  #引用包
  # library(ggalluvial)
  # library(ggplot2)
  # library(dplyr)
  
  
  m6aCluFile="m6aCluster.txt"         #m6A分型文件
  geneCluFile="geneCluster_2_0.txt"       #基因分型文件
  scoreFile="m6Ascore.group_PCA.txt"      #m6A打分的分组文件
  cliFile="clinical.txt"              #临床数据文件
  
  
  trait="Fustat"                      #临床性状
  # setwd("D:\\biowolf\\m6aTME\\39.ggalluvial")     #设置工作目录
  
  #读取输入文件
  m6aClu=read.table(m6aCluFile, header=T, sep="\t", check.names=F, row.names=1)
  
  
  geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
  
  # score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
  # DataTimeSurM6aScore
  
  #合并数据
  twoCluster=cbind(m6aClu, geneClu)
  rownames(twoCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(twoCluster))
  sameSample=intersect(row.names(twoCluster), row.names(DataTimeSurM6aScore))
  scoreClu=cbind(DataTimeSurM6aScore[sameSample,,drop=F], twoCluster[sameSample,,drop=F])
  sameSample=intersect(row.names(scoreClu), row.names(cli))
  ggalluvial_rt=cbind(scoreClu[sameSample,], cli[sameSample,])
  
  #准备桑基图输入文件
  ggalluvial_rt=ggalluvial_rt[,c("m6Acluster", "geneCluster", "group", trait)]
  colnames(ggalluvial_rt)=c("m6Acluster", "geneCluster", "m6Ascore", trait)
  corLodes=to_lodes_form(ggalluvial_rt, axes = 1:ncol(ggalluvial_rt), id = "Cohort")
  # method = "a"
  # method <- NULL
  #得到输出文件
  pdf(file=paste0(resultPath, "/","ggalluvial_",method,"_3cluster.pdf"), width=6, height=6)
  mycol=rep(c("#0066FF","#FF9900","#FF0000","#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
  ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
    scale_x_discrete(expand = c(0, 0)) +  
    #用aes.flow控制线条颜色，forward说明颜色和前面的柱状图一致，backward说明和后面的柱状图一致。
    geom_flow(width = 2/10,aes.flow = "forward") + 
    geom_stratum(alpha = .9,width = 2/10) +
    scale_fill_manual(values = mycol) +
    #size=3代表字体大小
    geom_text(stat = "stratum", size = 3,color="black") +
    xlab("") + ylab("") + theme_bw() + 
    theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #去掉坐标轴
    theme(panel.grid =element_blank()) + 
    theme(panel.border = element_blank()) + 
    ggtitle("") + guides(fill = FALSE)                            
  dev.off()
  
  
  # 原step40.00-----------------------------------------------------------------
  # 原step40.00-----------------------------------------------------------------
  # 原step40.00-----------------------------------------------------------------
  
  #install.packages("corrplot")
  
  
  #引用包
  library(corrplot)
  # scoreFile="m6Ascore.txt"         #m6A打分文件
  #读取m6A打分文件
  # score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
  
  immFile="ssGSEA.result.txt"      #ssGSEA结果文件                       # 可替换XXXX
  # setwd("D:\\biowolf\\m6aTME\\40.scoreImm")     #设置工作目录
  
  #读取免疫细胞文件
  immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
  immune=t(immune)
  
  #数据合并
  sameSample=intersect(row.names(m6Ascore), row.names(immune))
  data_cor=cbind(m6Ascore[sameSample,,drop=F], immune[sameSample,,drop=F])
  
  #相关性矩阵
  M=cor(data_cor)
  res1_cor=cor.mtest(data_cor, conf.level = 0.95)
  
  #绘制相关性图形
  pdf(file=paste0(resultPath,"/","cor_",method,".pdf"), width=6, height=6)
  # pdf(file="cor_PCA.pdf", width=8, height=8)
  corrplot(M,
           order="original",
           method = "circle",
           type = "upper",
           tl.cex=0.8, pch=T,
           p.mat = res1_cor$p,
           insig = "label_sig",
           pch.cex = 1.6,
           sig.level=0.05,
           number.cex = 1,
           col=colorRampPalette(c("blue", "white", "red"))(50),
           tl.col="black")
  dev.off()
  
  
  # 原step40.01-----------------------------------------------------------------
  # 原step40.01-----------------------------------------------------------------
  # 原step40.01-----------------------------------------------------------------
  
  
  
  # 原step41.00-----------------------------------------------------------------
  # 原step41.00-----------------------------------------------------------------
  # 原step41.00-----------------------------------------------------------------
  
  
  
  #引用包
  library(limma)
  library(ggpubr)
  
  
  # m6aCluFile="m6aCluster.txt"        #m6A分型文件
  # geneCluFile="geneCluster.txt"      #基因分型文件
  geneCluFile="geneCluster_2_0.txt"      #基因分型文件
  
  # scoreFile="m6Ascore.txt"           #m6A打分文件
  
  
  
  # setwd("D:\\biowolf\\m6aTME\\41.clusterScore")     #设置工作目录
  
  #读取输入文件
  # m6aClu=read.table(m6aCluFile, header=T, sep="\t", check.names=F, row.names=1)
  geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)
  
  # score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
  # m6Ascore
  
  
  
  #合并数据
  twoCluster=cbind(m6aClu, geneClu)
  sameSample=intersect(row.names(twoCluster), row.names(m6Ascore))
  
  
  data_m6Ascorem6AclugeneClu=cbind(m6Ascore[sameSample,,drop=F], twoCluster[sameSample,,drop=F])
  
  #######m6A分型与打分相关性########
  #设置比较组
  data_m6Ascorem6AclugeneClu$m6Acluster=factor(data_m6Ascorem6AclugeneClu$m6Acluster, levels=levels(factor(data_m6Ascorem6AclugeneClu$m6Acluster)))
  group=levels(factor(data_m6Ascorem6AclugeneClu$m6Acluster))
  comp=combn(group, 2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  #定义颜色
  bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  bioCol=bioCol[1:length(levels(factor(data_m6Ascorem6AclugeneClu$m6Acluster)))]
  
  #绘制boxplot
  boxplot=ggboxplot(data_m6Ascorem6AclugeneClu, x="m6Acluster", y="m6Ascore", color="m6Acluster",
                    xlab="m6Acluster",
                    ylab="m6Ascore",
                    legend.title="m6Acluster",
                    palette=bioCol,
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  
  #输出图片
  pdf(file=paste0(resultPath,"/","m6Acluster_",method,".pdf"), width=5, height=4.5)
  # pdf(file="m6Acluster.pdf", width=5, height=4.5)
  print(boxplot)
  dev.off()
  #######m6A分型与打分相关性########
  
  
  #######基因分型与打分相关性########
  #设置比较组
  data_m6Ascorem6AclugeneClu=cbind(m6Ascore[sameSample,,drop=F], twoCluster[sameSample,,drop=F])
  # data_m6Ascorem6AclugeneClu$geneCluster=factor(data_m6Ascorem6AclugeneClu$geneCluster, levels=levels(data_m6Ascorem6AclugeneClu$geneCluster))
  group=levels(factor(data_m6Ascorem6AclugeneClu$geneCluster))
  comp=combn(group, 2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  #定义颜色
  bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  bioCol=bioCol[1:length(levels(factor(data_m6Ascorem6AclugeneClu$geneCluster)))]
  
  #绘制boxplot
  boxplot=ggboxplot(data_m6Ascorem6AclugeneClu, x="geneCluster", y="m6Ascore", color="geneCluster",
                    xlab="geneCluster",
                    ylab="m6Ascore",
                    legend.title="geneCluster",
                    palette=bioCol,
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  
  #输出图片
  # pdf(file=paste(resultPath,"cor_PCA.pdf",sep = "/"), width=6, height=6)
  pdf(file=paste0(resultPath,"/","geneCluster_",method,".pdf"), width=5, height=4.5)
  # pdf(file="geneCluster.pdf", width=5, height=4.5)
  print(boxplot)
  dev.off()
  #######基因分型与打分相关性########
  
  
  
  # 原step41.01-----------------------------------------------------------------
  # 原step41.01-----------------------------------------------------------------
  # 原step41.01-----------------------------------------------------------------
  
  
  
  # 原step42.00-----------------------------------------------------------------
  # 原step42.00-----------------------------------------------------------------
  # 原step42.00-----------------------------------------------------------------
  
  
  #引用包
  library(ggpubr)
  library(reshape2)
  
  tmbFile="TMB.txt"     #肿瘤突变负荷文件
  # scoreFile="m6Ascore.group.txt"     #m6A打分的分组文件
  # DataTimeSurM6aScore
  
  
  # cluFile="geneCluster.txt"          #基因分型文件
  # geneClu
  
  # setwd("D:\\biowolf\\m6aTME\\42.scoreTMB")       #修改工作目录
  
  #读取输入文件
  tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)        #读取TMB数据文件
  
  # score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)    #读取m6A打分的分组文件
  # DataTimeSurM6aScore
  
  
  # clu=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)        #读取基因分型类文件
  
  #合并数据
  tmb=as.matrix(tmb)
  tmb[tmb>quantile(tmb,0.975)]=quantile(tmb,0.975)
  
  
  sameSample=intersect(row.names(tmb), row.names(DataTimeSurM6aScore))
  
  tmb=tmb[sameSample,,drop=F]
  
  DataTimeSurM6aScoreSS=DataTimeSurM6aScore[sameSample,,drop=F]
  rownames(geneClu)=gsub("(.*?)\\_(.*?)", "\\2", rownames(geneClu))
  geneCluSS=geneClu[sameSample,,drop=F]
  
  
  data_tmbgeneClu=cbind(DataTimeSurM6aScoreSS, tmb, geneCluSS)
  data_tmbgeneClu=data_tmbgeneClu[,c("m6Ascore", "group", "geneCluster", "TMB")]
  
  #设置比较组
  data_tmbgeneClu$group=factor(data_tmbgeneClu$group, levels=c("Low", "High"))
  group=levels(factor(data_tmbgeneClu$group))
  comp=combn(group, 2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  #设置颜色
  bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  bioCol=bioCol[1:length(group)]
  
  #绘制箱线图
  boxplot=ggboxplot(data_tmbgeneClu, x="group", y="TMB", fill="group",
                    xlab="",
                    ylab="Tumor Burden Mutation",
                    legend.title="m6AScore",
                    palette = bioCol )+ 
    stat_compare_means(comparisons = my_comparisons)
  
  # pdf(file=paste(resultPath,"cor_PCA.pdf",sep = "/"), width=6, height=6)
  pdf(file=paste0(resultPath,"/","tmbgeneClu_boxplot_",method,".pdf"),width=5,height=4.5)
  print(boxplot)
  dev.off()
  
  #相关性图形
  length=length(levels(factor(data_tmbgeneClu$geneCluster)))
  bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  p1=ggplot(data_tmbgeneClu, aes(m6Ascore, TMB)) + 
    xlab("m6Ascore")+ylab("Tumor Burden Mutation")+
    geom_point(aes(colour=geneCluster))+
    scale_color_manual(values=bioCol[1:length])+ 
    geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = 'spearman', aes(x =m6Ascore, y =TMB))
  
  
  #相关性图形
  # pdf(file="tmbgeneClu_cor_PCA.pdf", width=6, height=4.5)
  pdf(file=paste0(resultPath,"/","tmbgeneClu_cor_",method,".pdf"), width=6, height=4.5)
  print(p1)
  dev.off()
  
  
  
  
  # 原step42.01-----------------------------------------------------------------
  # 原step42.01-----------------------------------------------------------------
  # 原step42.01-----------------------------------------------------------------
  
  
  
  
  
  # 原step43.00-----------------------------------------------------------------
  # 原step43.00-----------------------------------------------------------------
  # 原step43.00-----------------------------------------------------------------
  
  #引用包
  library(survival)
  library(survminer)
  
  
  # tmbFile="TMB.txt"                  #肿瘤突变负荷文件
  # scoreFile="m6Ascore.group.txt"     #m6A打分的分组文件
  
  
  
  # setwd("D:\\biowolf\\m6aTME\\43.tmbSur")       #修改工作目录
  
  #读取输入文件
  # score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)    #读取m6A打分的分组文件
  # DataTimeSurM6aScore
  
  
  
  # tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)        #读取TMB数据文件
  # tmb
  
  
  
  #合并数据
  sameSample=intersect(row.names(tmb), row.names(DataTimeSurM6aScore))
  tmbSS=tmb[sameSample,,drop=F]
  # DataTimeSurM6aScoreSS=DataTimeSurM6aScore[sameSample,,drop=F]
  data_TimeSurM6aScoreTMB=cbind(DataTimeSurM6aScoreSS, tmbSS)
  
  # 第1种可能----------------------
  #43.10获取最优cutoff--------------------------------
  res.cut=surv_cutpoint(data_TimeSurM6aScoreTMB, time = "futime", event = "fustat", variables =c("TMB"))
  cutoff=as.numeric(res.cut$cutpoint[1])
  
  tmbType=ifelse(data_TimeSurM6aScoreTMB[,"TMB"]<=cutoff, "L-TMB", "H-TMB")
  scoreType=ifelse(data_TimeSurM6aScoreTMB$group=="Low", "L-m6Ascore", "H-m6Ascore")
  mergeType=paste0(tmbType, "+", scoreType)
  
  
  #生存曲线函数
  bioSurvival=function(surData=null, outFile=null){
    diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
    length=length(levels(factor(surData[,"group"])))
    pValue=1-pchisq(diff$chisq, df=length-1)
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("p=",sprintf("%.03f",pValue))
    }
    fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
    #print(surv_median(fit))
    
    #绘制生存曲线
    width=6.5
    height=5.5
    if(length(levels(factor(surData[,"group"])))>2){
      width=8
      height=6.5
    }
    bioCol=c("#FF0000","#0066FF","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
    bioCol=bioCol[1:length]
    surPlot=ggsurvplot(fit, 
                       data=surData,
                       conf.int=F,
                       pval=pValue,
                       pval.size=6,
                       legend.title="",
                       legend.labs=levels(factor(surData[,"group"])),
                       font.legend=10,
                       legend = c(0.8, 0.8),
                       xlab="Time(years)",
                       break.time.by = 1,
                       palette = bioCol,
                       surv.median.line = "hv",
                       risk.table=T,
                       cumevents=F,
                       risk.table.height=.25)
    #输出图形
    pdf(file=outFile, onefile = FALSE, width=width, height=height)
    print(surPlot)
    dev.off()
  }
  
  #绘制TMB的生存曲线
  data_TimeSurM6aScoreTMB$group=tmbType
  # pdf(file=paste(resultPath,"tmbgeneClu_cor_PCA.pdf",sep = "/"), width=6, height=4.5)
  
  bioSurvival(surData=data_TimeSurM6aScoreTMB, outFile=paste0(resultPath,"/","TMB.survival_",method,"_最佳cutoff.pdf"))
  
  #绘制TMB联合m6A打分的生存曲线
  data_TimeSurM6aScoreTMB$group=mergeType
  # bioSurvival(surData=data, outFile="TMB-score.survival.pdf")
  bioSurvival(surData=data_TimeSurM6aScoreTMB, outFile=paste0(resultPath,"/","TMB-score.survival_",method,"_最佳cutoff.pdf"))
  
  
  # #第2种可能--------------
  # #43.20获取中位cutoff无意义--------------------------------
  # # res.cut=surv_cutpoint(data_TimeSurM6aScoreTMB, time = "futime", event = "fustat", variables =c("TMB"))
  # # cutoff=as.numeric(res.cut$cutpoint[1])
  # 
  # tmbType=ifelse(data_TimeSurM6aScoreTMB[,"TMB"]<=median(data_TimeSurM6aScoreTMB[,"TMB"]), "L-TMB", "H-TMB")
  # scoreType=ifelse(data_TimeSurM6aScoreTMB$group=="Low", "L-m6Ascore", "H-m6Ascore")
  # mergeType=paste0(tmbType, "+", scoreType)
  # 
  # 
  # #生存曲线函数
  # bioSurvival=function(surData=null, outFile=null){
  #   diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
  #   length=length(levels(factor(surData[,"group"])))
  #   pValue=1-pchisq(diff$chisq, df=length-1)
  #   if(pValue<0.001){
  #     pValue="p<0.001"
  #   }else{
  #     pValue=paste0("p=",sprintf("%.03f",pValue))
  #   }
  #   fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
  #   #print(surv_median(fit))
  #   
  #   #绘制生存曲线
  #   width=6.5
  #   height=5.5
  #   if(length(levels(factor(surData[,"group"])))>2){
  #     width=8
  #     height=6.5
  #   }
  #   bioCol=c("#FF0000","#0066FF","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  #   bioCol=bioCol[1:length]
  #   surPlot=ggsurvplot(fit, 
  #                      data=surData,
  #                      conf.int=F,
  #                      pval=pValue,
  #                      pval.size=6,
  #                      legend.title="",
  #                      legend.labs=levels(factor(surData[,"group"])),
  #                      font.legend=10,
  #                      legend = c(0.8, 0.8),
  #                      xlab="Time(years)",
  #                      break.time.by = 1,
  #                      palette = bioCol,
  #                      surv.median.line = "hv",
  #                      risk.table=T,
  #                      cumevents=F,
  #                      risk.table.height=.25)
  #   #输出图形
  #   pdf(file=outFile, onefile = FALSE, width=width, height=height)
  #   print(surPlot)
  #   dev.off()
  # }
  # 
  # #绘制TMB的生存曲线
  # data_TimeSurM6aScoreTMB$group=tmbType
  # # pdf(file=paste(resultPath,"tmbgeneClu_cor_PCA.pdf",sep = "/"), width=6, height=4.5)
  # 
  # bioSurvival(surData=data_TimeSurM6aScoreTMB, outFile=paste0(resultPath,"/","TMB.survival_",method,"_中位cutoff.pdf"))
  # 
  # #绘制TMB联合m6A打分的生存曲线
  # data_TimeSurM6aScoreTMB$group=mergeType
  # # bioSurvival(surData=data, outFile="TMB-score.survival.pdf")
  # bioSurvival(surData=data_TimeSurM6aScoreTMB, outFile=paste0(resultPath,"/","TMB-score.survival_",method,"_中位cutoff.pdf"))
  # #43.21获取中位cutoff无意义--------------------------------
  
  
  
  # 原step43.01-----------------------------------------------------------------
  # 原step43.01-----------------------------------------------------------------
  # 原step43.01-----------------------------------------------------------------
  
  
  
  
  
  # 原step45.00-----------------------------------------------------------------
  # 原step45.00-----------------------------------------------------------------
  # 原step45.00-----------------------------------------------------------------
  
  
  
  library(maftools)       #引用包
  # setwd("D:\\biowolf\\m6aTME\\45.maftools")      #设置工作目录
  
  #读取m6A评分的分组文件
  # score=read.table("m6Ascore.group.txt", header=T, sep="\t", check.names=F)
  # DataTimeSurM6aScore
  table(DataTimeSurM6aScore$group)
  
  
  outTab_SampleBarcodeM6aScore=DataTimeSurM6aScore %>% 
    tibble::rownames_to_column(var = "Tumor_Sample_Barcode") %>% 
    dplyr::select(c("Tumor_Sample_Barcode", "m6Ascore"))
  
  # colnames(outTab)=c("Tumor_Sample_Barcode", "m6Ascore")
  
  write.table(outTab_SampleBarcodeM6aScore, file="ann.txt", sep="\t", quote=F, row.names=F)
  
  #读取基因突变文件
  geneNum=20
  geneMut=read.table("geneMut.txt", header=T, sep="\t", check.names=F, row.names=1)
  gene=row.names(geneMut)[1:geneNum]
  
  #颜色
  ann_colors=list()
  col=c("#0066FF","#FF0000")
  names(col)=c("Low", "High")
  ann_colors[["m6Ascore"]]=col
  
  #低评分组瀑布图
  pdf(file=paste0(resultPath,"/","low_低评分组瀑布图_",method,".pdf"), width=6, height=6)
  maf=read.maf(maf="low.maf", clinicalData="ann.txt")
  oncoplot(maf=maf, clinicalFeatures="m6Ascore", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
  dev.off()
  
  #高评分组瀑布图
  # pdf(file="high_高评分组瀑布图_PCA.pdf", width=6, height=6)
  pdf(file=paste0(resultPath,"/","high_高评分组瀑布图_",method,".pdf"), width=6, height=6)
  maf=read.maf(maf="high.maf", clinicalData="ann.txt")
  oncoplot(maf=maf, clinicalFeatures="m6Ascore", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
  dev.off()
  
  # 原step45.01-----------------------------------------------------------------
  # 原step45.01-----------------------------------------------------------------
  # 原step45.01-----------------------------------------------------------------
  
  
  
  
  # 原step46.00-----------------------------------------------------------------
  # 原step46.00-----------------------------------------------------------------
  # 原step46.00-----------------------------------------------------------------
  #install.packages("ggplot2")
  #install.packages("ggpubr")
  
  
  #引用包
  library(plyr)
  library(ggplot2)
  library(ggpubr)
  
  
  # scoreFile="m6Ascore.group.txt"    #m6A打分文件
  
  # cliFile="clinical.txt"            #临床数据文件
  # trait="Fustat"                    #临床性状
  # setwd("D:\\biowolf\\m6aTME\\46.scoreCli")     #设置工作目录
  
  #读取输入文件
  # score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
  # DataTimeSurM6aScore
  
  
  # cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
  # sameSample=intersect(row.names(score), row.names(cli))
  sameSample=intersect(row.names(DataTimeSurM6aScore), row.names(cli))
  
  
  data_TimeSurM6aScoreCli=cbind(DataTimeSurM6aScore[sameSample,,drop=F], cli[sameSample,,drop=F])
  
  # DataTimeSurM6aScore
  
  
  #定义临床性状的颜色
  bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  bioCol=bioCol[1:length(unique(data_TimeSurM6aScoreCli[,trait]))]
  
  #统计高低评分组病人数目
  data_TimeSurM6aScoreCliGroup=data_TimeSurM6aScoreCli[,c(trait, "group")]
  colnames(data_TimeSurM6aScoreCliGroup)=c("trait", "group")
  data_TimeSurM6aScoreCliGroup_df=as.data.frame(table(data_TimeSurM6aScoreCliGroup))
  
  df = data_TimeSurM6aScoreCliGroup_df
  #计算高低评分组的百分率
  df=ddply(df, .(group), transform, percent = Freq/sum(Freq) * 100)
  #百分比位置
  df=ddply(df, .(group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
  df$label=paste0(sprintf("%.0f", df$percent), "%")
  df$group=factor(df$group, levels=c("Low", "High"))
  #绘制百分率图
  p=ggplot(df, aes(x = factor(group), y = percent, fill = trait)) +
    geom_bar(position = position_stack(), stat = "identity", width = .7) +
    scale_fill_manual(values=bioCol)+
    xlab("m6Ascore")+ ylab("Percent weight")+  guides(fill=guide_legend(title=trait))+
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
    #coord_flip()+
    theme_bw()
  
  
  pdf(file=paste0(resultPath,"/","TimeSurM6aScoreCliGroup_barplot_group_",method,".pdf"), width=4, height=5)
  # pdf(file="TimeSurM6aScoreCliGroup_barplot_PCA.pdf", width=4, height=5)
  print(p)
  dev.off()
  
  #设置比较组
  data_TimeSurM6aScoreCli2=data_TimeSurM6aScoreCli[,c(trait, "m6Ascore")]
  colnames(data_TimeSurM6aScoreCli2)=c("trait", "m6Ascore")
  type=levels(factor(data_TimeSurM6aScoreCli2[,"trait"]))
  comp=combn(type, 2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  #绘制箱线图
  boxplot=ggboxplot(data_TimeSurM6aScoreCli2, x="trait", y="m6Ascore", fill="trait",
                    xlab=trait,
                    ylab="m6Ascore",
                    legend.title=trait,
                    palette=bioCol  )+ 
    stat_compare_means(comparisons=my_comparisons)
  
  
  pdf(file=paste0(resultPath,"/","TimeSurM6aScoreCliGroup_barplot_status_",method,".pdf"), width=4, height=4.5)
  # pdf(file="boxplot.pdf",width=4,height=4.5)
  print(boxplot)
  dev.off()
  
  
  # 原step46.01-----------------------------------------------------------------
  # 原step46.01-----------------------------------------------------------------
  # 原step46.01-----------------------------------------------------------------
  
  
  # 原step47.00-----------------------------------------------------------------
  # 原step47.00-----------------------------------------------------------------
  # 原step47.00-----------------------------------------------------------------
  
  #install.packages("survival")
  #install.packages("survminer")
  
  
  #引用包
  # library(survival)
  # library(survminer)
  # scoreFile="m6Ascore.group.txt"     #m6A打分分组文件
  
  cliFile="clinical_merge_binary.txt"             #临床数据文件
  #临床性状
  # setwd("D:\\biowolf\\m6aTME\\47.cliGroupSur")         #设置工作目录
  
  #读取输入文件
  # score=read.table(scoreFile, header=T, sep="\t",check.names=F, row.names=1)
  
  
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
  sameSample=intersect(row.names(cli), row.names(DataTimeSurM6aScore))
  DataTimeSurM6aScore=DataTimeSurM6aScore[sameSample,]
  # DataTimeSurM6aScore
  cli=cli[sameSample,]
  
  
  DataTimeSurM6aScoreclin2=cbind(futime=DataTimeSurM6aScore[,1], fustat=DataTimeSurM6aScore[,2],
                                 cli, group=DataTimeSurM6aScore[,"group"])
  
  data = DataTimeSurM6aScoreclin2
  
  
  survplot_cli_trait <- function(data, trait) {
    #提取临床数据
    # trait="Age" 
    rt=data[,c("futime", "fustat", trait, "group")]
    rt=rt[(rt[,trait]!="unknow"),]
    colnames(rt)=c("futime", "fustat", "clinical", "group")
    tab=table(rt[,"clinical"])
    tab=tab[tab!=0]
    
    #对每个临床信息里面的分类进行循环
    for(j in names(tab)){
      rt1=rt[(rt[,"clinical"]==j),]
      tab1=table(rt1[,"group"])
      tab1=tab1[tab1!=0]
      labels=names(tab1)
      if(length(labels)==2){
        titleName=j
        if((trait=="age") | (trait=="Age") | (trait=="AGE")){
          titleName=paste0("age", j)
        }
        diff=survdiff(Surv(futime, fustat) ~group,data = rt1)
        pValue=1-pchisq(diff$chisq, df=1)
        if(pValue<0.001){
          pValue="p<0.001"
        }else{
          pValue=paste0("p=",sprintf("%.03f", pValue))
        }
        fit <- survfit(Surv(futime, fustat) ~ group, data = rt1)
        bioCol=c("#FF0000","#0066FF","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
        bioCol=bioCol[1:length(levels(factor(rt1[,"group"])))]
        #绘制生存曲线
        surPlot=ggsurvplot(fit, 
                           data=rt1,
                           conf.int=F,
                           pval=pValue,
                           pval.size=6,
                           title=paste0("Patients with ",titleName),
                           legend.title="m6Ascore",
                           legend.labs=labels,
                           font.legend=12,
                           xlab="Time(years)",
                           break.time.by = 1,
                           palette=bioCol,
                           risk.table=TRUE,
                           risk.table.title="",
                           risk.table.col = "strata",
                           risk.table.height=.25)
        #输出图片
        j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
        pdf(file=paste0(resultPath,"/",trait,"_",j,"_",method,".pdf"), onefile = FALSE, width = 5.5, height =5)
        print(surPlot)
        dev.off()
      }
    }
  }
  
  # paste(resultPath,"TimeSurM6aScoreCliGroup_barplot_status_PCA.pdf",sep = "/")
  
  dput(colnames(DataTimeSurM6aScoreclin2))
  # c("futime", "fustat", "Fustat", "Age", "Gender", "Stage", "group")
  
  survplot_cli_trait(data = DataTimeSurM6aScoreclin2, trait = "Age")
  
  survplot_cli_trait(data = DataTimeSurM6aScoreclin2, trait = "Gender")
  
  survplot_cli_trait(data = DataTimeSurM6aScoreclin2, trait = "Stage")
  
  # DOvar = c("Age", "Gender", "Stage")
  # sapply(data, DOvar,survplot_cli_trait)
  
  
  # 原step47.01-----------------------------------------------------------------
  # 原step47.01-----------------------------------------------------------------
  # 原step47.01-----------------------------------------------------------------
  
  
  
  # 原step48.00-----------------------------------------------------------------
  # 原step48.00-----------------------------------------------------------------
  # 原step48.00-----------------------------------------------------------------
  
  #引用包
  library(limma)
  library(ggpubr)
  gene="CD274"            #基因的标准名字
  showName="PD-L1"        #图形里面显示的基因名称
  
  
  # data_expFile_merge="merge.txt"                #表达数据文件
  data_expFile_merge="merge_TCGA_GSE133057.txt"                #表达数据文件
  # scoreFile="m6Ascore.group.txt"     #m6A打分分组文件
  
  
  # setwd("D:\\biowolf\\m6aTME\\48.scoreGene")      #设置工作目录
  
  #读取表达数据文件
  data_expFile_merge_rt=read.table(data_expFile_merge, header=T, sep="\t", check.names=F)
  data_expFile_merge_rt=as.matrix(data_expFile_merge_rt)
  rownames(data_expFile_merge_rt)=data_expFile_merge_rt[,1]
  exp=data_expFile_merge_rt[,2:ncol(data_expFile_merge_rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data_expFile_merge_rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data_expFile_merge_rt=avereps(data_expFile_merge_rt)
  colnames(data_expFile_merge_rt)=gsub("(.*?)\\_(.*?)", "\\2", colnames(data_expFile_merge_rt))
  
  #提取目标基因表达量
  data=rbind(data_expFile_merge_rt, gene=data_expFile_merge_rt[gene,])
  exp=t(data[c("gene",gene),])
  exp=avereps(exp)
  
  #读取m6A打分分组文件
  # score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
  # DataTimeSurM6aScore
  
  #合并数据
  sameSample=intersect(row.names(exp), row.names(DataTimeSurM6aScore))
  exp=exp[sameSample,]
  
  exp[exp>quantile(exp,0.975)]=quantile(exp,0.975)
  DataTimeSurM6aScoreSS=DataTimeSurM6aScore[sameSample,]
  data=cbind(as.data.frame(exp), as.data.frame(DataTimeSurM6aScoreSS))
  
  #设置比较组
  data$group=factor(data$group, levels=c("Low", "High"))
  group=levels(factor(data$group))
  comp=combn(group,2)
  
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  #绘制boxplot
  boxplot=ggboxplot(data, x="group", y="gene", fill="group",
                    xlab="m6Ascore",
                    ylab=paste(showName, "expression"),
                    legend.title="m6Ascore",
                    palette=c("#0066FF","#FF0000"))+ 
    stat_compare_means(comparisons = my_comparisons)
  boxplot
  
  
  #输出图片
  # pdf(file=paste0(gene, ".pdf"), width=5, height=4.5)
  pdf(file=paste0(resultPath,"/",gene, "_",method,".pdf"), width = 5, height =4.5)
  print(boxplot)
  dev.off()
  
  
  
  
  # 原step48.01-----------------------------------------------------------------
  # 原step48.01-----------------------------------------------------------------
  # 原step48.01-----------------------------------------------------------------
  
  # 原step49.00-----------------------------------------------------------------
  # 原step49.00-----------------------------------------------------------------
  # 原step49.00-----------------------------------------------------------------
  
  library(ggpubr)                    #引用包
  tciaFile="TCIA.txt"                #免疫治疗打分文件
  # scoreFile="m6Ascore.group.txt"     #m6A打分分组文件
  # setwd("D:\\biowolf\\m6aTME\\49.IPS")     #修改工作目录
  
  #读取免疫治疗打分文件
  ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)
  
  #读取m6A打分分组文件
  # score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
  # DataTimeSurM6aScore
  
  
  
  #合并数据
  sameSample=intersect(row.names(ips), row.names(DataTimeSurM6aScore))
  ips=ips[sameSample, , drop=F]
  
  DataTimeSurM6aScore_group=DataTimeSurM6aScore[sameSample, "group", drop=F]
  
  ipsFun <- function(ips, DataTimeSurM6aScore_group, i, resultPath) {
    data=cbind(ips, DataTimeSurM6aScore_group)
    #设置比较组
    data$group=factor(data$group, levels=c("Low", "High"))
    group=levels(factor(data$group))
    comp=combn(group, 2)
    my_comparisons=list()
    for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
    
    #对免疫治疗打分进行循环,分别绘制小提琴图
    for(i in colnames(data)[1:(ncol(data)-1)]){
      rt=data[,c(i, "group")]
      colnames(rt)=c("IPS", "group")
      gg1=ggviolin(rt, x="group", y="IPS", fill = "group", 
                   xlab="m6Ascore", ylab=i,
                   legend.title="m6Ascore",
                   palette=c("#0066FF", "#FF0000"),
                   add = "boxplot", add.params = list(fill="white"))+ 
        stat_compare_means(comparisons = my_comparisons)
      #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
      
      pdf(file=paste0(resultPath,"/",i, "_",method,".pdf"), width=6, height=5)
      print(gg1)
      dev.off()
    }
    
  }
  
  ipsFun(ips, DataTimeSurM6aScore_group, i, resultPath)
  
  # 原step49.01-----------------------------------------------------------------
  # 原step49.01-----------------------------------------------------------------
  # 原step49.01-----------------------------------------------------------------
  
  # 原step50.00-----------------------------------------------------------------
  # 原step50.00-----------------------------------------------------------------
  # 原step50.00-----------------------------------------------------------------
  
  #引用包
  library(plyr)
  library(ggplot2)
  library(ggpubr)
  
  # scoreFile="m6Ascore.group.txt"    #m6A评分的分组文件
  cliFile_MSI="MSI.txt"                 #微卫星不稳定性文件
  trait="MSI"                       #MSI
  # setwd("D:\\biowolf\\m6aTME\\50.MSI")     #设置工作目录
  
  #读取输入文件
  # score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
  # DataTimeSurM6aScore
  
  data_cliFile_MSI=read.table(cliFile_MSI, header=T, sep="\t", check.names=F, row.names=1)
  sameSample=intersect(row.names(DataTimeSurM6aScore), row.names(data_cliFile_MSI))
  
  DataTimeSurM6aScoreMSI=cbind(DataTimeSurM6aScore[sameSample,,drop=F], data_cliFile_MSI[sameSample,,drop=F])
  
  
  MSIFun <- function(rt,trait) {
    rt$MSI=factor(rt$MSI, levels=c("MSS", "MSI-L", "MSI-H"))
    rt$group=factor(rt$group, levels=c("Low", "High"))
    bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
    bioCol=bioCol[1:length(unique(rt[,trait]))]
    
    #统计高低评分组病人数目
    rt1=rt[,c(trait, "group")]
    colnames(rt1)=c("trait", "group")
    df=as.data.frame(table(rt1))
    #计算高低评分组的百分率
    df=ddply(df, .(group), transform, percent = Freq/sum(Freq) * 100)
    #百分比位置
    df=ddply(df, .(group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
    df$label=paste0(sprintf("%.0f", df$percent), "%")
    #绘制百分率图
    p=ggplot(df, aes(x = factor(group), y = percent, fill = trait)) +
      geom_bar(position = position_stack(), stat = "identity", width = .7) +
      scale_fill_manual(values=bioCol)+
      xlab("m6Ascore")+ ylab("Percent weight")+  guides(fill=guide_legend(title=trait))+
      geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
      #coord_flip()+
      theme_bw()
    
    pdf(file=paste0(resultPath,"/",method,"_MSI_barplot.pdf"), width=4, height=5)
    print(p)
    dev.off()
    
    #设置比较组
    rt2=rt[,c(trait, "m6Ascore")]
    colnames(rt2)=c("trait", "m6Ascore")
    type=levels(factor(rt2[,"trait"]))
    comp=combn(type, 2)
    my_comparisons=list()
    for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
    #绘制箱线图
    boxplot=ggboxplot(rt2, x="trait", y="m6Ascore", fill="trait",
                      xlab="",
                      ylab="m6Ascore",
                      legend.title=trait,
                      palette=bioCol)+ 
      stat_compare_means(comparisons=my_comparisons)
    
    pdf(file=paste0(resultPath,"/",method,"_MSI_boxplot.pdf"), width=4, height=5)
    print(boxplot)
    dev.off()
  }
  
  
  
  MSIFun(rt = DataTimeSurM6aScoreMSI,trait = "MSI"  )
  
  # 原step50.01-----------------------------------------------------------------
  # 原step50.01-----------------------------------------------------------------
  # 原step50.01-----------------------------------------------------------------
  
  save.image(file = paste0("20220125_0223_",method,".RData"))
  
   
  
}



# end ---------------------------------------------------------------------
# end ---------------------------------------------------------------------
# end ---------------------------------------------------------------------


