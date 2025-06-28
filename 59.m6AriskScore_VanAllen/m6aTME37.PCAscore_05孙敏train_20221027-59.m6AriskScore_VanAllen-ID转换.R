





# GB2312



# O:\m6a-esca\57.m6AriskScoreIMvigor210_immunotherapy


rm(list = ls())

# setwd("O:\\m6a-esca\\57.m6AriskScoreIMvigor210_immunotherapy")     #设置工作目录

setwd("I:/m6a-esca/59.m6AriskScore_VanAllen")

Sys.time()
# [1] "2022-10-26 22:19:19 CST"




expFile="uniSigGeneExp.txt"      #表达输入文件
# setwd("C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\37.PCAscore")     #设置工作目录

#读取输入文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
# sampleGeneMatrix

m6a_escaGene = colnames(data)

#
#PCA分析
# pca=prcomp(data, scale=TRUE)
# value=predict(pca)
# m6Ascore=value[,1]+value[,2]
# m6Ascore=as.data.frame(m6Ascore)


# scoreOut=rbind(id=colnames(m6Ascore), m6Ascore)
# write.table(scoreOut, file="m6Ascore.txt", sep="\t", quote=F, col.names=F)

#引用包
library(survival)
library(survminer)


#引用包
library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)





# 以下不运行00------------------------- ----------------------------------------------






scoreFile="m6Ascore.txt"     #m6A打分文件
cliFile="time.txt"           #生存数据文件
# setwd("C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\38.scoreSur")      #设置工作目录

#读取输入文件
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
sampleType=gsub("(.*?)\\_.*", "\\1", row.names(score))
score=cbind(score, sampleType)
rownames(score)=gsub("(.*?)\\_(.*?)", "\\2", rownames(score))
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

#数据合并
sameSample=intersect(row.names(score), row.names(cli))
data=cbind(cli[sameSample,], score[sameSample,])

#获取最优cutoff
res.cut=surv_cutpoint(data, time="futime", event="fustat", variables=c("m6Ascore"))
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)
Type=ifelse(data[,"m6Ascore"]<=cutoff, "Low", "High")
data$group=Type
outTab=rbind(id=colnames(data), data)
write.table(outTab, file="m6Ascore.group.txt", sep="\t", quote=F, col.names=F)

#计算高低风险组生存差异
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = data)
#print(surv_median(fit))

#绘制生存曲线
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data=data,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="m6Ascore",
                   legend.labs=levels(factor(data[,"group"])),
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
pdf(file="survival-2.pdf", onefile = FALSE, width=5.5, height=5.5)
print(surPlot)
dev.off()


#保存图片
pdf(file="survival-3.pdf", onefile = FALSE, width=5, height=5.5)
print(surPlot)
dev.off()




# 以上不运行01------------------------- ----------------------------------------------




#读取IMvigor210表达数据
# vigor=read.table("exp.txt", header=T, sep="\t", check.names=F, row.names=1)

library(data.table)
# VanAllen = data.table::fread("rld.BMS038.20171011_geneSymbol.csv", header=T, sep=",", check.names=F)

# VanAllen = data.table::fread("20160304_MEL-TPM_noDups.txt", header=T, sep="\t", check.names=F)
VanAllen = data.table::fread("20160304_MEL-TPM_noDups.txt", header=T, sep="\t", check.names=F)
VanAllen = VanAllen[order(gene_id),]
dim(VanAllen)
# [1] 57820    52


VanAllen  = VanAllen %>% tidyr::separate(gene_id,into = c("gene_id","drop"),sep="\\.") %>% 
  dplyr::select(-c("drop"))

# class(VanAllen )
VanAllen = as.data.frame(VanAllen)

# expr_df_nopoint <- expr_df %>% 
#   tidyr::separate(gene_id,into = c("gene_id"),sep="\\.") 

# row.names(VanAllen) = as.character(VanAllen$gene_id)
row.names(VanAllen) = as.character(VanAllen$gene_id)

# exprSet = VanAllen

exprSet = VanAllen[,-1]
row.names(exprSet)


ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
LogC

if (LogC) { ex[which(ex <= 0)] <- NaN
  exprSet <- log2(ex)
  print("log2 transform finished")}else{print("log2 transform not needed")}



library(dplyr)
library(tidyr)
library(limma)

# annotateFile = data.table::fread("ENSG.数字-探针转换为geneSymbol_ann.txt", header=T, sep="\t", check.names=F)
# annotateFile = data.table::fread("ENSG_2_geneSymbolName.txt", header=T, sep="\t", check.names=F)

##加载gtf文件---------------------------
load(file = "gtf_df.Rda")
# test <- gtf_df[1:100,]

## mRNA
# mRNA_exprSet <- gtf_df %>% 
#   dplyr::filter(type=="gene",gene_biotype=="protein_coding") %>% #筛选gene,和编码指标
#   dplyr::select(c(gene_name,gene_id,gene_biotype)) %>% 
#   # dplyr::inner_join(expr_df_nopoint,by ="gene_id") %>% 
#   dplyr::inner_join(exprSet,by ="gene_id")

annotateFile = gtf_df %>% 
  dplyr::select(c(gene_name,gene_id)) 


# annotateFile = annotateFile[order(ID),]
# probe2symbol = annotateFile  %>% dplyr::distinct(gene_id)
# probe2symbol = data.frame(probeset,symbol,stringsAsFactors = F)


probe2symbol = annotateFile
# colnames(probe2symbol) = c("probeset","symbol")
colnames(probe2symbol) = c("symbol","probeset")

# 6.探针转换与基因去重00------- ---------------------


library(dplyr)
library(tibble)

# intersect(probe2symbol$probeset,rownames(exprSet))


exprSet2 <- exprSet %>% 
  rownames_to_column(var="probeset") %>% 
  #合并探针的信息
  inner_join(probe2symbol,by="probeset") %>% 
  #去掉多余信息
  select(-probeset) %>% 
  #重新排列
  select(symbol,everything()) %>% 
  #求出平均数(这边的点号代表上一步产出的数据)
  mutate(rowMean =rowMeans(.[grep("MEL", names(.))])) %>% 
  #去除symbol中的NA
  filter(symbol != "NA") %>% 
  #把表达量的平均值按从大到小排序
  arrange(desc(rowMean)) %>% 
  # symbol留下第一个
  distinct(symbol,.keep_all = T) %>% 
  #反向选择去除rowMean这一列
  select(-rowMean) %>% 
  # 列名变成行名
  column_to_rownames(var = "symbol")

# exprSet2SampleName = as.data.frame(colnames(exprSet2))
# colnames(exprSet2SampleName) = sample

# exprSet2SampleName= strsplit(exprSet2SampleName,"\\_")[[1]]

exprSet2SampleName = colnames(exprSet2)

Sampleid = c()

for (i in colnames(exprSet2)){
  
  test=strsplit(i,split="_")
  
  Sampleid = c(Sampleid,unlist(test)[2])
  
}


SampleName = c()

for (i in Sampleid){
  
  test=strsplit(i,split="-")
  
  SampleName = c(SampleName,unlist(test)[1])
  
}

colnames(exprSet2) = SampleName 



# for (i in colnames(exprSet2)){
#   
#   test=strsplit(i,split="_")
#   
#   Sampleid = c(Sampleid,print(paste(unlist(test)[1:2],collapse="\\-")))
#   
# }






# library(stringi)
# T <- as.data.frame(do.call(rbind,  stri_split_fixed(exprSet2SampleName, "\\_", 1)))







# 读取表达数据00 ----------------------------------------------------------------


# VanAllen1 = VanAllen %>% dplyr::select(V109, everything()) %>% .[-c(1:27),-2]

# VanAllen2  = avereps(VanAllen1 ,ID = VanAllen1$V109)

# VanAllen2 <- data.frame(VanAllen2)

# VanAllen3 =  VanAllen2 %>% tidyr::drop_na(V109)
# rownames(VanAllen3 ) <- VanAllen3$V109
# VanAllen3<- VanAllen3[,-1]
# head(VanAllen3)[,1:6]
# VanAllen4 <- apply(VanAllen3,2,as.numeric) %>% as.data.frame()
# rownames(VanAllen4 )  = rownames(VanAllen3)
# 读取表达数据01 ----------------------------------------------------------------



# 读取临床数据00 ----------------------------------------------------------------
library(openxlsx)

clinic = openxlsx::read.xlsx("tables2.clinical_and_genome_characteristics_each_patient.xlsx",sheet = 2)


# clinic = data.table::fread("bms038_clinical_data.csv", header=T, sep=",", check.names=F)
table(duplicated(clinic$patient))
# FALSE 
#    42


exprSet3 = exprSet2[,clinic$patient]


# colnames(VanAllen4 )  = tolower(colnames(VanAllen4 ))

# clinic$Sample  = tolower(clinic$Sample)

# clinic$Sample %in% colnames(VanAllen4)


# VanAllen5 = VanAllen4 [,colnames(VanAllen4) %in% clinic$Sample]

# clinic2 = clinic[match(colnames(VanAllen5),clinic$Sample),]


m6a_escaGene %in% rownames(exprSet3)


table(m6a_escaGene %in% rownames(exprSet3))

# VanAllen6 = VanAllen5[match(m6a_escaGene,rownames(VanAllen5)),] %>% na.omit()

exprSet4 = exprSet3[m6a_escaGene,] 


# 读取临床数据01 ----------------------------------------------------------------

# rownames(VanAllen2 ) <- VanAllen2$V109
# exprset <- exprset[,-c(1:2)]



# row.names(VanAllen) = VanAllen$V109

# vigor=t(vigor[coxGene,])

# vigor=t(vigor[m6a_escaGene,])
data = t(exprSet4  )  #sampleGeneMatrix



#PCA分析
pca=prcomp(data, scale=TRUE)
value=predict(pca)
# m6Ascore=value[,1]+value[,2]
m6Ascore=value[,1]-value[,2]
m6Ascore=as.data.frame(m6Ascore)










#读取生存数据文件
# cli=read.table("time.txt", header=T, sep="\t", check.names=F, row.names=1)

cli=clinic[,c("patient","overall_survival","dead")]

cli= cli %>% tibble::column_to_rownames("patient")
cli$overall_survival=cli$overall_survival/30   #month




colnames(cli)=c("futime", "fustat")

#合并IMvigor210数据
# sameSample=intersect(row.names(vigor), row.names(cli))
# vigorTime=cbind(cli[sameSample,,drop=F], vigor[sameSample,,drop=F])
#vigorTime[,3:ncol(vigorTime)]=vigorTime[,3:ncol(vigorTime)]*median(as.matrix(rt[,coxGene]))/median(as.matrix(vigorTime[,3:ncol(vigorTime)]))

#输出IMvigor210的风险值
# vigorScore=predict(multiCox, type="risk", newdata=vigorTime)

range(cli$futime)
# [1]  1.133333 54.400000


# 代码复用00 -------------------------- ------------------------------------------
# 代码复用00 -------------------------- ------------------------------------------
# 代码复用00 -------------------------- ------------------------------------------
# 代码复用00 -------------------------- ------------------------------------------
# 代码复用00 -------------------------- ------------------------------------------


vigorTime = cli

vigorScore = m6Ascore

# Risk=as.vector(ifelse(vigorScore>median(vigorScore), "high", "low"))
Risk=as.vector(ifelse(vigorScore>median(vigorScore$m6Ascore), "high", "low"))
vigorRiskOut=cbind(vigorTime, riskScore=as.vector(vigorScore$m6Ascore), Risk)
vigorRiskOut2=cbind(id=rownames(vigorRiskOut), vigorRiskOut)


write.table(vigorRiskOut2,file="risk.IMvigor.txt",sep="\t",quote=F,row.names=F)


#读取输入文件
# rt=utils::read.table(inputFile)
# rt=utils::read.table(inputFile, header = T)
# rt=utils::read.table(inputFile,header = T,sep = "\t", check.names=F)
# rt=read.table(inputFile, header=T, sep="\t", check.names=F)
#比较高低风险组生存差异，得到显著性p值


rt = vigorRiskOut
diff=survdiff(Surv(futime, fustat) ~Risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}

pValue

fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt)

#绘制生存曲线
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="Risk",
                   legend.labs=c("High risk", "Low risk"),
                   xlab="Time(Months)",
                   break.time.by = 2,
                   palette=c("red", "blue"),
                   risk.table=F,
                   risk.table.title="",
                   risk.table.col = "strata",
                   risk.table.height=.25)

print(surPlot)

outFile="sur.IMvigor-00.pdf"
pdf(file=outFile,onefile = FALSE,width = 6,height =5)
print(surPlot)
dev.off()




#绘制生存曲线
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="Risk",
                   legend.labs=c("High risk", "Low risk"),
                   xlab="Time(Months)",
                   break.time.by = 2,
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.title="",
                   risk.table.col = "strata",
                   risk.table.height=.25)

print(surPlot)


outFile="sur.IMvigor-01.pdf"
pdf(file=outFile,onefile = FALSE,width = 5,height =5)
print(surPlot)
dev.off()




#绘制生存曲线
surPlot2=ggsurvplot(fit,
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="Risk",
                   legend.labs=c("High risk", "Low risk"),
                   xlab="Time(Months)",
                   break.time.by = 6,
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.title="",
                   risk.table.col = "strata",
                   risk.table.height=.25)

print(surPlot2)

outFile="sur.IMvigor-02.pdf"
pdf(file=outFile,onefile = FALSE,width = 6,height =5)
print(surPlot2)
dev.off()






#绘制生存曲线函数00----以下不运行00--------------------
inputFile="risk.IMvigor.txt"
outFile="sur.IMvigor.pdf"


bioSurvival=function(inputFile=null, outFile=null){
  #读取输入文件
  # rt=utils::read.table(inputFile)
  # rt=utils::read.table(inputFile, header = T)
  rt=utils::read.table(inputFile,header = T,sep = "\t", check.names=F)
  # rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  #比较高低风险组生存差异，得到显著性p值
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="Risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(Months)",
                     break.time.by = 2,
                     palette=c("red", "blue"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  pdf(file=outFile,onefile = FALSE,width = 6,height =5)
  print(surPlot)
  dev.off()
}




bioSurvival(inputFile="risk.IMvigor.txt", outFile="sur.IMvigor-2.pdf")





range(cli$futime)

#定义绘制ROC曲线函数
# 1yearsurvival ROC
bioROC=function(inputFile=null, outFile=null){
  #读取输入文件
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  #ROC曲线
  # ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
  #                marker=rt$riskScore,cause=1,
  #                weighting='aalen',
  #                times=c(1*12,3*12,5*12),ROC=TRUE)
  
  ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
                 marker=rt$riskScore,cause=1,
                 weighting='aalen',
                 times=c(1*12),ROC=TRUE)
  
  
  pdf(file=outFile, width=5, height=5)
  plot(ROC_rt,time=1*12,col='green',title=FALSE,lwd=2)
  # plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
  # plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
  # legend('bottomright',
  #        c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
  #          paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
  #          paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
  #        col=c("green",'blue','red'),lwd=2,bty = 'n')
  
  legend('bottomright',
         # c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
          
          c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[2]))),
         col=c("green"),lwd=2,bty = 'n')
  
  dev.off()
}
#绘制生存曲线函数01------------------------



#调用函数，绘制生存曲线和ROC曲线00---------------- -----

bioROC(inputFile="risk.IMvigor.txt", outFile="ROC.IMvigor.pdf")


#绘制生存曲线函数00----以上不运行00--------------------


# -------------------------------------------------------------------------


# # 第2种可能最优cutoff ----------------------------- --------------------------------------
#获取最优cutoff
# res.cut=surv_cutpoint(vigorRiskOut, time="futime", event="fustat", variables=c("riskScore"),minprop = 0.3)

# vigorRiskOut = vigorRiskOut[order(vigorRiskOut$futime),]
# vigorRiskOut = vigorRiskOut[-c(41,42),]




res.cut=surv_cutpoint(vigorRiskOut, time="futime", event="fustat", variables=c("riskScore"),minprop = 0.5)
# res.cut=surv_cutpoint(vigorRiskOut, time="futime", event="fustat", variables=c("riskScore"),minprop = 0.35)
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)

# cutoff = median(riskScore)


# Type=ifelse(vigorRiskOut[,"riskScore"]>cutoff, "Low", "High")

Type=ifelse(vigorRiskOut[,"riskScore"]>cutoff, "High", "Low")
vigorRiskOut$group=Type

# outTab=rbind(id=colnames(DataTimeSurM6aScore), DataTimeSurM6aScore)
# write.table(outTab, file=paste0("m6Ascore.group_",method,".txt"), sep="\t", quote=F, col.names=F)

#计算高低风险组生存差异
vigorRiskOut$group=factor(vigorRiskOut$group, levels=c("Low", "High"))

diff=survdiff(Surv(futime, fustat) ~ group, data = vigorRiskOut)

length=length(levels(factor(vigorRiskOut[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
pValue



fit <- survfit(Surv(futime, fustat) ~ group, data = vigorRiskOut)
fit


#print(surv_median(fit))

# sdf = diff
# HR = (sdf$obs[1]/sdf$exp[1])/(sdf$obs[2]/sdf$exp[2])  #高风险生存率低yes
# up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdf$exp[1]+1/sdf$exp[2]))
# low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdf$exp[1]+1/sdf$exp[2]))
# HR <- paste("Hazard Ratio = ", round(HR ,3), sep = "")
# CI <- paste("95% CI: ", paste(round(low95,3), round(up95,3), sep = " - "), sep = "")
# HR
# CI

sdf = diff
HR = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])  #高风险生存率低yes
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
HR <- paste("Hazard Ratio = ", round(HR ,3), sep = "")
CI <- paste("95% CI: ", paste(round(low95,3), round(up95,3), sep = " - "), sep = "")
HR
CI



#绘制生存曲线
# bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
# bioCol=bioCol[1:length]


# 调用调色板
library("RColorBrewer")
display.brewer.all()

cols <- brewer.pal(14, "Set2")
length(cols)
bioCol = cols[1:length] 

fit

surPlot=ggsurvplot(fit,
                   data=vigorRiskOut,
                   conf.int=F,
                   # pval=pValue,
                   pval=paste(pValue,HR, CI, sep = "\n"),
                   pval.size=6,
                   legend.title="m6Ascore",
                   legend.labs=levels(factor(vigorRiskOut[,"group"])),
                   legend = c(0.8, 0.8),
                   font.legend=8,
                   xlab="Time(Months)",
                   break.time.by = 6,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)


print(surPlot)






# 代码复用01-------------------------- ------------------------------------------
# 代码复用01-------------------------- ------------------------------------------
# 代码复用01-------------------------- ------------------------------------------
# 代码复用01-------------------------- ------------------------------------------
# 代码复用01-------------------------- ------------------------------------------













#保存图片
# pdf(file=paste0(resultPath,"/","m6AscoreSurvival_",method,"_最优cutoff.pdf"), onefile = FALSE, width=7, height=5.5)
# print(surPlot)
# dev.off()

pdf(file=paste0("m6AscoreSurvival_","_最优cutoff-2.pdf"), onefile = FALSE, width=5.5, height=5.5)
# pdf(file=paste0("/","m6AscoreSurvival_","_最优cutoff.pdf"), onefile = FALSE, width=7, height=5.5)
print(surPlot)
dev.off()



# 临床特点的关系00 -----------------------------------------------------------------





#install.packages("ggpubr")


#引用包
library(limma)
library(ggpubr)



riskFile="risk.IMvigor.txt"      #风险文件
# cliFile="clinical.txt"           #临床数据文件
# setwd("C:\\biowolf\\NMF\\38.IMvigorCliCor")     #设置工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
# risk$riskScore[risk$riskScore>quantile(risk$riskScore,0.99)]=quantile(risk$riskScore,0.99)

#读取临床数据文件
# cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)






#合并数据
# samSample=intersect(row.names(risk), row.names(cli))
# risk=risk[samSample,"riskScore",drop=F]
# cli=cli[samSample,,drop=F]
# rt=cbind(risk, cli)

risk=risk["riskScore",drop=F]


dput(colnames(clinic))
# c("PatientID", "Sample", "SampleType", "Cohort", "SubtypeEZ", 
#   "TRTGRP", "BOR", "myBOR", "PFS_SOR", "OS_SOR", "OS", "OSWK", 
#   "IBOR", "PFS", "PFSWK", "myBOR2", "myBOR3")


c("patient_number", "patient", "gct", "RECIST", "overall_survival", 
  "progression_free", "progression", "primary", "survival_1yr", 
  "survival_2yr", "response_6mo", "response_1yr", "response_2yr", 
  "response", "nonresponse", "long.survival", "group", "histology", 
  "stage", "M", "gender", "date_birth", "LDH", "pre_therapies", 
  "pre_BRAF", "post_BRAF", "therapy_start", "therapy_end", "date_death", 
  "dead", "date_progression", "age_start")

clinic2=clinic [,c("patient","RECIST")]

clinic2$BOR = ifelse(clinic$RECIST == "PD","PD","CRSD")

# table(rt$clinic2)
# # CR PD PR SD  X 
# # 2 27  5  7  1 
clinic2$myBOR <- ifelse(
  clinic$RECIST == "X", "PDSD",
     ifelse(
       clinic$RECIST == "SD", "PDSD",
         ifelse(
           clinic$RECIST == "CR", "CRPR",
             ifelse(
               clinic$RECIST == "PR", "CRPR",
                  "PDSD"
                   )
               )
           )
       )



# clinic3 =clinic2 [,c("SubtypeEZ", "BOR", "myBOR", "IBOR", "myBOR2", "myBOR3")]

rt=cbind(risk, clinic2)

# table(rt$clinic2)
# CR PD PR SD  X 
# 2 27  5  7  1 

rt$patient <- NULL


#临床相关性分析，输出图形结果
for(clinical in colnames(rt)[2:ncol(rt)]){
  data=rt[c("riskScore", clinical)]
  colnames(data)=c("riskScore", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  # data=data[(data[,"clinical"] != "NA"),]
  data = na.omit(data)
  #设置比较组
  group=levels(factor(data$clinical))
  data$clinical=factor(data$clinical, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  #绘制箱线图
  boxplot=ggboxplot(data, x="clinical", y="riskScore", color="clinical",
                    xlab="",
                    ylab="Risk score",
                    legend.title=clinical,
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  #输出图片
  pdf(file=paste0("cliCor.", clinical, ".pdf"), width=6, height=5)
  print(boxplot)
  dev.off()
}



# 临床特点的关系01-----------------------------------------------------------------



# bar图00-----------------------------------------------------------------
library(plyr)
library(ggplot2)
library(ggpubr)



dput(colnames(rt))
# c("riskScore", "Best Confirmed Overall Response", "binaryResponse", 
# "Enrollment IC", "IC Level", "TC Level", "Immune phenotype", 
# "FMOne mutation burden per MB", "Sex", "Race", "Intravesical BCG administered", 
# "Baseline ECOG Score", "Tobacco Use History", "Met Disease Status", 
# "Sample age", "Tissue", "Received platinum", "Sample collected pre-platinum", 
# "Neoantigen burden per MB", "sizeFactor", "ANONPT_ID", "os", 
# "censOS", "Lund", "Lund2", "TCGA Subtype")

Type=ifelse(rt[,"riskScore"]<cutoff, "Low", "High")
rt$group=Type


# -----1------------------------------------------------------------------

# trait="binaryResponse" 
# trait="myBOR3"

trait="BOR"



# 调用调色板
library("RColorBrewer")
display.brewer.all()

cols <- brewer.pal(14, "Set1")
length(cols)
# 调用wgcna的标准颜色---------------------
# library(WGCNA)
# WGCNA::standardColors()

# color14<-WGCNA::standardColors()[1:14]
# color14



# bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
# bioCol = WGCNA::standardColors()
bioCol <- brewer.pal(14, "Set2")

bioCol=bioCol[1:length(unique(rt[,trait]))]

rt1=rt[,c(trait, "group")]


colnames(rt1)=c("trait", "group")
df=as.data.frame(table(rt1))


df=ddply(df, .(group), transform, percent = Freq/sum(Freq) * 100)

df=ddply(df, .(group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
df$group=factor(df$group, levels=c("Low", "High"))

p=ggplot(df, aes(x = factor(group), y = percent, fill = trait)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("m6Ascore")+ ylab("Percent weight")+  guides(fill=guide_legend(title=trait))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()

pdf(file="barplot.pdf", width=4, height=5)
print(p)
dev.off()




# -----2以下分组无意义00 -------------  -----------------------------------------------------


# trait="binaryResponse" 
# trait="myBOR3"

trait="myBOR"



# 调用调色板
library("RColorBrewer")
display.brewer.all()

cols <- brewer.pal(14, "Set1")
length(cols)
# 调用wgcna的标准颜色---------------------
# library(WGCNA)
# WGCNA::standardColors()

# color14<-WGCNA::standardColors()[1:14]
# color14



# bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
# bioCol = WGCNA::standardColors()
bioCol <- brewer.pal(14, "Set2")

bioCol=bioCol[1:length(unique(rt[,trait]))]

rt1=rt[,c(trait, "group")]


colnames(rt1)=c("trait", "group")
df=as.data.frame(table(rt1))


df=ddply(df, .(group), transform, percent = Freq/sum(Freq) * 100)

df=ddply(df, .(group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
df$group=factor(df$group, levels=c("Low", "High"))

p=ggplot(df, aes(x = factor(group), y = percent, fill = trait)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("m6Ascore")+ ylab("Percent weight")+  guides(fill=guide_legend(title=trait))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()

pdf(file="barplot-myBOR2.pdf", width=4, height=5)
print(p)
dev.off()


# -----2以上分组无意义01----------------- -------------------------------------------------




save.image(file = "01.RData")



