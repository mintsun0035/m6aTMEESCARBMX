
# GB2312



# O:\m6a-esca\57.m6AriskScoreIMvigor210_immunotherapy


rm(list = ls())

# setwd("O:\\m6a-esca\\57.m6AriskScoreIMvigor210_immunotherapy")     #设置工作目录

setwd("I:/m6a-esca/58.m6AriskScore_Riaz")


Sys.time()
# [1] "2022-10-26 22:19:19 CST"

setwd("I:/m6a-esca/58.m6AriskScore_Riaz")


expFile="uniSigGeneExp.txt"      #表达输入文件
# setwd("C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\37.PCAscore")     #设置工作目录

#读取输入文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
# sampleGeneMatrix

m6a_escaGene = colnames(data)

#
#PCA分析
pca=prcomp(data, scale=TRUE)
value=predict(pca)
m6Ascore=value[,1]+value[,2]
m6Ascore=as.data.frame(m6Ascore)


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
riaz = data.table::fread("rld.BMS038.20171011_geneSymbol.csv", header=T, sep=",", check.names=F)
riaz = riaz[order(riaz$V109),]

library(dplyr)
library(tidyr)
library(limma)


# 读取表达数据00 ----------------------------------------------------------------


riaz1 = riaz %>% dplyr::select(V109, everything()) %>% .[-c(1:27),-2]

riaz2  = avereps(riaz1 ,ID = riaz1$V109)

riaz2 <- data.frame(riaz2)

riaz3 =  riaz2 %>% tidyr::drop_na(V109)
rownames(riaz3 ) <- riaz3$V109
riaz3<- riaz3[,-1]
head(riaz3)[,1:6]
riaz4 <- apply(riaz3,2,as.numeric) %>% as.data.frame()
rownames(riaz4 )  = rownames(riaz3)
# 读取表达数据01 ----------------------------------------------------------------



# 读取临床数据00 ----------------------------------------------------------------
clinic = data.table::fread("bms038_clinical_data.csv", header=T, sep=",", check.names=F)





colnames(riaz4 )  = tolower(colnames(riaz4 ))

clinic$Sample  = tolower(clinic$Sample)

clinic$Sample %in% colnames(riaz4)


riaz5 = riaz4 [,colnames(riaz4) %in% clinic$Sample]

clinic2 = clinic[match(colnames(riaz5),clinic$Sample),]


m6a_escaGene %in% rownames(riaz5)
table(m6a_escaGene %in% rownames(riaz5))

riaz6 = riaz5[match(m6a_escaGene,rownames(riaz5)),] %>% na.omit()



# 读取临床数据01 ----------------------------------------------------------------

# rownames(riaz2 ) <- riaz2$V109
# exprset <- exprset[,-c(1:2)]



# row.names(riaz) = riaz$V109

# vigor=t(vigor[coxGene,])

# vigor=t(vigor[m6a_escaGene,])
data = t(riaz6 )  #sampleGeneMatrix



#PCA分析
pca=prcomp(data, scale=TRUE)
value=predict(pca)
m6Ascore=value[,1]-value[,2]
m6Ascore=as.data.frame(m6Ascore)







#读取生存数据文件
# cli=read.table("time.txt", header=T, sep="\t", check.names=F, row.names=1)

cli=clinic2[,c("PatientID","OS","OS_SOR")]

cli= cli %>% tibble::column_to_rownames("PatientID")
cli$OS=cli$OS/30




colnames(cli)=c("futime", "fustat")

#合并IMvigor210数据
# sameSample=intersect(row.names(vigor), row.names(cli))
# vigorTime=cbind(cli[sameSample,,drop=F], vigor[sameSample,,drop=F])
#vigorTime[,3:ncol(vigorTime)]=vigorTime[,3:ncol(vigorTime)]*median(as.matrix(rt[,coxGene]))/median(as.matrix(vigorTime[,3:ncol(vigorTime)]))

#输出IMvigor210的风险值
# vigorScore=predict(multiCox, type="risk", newdata=vigorTime)

range(cli$futime)
# [1]  0.3333333 38.1333333

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
vigorRiskOut=cbind(id=rownames(vigorRiskOut), vigorRiskOut)


write.table(vigorRiskOut,file="risk.IMvigor.txt",sep="\t",quote=F,row.names=F)




#绘制生存曲线函数00------------------------
bioSurvival=function(inputFile=null, outFile=null){
  #读取输入文件
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
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
bioSurvival(inputFile="risk.IMvigor.txt", outFile="sur.IMvigor.pdf")
bioROC(inputFile="risk.IMvigor.txt", outFile="ROC.IMvigor.pdf")




# -------------------------------------------------------------------------


# # 第2种可能最优cutoff ----------------------------- --------------------------------------
#获取最优cutoff
res.cut=surv_cutpoint(vigorRiskOut, time="futime", event="fustat", variables=c("riskScore"),minprop = 0.3)
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)

Type=ifelse(vigorRiskOut[,"riskScore"]>cutoff, "Low", "High")
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

#print(surv_median(fit))

sdf = diff
HR = (sdf$obs[1]/sdf$exp[1])/(sdf$obs[2]/sdf$exp[2])  #高风险生存率低yes
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdf$exp[1]+1/sdf$exp[2]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdf$exp[1]+1/sdf$exp[2]))

HR <- paste("Hazard Ratio = ", round(HR ,3), sep = "")
CI <- paste("95% CI: ", paste(round(low95,3), round(up95,3), sep = " - "), sep = "")

HR
CI

sdf = diff
HR = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])  #高风险生存率低yes
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))



HR <- paste("Hazard Ratio = ", round(HR ,3), sep = "")
CI <- paste("95% CI: ", paste(round(low95,3), round(up95,3), sep = " - "), sep = "")
HR
CI



#绘制生存曲线
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]


# 调用调色板
library("RColorBrewer")
display.brewer.all()

cols <- brewer.pal(14, "Set2")
length(cols)
bioCol = cols 

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
                   break.time.by = 2,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25                   )


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

pdf(file=paste0("m6AscoreSurvival_","_最优cutoff-1.pdf"), onefile = FALSE, width=5.5, height=5.5)
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


dput(colnames(clinic2))
c("PatientID", "Sample", "SampleType", "Cohort", "SubtypeEZ", 
  "TRTGRP", "BOR", "myBOR", "PFS_SOR", "OS_SOR", "OS", "OSWK", 
  "IBOR", "PFS", "PFSWK", "myBOR2", "myBOR3")

clinic3 =clinic2 [,c("SubtypeEZ", "BOR", "myBOR", "IBOR", "myBOR2", "myBOR3")]

rt=cbind(risk, clinic3)


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


# trait="binaryResponse" 
trait="myBOR3"


# 调用调色板
library("RColorBrewer")
display.brewer.all()

cols <- brewer.pal(14, "Set1")
length(cols)
# 调用wgcna的标准颜色---------------------
library(WGCNA)
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



save.image(file = "01.RData")

