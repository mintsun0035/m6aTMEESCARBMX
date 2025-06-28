#install.packages("survivalROC")


rm(list = ls())

setwd("O:\\m6a-esca\\46.scoreCli_ROC")




library(survivalROC)
library(dplyr)
library(timeROC)


# rt=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)
rt=read.table("m6Ascore.group.txt",header=T,sep="\t",check.names=F,row.names=1)
rt = dplyr::rename (rt,riskScore = m6Ascore)




#1年ROC
pdf(file="ROC-1.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
      predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
  xlab="1-Specificity", ylab="Sensitivity",
  main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#3年ROC
pdf(file="ROC-3.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                predict.time =3, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#5年ROC
pdf(file="ROC-5.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                predict.time =5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()


par()

#整合1，3，5年ROC
pdf(file="ROC.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)

roc1=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                 predict.time =1, method="KM")
plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
text(locator(1), paste("1 year",round(roc1$AUC,3),sep=":"),col="red")


roc2=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                 predict.time =3, method="KM")   #在此更改时间，单位为年
lines(roc2$FP,roc2$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="blue",lwd=2)
text(locator(1), paste("3 year",round(roc2$AUC,3),sep=":"),col="blue")

roc3=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                 predict.time =2, method="KM")   #在此更改时间，单位为年
lines(roc3$FP,roc3$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="green",lwd=2)
text(locator(1), paste("3.1 year",round(roc3$AUC,3),sep=":"),col="green")

# legend("bottomright", 
#        c("1-year AUC:0.780","3-year AUC:0.748","5-year AUC:0.706"),
#        lwd=2,
#        col=c("red","blue","green"))
dev.off()





#引用包
library(survival)
library(survminer)
library(timeROC)
library(stringr)



risk = rt

#定义颜色
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)





######绘制1 2 3年的ROC曲线######
ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
               marker=risk$riskScore,cause=1,
               weighting='aalen',
               times=c(1,2,3),ROC=TRUE)
pdf(file="ROC123.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=2,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=bioCol[1:3], lwd=2, bty = 'n')
dev.off()



# filter_all(any_vars(str_detect(., pattern = "Ca")))
# filter_at(mtcars, vars(starts_with("d")), any_vars((. %% 2) == 0))
risk_TCGA = risk %>%  tibble::rownames_to_column("id") %>% 
  dplyr:: filter_all(any_vars(str_detect(., pattern = "TCGA")))


######绘制1 2 3年的ROC曲线######
ROC_rt_TCGA=timeROC(T=risk_TCGA$futime,delta=risk_TCGA$fustat,
               marker=risk_TCGA$riskScore,cause=1,
               weighting='aalen',
               times=c(1,2,3),ROC=TRUE)
pdf(file="ROC123_rt_TCGA.pdf", width=5, height=5)
plot(ROC_rt_TCGA,time=1,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt_TCGA,time=2,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt_TCGA,time=3,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt_TCGA$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt_TCGA$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt_TCGA$AUC[3]))),
       col=bioCol[1:3], lwd=2, bty = 'n')
dev.off()














# ----------- -------------------------------------------------------------












#引用包
library(survival)
library(survminer)
library(timeROC)
riskFile="risk.txt"         #风险输入文件
cliFile="clinical.txt"      #临床数据文件
# setwd("D:\\biowolf\\FerrLnc\\19.ROC")     #修改工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#定义颜色
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)


######绘制1 2 3年的ROC曲线######
ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
               marker=risk$riskScore,cause=1,
               weighting='aalen',
               times=c(1,2,3),ROC=TRUE)
pdf(file="ROC.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=2,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=bioCol[1:3], lwd=2, bty = 'n')
dev.off()


######绘制临床的ROC曲线######
predictTime=1     #定义预测年限
aucText=c()
pdf(file="cliROC.pdf", width=6, height=6)
#绘制风险得分的ROC曲线
i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)
#对临床数据进行循环，绘制临床数据的ROC曲线
for(i in 4:ncol(rt)){
  ROC_rt=timeROC(T=rt$futime,
                 delta=rt$fustat,
                 marker=rt[,i], cause=1,
                 weighting='aalen',
                 times=c(predictTime),ROC=TRUE)
  plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2, add=TRUE)
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}
#绘制图例，得到ROC曲线下的面积
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()














