
# uft-8-------------


rm(list = ls())
setwd("O:/m6a-esca/56.uniCoxMultliCoxNomogram")

library(survival)     

library(data.table)
library(openxlsx)
library(dplyr)
library(tidyr)
library(tidyverse)
library(tibble)
library(survivalROC)         #引用包



rt3  = fread(file = "m6Ascore.group.txt")



# tcga00.00 --------------------------------------------------------------------

# rt1 =  fread(file = "clinical_TCGA.xls") %>% filter(Id %in% rt3$id)  
# rt1 =  fread(file = "TCGA-ESCA_clinical剔除生存数据不明确.xls") %>% filter(Id %in% rt3$id)  

# rt1 =  fread(file = "食管癌.txt") %>% filter(Id %in% rt3$id)  
rt1 =  fread(file = "食管癌.txt") %>% filter(id %in% rt3$id)  


# rt1 = rt1 %>% dplyr::select("Id","Age","Gender","Stage") %>% as.data.frame()
dput(colnames(rt1))
c("id", "age", "sex", "race", "smoking", "alcohol", "radiation", 
"pharmaceutical", "tumor_grade", "pathologic_stage", "stage_T", 
"stage_M", "stage_N", "status", "survival_time")

# rt1$pathologic_stage

rt1$pathologic_stage =   str_match(rt1$pathologic_stage, "Stage ([^ABC]*)")[,2]


rt1b = rt1


# rt1$Stage =   str_match(rt1$Stage, "Stage ([^ABC]*)")[,2]
table(rt1$pathologic_stage)



rt1$pathologic_stage = tolower(rt1$pathologic_stage)
table(rt1$pathologic_stage)




# rt1$Stage<-gsub(pattern = "stage iva", replacement = "4", rt1$Stage)
# rt1$Stage<-gsub(pattern = "stage iv", replacement = "4", rt1$Stage) 
# rt1$Stage<-gsub(pattern = "stage iiic", replacement = "3", rt1$Stage) 
# rt1$Stage<-gsub(pattern = "stage iiib", replacement = "3", rt1$Stage) 
# rt1$Stage<-gsub(pattern = "stage iiia", replacement = "3", rt1$Stage) 
# rt1$Stage<-gsub(pattern = "stage iii", replacement = "3", rt1$Stage)
# rt1$Stage<-gsub(pattern = "stage iic", replacement = "2", rt1$Stage) 
# rt1$Stage<-gsub(pattern = "stage iib", replacement = "2", rt1$Stage) 
# rt1$Stage<-gsub(pattern = "stage iia", replacement = "2", rt1$Stage) 
# rt1$Stage<-gsub(pattern = "stage ii", replacement = "2", rt1$Stage)
# rt1$Stage<-gsub(pattern = "stage i", replacement = "1", rt1$Stage)


rt1$pathologic_stage=gsub(pattern = "iv", replacement = "4", rt1$pathologic_stage)
rt1$pathologic_stage=gsub(pattern = "iii", replacement = "3", rt1$pathologic_stage)
rt1$pathologic_stage=gsub(pattern = "ii", replacement = "2", rt1$pathologic_stage)
rt1$pathologic_stage=gsub(pattern = "i", replacement = "1", rt1$pathologic_stage)
table(rt1$pathologic_stage)

rt1$Stage <- NULL
rt1$Gender <- NULL


# rt1$Stage<-gsub(pattern = "unknow", replacement = "1", rt1$Stage)
# table(rt1$Stage)

rt1 = dplyr::rename(rt1, gender = sex)
rt1$gender = ifelse(rt1$gender == "FEMALE", "0", "1" )

table(rt1$race)
rt1$race <- NULL
table(rt1$smoke)



rt1_filter = rt1 %>% dplyr::select("id","age","gender","tumor_grade", "pathologic_stage", "stage_T", 
                            "stage_M", "stage_N") %>% as.data.frame()


rt1_filter = rt1_filter %>% dplyr::rename(grade = tumor_grade) %>% 
  rename(stage = pathologic_stage)


table(rt1_filter$grade)
# rt1_filter$grade = ifelse()

rt1_filter$grade =   gsub(pattern = "G", replacement = "", rt1_filter$grade)
rt1_filter$grade =   gsub(pattern = "X", replacement = "4", rt1_filter$grade)

table(rt1_filter$stage_T)
rt1_filter$stage_T =   gsub(pattern = "a", replacement = "", rt1_filter$stage_T)
rt1_filter$stage_T =   gsub(pattern = "T", replacement = "", rt1_filter$stage_T)

table(rt1_filter$stage_M)
rt1_filter$stage_M =   gsub(pattern = "a", replacement = "", rt1_filter$stage_M)
rt1_filter$stage_M =   gsub(pattern = "M", replacement = "", rt1_filter$stage_M)
rt1_filter$stage_M =   gsub(pattern = "X", replacement = "1", rt1_filter$stage_M)
table(rt1_filter$stage_M)


table(rt1_filter$stage_N)
rt1_filter$stage_N =   gsub(pattern = "N", replacement = "", rt1_filter$stage_N)
rt1_filter$stage_N =   gsub(pattern = "X", replacement = "0", rt1_filter$stage_N)


# rt1_filter = apply(rt1_filter, 2:ncol(rt1_filter),as.numeric)



#merge riskdata and clinic data------------------------------------------------------------------
# samSample=intersect(row.names(risk), row.names(cli))
# risk1=risk[samSample,,drop=F]
# cli=cli[samSample,,drop=F]
# rt=cbind(risk1, cli)

cliRiskMerge = merge(rt1_filter, rt3,by = "id")
cliRiskMerge$sampleType <- NULL
cliRiskMerge$group<- NULL

cliRiskMergeNOid = cliRiskMerge %>% tibble::column_to_rownames(var = "id") %>% 
  dplyr::select(fustat, everything()) %>% dplyr::select(futime, everything()) 
cliRiskMergeNOid = apply(cliRiskMergeNOid,2,as.numeric) %>% as.data.frame()

rownames(cliRiskMergeNOid) = cliRiskMerge$id






# colnames(rt1) =tolower(colnames(rt1))
# tcga00.01--------------------------------------------------------------------






library(survival)       


########################
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
	#
	rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrLow[hrLow<0.001]=0.001
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
	#
	pdf(file=forestFile, width=6.5, height=4.8)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	#
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
	text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
		
	#
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	LOGindex=2 
	hrLow = log(as.numeric(hrLow),LOGindex)
	hrHigh = log(as.numeric(hrHigh),LOGindex)
	hr = log(as.numeric(hr),LOGindex)
	xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
	abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
	boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), forestCol, forestCol)
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=2)
	a1 = axis(1,labels=F,tick=F)
	axis(1,a1,LOGindex^a1)
	dev.off()
}
########################

#
# indep=function(riskFile=null,cliFile=null,uniOutFile=null,multiOutFile=null,uniForest=null,multiForest=null){
# 	risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #
# 	cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #
# 	
# 	#
# 	sameSample=intersect(row.names(cli),row.names(risk))
# 	risk=risk[sameSample,]
# 	cli=cli[sameSample,]
# 	rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
# 	
# 	#
# 	uniTab=data.frame()
# 	for(i in colnames(rt[,3:ncol(rt)])){
# 		 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
# 		 coxSummary = summary(cox)
# 		 uniTab=rbind(uniTab,
# 		              cbind(id=i,
# 		              HR=coxSummary$conf.int[,"exp(coef)"],
# 		              HR.95L=coxSummary$conf.int[,"lower .95"],
# 		              HR.95H=coxSummary$conf.int[,"upper .95"],
# 		              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
# 		              )
# 	}
# 	write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
# 	bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="green")
# 
# 	#
# 	uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
# 	rt1=rt[,c("futime", "fustat", as.vector(uniTab[,"id"]))]
# 	multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
# 	multiCoxSum=summary(multiCox)
# 	multiTab=data.frame()
# 	multiTab=cbind(
# 	             HR=multiCoxSum$conf.int[,"exp(coef)"],
# 	             HR.95L=multiCoxSum$conf.int[,"lower .95"],
# 	             HR.95H=multiCoxSum$conf.int[,"upper .95"],
# 	             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
# 	multiTab=cbind(id=row.names(multiTab),multiTab)
# 	write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
# 	bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="red")
# }

indep2=function(rt,uniOutFile=null,multiOutFile=null,uniForest=null,multiForest=null){
  # risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #
  # cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #
  # #????????
  # sameSample=intersect(row.names(cli),row.names(risk))
  # risk=risk[sameSample,]
  # cli=cli[sameSample,]
  # rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
  
  #
  uniTab=data.frame()
  for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    uniTab=rbind(uniTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
  write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="green")
  
  #
  uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
  rt1=rt[,c("futime", "fustat", as.vector(uniTab[,"id"]))]
  multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
  multiCoxSum=summary(multiCox)
  multiTab=data.frame()
  multiTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  multiTab=cbind(id=row.names(multiTab),multiTab)
  write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="red")
}






# indep2(rt = rt6,
#       uniOutFile="uniCox.txt",
#       multiOutFile="multiCox.txt",
#       uniForest="uniForest.pdf",
#       multiForest="multiForest.pdf")


indep2(rt = cliRiskMergeNOid,
       uniOutFile="uniCox.txt",
       multiOutFile="multiCox.txt",
       uniForest="uniForest.pdf",
       multiForest="multiForest.pdf")


# ----- -------------------------------------------------------------------



# utf-8----
#绘制risk打分的ROC曲线00--------------------------------------------------
library(survivalROC)         #引用包

rt_cliRiskMergeNOid = cliRiskMergeNOid

rocCol=rainbow(ncol(rt_cliRiskMergeNOid)-2)
aucText=c()
pdf(file="cliROC.pdf", width=6, height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
# roc=survivalROC(Stime=risk$futime, status=risk$fustat, marker=risk$riskScore, predict.time=1, method="KM")
roc=survivalROC(Stime=rt_cliRiskMergeNOid$futime, status=rt_cliRiskMergeNOid$fustat, marker = rt_cliRiskMergeNOid$age, predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("age"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)

#绘制其他临床性状的ROC曲线
j=1
for(i in colnames(rt_cliRiskMergeNOid[,4:ncol(rt_cliRiskMergeNOid)])){
  roc=survivalROC(Stime=rt_cliRiskMergeNOid$futime, status=rt_cliRiskMergeNOid$fustat, marker = rt_cliRiskMergeNOid[,i], predict.time =1, method="KM")
  j=j+1
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
  aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
}
legend("bottomright", aucText, lwd=2, bty="n", col=rocCol)
dev.off()



#绘制risk打分的ROC曲线01--------------------------------------------------



# nomogram00 --------------------------------------------------------------
dput(colnames(rt1_filter))
# c("id", "age", "gender", "grade", "stage", "stage_T", "stage_M",
#   "stage_N")

rt4 = rt1[,c("id", "age", "gender", "tumor_grade", "pathologic_stage", "stage_T", "stage_M",
             "stage_N")]
rt4riskScoreClinic = merge(rt4, rt3,by = "id")
rt4riskScoreClinic$sampleType <- NULL
rt4riskScoreClinic$group<- NULL



library(rms)
library(foreign)
library(survival)
library(regplot)
library(mstate)

# rt<-read.table("cox_input.txt",header=T,sep="\t")

rt = rt4riskScoreClinic

# rt$age <- factor(rt$age,labels=c())
rt$gender = ifelse(rt$gender == "0", "Female", "Male" )


rt$gender <- factor(rt$gender,labels=c("male", "female"))
# rt$child_pugh_classification_grade <- factor(rt$child_pugh_classification_grade,labels=c("A", "B", "C", "N/A"))
# rt$history_of_neoadjuvant_treatment <- factor(rt$history_of_neoadjuvant_treatment,labels=c("No", "Yes"))
# rt$grade <- factor(rt$grade,labels=c("G1", "G2", "G3", "G4", "N/A"))
# rt$T <- factor(rt$T,labels=c("T1", "T2", "T3", "T4", "TX"))
# rt$N <- factor(rt$N,labels=c("N0", "N1", "NX"))
# rt$M <- factor(rt$M,labels=c("M0", "M1", "MX"))
# rt$stage <- factor(rt$stage,labels=c("I", "II", "III", "IV", "N/A"))
# rt$risk_score <- factor(rt$risk_score,labels=c("high", "low"))

rt$id<- NULL
rt$tumor_grade = factor(rt$tumor_grade,labels=c("G1","G2","G3","GX"))
table(rt$tumor_grade)
# G1 G2 G3 GX 
# 13 56 35 24 

rt$pathologic_stage = factor(rt$pathologic_stage,labels=c("1","2","3","4"))

rt$stage_T =   gsub(pattern = "a", replacement = "", rt$stage_T)

table(rt$stage_T)

rt$stage_T = factor(rt$stage_T,labels=c("T1","T2","T3","T4" ))

table(rt$stage_M)

rt$stage_M =   gsub(pattern = "a", replacement = "", rt$stage_M)
rt$stage_M =   gsub(pattern = "X", replacement = "1", rt$stage_M)
table(rt$stage_M)
rt$stage_M = factor(rt$stage_M,labels=c("M0","M1 " ))

table(rt$stage_N)

# rt$stage_N =   gsub(pattern = "N", replacement = "", rt$stage_N)
rt$stage_N =   gsub(pattern = "X", replacement = "0", rt$stage_N)
rt$stage_N<-factor(rt$stage_N,labels=c("N0","N1","N2","N3"))



ddist <- datadist(rt)
options(datadist='ddist')

#传统Nomogram
# cox <- cph(Surv(futime,fustat) ~ child_pugh_classification_grade+T+M+risk_score,surv=T,x=T, y=T,data=rt) 

paste(colnames(rt),collapse = "+")
# "age+gender+tumor_grade+pathologic_stage+stage_T+stage_M+stage_N+futime+fustat+m6Ascore"


cox <- cph(Surv(futime,fustat) ~ age+gender+tumor_grade+pathologic_stage+stage_T+
             stage_M+stage_N+m6Ascore,surv=T,x=T, y=T,data=rt) 

surv <- Survival(cox)
# sur_1_year<-function(x)surv(1*365,lp=x)
# sur_3_year<-function(x)surv(1*365*3,lp=x)
# sur_5_year<-function(x)surv(1*365*5,lp=x)

sur_1_year<-function(x)surv(1,lp=x)
sur_3_year<-function(x)surv(3,lp=x)
sur_5_year<-function(x)surv(5,lp=x)



nom_sur <- nomogram(cox,fun=list(sur_1_year,sur_3_year,sur_5_year),lp= F,funlabel=c('1-Year survival','3-Year survival','5-Year survival'),maxscale=100,fun.at= c('1.0','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0'))
pdf("nomogram.pdf")
plot(nom_sur,xfrac=0.4)
dev.off()


nom_sur <- nomogram(cox,fun=list(sur_1_year,sur_3_year,sur_5_year),lp= F,funlabel=c('1-Year survival','3-Year survival','5-Year survival'),maxscale=100,fun.at=c('0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1'))


pdf("nom.pdf",15,10)
plot(nom_sur,xfrac=0.25)
dev.off()


pdf("nom3.pdf",10,8)
plot(nom_sur,xfrac=0.25)
dev.off()



# nomogram01 --------------------------------------------------------------




# end ---------------------------------------------------------------------
# end ---------------------------------------------------------------------
# end ---------------------------------------------------------------------
# end ---------------------------------------------------------------------
# end ---------------------------------------------------------------------





#SurvivalROC for time to event data---------------------


#Example1: Validate the reported prognosis model in your own patients. The model produces high points. You want to see whether higher point is correlated with longer or shorter survival using your survival data. 
#Example2: You can also validate your own survival nomogram as well
library(survivalROC)

f2 = cox


# total_point<-f$linear.predictors
total_point<-f2$linear.predictors

# data_surv2$total_point<-f2$linear.predictors


# data_surv3= data_surv

# data_surv = data_surv2

data_surv = rt4riskScoreClinic[,c("futime","fustat")]
colnames(data_surv) = c("time","status")


summary(data_surv$time)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.120   1.863   3.983   3.823   4.997  10.052

data_surv$time = data_surv$time*12
summary(data_surv$time)

data_surv$total_point<-f2$linear.predictors



# survivalROC(data_surv$time,data_surv$status,data_surv$total_point,predict.time=60,method ="KM")



# test,start --------------------------------------------------------------------

dir.create("clinicalRiskScore_ROC_timeDependent")

# Stime=rt$futime

# cutoffTime<- c(1*12,2*12,3*12,4*12,5*12,6*12,7*12)

cutoffTime<- c(1*12,3*12,5*12)

# cutoffTime<- as.vector(seq(1:7))

roc_cutoffTime_km <- function(cutoffTime) {
  roc=survivalROC(data_surv$time,data_surv$status,data_surv$total_point,
                  predict.time =cutoffTime, method="KM")
  par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
  pdf(file=paste("clinicalRiskScore_ROC_timeDependent/","0102AA_ROC_",cutoffTime,".pdf",sep = ""))
  plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
       xlab="False positive rate", ylab="True positive rate",
       main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"," \n Year=",cutoffTime),
       lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  abline(0,1)
  text(0.89, 0.14, paste("AUC = ", format(roc$AUC, digits = 3)))
  dev.off()
}

sapply(cutoffTime,roc_cutoffTime_km)



# cutoffTime<- c(1*12,2*12,3*12,4*12,5*12,6*12)
roc_cutoffTime_KM <- function(cutoffTime) {
  roc=survivalROC(data_surv$time,data_surv$status,data_surv$total_point,
                  predict.time =cutoffTime, method="KM")
  par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
  pdf(file=paste("clinicalRiskScore_ROC_timeDependent/","0103AA_ROC_",cutoffTime,".pdf",sep = ""))
  plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
       xlab="False positive rate", ylab="True positive rate",
       main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")","Year=",cutoffTime),
       lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  abline(0,1)
  text(0.89, 0.14, paste("AUC = ", format(roc$AUC, digits = 3)))
  dev.off()
}
sapply(cutoffTime,roc_cutoffTime_KM)





# test,end --------------------------------------------------------------------
# cutoffTime <-seq(1,20,1),start-----------------------------------------------
# cutoffTime <-seq(1,6,1)

# cutoffTime<- c(1*12,2*12,3*12,4*12,5*12,6*12)


# cutoffTime<- c(1,2,3)
for (i in c(1:length(cutoffTime))){ 
  roc=survivalROC(data_surv$time,data_surv$status,data_surv$total_point,
                  predict.time =cutoffTime[i], method="KM")
  par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
  pdf(file=paste("clinicalRiskScore_ROC_timeDependent/","0101AB_ROC_",cutoffTime[i],".pdf",sep = ""))
  plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
       xlab="False positive rate(1-Specificity)", ylab="True positive rate(Sensitivity)",
       main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"," \n Year=",cutoffTime[i]),
       lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  abline(0,1)
  text(0.89, 0.14, paste("AUC = ", format(roc$AUC, digits = 3)))
  dev.off()
}



# cutoffTime <-seq(1,20,1),end-----------------------------------------------





# test,start-------------------------
# cutoffTime_One2Ten<-seq(1,7,1)
# cutoffTime_One2Ten<- c(1*12,2*12,3*12,4*12,5*12,6*12,7*12)
cutoffTime_One2Ten<- c(1*12,2*12,3*12,4*12,5*12,6*12,7*12)

sroc <- lapply(cutoffTime_One2Ten, function(x) {
  stroc=survivalROC(Stime=data_surv$time, status=data_surv$status, marker = data_surv$total_point, 
                    predict.time =x, method="KM", 
                    span = .25 * 350^(-.2))
  data.frame(TPF = stroc[["TP"]], FPF = stroc[["FP"]], 
             c = stroc[["cut.values"]], 
             time = rep(stroc[["predict.time"]], length(stroc[["FP"]])),
             AUC = rep(stroc[["AUC"]], length(stroc[["FP"]])) )
  # stroc$AUC
})
sroclong <- do.call(rbind, sroc)


auc_3<-unique(sroclong$AUC)
auc_3

# 0.8936692 0.8901178 0.8719516 0.8641216 0.8967544





sroc_AUC <- lapply(cutoffTime_One2Ten, function(x) {
  stroc=survivalROC(Stime=data_surv$time, status=data_surv$status, marker = data_surv$total_point, 
                    predict.time =x, method="KM"      ) # span = .25 * 350^(-.2))
  data.frame(TPF = stroc[["TP"]], FPF = stroc[["FP"]], 
             c = stroc[["cut.values"]], 
             time = rep(stroc[["predict.time"]], length(stroc[["FP"]])),
             AUC = rep(stroc[["AUC"]], length(stroc[["FP"]])) )
  # stroc$AUC
})



sroclong2 <- do.call(rbind, sroc_AUC)
auc_sroclong2 <-unique(sroclong2$AUC)
auc_sroclong2 




# ggplot(sroclong, aes(x = FPF, y = TPF, label = c, color = time)) + 
#   geom_roc(labels = FALSE, stat = "identity") + style_roc(theme = theme_gray)

# install.packages("gridSVG")

# library(devtools)
library(ggplot2)

# devtools::install_github("sachsmc/plotROC")


library(plotROC)
basicplot<-ggplot(sroclong, aes(x = FPF, y = TPF, label = c, color = time)) + 
  geom_roc(labels = FALSE, stat = "identity") + style_roc(theme = theme_gray)

ggplot(sroclong, aes(x = FPF, y = TPF, label = c, color = time)) + 
  geom_roc(labels = FALSE, stat = "identity") + style_roc(theme = theme_gray)

ggplot(sroclong, aes(x = FPF, y = TPF, label = c, color = time)) + 
  geom_roc(labels = TRUE, stat = "identity") + style_roc(theme = theme_bw())

ggplot(sroclong, aes(x = FPF, y = TPF, label = c, color = time)) + 
  geom_roc(labels = FALSE, stat = "identity") + style_roc()

ggplot(sroclong, aes(x = FPF, y = TPF, label = c, color = as.factor(time))) + 
  geom_roc(labels = FALSE, stat = "identity") + style_roc()

auc_3<-unique(sroclong$AUC)
auc_3

ggplot(sroclong, aes(x = FPF, y = TPF, label = c, color = as.factor(time))) + 
  geom_roc(labels = FALSE, stat = "identity") + style_roc()+ 
  annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot)$AUC, 3)))

text1=paste("AUC at 1 year : ",format(auc_sroclong2[1],digits=3),sep="")
text2=paste("AUC at 2 years : ",format(auc_sroclong2[2],digits=3),sep="")
text3=paste("AUC at 3 years : ",format(auc_sroclong2[3],digits=3),sep="")
text4=paste("AUC at 4 years : ",format(auc_sroclong2[4],digits=3),sep="")
text5=paste("AUC at 5 years : ",format(auc_sroclong2[5],digits=3),sep="")
text6=paste("AUC at 6 year : ",format(auc_sroclong2[6],digits=3),sep="")
text7=paste("AUC at 7 years : ",format(auc_sroclong2[7],digits=3),sep="")
# text8=paste("AUC at 8 years : ",format(auc_sroclong2[8],digits=3),sep="")
# text9=paste("AUC at 9 years : ",format(auc_sroclong2[9],digits=3),sep="")
# text10=paste("AUC at 10 years : ",format(auc_sroclong2[10],digits=3),sep="")


value<-0.6



# p<- ggplot(sroclong, aes(x = FPF, y = TPF, label = c, color = as.factor(time))) + 
#   geom_roc(labels = FALSE, stat = "identity") + style_roc(theme = theme_bw())
p<- ggplot(sroclong, aes(x = FPF, y = TPF, label = c, color = as.factor(time))) + 
  geom_roc(labels = FALSE, stat = "identity") + style_roc()
# p<- p + annotate("text", x=0.8, y = c(value*0.95,value*0.90,value*0.85,value*0.80,value*0.75,
#                                       value*0.70,value*0.65,value*0.60,value*0.55,value*0.50), 
#                  label = c(text1,text2,text3,text4,text5,text6,text7,text8,text9,text10),size=4.0,col="violetred1")

p<- p + annotate("text", x=0.8, y = c(value*0.95,value*0.90,value*0.85,value*0.80,value*0.75,
                                      value*0.70,value*0.65), 
                 label = c(text1,text2,text3,text4,text5,text6,text7),size=4.0,col="black")

p<- p + theme_set(theme_classic())+theme(panel.grid.major=element_line(colour=NA))
print(p)


# ggsave( file = "rocOne2Ten_ggplot_01.png", width = 8, height = 8, type = "cairo", dpi = 300) 
# ggsave( file = "rocOne2Ten_ggplot_02.pdf", width = 8, height = 8)
# 
# ggsave( file = "rocOne2Ten_ggplot_03.jpeg", width = 10, height = 10, units="in", dpi = 300)
# 
# ggsave( file = "rocOne2Ten_ggplot_04.pdf", width = 5, height = 5)
# ggsave( file = "rocOne2Ten_ggplot_05.pdf", width = 10, height = 10)
# 
# ggsave( file = "rocOne2Ten_ggplot_06.pdf")
# 
# ggsave( file = "rocOne2Ten_ggplot_07.pdf")





# 重新做(---------------------------------------------------------------------
# 重新做( ---------------------------------------------------------------------
# 重新做( ---------------------------------------------------------------------


rm(p)

dir.create("ROC_timeDependent")

sroc <- lapply(cutoffTime_One2Ten, function(x) {
  stroc=survivalROC(Stime=data_surv$time, status=data_surv$status, marker = data_surv$total_point, 
                    predict.time =x, method="NNE", 
                    span = .25 * 350^(-.2))
  data.frame(TPF = stroc[["TP"]], FPF = stroc[["FP"]], 
             c = stroc[["cut.values"]], 
             time = rep(stroc[["predict.time"]], length(stroc[["FP"]])),
             AUC = rep(stroc[["AUC"]], length(stroc[["FP"]])) )
  # stroc$AUC
})
sroclong <- do.call(rbind, sroc)
sroclong$AUC

auc_3<-unique(sroclong$AUC)
auc_3





sroc_AUC <- lapply(cutoffTime_One2Ten, function(x) {
  stroc=survivalROC(Stime=data_surv$time, status=data_surv$status, marker = data_surv$total_point, 
                    predict.time =x, method="KM"      ) # span = .25 * 350^(-.2))
  data.frame(TPF = stroc[["TP"]], FPF = stroc[["FP"]], 
             c = stroc[["cut.values"]], 
             time = rep(stroc[["predict.time"]], length(stroc[["FP"]])),
             AUC = rep(stroc[["AUC"]], length(stroc[["FP"]])) )
  # stroc$AUC
})


sroclong2 <- do.call(rbind, sroc_AUC)
auc_sroclong2 <-unique(sroclong2$AUC)
auc_sroclong2

# sroclong2 <- do.call(rbind, sroc_AUC)
# auc_sroclong2 <-sroclong2$AUC[sroclong2$time==cutoffTime_One2Ten]
# auc_sroclong2 
# length(auc_sroclong2)



# maxYear <- 7
maxYear <- floor(summary(data_surv$time)[6])
maxYear <- as.numeric(maxYear)
maxYear


maxYear <- 6
maxYear

maxYear <- 6


length(auc_sroclong2)/maxYear
# [1] 43.28571


auc_sroclong2<-auc_sroclong2[seq(1,length(auc_sroclong2),length(auc_sroclong2)/maxYear*12)]  #16要修改
auc_sroclong2
# [1] 0.8936692 0.8936692 0.8901178 0.8719516 0.8641216 0.8967544 0.8967544

# ------- -----------------------------------------------------------------
rm(p)

text1=paste("AUC at 1 year : ",format(auc_sroclong2[1],digits=3),sep="")
text2=paste("AUC at 2 years : ",format(auc_sroclong2[2],digits=3),sep="")
text3=paste("AUC at 3 years : ",format(auc_sroclong2[3],digits=3),sep="")
text4=paste("AUC at 4 years : ",format(auc_sroclong2[4],digits=3),sep="")
text5=paste("AUC at 5 years : ",format(auc_sroclong2[5],digits=3),sep="")
# text6=paste("AUC at 6 year : ",format(auc_sroclong2[6],digits=3),sep="")
# text7=paste("AUC at 7 years : ",format(auc_sroclong2[7],digits=3),sep="")
# text8=paste("AUC at 8 years : ",format(auc_sroclong2[8],digits=3),sep="")
# text9=paste("AUC at 9 years : ",format(auc_sroclong2[9],digits=3),sep="")
# text10=paste("AUC at 10 years : ",format(auc_sroclong2[10],digits=3),sep="")


value<-0.6


# p<- ggplot(sroclong2, aes(x = FPF, y = TPF, label = c, color = as.factor(time))) +
#   geom_roc(labels = FALSE, stat = "identity") + style_roc(theme = theme_bw())
p<- ggplot(sroclong, aes(x = FPF, y = TPF, label = c, color = as.factor(time))) + 
  geom_roc(labels = FALSE, stat = "identity") + style_roc()
# p<- p + annotate("text", x=0.8, y = c(value*0.95,value*0.90,value*0.85,value*0.80,value*0.75,
#                                       value*0.70,value*0.65,value*0.60,value*0.55,value*0.50), 
#                  label = c(text1,text2,text3,text4,text5,text6,text7,text8,text9,text10),size=4.0,col="violetred1")

p<- p + annotate("text", x=0.9, y = c(value*0.95,value*0.90,value*0.85,value*0.80,value*0.75,
                                      value*0.70,value*0.65), 
                 label = c(text1,text2,text3,text4,text5,text6,text7)
                 ,size=4.0,col="black")

# p<- p + theme_set(theme_classic())+theme(panel.grid.major=element_line(colour=NA))
# p<- p + theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA))
# P <- p+geom_abline(intercept = 0,slope = 1,linetype = "dashed")
# p<- p +theme(panel.grid.major=element_line(colour=NA))


p <-p + geom_abline(intercept = 0, slope = 1, linetype = "dashed")+theme_classic() 

p <-p +theme(legend.position = c(0.8,0.2))

# p <-p +geom_path(aes(color= AUC))
print(p)

dir.create("clinicRiskScore_rocOne2Ten")
ggsave( file = "clinicRiskScore_rocOne2Ten/rocOne3Ten_ggplot_01.png", width = 8, height = 8, type = "cairo", dpi = 300) 
ggsave( file = "clinicRiskScore_rocOne2Ten/rocOne3Ten_ggplot_02.pdf", width = 8, height = 8)

ggsave( file = "clinicRiskScore_rocOne2Ten/rocOne3Ten_ggplot_03.jpeg", width = 10, height = 10, units="in", dpi = 300)

ggsave( file = "clinicRiskScore_rocOne2Ten/rocOne3Ten_ggplot_04.pdf", width = 5, height = 5)
ggsave( file = "clinicRiskScore_rocOne2Ten/rocOne3Ten_ggplot_05.pdf", width = 10, height = 10)
ggsave( file = "clinicRiskScore_rocOne2Ten/rocOne3Ten_ggplot_06.pdf", width = 8, height = 8)


# 重新做)---------------------------------------------------------------------
# 重新做) ---------------------------------------------------------------------
# 重新做) ---------------------------------------------------------------------



# 后期添加add-start -----------------------------------------------------------

fp <- predict(f)
# summary(f,data=data_surv)
library(Hmisc)
options(scipen=200)
with(data_surv,rcorr.cens(fp,Surv(time, status)  ))

#     C Index            Dxy           S.D.              n        missing     uncensored Relevant Pairs 
#    0.16474385    -0.67051229     0.03879526   302.00000000     0.00000000    65.00000000 29124.00000000 
#    Concordant      Uncertain 
# 4798.00000000 61778.00000000 

# 后期添加add-end -----------------------------------------------------------




# ROC_nomogram_totol_clinicalRisk00 --------------------------------------------------
library(timeROC)
#定义颜色
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)

risk = data_surv
colnames(risk) = c("futime","fustat","riskScore") 
risk$futime = risk$futime/12

######绘制1 2 3年的ROC曲线######
ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
               marker=risk$riskScore,cause=1,
               weighting='aalen',
               times=c(1,2,3),ROC=TRUE)
pdf(file="ROC_nomogram_totol_clinicalRisk.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=2,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=bioCol[1:3], lwd=2, bty = 'n')
dev.off()


# ROC_nomogram_totol_clinicalRisk01 --------------------------------------------------


# end ---------------------------------------------------------------------
# end ---------------------------------------------------------------------
# end ---------------------------------------------------------------------
# end ---------------------------------------------------------------------
# end ---------------------------------------------------------------------












############??????????????############
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  #????????????
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrLow[hrLow<0.001]=0.001
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #????????
  pdf(file=forestFile, width=6.5, height=4.8)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #????????????????????????
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
  
  #??????????
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  LOGindex=2 
  hrLow = log(as.numeric(hrLow),LOGindex)
  hrHigh = log(as.numeric(hrHigh),LOGindex)
  hr = log(as.numeric(hr),LOGindex)
  xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=2)
  a1 = axis(1,labels=F,tick=F)
  axis(1,a1,LOGindex^a1)
  dev.off()
}
############??????????????############

indep3=function(riskFile=null,uniOutFile=null,multiOutFile=null,uniForest=null,multiForest=null){
  rt=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #????????????
  # cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #????????????
  # #????????
  # sameSample=intersect(row.names(cli),row.names(risk))
  # risk=risk[sameSample,]
  # cli=cli[sameSample,]
  # rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
  
  #??????????????????????
  uniTab=data.frame()
  for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    uniTab=rbind(uniTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
  write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="green")
  
  #??????????????????????
  uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
  rt1=rt[,c("futime", "fustat", as.vector(uniTab[,"id"]))]
  multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
  multiCoxSum=summary(multiCox)
  multiTab=data.frame()
  multiTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  multiTab=cbind(id=row.names(multiTab),multiTab)
  write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="red")

}




indep3(riskFile="input.txt",
       uniOutFile="uniCox.txt",
       multiOutFile="multiCox.txt",
       uniForest="uniForest.pdf",
       multiForest="multiForest.pdf")











# geo??????????00.00 -----------------------------------------------------------


rt2 = openxlsx::read.xlsx(xlsxFile = "geoClinical_GSE133057.xlsx",sheet = 3)
rt2 = rt2 %>% filter(sample %in% rt3$id) %>% 
  dplyr::select("sample","age","gender","stage")
rt2 = rt2 %>% dplyr::rename( id = sample) 
table(rt2$gender)
rt2$gender = ifelse(tolower(rt2$gender) == "female", "0", "1" )
table(rt2$gender)
# geo??????????00.01 -----------------------------------------------------------



rt4 = rbind(rt1,rt2)

rt5 = rt3 %>% inner_join(rt4,by = "id") %>% dplyr::select(-sampleType,-group) %>% as.data.frame()
# class(rt5)
table(rt5$gender)

# rt5$gender = ifelse(tolower(rt5$gender) == "female", "0", "1" )
# rt5$gender = ifelse(rt5$gender == "female"|"FEMALE", "0", "1" )

# table(rt5$gender)

table(rt5$stage)


# rt5$stage =   str_match(rt5$stage, "Stage ([^ABC]*)")[,2]


#install.packages('survival')

rt5$stage <- as.numeric(rt5$stage)
rt5$gender <- as.numeric(rt5$gender )

rt6 = rt5 %>% tibble::column_to_rownames(var = "id") %>% 
  mutate(riskScore = m6Ascore)  %>% select(-m6Ascore)


head(rt6)

write.table(rt6,file = "input.txt",sep="\t",row.names=T,col.names = NA)

sessionInfo()
# China.utf8
