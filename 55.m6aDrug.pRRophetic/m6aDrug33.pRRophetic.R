
rm(list = ls())

setwd("C:\\Users\\sun\\Desktop\\m6a-esca\\55.m6aDrug.pRRophetic")




#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("car", "ridge", "preprocessCore", "genefilter", "sva"))
#BiocManager::install(c("pRRophetic"))

#install.packages("ggpubr")

#install.packages("C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\55.m6aDrug.pRRophetic", repos = NULL, type = "source")

#???ð?

library(ridge)

library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
set.seed(12345)



#ҩ???б?
# expFile="symbol.txt"         #?????????ļ?
# riskFile="risk.all.txt"      #?????????ļ?

# ?õ???ԭʼ???ݣ?????ɸѡ??ѡ??????????????-------
expFile="merge.txt"         #?????????ļ?


riskFile="m6Ascore.group.txt"      #?????????ļ?


# setwd("C:\\biowolf\\m6aDrug\\33.pRRophetic")     #???ù???Ŀ¼



allDrugs=c("A.443654", "A.770041", "ABT.263", "ABT.888", "AG.014699", "AICAR", "AKT.inhibitor.VIII", "AMG.706", 
           "AP.24534", "AS601245", "ATRA", "AUY922", "Axitinib", "AZ628", "AZD.0530", "AZD2281", 
           "AZD6244", "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", "Bexarotene", "BI.2536", "BIBW2992", 
           "Bicalutamide", "BI.D1870", "BIRB.0796", "Bleomycin", "BMS.509744", "BMS.536924", 
           "BMS.708163", "BMS.754807", "Bortezomib", "Bosutinib", "Bryostatin.1", "BX.795", 
           "Camptothecin", "CCT007093", "CCT018159", "CEP.701", "CGP.082996", "CGP.60474", 
           "CHIR.99021", "CI.1040", "Cisplatin", "CMK", "Cyclopamine", "Cytarabine", 
           "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin", "EHT.1864", "Elesclomol", "Embelin", "Epothilone.B", 
           "Erlotinib", "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib", "Gemcitabine", 
           "GNF.2", "GSK269962A", "GSK.650394", "GW.441756", "GW843682X", "Imatinib", "IPA.3", "JNJ.26854165", 
           "JNK.9L", "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", "KU.55933", "Lapatinib", "Lenalidomide", 
           "LFM.A13", "Metformin", "Methotrexate", "MG.132", "Midostaurin", "Mitomycin.C", "MK.2206", "MS.275", 
           "Nilotinib", "NSC.87877", "NU.7441", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", 
           "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide", "Pazopanib", "PD.0325901", "PD.0332991", "PD.173074", 
           "PF.02341066", "PF.4708671", "PF.562271", "PHA.665752", "PLX4720", "Pyrimethamine", "QS11", "Rapamycin", 
           "RDEA119", "RO.3306", "Roscovitine", "Salubrinal", "SB.216763", "SB590885", "Shikonin", "SL.0101.1", 
           "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", 
           "Vinblastine", "Vinorelbine", "Vorinostat", "VX680", "VX.702", "WH.4.023", "WO2009093972", 
           "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", "Z.LLNle.CHO", "ZM.447439")


# A.443654, A.770041, ABT.263, ABT.888, AG.014699, AICAR, AKT.inhibitor.VIII, AMG.706, AP.24534, AS601245,
# ATRA, AUY922, Axitinib, AZ628, AZD.0530, AZD.2281, AZD6244, AZD6482, AZD7762, AZD8055, BAY.61.3606, Bexarotene,
# BI.2536, BIBW2992, Bicalutamide, BI.D1870, BIRB.0796, Bleomycin, BMS.509744, BMS.536924, BMS.708163, BMS.754807,
# Bortezomib, Bosutinib, Bryostatin.1, BX.795, Camptothecin, CCT007093, CCT018159, CEP.701, CGP.082996, CGP.60474,
# CHIR.99021, CI.1040, Cisplatin, CMK, Cyclopamine, Cytarabine, Dasatinib, DMOG, Docetaxel, Doxorubicin, EHT.1864,
# Elesclomol, Embelin, Epothilone.B, Erlotinib, Etoposide, FH535, FTI.277, GDC.0449, GDC0941, Gefitinib,
# Gemcitabine, GNF.2, GSK269962A, GSK.650394, GW.441756, GW843682X, Imatinib, IPA.3, JNJ.26854165, JNK.9L,
# JNK.Inhibitor.VIII, JW.7.52.1, KIN001.135, KU.55933, Lapatinib, Lenalidomide, LFM.A13,
# Metformin, Methotrexate, MG.132, Midostaur



# possibleDrugs2016 <- c("Erlotinib", "Rapamycin", "Sunitinib", "PHA-665752", "MG-132", "Paclitaxel", 
#                        "Cyclopamine", "AZ628", "Sorafenib", "VX-680", "Imatinib", "TAE684", "Crizotinib", 
#                        "Saracatinib", "S-Trityl-L-cysteine", "Z-LLNle-CHO", "Dasatinib", "GNF-2", "CGP-60474", 
#                        "CGP-082996", "A-770041", "WH-4-023", "WZ-1-84", "BI-2536", "BMS-536924", "BMS-509744", 
#                        "CMK", "Pyrimethamine", "JW-7-52-1", "A-443654", "GW843682X", "MS-275", "Parthenolide", 
#                        "KIN001-135", "TGX221", "Bortezomib", "XMD8-85", "Roscovitine", "Salubrinal", "Lapatinib", 
#                        "GSK269962A", "Doxorubicin", "Etoposide", "Gemcitabine", "Mitomycin C", "Vinorelbine", 
#                        "NSC-87877", "Bicalutamide", "QS11", "CP466722", "Midostaurin", "CHIR-99021", "AP-24534", 
#                        "AZD6482", "JNK-9L", "PF-562271", "HG-6-64-1", "JQ1", "JQ12", "DMOG", "FTI-277", 
#                        "OSU-03012", "Shikonin", "AKT inhibitor VIII", "Embelin", "FH535", "PAC-1", "IPA-3", 
#                        "GSK-650394", "BAY 61-3606", "5-Fluorouracil", "Thapsigargin", "Obatoclax Mesylate", 
#                        "BMS-754807", "Lisitinib", "Bexarotene", "Bleomycin", "LFM-A13", "GW-2580", "AUY922", 
#                        "Phenformin", "Bryostatin 1", "Pazopanib", "LAQ824", "Epothilone B", "GSK1904529A", 
#                        "BMS345541", "Tipifarnib", "BMS-708163", "Ruxolitinib", "AS601245", "Ispinesib Mesylate", 
#                        "TL-2-105", "AT-7519", "TAK-715", "BX-912", "ZSTK474", "AS605240", "Genentech Cpd 10", 
#                        "GSK1070916", "KIN001-102", "LY317615", "GSK429286A", "FMK", "QL-XII-47", "CAL-101", 
#                        "UNC0638", "XL-184", "WZ3105", "XMD14-99", "AC220", "CP724714", "JW-7-24-1", 
#                        "NPK76-II-72-1", "STF-62247", "NG-25", "TL-1-85", "VX-11e", "FR-180204", "Tubastatin A", 
#                        "Zibotentan", "YM155", "NSC-207895", "VNLG/124", "AR-42", "CUDC-101", 
#                        "Belinostat", "I-BET-762", "CAY10603", "Linifanib ", "BIX02189", "CH5424802", 
#                        "EKB-569", "GSK2126458", "KIN001-236", "KIN001-244", "KIN001-055", "KIN001-260", 
#                        "KIN001-266", "Masitinib", "MP470", "MPS-1-IN-1", "BHG712", "OSI-930", 
#                        "OSI-027", "CX-5461", "PHA-793887", "PI-103", "PIK-93", "SB52334", "TPCA-1", 
#                        "TG101348", "Foretinib", "Y-39983", "YM201636", "Tivozanib", "GSK690693", 
#                        "SNX-2112", "QL-XI-92", "XMD13-2", "QL-X-138", "XMD15-27", "T0901317", "EX-527", 
#                        "THZ-2-49", "KIN001-270", "THZ-2-102-1", "AICAR", "Camptothecin", "Vinblastine", 
#                        "Cisplatin", "Cytarabine", "Docetaxel", "Methotrexate", "ATRA", "Gefitinib", "Navitoclax", 
#                        "Vorinostat", "Nilotinib", "RDEA119", "CI-1040", "Temsirolimus", "Olaparib", "Veliparib", 
#                        "Bosutinib", "Lenalidomide", "Axitinib", "AZD7762", "GW 441756", "CEP-701", "SB 216763", 
#                        "17-AAG", "VX-702", "AMG-706", "KU-55933", "Elesclomol", "Afatinib", "GDC0449", "PLX4720", 
#                        "BX-795", "NU-7441", "SL 0101-1", "BIRB 0796", "JNK Inhibitor VIII", "681640", "Nutlin-3a (-)", 
#                        "PD-173074", "ZM-447439", "RO-3306", "MK-2206", "PD-0332991", "BEZ235", "GDC0941", "AZD8055", 
#                        "PD-0325901", "SB590885", "selumetinib", "CCT007093", "EHT 1864", "Cetuximab", "PF-4708671", 
#                        "JNJ-26854165", "HG-5-113-01", "HG-5-88-01", "TW 37", "XMD11-85h", "ZG-10", "XMD8-92", "QL-VIII-58", 
#                        "CCT018159", "AG-014699", "SB 505124", "Tamoxifen", "QL-XII-61", "PFI-1", "IOX2", "YK 4-279", 
#                        "(5Z)-7-Oxozeaenol", "piperlongumine", "FK866", "Talazoparib", "rTRAIL", "UNC1215", "SGC0946", "XAV939", 
#                        "Trametinib", "Dabrafenib", "Temozolomide", "Bleomycin (50 uM)", "SN-38", "MLN4924")

#ҩ??????
GCP.drug <- read.table("drug.txt") #????Ҫ???ĵ?��??ҩ??ͻ???drug_eg.txt
GCP.drug <- GCP.drug$V1

allDrugs = GCP.drug





# ------------------- -----------------------------------------------------




#??ȡ?????????ļ?,???????ݽ??д???
rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
# data=data[rowMeans(data)>0.5,]
data=data[rowMeans(data)>0,]



#ɾ????????Ʒ
# group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
# group=sapply(strsplit(group,""), "[", 1)
# group=gsub("2","1",group)
# data=data[,group==0]
data=t(data)
# rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))

rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

data=avereps(data)
data=t(data)







#??ȡ?????????ļ?
# ?????
riskRT=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(riskRT)[colnames(riskRT) == "group"] <- "risk" 

data=data[,rownames(riskRT)]



for(drug in allDrugs){
  
  
  # drug=c("A.443654")
  
	#Ԥ??ҩ?
	senstivity=pRRopheticPredict(data, drug, selection=1)
	senstivity=senstivity[senstivity!="NaN"]
	#senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)
	
	#?????ļ???ҩ???????Խ????ϲ?
	sameSample=intersect(row.names(riskRT), names(senstivity))
	# risk=riskRT[sameSample, "risk",drop=F]
	risk=riskRT[sameSample, "risk",drop=F]
	senstivity=senstivity[sameSample]
	rt=cbind(risk, senstivity)
	
	#???ñȽ???
	# rt$risk=factor(rt$risk, levels=c("low", "high"))
	rt$risk=factor(rt$risk, levels=c("Low", "High"))
	type=levels(factor(rt[,"risk"]))
	
	
	comp=combn(type, 2)
	my_comparisons=list()
	for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	
	#??ȡ?ߵͷ?????????pvalue
	test=wilcox.test(senstivity~risk, data=rt)
	
	if(test$p.value<0.05){
		#????????ͼ
		boxplot=ggboxplot(rt, x="risk", y="senstivity", fill="risk",
					      xlab="Risk",
					      ylab=paste0(drug, " senstivity (IC50)"),
					      legend.title="Risk",
					      palette=c("#0066FF","#FF0000")
					     )+ 
			stat_compare_means(comparisons=my_comparisons)
		pdf(file=paste0("durgSenstivity.", drug, ".pdf"), width=5, height=4.5)
		print(boxplot)
		dev.off()
	}
}







# ???????޸?drugs???ٴ?????00 -------------------------------------------------------
allDrugs=c("Rapamycin", 
           "RDEA119", "RO.3306", "Roscovitine", "Salubrinal", "SB.216763", "SB590885", "Shikonin", "SL.0101.1", 
           "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", 
           "Vinblastine", "Vinorelbine", "Vorinostat", "VX.702", "WH.4.023", "WO2009093972", 
           "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", "Z.LLNle.CHO", "ZM.447439")




# data, drug

drug <-NULL

Gcp138drugs <- function(drug) {
  cat(drug," start\n") # ??ʾ??ǰҩ?￪ʼ
  senstivity=pRRopheticPredict(data, drug, selection=1)
  
  if(!all(names(senstivity)==rownames(risk))) {stop("Name mismatched!\n")} # ?????ֲ?ƥ???򱨴??˳?
  
  senstivity=senstivity[senstivity!="NaN"]
  #senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)
  
  #?????ļ???ҩ???????Խ????ϲ?
  sameSample=intersect(row.names(riskRT), names(senstivity))
  # risk=riskRT[sameSample, "risk",drop=F]
  risk=riskRT[sameSample, "risk",drop=F]
  senstivity=senstivity[sameSample]
  rt=cbind(risk, senstivity)
  
  #???ñȽ???
  # rt$risk=factor(rt$risk, levels=c("low", "high"))
  rt$risk=factor(rt$risk, levels=c("Low", "High"))
  type=levels(factor(rt[,"risk"]))
  
  
  comp=combn(type, 2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  #??ȡ?ߵͷ?????????pvalue
  test=wilcox.test(senstivity~risk, data=rt)
  print(test$p.value)
  
  if(test$p.value<0.05){
    #????????ͼ
    boxplot=ggboxplot(rt, x="risk", y="senstivity", fill="risk",
                      xlab="Risk",
                      ylab=paste0(drug, " senstivity (IC50)"),
                      legend.title="Risk",
                      palette=c("#0066FF","#FF0000")
    )+ stat_compare_means(comparisons=my_comparisons)
    # pdf(file=paste0("durgSenstivity.", drug, ".pdf"), width=5, height=4.5)
    # print(boxplot)
    # dev.off()
    
    plotp[[drug]] <- boxplot # ???????б??﹩?ϲ?ͼƬ??
  }
  cat(drug," has been finished!\n") # ??ʾ??ǰҩ???ѷ???????
}


plotp <- list()
# Gcp138drugs ?ǹ??ܰ?
plotp = lapply(allDrugs, Gcp138drugs)
# plotp ???ǿ?ֵ

# plotp ??????????????
# ?ʺ?չʾ??????
#原代码p2 <- plot_grid(plotlist=plotp, ncol=6)
p2 <- plot(plotlist=plotp, ncol=6)
ggsave("boxplot of predicted IC50_multiple.pdf", width = 12, height = 36)



