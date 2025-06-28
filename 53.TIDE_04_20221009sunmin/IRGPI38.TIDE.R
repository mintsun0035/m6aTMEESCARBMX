#gb2312



# setwd("C:\\Users\\sun\\Desktop\\m6a-esca\\53.TIDE_02")



rm(list = ls())
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#install.packages("ggpubr")
#???ð?
library(limma)
library(ggpubr)
# tideFile="TIDE.txt"          #TIDE?ļ?
tideFile="TIDE_output_ESCA.csv"          #TIDE?ļ?



# riskFile="risk.TCGA.txt"     #?????ļ?
# riskFile="m6Ascore.group_PCA.txt"     #?????ļ?

riskFile="m6Ascore.group.txt"     #?????ļ?


#??ȡTIDE????
# tide=read.table(tideFile, header=T, sep="\t", check.names=F, row.names=1)

tide=read.table(tideFile, header=T,sep=",",row.names=1)

# group=sapply(strsplit(row.names(tide),"\\-"), "[", 4)
# group=sapply(strsplit(group,""), "[", 1)
# group=gsub("2", "1", group)
# tide=tide[group==0,]
# row.names(tide)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))
# tide=avereps(tide)
row.names(tide)=gsub("(.*?)\\_(.*?)", "\\2", row.names(tide))



#??ȡ?????????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
#?ϲ?????    row.names??ʾ?е?????
sameSample=intersect(row.names(tide), row.names(risk))
# intersectR?????еĺ??????ڲ???��???????Ľ??????˺?????��??????(?? Vectors??dataframes ??)??Ϊ???????????ɾ???��???????Ĺ??????ݵĵ???????????
tide=tide[sameSample, , drop=F]
# drop=T?ǽ?ά????˼  ʹ?ú???colnames()??rownames()?ֱ??????????к???????
colnames(risk)[colnames(risk)=="group"] = "Risk"
risk=risk[sameSample, "Risk", drop=F] 
data=cbind(tide, risk)
#???ñȽ???
data$Risk=ifelse(data$Risk=="High", "High", "Low")
# ????????data??RiskΪ"High"????????"High-m6Ascore"??????????"Low-m6Ascore"
group=levels(factor(data$Risk))
# R????ʹ??factor???????ַ?????��ת??Ϊ??????��??ʹ??levels?????�????ӵ?ˮƽ
data$Risk=factor(data$Risk, levels=c("Low", "High"))
# group=levels(factor(data$Risk))
 comp=combn(group,2)
# combnֻ?Ǹ??????п??ܵ?????????????ԭ???ݵ?˳?????е?
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#??tide???ֽ???ѭ??,?ֱ?????С????ͼ
for(i in colnames(data)[1:(ncol(data)-1)]){
	gg1=ggviolin(data, x="Risk", y=i, fill = "Risk", 
	         xlab="", ylab=i,
	         palette=c("#0066FF","#FF0000"),
	         legend.title="m6Ascore",
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	
	pdf(file=paste0("violin.", i, ".pdf"), width=6, height=5)
	print(gg1)
	dev.off()
}




# sessionInfo()




library(reshape2)
library(dplyr)
library(tidyr)

data2 = dplyr::select(data, -No.benefits,-CTL.flag,-Responder) %>% dplyr::select(Risk,everything())


# melt(data,id.vars,measure.vars,variable.name=??variable??,??,na.rm=FALSE,value.name=??value??,
#      factorsAsStrings=TRUE)

# exp=cbind(exp, Type=sampleType)
# exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
data3=melt(data = data2, id.vars=c("Risk"))

colnames(data3)=c("Risk", "TIDEtype", "Expression_value")


p=ggboxplot(data3, x="TIDEtype", y="Expression_value", color = "Risk", 
            ylab="TIDEtype expression",
            xlab="",
            legend.title="Risk",
            palette = c("blue", "red"),
            width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Risk),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

# print(p1)




pdf(file="boxplot-01.pdf", width=7, height=5)
print(p1)
dev.off()





p2=p+stat_compare_means(aes(group=Risk),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.format")


# print(p2)




pdf(file="boxplot-02.pdf", width=7, height=5)
print(p2)
dev.off()





# geo data---------------------------------------------------------------------


data_geo = data[grep(pattern = "^GSM", x = row.names(data)),]

#??tide???ֽ???ѭ??,?ֱ?????С????ͼ
for(i in colnames(data_geo)[1:(ncol(data_geo)-1)]){
  gg1=ggviolin(data_geo, x="Risk", y=i, fill = "Risk", 
               xlab="", ylab=i,
               palette=c("#0066FF","#FF0000"),
               legend.title="m6Ascore",
               add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  
  pdf(file=paste0("violin.", i, "_geo.pdf"), width=6, height=5)
  print(gg1)
  dev.off()
}








library(reshape2)
library(dplyr)
library(tidyr)

# data4 = dplyr::select(data, -Responder) %>% dplyr::select(Risk,everything())
data4 = data

trait="Responder"


rt1=data4[,c(trait, "Risk")]
colnames(rt1)=c("trait", "group")
df=as.data.frame(table(rt1))


library(plyr)

df=ddply(df, .(group), transform, percent = Freq/sum(Freq) * 100)

df=ddply(df, .(group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
df$group=factor(df$group, levels=c("Low", "High"))


bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$Risk)))]


p=ggplot(df, aes(x = factor(group), y = percent, fill = trait)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("m6Ascore")+ ylab("Percent weight")+  guides(fill=guide_legend(title=trait))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  theme_bw()

print(p)


pdf(file="barplot.pdf", width=4, height=5)
print(p)
dev.off()









# rt2=rt[,c(trait, "m6Ascore")]
# colnames(rt2)=c("trait", "m6Ascore")
# type=levels(factor(rt2[,"trait"]))
# comp=combn(type, 2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}



# boxplot=ggboxplot(rt2, x="trait", y="m6Ascore", fill="trait",
#                   xlab=trait,
#                   ylab="m6Ascore",
#                   legend.title=trait,
#                   palette=bioCol
# )+ 
#   stat_compare_means(comparisons=my_comparisons)



# pdf(file="boxplot.pdf",width=4,height=4.5)
# print(boxplot)
# dev.off()
