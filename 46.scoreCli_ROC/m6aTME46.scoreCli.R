

#install.packages("ggplot2")
#install.packages("ggpubr")


#???ð?
library(plyr)
library(ggplot2)
library(ggpubr)
scoreFile="m6Ascore.group.txt"    #m6A?????ļ?
cliFile="clinical.txt"            #?ٴ??????ļ?
trait="Fustat"                    #?ٴ???״?ɸ???ԭΪFustat
setwd("O:\\m6a-esca\\46.scoreCli")     #???ù???Ŀ¼

#??ȡ?????ļ?
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(score), row.names(cli))
rt=cbind(score[sameSample,,drop=F], cli[sameSample,,drop=F])

#?????ٴ???״????ɫ
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt[,trait]))]

#ͳ?Ƹߵ??��??鲡????Ŀ
rt1=rt[,c(trait, "group")]
colnames(rt1)=c("trait", "group")
df=as.data.frame(table(rt1))
#?????ߵ??��????İٷ???
df=ddply(df, .(group), transform, percent = Freq/sum(Freq) * 100)
#?ٷֱ?λ??
df=ddply(df, .(group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
df$group=factor(df$group, levels=c("Low", "High"))
#???ưٷ???ͼ
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

#???ñȽ???
rt2=rt[,c(trait, "m6Ascore")]
colnames(rt2)=c("trait", "m6Ascore")
type=levels(factor(rt2[,"trait"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
#????????ͼ
boxplot=ggboxplot(rt2, x="trait", y="m6Ascore", fill="trait",
		          xlab=trait,
		          ylab="m6Ascore",
		          legend.title=trait,
		          palette=bioCol
		          )+ 
	    stat_compare_means(comparisons=my_comparisons)
pdf(file="boxplot.pdf",width=4,height=4.5)
print(boxplot)
dev.off()


######Vi