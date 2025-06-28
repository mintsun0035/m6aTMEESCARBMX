######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

expFile="uniSigGeneExp.txt"      #表达输入文件
setwd("C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\37.PCAscore")     #设置工作目录

#读取输入文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

#PCA分析
pca=prcomp(data, scale=TRUE)
value=predict(pca)
m6Ascore=value[,1]+value[,2]
m6Ascore=as.data.frame(m6Ascore)
scoreOut=rbind(id=colnames(m6Ascore), m6Ascore)
write.table(scoreOut, file="m6Ascore.txt", sep="\t", quote=F, col.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
