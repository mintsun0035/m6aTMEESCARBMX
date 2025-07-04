######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")


library(ConsensusClusterPlus)       #引用包
expFile="m6aGeneExp.txt"            #表达输入文件
workDir="C:\\Users\\TLSP\\Desktop\\m6a-ESCA\\23.m6aCluster"     #工作目录
setwd(workDir)       #设置工作目录

#读取输入文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

#聚类
maxK=9
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="pam",
              distance="euclidean",
              seed=123456,
              plot="png")


#输出分型结果
clusterNum=3        #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("m6Acluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$m6Acluster))
cluster$m6Acluster=letter[match(cluster$m6Acluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="m6aCluster.txt", sep="\t", quote=F, col.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
