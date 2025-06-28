

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("sva")
# 

#???ð?
rm(list = ls())


library(limma)
library(sva)
library(data.table)


# setwd("D:\\biowolf\\m6aTME\\16.merge")     #???ù???Ŀ¼
setwd("C:\\Users\\sun\\Desktop\\16.merge")     #???ù???Ŀ¼

# data.table::fread()


 files=c("symbol.txt", "GSE13898.txt")                  #?????ļ?????

# files=dir(pattern = ".txt$")

gc()  #?????ڴ?

#??ȡ????????
geneList=list()


for(i in 1:length(files)){
  
    # i = 4
    inputFile=files[i]
    # rt=read.table(inputFile, header=T, sep="\t",check.names=F)   #??1?в???Ϊ????
    # rt=data.table::fread(inputFile, header=T, sep="\t",check.names=F)   #??1?в???Ϊ????
    rt=data.table::fread(inputFile)   #??1?в???Ϊ????
    header=unlist(strsplit(inputFile, "\\.|\\-"))  #???ݷָ?
    geneList[[header[1]]]=unlist(as.vector(rt[,1]))    # ?????޸Ĵ˴?
    
    
}


intersectGenes=Reduce(intersect, geneList)

#??��??ȡ??????lapply??????????




# 222222 ------------------------------------------------------------------
# 222222 ------------------------------------------------------------------
# 222222 ------------------------------------------------------------------
# 222222 ------------------------------------------------------------------
# 222222 ------------------------------------------------------------------



#???ݺϲ?
allTab=data.frame()
batchType=c()



for(i in 1:length(files)){
  
  
  
    inputFile=files[i]
    header=unlist(strsplit(inputFile, "\\.|\\-"))
    #??ȡ?????ļ????????????ļ?????????
    rt=data.table::fread(inputFile, header=T, sep="\t", check.names=F)
    rt=as.matrix(rt)
    rownames(rt)=rt[,1]
    exp=rt[,2:ncol(rt)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
    rt=avereps(data)
    colnames(rt)=paste0(header[1], "_", colnames(rt))
    #??TCGAɾ????????Ʒ
    if(header[1] == "TCGA"){
		group=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
		group=sapply(strsplit(group,""), "[", 1)
		rt=rt[,group==0]
		rt=t(rt)
		row.names(rt)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(rt))
		rt=avereps(rt)
		rt=t(rt)
    }
    #????ֵ????????ȡlog2?????���
    qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
    if(LogC){
    	rt[rt<0]=0
        rt=log2(rt+1)}
    if(header[1] != "TCGA"){
    	rt=normalizeBetweenArrays(rt)
    }
     #???ݺϲ?
    if(i==1){
    	allTab=rt[intersectGenes,]
    }else{
    	allTab=cbind(allTab, rt[intersectGenes,])
    }
    batchType=c(batchType, rep(i,ncol(rt)))
}




#?????ݽ??н????????????????Ľ???
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge_GSE13898.txt", sep="\t", quote=F, col.names=F)

save.image(file = "m6aTME16.merge.RData")

