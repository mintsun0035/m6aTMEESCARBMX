######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

#install.packages("RCircos")


library("RCircos")       #???ð?
setwd("C:\\Users\\sun\\Desktop\\m6a-esca\\14.Rcircos")    #???ù???Ŀ¼

#??ʼ??Ȧͼ
cytoBandIdeogram=read.table("refer.txt", header=T, sep="\t")
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 5
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

#????Ȧͼ????
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size=1
rcircos.params$point.size=5
RCircos.Reset.Plot.Parameters(rcircos.params)

#?????ļ?
pdf(file="RCircos.pdf", width=8, height=8)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

#ɢ??ͼ
RCircos.Scatter.Data=read.table("Rcircos.scatter.txt", header=T, sep="\t", check.names=F)
data.col <- 4
track.num <- 1
side <- "in"
RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col, track.num, side, by.fold=0.1)

#???ϻ???????
RCircos.Gene.Label.Data=read.table("Rcircos.geneLabel.txt", header=T, sep="\t", check.names=F)
name.col <- 4
side <- "in"
track.num <- 2
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
track.num <- 3
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056
