tem <- read.table(file="D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_7_19/lncpcgMir.TCGA24diff.table.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
matrix.padj <- tem[,colnames(tem)[grep("padj",colnames(tem))]]
any(na.omit(as.logical(matrix.padj == 0))) 

format.tem <- format(tem[,-c(1,2)],digits=3,scientific=T)
format.tem <- cbind.data.frame(tem[,c(1,2)],format.tem)
write.table(format.tem,file="D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_7_19/lncpcgMir.TCGA24diff.table.sientific.format.txt",sep="\t",quote=F,row.names=F,col.names=T)




BarCode <- read.table(file="D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_7_19/lncpcgMir.TCGA33cancerBarcode.txt",sep="\t",header=F,fill=T,quote=NULL,stringsAsFactors=F)
head(BarCode)
BarCode
colnames(BarCode) <- c("SampleID","CancerType")
library(stringr)
tt <- cbind.data.frame(BarCode$SampleID,str_split(BarCode$CancerType,"_",simplify = T))
colnames(tt) <- c("SampleID","CancerType","NorTur")
tempList <- split(tt,tt$CancerType)
lapply(tempList,function(x){
  table(x$NorTur)
})
table(tempList[[1]]$NorTur)
table(tempList[[2]]$NorTur)
grep("normal",as.character(tt$NorTur))

