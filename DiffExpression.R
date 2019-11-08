#将差异表达的平均表达值，FDR，log2foldchange全部取小数点后4位展示
DiffExpression <-  read.table("D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_7_19/lncpcgMir.TCGA24diff.table.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
tempExpression <- as.matrix(DiffExpression[,c(-1,-2)])
tempExpression <- round(tempExpression,digits = 4) 
DiffExpression <- cbind.data.frame(DiffExpression[,1:2],tempExpression)
write.table(DiffExpression,file="D:/Rsources/Project/DataBase/HDfregulatory/processingData2019_7_19/lncpcgMir.TCGA24diff.table.txt",sep="\t",quote=F,row.names=F,col.names=T,na = "NA",fileEncoding="UTF-8")

mart_export <- read.table("D:/Rsources/Project/Ԥ??????/??¡Ԥ??????/??ֱ??????¡Ԥ??????/colorectal.expression.data/mart_export.txt",
                          sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
library(dplyr)
left_join(DiffExpression,mart_export,by=c("Gene" = "Gene.stable.ID"))

