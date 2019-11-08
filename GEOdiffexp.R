#'-------------------------------------
#'
#'   GEO 差异表达的表格
#'
#'-------------------------------------

## 读入差异表达数据
targetPairs.diffExp <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_9_26/targetPairs.diffExp.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
head(targetPairs.diffExp)
#id rid genename type gseid meanexp log2foldchange padj control cases

#id rid target gseid rmeanexp rlog2foldchange rpadj  regulator tmeanexp tlog2foldchange tpadj control cases
#1 RID00914 LINC00312 GSE10072 8.99 -1.4000 2.91e-26 HOXA5 8.99 -1.4000 2.91e-26 49 58
tt <- targetPairs.diffExp[targetPairs.diffExp$RID == "RID00914",]
tt$type

LncTarD <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/2019_9_25/LncTarD4.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
Evidence <- gsub('"','',LncTarD$Evidence)
LncTarD$Evidence <- Evidence
dim(LncTarD)
LncTarD[,c(31,32)]
write.table(LncTarD[,c(-31,-32)],file="D:/Rsources/Project/DataBase/HDfregulatory/2019_9_25/LncTarD5.txt",sep="\t",quote=F,row.names=F,col.names=T)

LncTarD <- LncTarD[,c(-31,-32)]
head(LncTarD)

