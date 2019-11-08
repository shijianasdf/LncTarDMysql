#'---------------------------
#'   
#'   制作下载GSE列表的数据
#'   
#'---------------------------
## 读入数据
LncTarD <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/LncTarD.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) 
head(LncTarD)
colnames(LncTarD)
cbind.data.frame(LncTarD$RID,LncTarD$Data.accession) 

