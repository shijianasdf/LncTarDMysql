#'-----------------------------------
#'       生成下载数据的大表
#'-----------------------------------

LncTarD <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/2019_9_20/LncTarD2.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) 
head(LncTarD)
data.accession <- cbind.data.frame(LncTarD$RID,LncTarD$DiseaseName,LncTarD$DiseaseName2,LncTarD$diseaseCategory,LncTarD$Experimental.method.for.lncRNA.expression,LncTarD$PubMedID,LncTarD$Data.accession,stringsAsFactors=F)
head(data.accession)

# 去除实验证实的差异表达，只保留GSE以及TCGA数据库的数据 
data.accession <- data.accession[which(!is.na(LncTarD$Data.accession)),]
data.accession <- data.accession[grep("GSE|TCGA",data.accession$`LncTarD$Data.accession`),]
colnames(data.accession) <- c("RID","DiseaseName","DiseaseName2","diseaseCategory","Experimental.method.for.lncRNA.expression","PubmedID","Data.accession")
head(data.accession)
library(stringr)
tt <- str_split(data.accession$Data.accession,";",simplify = T)
ttt <- str_split(data.accession$Data.accession,";")
head(tt)
head(ttt)
?str_split
data.accession.new <- data.frame(RID=character(0),DiseaseName=character(0),DiseaseName2=character(0),diseaseCategory=character(0),Experimental.method.for.lncRNA.expression=character(0),pubmedID=character(0),Data.accession=character(0))
for( i in 1:length(data.accession$RID) ){
  for( j in 1:length(ttt[[i]]) ){
    temp <- cbind.data.frame(data.accession[i,1:6],ttt[[i]][j])
    data.accession.new <- rbind.data.frame(data.accession.new,temp)
  }
}
head(data.accession.new)
# RID         DiseaseName      DiseaseName2 diseaseCategory Experimental.method.for.lncRNA.expression PubmedID ttt[[i]][j]
# 5  RID00005      Gastric cancer    Gastric cancer          Cancer                qPCR;microarray;sequencing 29985481     GSE1279
# 51 RID00005      Gastric cancer    Gastric cancer          Cancer                qPCR;microarray;sequencing 29985481        TCGA
# 27 RID00027 Lung adenocarcinoma       Lung cancer          Cancer                        qRT-PCR;sequencing 29857296        TCGA
# 35 RID00035            Melanoma          Melanoma          Cancer                        microarray;qRT-PCR 29956757     GSE3189
# 45 RID00045   Colorectal cancer Colorectal cancer          Cancer                                sequencing 29729381        TCGA
# 60 RID00060      Gastric cancer    Gastric cancer          Cancer                         microarray;RT-PCR 29743591    GSE95667
colnames(data.accession.new) <- c("RID","DiseaseName","DiseaseName2","diseaseCategory","Experimental.method.for.lncRNA.expression","PubmedID","Data.accession")

write.table(data.accession.new,file="D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/data.accession.new.txt",sep="\t",quote=F,row.names=F,col.names=T)















