Dat <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/example/合并9.0--HDFRegulatory_lncRNA.2017.04-2018.12.OK.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
Dat <- read.table("D:/LncTarD7_XAM.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
Dat <- read.table("D:/LncTarD7_XAMtoTree.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
head(Dat)
#查看entrezID是否为数字类???
class(Dat$RegulatorEntrezID)
class(Dat$TargetEntrezID)

#将evidence字段的双引号去掉
Evidence <- gsub('"','',Dat$Evidence)
Dat$Evidence <- Evidence
Dat$lncRegulatorPosition  
Dat$lncTargetPosition
write.table(Dat,file="D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_7_19/LncTarD2019-10-11.txt",sep="\t",quote=F,row.names=F,col.names=T)

#checkbox searchregulatorymechanism options
library(stringr)
regulatoryCategory <- str_split(unique(Dat$SearchregulatoryMechanism),";",simplify = T)
unique(regulatoryCategory[which(regulatoryCategory != "")])

#select drugs options
drugs <- str_split(na.omit(Dat$Drugs),";",simplify = T)
unique(drugs[which(drugs != "")])
drugs <- data.frame(drugs=unique(drugs[which(drugs != "")]))
write.table(drugs,file="D:/Rsources/Project/DataBase/HDfregulatory/drugs_select.txt",sep="\t",quote=F,row.names=F,col.names=F)
drugsOriginal <- read.table(file="D:/Rsources/Project/DataBase/HDfregulatory/drugs_select.txt",sep="\t",header=F,fill=T,quote=NULL,stringsAsFactors=F)
class(drugsOriginal)
class(drugs)
drugsOriginal[!drugsOriginal[,1] %in% drugs[,1]]
#xml function.xml
influencedFunction <- str_split(Dat$influencedFunction,";",simplify = T)
sink("D:/Rsources/Project/DataBase/HDfregulatory/xml/functions.txt")
paste(unique(influencedFunction[which(influencedFunction != "")]),collapse = ",")
sink()
#xml cancer.xml
sink("D:/Rsources/Project/DataBase/HDfregulatory/xml/cancer.txt")
paste(unique(Dat$DiseaseName),collapse=",")
sink()
#xml gene.xml
sink("D:/Rsources/Project/DataBase/HDfregulatory/xml/genes.txt")
paste(c(unique(na.omit(c(Dat$Regulator,Dat$Target))),unique(na.omit(c(Dat$TargetEntrezID,Dat$RegulatorEntrezID))),unique(na.omit(c(Dat$RegulatorEnsembleID,Dat$TargetEnsembleID)))),collapse = ",")
sink()


