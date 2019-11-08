#'-------------------------------------
#'
#'       生成hdfnetwork
#'-------------------------------------

## 读入数据
example <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/2019_9_25/LncTarD5.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) 
example <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_7_19/LncTarD2019-10-11.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) 
head(example)
#清洗数据，将regulatoryMechanism中的NA根据RegulationDiretion中的positively和negativelay替换成upregulated和downregulated
table(example$regulatoryMechanism,useNA = "ifany")
pos <- which(is.na(example$regulatoryMechanism))
length(pos)
example$regulatoryMechanism[pos] <- example$RegulationDiretion[pos]
table(example$regulatoryMechanism[pos])
example$regulatoryMechanism[pos][grep("positively",example$regulatoryMechanism[pos])] <- "upregulated"
example$regulatoryMechanism[pos][grep("negatively",example$regulatoryMechanism[pos])] <- "downregulated"
netexample <- cbind.data.frame(example$DiseaseName,example$Regulator,example$RegulatorEntrezID,example$RegulatorEnsembleID,example$RegulatorAliases, example$Target,example$TargetEntrezID,example$TargetEnsembleID,example$TargetAliases,example$regulatoryMechanism,
                               example$influencedFunction,example$Drugs,example$RID,example$SearchregulatoryMechanism,example$diseaseCategory,example$DiseaseName2,stringsAsFactors=F)
head(netexample)
write.table(netexample,file="D:/Rsources/Project/DataBase/HDfregulatory/2019_9_25/netEcharts2019_10_11.txt",sep="\t",quote=F,row.names=F,col.names=T)

