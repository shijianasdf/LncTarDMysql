#----------------------------------------
#'      生成echart tree表格
#----------------------------------------
example <- read.table(file="D:/Rsources/Project/DataBase/HDfregulatory/example/合并9.0--HDFRegulatory_lncRNA.2017.04-2018.12.OK.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F,encoding = "UTF-8")
head(example)

#evidence数据清洗
Evidence <- gsub('"','',example$Evidence)
example$Evidence <- Evidence
#Data Accession如果为""变为NA
sort(table(example$Data.accession,useNA = "ifany"))
#table(example$Data.accession == "")
#example$Data.accession[example$Data.accession == ""] <- NA
#将主表输出
write.table(example,file="D:/Rsources/Project/DataBase/HDfregulatory/example/example.txt",sep="\t",quote=F,row.names=F,col.names=T)


#清洗数据，将regulatoryMechanism中的NA 根据RegulationDiretion中的positively和negativelay替换成upregulated和downregulated
table(example$regulatoryMechanism,useNA = "ifany")
pos <- which(is.na(example$regulatoryMechanism))
length(pos)
example$regulatoryMechanism[pos] <- example$RegulationDiretion[pos]
table(example$regulatoryMechanism[pos])
example$regulatoryMechanism[pos][grep("positively",example$regulatoryMechanism[pos])] <- "upregulated"
example$regulatoryMechanism[pos][grep("negatively",example$regulatoryMechanism[pos])] <- "downregulated"
#清洗diseasecategory和diseasename2,将字符串中的前后空格去掉
example$DiseaseName2 <- gsub("^\\s+|\\s+$", "",example$DiseaseName2)
example$diseaseCategory <- gsub("^\\s+|\\s+$", "",example$diseaseCategory)
example$DiseaseName <- gsub("^\\s+|\\s+$", "",example$DiseaseName)

#生成echart网络数据
netexample <- cbind.data.frame(example$DiseaseName,example$Regulator,example$RegulatorEntrezID,example$RegulatorEnsembleID,example$RegulatorAliases, example$Target,example$TargetEntrezID,example$TargetEnsembleID,example$TargetAliases,example$regulatoryMechanism,
                               example$influencedFunction,example$Drugs,example$RID,example$SearchregulatoryMechanism,example$diseaseCategory,example$DiseaseName2,stringsAsFactors=F)
head(netexample)
class(netexample$`example$DiseaseName2`)
netexample$`example$DiseaseName2` <- gsub("^\\s+|\\s+$", "",netexample$`example$DiseaseName2`) 
netexample$`example$diseaseCategory` <- gsub("^\\s+|\\s+$", "",netexample$`example$diseaseCategory`) 
write.table(netexample,file="D:/Rsources/Project/DataBase/HDfregulatory/example/netexample.txt",sep="\t",quote=F,row.names=F,col.names=F)



#生成echart网络数据
example <- read.table(file="D:/Rsources/Project/DataBase/HDfregulatory/2019_9_20/LncTarD2.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F,encoding = "UTF-8")
head(example)
colnames(example)
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
write.table(netexample,file="D:/Rsources/Project/DataBase/HDfregulatory/2019_9_20/netEcharts1.txt",sep="\t",quote=F,row.names=F,col.names=F)











