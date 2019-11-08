#'-----------------------------------------------
#'              统计数据基本信息
#'      
#'-----------------------------------------------

## 读入数据
LncTarD <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/2019_9_25/LncTarD5.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
LncTarD <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_7_19/LncTarD2019-10-11.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
head(LncTarD)
length(unique(na.omit(LncTarD$Target)))
table(LncTarD$SearchregulatoryMechanism)
length(LncTarD$SearchregulatoryMechanism[LncTarD$SearchregulatoryMechanism != "expression association"])
##  lncRNA-target	biological functions	human diseases	lncRNAs	miRNAs	PCGs
length(unique(LncTarD$RID)) #lncRNA-target 2822
length(unique(LncTarD$DiseaseName)) #human diseases 178
length(unique(LncTarD$DiseaseName[LncTarD$Regulator=="HOTAIR"]))
length(unique(LncTarD$DiseaseName2[LncTarD$Regulator=="HOTAIR" | LncTarD$Target =="HOTAIR"]))
length(unique(LncTarD$DiseaseName[LncTarD$Regulator=="TP53" | LncTarD$Target =="TP53"]))

length(unique(LncTarD$Regulator)) ##多少个regulaotr
length(unique(LncTarD$Target))
tempRegulator <- c()
for(i in 1:length(unique(LncTarD$Regulator))){
  tempRegulator[i] <- length(unique(LncTarD$DiseaseName[LncTarD$Regulator == unique(LncTarD$Regulator)[i]]))
}
length(tempRegulator[tempRegulator != 1])
tempTarget <- c()
for(i in 1:length(unique(LncTarD$Target))){
  tempTarget[i] <- length(unique(LncTarD$DiseaseName[LncTarD$Target == unique(LncTarD$Target)[i]]))
}
length(tempTarget[tempTarget != 1])

unique(LncTarD$DiseaseName[!is.na(LncTarD$Drugs)])

newLncTarD <- cbind.data.frame(LncTarD, paste(LncTarD$Regulator,LncTarD$Target) )
temp <- c()
for( i in 1:length( unique(newLncTarD$`paste(LncTarD$Regulator, LncTarD$Target)`))){
  temp[i] <- length( unique(newLncTarD$DiseaseName[ which( newLncTarD$`paste(LncTarD$Regulator, LncTarD$Target)` == unique(newLncTarD$`paste(LncTarD$Regulator, LncTarD$Target)`)[i] ) ] ) )
}
length(temp)
length(temp[temp != 1])
length(temp[temp == 1])
dim(newLncTarD)

length(unique(c(LncTarD$regulatoryType,LncTarD$TargetType)))
table(LncTarD$RegulatorType)
table(LncTarD$TargetType)
unique(c(LncTarD$Regulator[which(LncTarD$RegulatorType == "lncRNA")],LncTarD$Target[which(LncTarD$TargetType == "lncRNA")])) #lncRNA 475
unique(c(LncTarD$Regulator[which(LncTarD$RegulatorType == "miRNA")],LncTarD$Target[which(LncTarD$TargetType == "miRNA")])) #miRNA 222

temp <- unlist(strsplit(LncTarD$regulatoryMechanism[grep("miR|let",LncTarD$regulatoryMechanism)],"\\("))[grep("miR|let",unlist(strsplit(LncTarD$regulatoryMechanism[grep("miR|let",LncTarD$regulatoryMechanism)],"\\(")))]
unique(c(LncTarD$Regulator[which(LncTarD$RegulatorType == "miRNA")],LncTarD$Target[which(LncTarD$TargetType == "miRNA")],unlist(strsplit(unlist(strsplit(temp,";")),")")))) #miRNA 391

length(unique(c(LncTarD$Regulator[which(LncTarD$RegulatorType == "PCG" | LncTarD$RegulatorType == "TF")],LncTarD$Target[which(LncTarD$TargetType == "PCG" | LncTarD$TargetType == "TF")]))) #PCG 783

## drugs
drugs <- str_split(na.omit(LncTarD$Drugs),";",simplify = T)
unique(drugs[which(drugs != "")])

## 统计影响功能数
library(stringr)
class(unlist(str_split(unique(LncTarD$influencedFunction),";")))
unlist(str_split(unique(LncTarD$influencedFunction),";"))
unique(substring(unlist(str_split(unique(LncTarD$influencedFunction),";")),1,
          nchar(unlist(str_split(unique(LncTarD$influencedFunction),";")))-3))



## ceRNA/sponge	epigenetic regulation	transcriptional regulation	interact with protein	interact with mRNA	chromatin looping
table(LncTarD$SearchregulatoryMechanism)
unique(unlist(strsplit(unique(LncTarD$SearchregulatoryMechanism),";")))
table(unlist(strsplit(LncTarD$SearchregulatoryMechanism,";")))
table(LncTarD$SearchregulatoryMechanism)


## 寻找expression association例子
head(LncTarD)
exp.LncTarD <- LncTarD[which(LncTarD$SearchregulatoryMechanism == "expression association"),]
head(exp.LncTarD)
table(exp.LncTarD$Data.accession2,useNA="ifany")
?table

exp.dataAccession.LncTarD <- exp.LncTarD[!is.na(exp.LncTarD$Data.accession2),]
head(exp.dataAccession.LncTarD)
table(paste(exp.dataAccession.LncTarD$Regulator,exp.dataAccession.LncTarD$Target,sep="_"))

## 统计前20种癌症调控个数
head(LncTarD)
rev(sort(table(LncTarD$DiseaseName2)))[1:20]

## 统计前20个lncRNA的调控个数以及疾病数
lnc <- c(LncTarD$Regulator[LncTarD$RegulatorType == "lncRNA"],LncTarD$Target[LncTarD$TargetType == "lncRNA"])   
rev(sort(table(lnc)))[1:20]
top20genes <- c("MALAT1",	"HOTAIR",	"H19",	"MEG3",	"UCA1",	"PVT1",	"CDKN2B-AS1",	"GAS5",	"NEAT1",	"TUG1",	"XIST",	"CCAT1",	"HOTTIP",	"HULC",	"SNHG1",	"LINC-ROR",	"HOXA11-AS",	"CRNDE",	"HNF1A-AS1",	"CASC2")

length(unique(LncTarD$DiseaseName[which(LncTarD$Regulator == "MALAT1" | LncTarD$Target == "MALAT1")]))
for(i in 1:length(top20genes)){
  print(length(unique(LncTarD$DiseaseName[which(LncTarD$Regulator == top20genes[i] | LncTarD$Target == top20genes[i])])))
}

## top 20个功能的调控个数以及疾病个数
Dat <- read.table("D:/LncTarD7_XAMtoTree.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
library(stringr)
#Dat$influencedFunction
functions <- str_split(Dat$influencedFunction,";",simplify = T)
sort(unique(as.character(functions)[as.character(functions) != ""])) 
class(functions)
influencedFunction <- as.character(functions)[which(as.character(functions)!="")]
rev(sort(table(influencedFunction)))[1:20]

top20influencedfunction <- names(rev(sort(table(influencedFunction)))[1:20])
for(i in 1:length(top20influencedfunction)){
  #which(LncTarD$Drugs == top20drugs[i])
  print( length( unique(LncTarD$DiseaseName[grep(top20influencedfunction[i],LncTarD$influencedFunction)]) ) )
}
tt <- c()
for(i in 1:length(top20influencedfunction)){
  #which(LncTarD$Drugs == top20drugs[i])
  tt[i] <- length( unique(LncTarD$DiseaseName[grep(top20influencedfunction[i],LncTarD$influencedFunction)]) ) 
}
paste(tt,collapse = ",")
## top 20 drug的调控个数以及对应人类疾病个数
drugs <- str_split(na.omit(LncTarD$Drugs),";",simplify = T)
uniquedrug <- as.character(drugs)[as.character(drugs) != ""] 
rev(sort(table(uniquedrug)))[1:22]

top20drugs <- c("cisplatin","5-fluorouracil","paclitaxel","temozolomide","doxorubicin","adriamycin","imatinib","oxaliplatin","platinum","gemcitabine","trastuzumab","docetaxel","tamoxifen","gefitinib","sorafenib","cetuximab","etoposide","lapatinib","sulforaphane","sunitinib")
for(i in 1:length(top20drugs)){
  #which(LncTarD$Drugs == top20drugs[i])
  print( length( unique(LncTarD$DiseaseName[grep(top20drugs[i],LncTarD$Drugs)]) ) )
}

## PRC2 EZH2 SUG12 EUD

## GAS5
dim(LncTarD[LncTarD$Regulator == "GAS5",])
GAS5_LncTarD <- LncTarD[LncTarD$Regulator == "GAS5",]
GAS5_LncTarD


