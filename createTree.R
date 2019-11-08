#Dat <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/example/合并9.0--HDFRegulatory_lncRNA.2017.04-2018.12.OK.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
## 读入LncTarD5.txt
Dat <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/2019_9_25/LncTarD5.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
Dat <- read.table("D:LncTarDtoTree.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
Dat <- read.table("D:/LncTarD7_XAMtoTree.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
head(Dat)

Tree <- read.table("D:/tree.txt",sep="\t",header=F,fill=T,quote=NULL,stringsAsFactors=F)
head(Tree)
Dat$Regulator
all(unique(Tree[which(Tree[,2] == "lncRNA" | Tree[,2]=="miRNA" | Tree[,2]=="others" | Tree[,2]=="PCG" | Tree[,2]=="TF"),][,1]) %in% c(Dat$Regulator,Dat$Target))
pos <- which(!unique(Tree[which(Tree[,2] == "lncRNA" | Tree[,2]=="miRNA" | Tree[,2]=="others" | Tree[,2]=="PCG" | Tree[,2]=="TF"),][,1]) %in% c(Dat$Regulator,Dat$Target))
Tree[which(Tree[,2] == "lncRNA" | Tree[,2]=="miRNA" | Tree[,2]=="others" | Tree[,2]=="PCG" | Tree[,2]=="TF"),][pos,]
## 900 miR-140-6p miRNA miR-140-6p _self FALSE FALSE
c(Dat$Regulator,Dat$Target)[which(!c(Dat$Regulator,Dat$Target) %in% unique(Tree[which(Tree[,2] == "lncRNA" | Tree[,2]=="miRNA" | Tree[,2]=="others" | Tree[,2]=="PCG" | Tree[,2]=="TF"),][,1]))]


#Function_tree
{
  library(stringr)
  #Dat$influencedFunction
  functions <- str_split(Dat$influencedFunction,";",simplify = T)
  class(functions)
  influencedFunction <- unique(as.character(functions)[which(as.character(functions)!="")])
  function.tree<- data.frame(id=character(0),pId=character(0),name=character(0),target=character(0),open=character(0),icon=character(0),stringsAsFactors = F)
  functionid <- influencedFunction
  functionpId <- rep("Functions",length(functionid))
  functionName <- functionid
  functiontarget <- rep("_self",length(functionid))
  functionopen <- rep("FALSE",length(functionid))
  functionicon <- rep("FALSE",length(functionid))
  function.tree <- cbind.data.frame(functionid,functionpId,functionName,functiontarget,functionopen,functionicon,stringsAsFactors = F)
  colnames(function.tree) <- c("id","pId","name","target","open","icon")
  #Functions root
  functionRoot.tree <- data.frame(id="Functions",pId="root",name="Functions",target="_self",open="FALSE",icon="FALSE",stringsAsFactors = F)
  #Function.tree
  function.trees <- rbind.data.frame(function.tree,functionRoot.tree,stringsAsFactors = F)
  function.trees <- unique(function.trees)
}
unique(function.trees[order(function.trees$id),]$id)
#Drug_tree
{
  library(stringr)
  drugs <- unique(na.omit(Dat$Drugs))
  drugs <- str_split(drugs,";",simplify = T)
  as.character(drugs)
  drugs <- unique(as.character(drugs)[which(as.character(drugs)!="")])
  Drug.tree<- data.frame(id=character(0),pId=character(0),name=character(0),target=character(0),open=character(0),icon=character(0),stringsAsFactors = F)
  Drugid <- drugs
  DrugpId <- rep("Drugs",length(Drugid))
  DrugName <- Drugid
  Drugtarget <- rep("_self",length(Drugid))
  Drugopen <- rep("FALSE",length(Drugid))
  Drugicon <- rep("FALSE",length(Drugid))
  Drug.tree <- cbind.data.frame(Drugid,DrugpId,DrugName,Drugtarget,Drugopen,Drugicon,stringsAsFactors = F)
  colnames(Drug.tree) <- c("id","pId","name","target","open","icon")
  #????root
  DrugsRoot.tree <- data.frame(id="Drugs",pId="root",name="Drugs",target="_self",open="FALSE",icon="FALSE",stringsAsFactors = F)
  #?ϲ?��??tree
  Drugs.trees <- rbind.data.frame(Drug.tree,DrugsRoot.tree,stringsAsFactors = F)
  Drugs.trees <- unique(Drugs.trees)
}

#regulatoryMechanism tree
{
  unique(Dat$SearchregulatoryMechanism)
  library(stringr)
  regulatoryCategory <- str_split(unique(Dat$SearchregulatoryMechanism),";",simplify = T)
  regulatoryCategory <- unique(as.character(regulatoryCategory)[which(as.character(regulatoryCategory) != "")])
  regulatoryCategory.tree <- data.frame(id=character(0),pId=character(0),name=character(0),target=character(0),open=character(0),icon=character(0),stringsAsFactors = F)
  regulatoryCategoryid <- regulatoryCategory
  regulatoryCategorypId <- rep("RegulatoryMechanism",length(regulatoryCategoryid))
  regulatoryCategoryName <- regulatoryCategoryid
  regulatoryCategorytarget <- rep("_self",length(regulatoryCategoryid))
  regulatoryCategoryopen <- rep("FALSE",length(regulatoryCategoryid))
  regulatoryCategoryicon <- rep("FALSE",length(regulatoryCategoryid))
  regulatoryCategory.tree <- cbind.data.frame(regulatoryCategoryid,regulatoryCategorypId,regulatoryCategoryName,regulatoryCategorytarget,regulatoryCategoryopen,regulatoryCategoryicon,stringsAsFactors = F)
  colnames(regulatoryCategory.tree) <- c("id","pId","name","target","open","icon")
  
  regulatoryCategoryRoot.tree <- data.frame(id="RegulatoryMechanism",pId="root",name="RegulatoryMechanism",target="_self",open="FALSE",icon="FALSE",stringsAsFactors = F)
  regulatoryCategory.trees <- rbind.data.frame(regulatoryCategory.tree,regulatoryCategoryRoot.tree,stringsAsFactors = F)
  regulatoryCategory.trees <- unique(regulatoryCategory.trees)
}

#Diseases_tree
{
  diseases.category <- cbind.data.frame(Dat$DiseaseName,Dat$diseaseCategory,Dat$DiseaseName2,stringsAsFactors=F)
  colnames(diseases.category) <- c("diseasename","diseasecategoryRoot","diseasecategory1")
  diseases.category <- unique(diseases.category)
  #??ȡdiseasename??diseasenamecategory1??ͬ???У????д??? ?õ?simpleDisease.tree
  Simple.diseases <- diseases.category[which(diseases.category$diseasename == diseases.category$diseasecategory1),]
  simpleDisease.tree <- data.frame(id=character(0),pId=character(0),name=character(0),target=character(0),open=character(0),icon=character(0),stringsAsFactors = F)
  simpleDiseaseTreeid <- Simple.diseases$diseasename
  simpleDiseaseTreepId <- Simple.diseases$diseasecategoryRoot
  simpleDiseaseTreeName <- simpleDiseaseTreeid
  simpleDiseaseTreetarget <- rep("_self",length(simpleDiseaseTreeid))
  simpleDiseaseTreeopen <- rep("FALSE",length(simpleDiseaseTreeid))
  simpleDiseaseTreeicon <- rep("FALSE",length(simpleDiseaseTreeid))
  simpleDisease.tree <- cbind.data.frame(simpleDiseaseTreeid,simpleDiseaseTreepId,simpleDiseaseTreeName,simpleDiseaseTreetarget,simpleDiseaseTreeopen,simpleDiseaseTreeicon,stringsAsFactors=F)
  colnames(simpleDisease.tree) <- c("id","pId","name","target","open","icon")
  simpleDisease.tree <- unique(simpleDisease.tree)
  #??ȡdiseasename??diseasenamecategory1????ͬ???У????ж??????????ദ?? ?õ?Complex.diseases
  Complex.diseases <- diseases.category[which(diseases.category$diseasename != diseases.category$diseasecategory1),]
  ComplexDisease.tree <- data.frame(id=character(0),pId=character(0),name=character(0),target=character(0),open=character(0),icon=character(0),stringsAsFactors = F)
  ComplexDiseaseTreeid <- c(Complex.diseases$diseasename,Complex.diseases$diseasecategory1)
  ComplexDiseaseTreepId <- c(Complex.diseases$diseasecategory1,Complex.diseases$diseasecategoryRoot)
  ComplexDiseaseTreeName <- ComplexDiseaseTreeid
  ComplexDiseaseTreetarget <- rep("_self",length(ComplexDiseaseTreeid))
  ComplexDiseaseTreeopen <- rep("FALSE",length(ComplexDiseaseTreeid))
  ComplexDiseaseTreeicon <- rep("FALSE",length(ComplexDiseaseTreeid))
  ComplexDisease.tree <- cbind.data.frame(ComplexDiseaseTreeid,ComplexDiseaseTreepId,ComplexDiseaseTreeName,ComplexDiseaseTreetarget,ComplexDiseaseTreeopen,ComplexDiseaseTreeicon,stringsAsFactors = F)
  colnames(ComplexDisease.tree) <- c("id","pId","name","target","open","icon")
  ComplexDisease.tree <- unique(ComplexDisease.tree)
  #????diseaseCategory????Ŀ¼?Ĳ㼶?ṹ???õ?DiseaseCategory.tree
  DiseaseCategory.tree <- data.frame(id=character(0),pId=character(0),name=character(0),target=character(0),open=character(0),icon=character(0),stringsAsFactors = F)
  DiseasesCategoryid <- diseases.category$diseasecategoryRoot
  DiseaseCategorypId <- rep("DiseaseCategory",length(DiseasesCategoryid))
  DiseaseCategoryName <- DiseasesCategoryid
  DiseaseCategorytarget <- rep("_self",length(DiseasesCategoryid))
  DiseaseCategoryopen <- rep("FALSE",length(DiseasesCategoryid))
  DiseaseCategoryicon <- rep("FALSE",length(DiseasesCategoryid))
  DiseaseCategory.tree <- cbind.data.frame(DiseasesCategoryid,DiseaseCategorypId,DiseaseCategoryName,DiseaseCategorytarget,DiseaseCategoryopen,DiseaseCategoryicon,stringsAsFactors=F)
  colnames(DiseaseCategory.tree) <- c("id","pId","name","target","open","icon")
  DiseaseCategory.tree <- unique(DiseaseCategory.tree)
  #??Ϊroot???ڵ?ΪDiseaseCategory??tree
  DiseaseRoot.tree <- data.frame(id="DiseaseCategory",pId="root",name="DiseaseCategory",target="_self",open="FALSE",icon="FALSE",stringsAsFactors = F)
  #????4??tree
  Disease.trees <- rbind.data.frame(simpleDisease.tree,ComplexDisease.tree,DiseaseCategory.tree,DiseaseRoot.tree,stringsAsFactors = F)
  Disease.trees <- unique(Disease.trees)
  #1)ɾ???????ֶζ?Ϊcancer????
  #2)??DiseaseCategory??Other?ڵ??ŵ?DiseaseCategory??????չʾ
  #3)????DiseaseCategory?нڵ?Cancer?ڵ??½ڵ???˳?????ø??ڵ??????棬û???ӽڵ??Ľڵ???????
  Disease.trees <- Disease.trees[-which(Disease.trees$id == "Cancer" & Disease.trees$name == "Cancer" & Disease.trees$pId == "Cancer"),]
}

#Genes_tree
#?鿴regulator??target?????????ͣ???circRNA??snoRNA??Ϊothers
unique(c(Dat$RegulatorType,Dat$TargetType)) 

#"lncRNA"  "TF"      "PCG"     "miRNA"   "circRNA"  "snoRNA"
{
  regulator.type <- cbind.data.frame(Dat$Regulator,Dat$RegulatorType)
  colnames(regulator.type) <- c("geneSymbol","geneType")
  target.type <- cbind.data.frame(Dat$Target,Dat$TargetType)
  colnames(target.type) <- c("geneSymbol","geneType")
  genes.type <- rbind.data.frame(regulator.type,target.type) 
  genes.type <- unique(genes.type) 
  genes.type$geneSymbol <- as.character(genes.type$geneSymbol)
  genes.type$geneType <- as.character(genes.type$geneType)
  pos <- which(genes.type$geneType == "circRNA" | genes.type$geneType == "snoRNA")
  genes.type$geneType[pos] <- "others"

  Genes.tree <- data.frame(id=character(0),pId=character(0),name=character(0),target=character(0),open=character(0),icon=character(0),stringsAsFactors = F)
  id <- genes.type$geneSymbol
  pId <- genes.type$geneType
  name <- id
  target <- rep("_self",length(id))
  open <- rep("FALSE",length(id))
  icon <- rep("FALSE",length(id))
  Genes.tree <- cbind.data.frame(id,pId,name,target,open,icon,stringsAsFactors = F)
  #????lncRNA,TF,PCG,miRNA,others?Լ?root?ڵ?
  GenesRoot.tree <- data.frame(id="Genes",pId="root",name="Genes",target="_self",open="FALSE",icon="FALSE",stringsAsFactors = F)
  GenesLncRNA.tree <- data.frame(id="lncRNA",pId="Genes",name="lncRNA",target="_self",open="FALSE",icon="FALSE",stringsAsFactors = F)
  GenesmiRNA.tree <- data.frame(id="miRNA",pId="Genes",name="miRNA",target="_self",open="FALSE",icon="FALSE",stringsAsFactors = F)
  GenesPCG.tree <- data.frame(id="PCG",pId="Genes",name="PCG",target="_self",open="FALSE",icon="FALSE",stringsAsFactors = F)
  GenesTF.tree <- data.frame(id="TF",pId="Genes",name="TF",target="_self",open="FALSE",icon="FALSE",stringsAsFactors = F)
  Genesothers.tree <- data.frame(id="others",pId="Genes",name="others",target="_self",open="FALSE",icon="FALSE",stringsAsFactors = F)
  #?ϲ?tree
  Genes.trees <- rbind.data.frame(Genes.tree,GenesRoot.tree,GenesLncRNA.tree,GenesmiRNA.tree,GenesPCG.tree,GenesTF.tree,Genesothers.tree,stringsAsFactors = F)
  Genes.trees <- unique(Genes.trees)
}

#最终tree
tree.new <- rbind.data.frame(Genes.trees,Disease.trees,function.trees,regulatoryCategory.trees,Drugs.trees,stringsAsFactors = F)
#tree.new[tree.new$pId == "Functions",]$id[1] <- gsub(" ","",tree.new[tree.new$pId == "Functions",]$id[1])
newTree <- tree.new[order(tree.new$pId,tree.new$id),]


write.table(newTree,file="D:/Rsources/Project/DataBase/HDfregulatory/processingData2019_7_19/tree2019_10_11.txt",sep="\t",quote=F,row.names=F,col.names=F)























