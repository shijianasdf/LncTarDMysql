#'--------------------------------------------
#'   
#'             tree按字母排序
#'   
#'--------------------------------------------
## 读入数据
tree <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/tree.txt",sep="\t",header=F,fill=T,quote=NULL,stringsAsFactors=F) 
head(tree)
colnames(tree) <- c("id","pId","name","target","open","icon")
##写了这么多程序，其实一个函数就搞定了
{
  ## root tree
  {
    ##root tree
    root.tree <- tree[tree$pId == "root",]
  }
  ## gene tree排序
  {
    ##genes tree
    Genes.tree <- tree[tree$pId == "Genes",]
    
    ##genes lncRNA内部排序
    lncRNA.tree <- tree[which(tree$pId == "lncRNA"),]
    dim(lncRNA.tree)
    head(lncRNA.tree)
    head(lncRNA.tree[order(lncRNA.tree$id),])
    dim(lncRNA.tree[order(lncRNA.tree$id),])
    lncRNA.tree <- lncRNA.tree[order(lncRNA.tree$id),]
    #tree <- tree[which(tree$pId == "lncRNA"),][order(tree[which(tree$pId == "lncRNA"),]$id),]
    
    ##genes miRNA内部排序
    miRNA.tree <- tree[which(tree$pId == "miRNA"),]
    dim(miRNA.tree)
    head(miRNA.tree)
    head(miRNA.tree[order(miRNA.tree$id),])
    dim(miRNA.tree[order(miRNA.tree$id),])
    miRNA.tree <- miRNA.tree[order(miRNA.tree$id),]
    #tree <- tree[which(tree$pId == "miRNA"),][order(tree[which(tree$pId == "miRNA"),]$id),]
    
    ##genes PCG内部排序
    PCG.tree <- tree[which(tree$pId == "PCG"),]
    dim(PCG.tree)
    head(PCG.tree)
    head(PCG.tree[order(PCG.tree$id),])
    dim(PCG.tree[order(PCG.tree$id),])
    PCG.tree <- PCG.tree[order(PCG.tree$id),]
    
    ##genes TF内部排序
    TF.tree <- tree[which(tree$pId == "TF"),]
    dim(TF.tree)
    head(TF.tree)
    head(TF.tree[order(TF.tree$id),])
    dim(TF.tree[order(TF.tree$id),])
    TF.tree <- TF.tree[order(TF.tree$id),]
    
    ##genes others内部排序
    others.tree <- tree[which(tree$pId == "others"),]
    dim(others.tree)
    head(others.tree)
    head(others.tree[order(others.tree$id),])
    dim(others.tree[order(others.tree$id),])
    others.tree <- others.tree[order(others.tree$id),]
  }
  ## Functions tree排序
  {
    Functions.tree <- tree[tree$pId == "Functions",]
    dim(Functions.tree)
    head(Functions.tree)
    head(Functions.tree[order(Functions.tree$id),])
    dim(Functions.tree[order(Functions.tree$id),])
    Functions.tree <- Functions.tree[order(Functions.tree$id),]
  }
  ## RegulatoryMechanism tree排序
  {
    RegulatoryMechanism.tree <- tree[tree$pId == "RegulatoryMechanism",]
    # dim(RegulatoryMechanism.tree)
    # head(RegulatoryMechanism.tree)
    # head(RegulatoryMechanism.tree[order(RegulatoryMechanism.tree$id),])
    # dim(RegulatoryMechanism.tree[order(RegulatoryMechanism.tree$id),])
    # RegulatoryMechanism.tree <- RegulatoryMechanism.tree[order(RegulatoryMechanism.tree$id),]
  }
  ## Drugs tree排序
  {
    Drugs.tree <- tree[tree$pId == "Drugs",]
    dim(Drugs.tree)
    head(Drugs.tree)
    head(Drugs.tree[order(Drugs.tree$id),])
    dim(Drugs.tree[order(Drugs.tree$id),])
    Drugs.tree <- Drugs.tree[order(Drugs.tree$id),]
  }
  ## disease tree排序
  {
    ##--------------------disease tree内部排序------------------------
    ##DiseaseCategory tree排序
    DiseaseCategory.tree <- tree[which(tree$pId == "DiseaseCategory"),]
    dim(DiseaseCategory.tree)
    head(DiseaseCategory.tree)
    head(DiseaseCategory.tree[order(DiseaseCategory.tree$id),])
    dim(DiseaseCategory.tree[order(DiseaseCategory.tree$id),])
    DiseaseCategory.tree <- DiseaseCategory.tree[order(DiseaseCategory.tree$id),]
    
    ##DiseaseCategory 中第一个元素cancer tree排序
    cancer.tree <- tree[which(tree$pId == DiseaseCategory.tree$id[1]),]
    cancer.tree <- cancer.tree[order(cancer.tree$id),]
    head(cancer.tree)
  }
  
}
head(tree[order(tree$pId,tree$id),],50)
tree[tree$pId == "Functions",]$id[1] <- gsub(" ","",tree[tree$pId == "Functions",]$id[1])
newTree <- tree[order(tree$pId,tree$id),]
newTree[which(newTree$pId == "RegulatoryMechanism"),]$id
write.table(newTree,file="D:/Rsources/Project/DataBase/HDfregulatory/tree1.txt",sep="\t",quote=F,row.names=F,col.names=F)
newTree[which(newTree$name == "other"),]
tree[tree$pId == "Functions",]

## 读入数据,将主表的SearchregulatoryMechanism字段的"other"转变为"expression association"
LncTarD <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/LncTarD.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) 
head(LncTarD)
LncTarD$SearchregulatoryMechanism[which(LncTarD$SearchregulatoryMechanism == "other")] <- "expression association"
write.table(LncTarD,file="D:/Rsources/Project/DataBase/HDfregulatory/2019_9_20/LncTarD1.txt",sep="\t",quote=F,row.names=F,col.names=T)

LncTarD <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/2019_9_20/LncTarD1.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) 
LncTarD$SearchregulatoryMechanism[which(LncTarD$SearchregulatoryMechanism == "Expression association")] <- "expression association"
write.table(LncTarD,file="D:/Rsources/Project/DataBase/HDfregulatory/2019_9_20/LncTarD2.txt",sep="\t",quote=F,row.names=F,col.names=T)








