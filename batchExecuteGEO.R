#'------------------------------------------------
#'
#'        @author : shijian
#'
#'        @description : 批量处理表达谱数据
#'
#'------------------------------------------------

## 读入批量处理的大表
Data.accession <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/data.accession.new1.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) 
LncTarD <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/LncTarD2.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
filepath <- "D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSETCGA/microarray"
list.files(filepath,recursive=T,full.names = T)
?list.files



exprSet <- read.table(file="D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSETCGA/microarray/GSE2658/RID00905-GSE2658-expression.MATRIX.txt",
                      row.names=NULL,sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
probeAnnotation <- read.table(file="D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSETCGA/microarray/GSE2658/RID00905-GSE2658-ProbeAnnotation.txt",
                              row.names=NULL,sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
head(exprSet)
head(exprSet[,-1])
exprSet[1:4,1:4]
colnames(exprSet)
head(probeAnnotation)
library(GEOquery)
gse <- getGEO("GSE2658",GSEMatrix = F) #获取soft信息
GPLid <- Meta(gse)$platform_id  
GPLid   ##GPL570 平台注释信息

colnames(probeAnnotation)
ids <- probeAnnotation[ ,c( "ID", "Gene.Symbol" ) ]
ids <- ids[which(!(ids$'Gene.Symbol' == "")),]  #过滤掉没有基因注释的探针
a <- strsplit( as.character(ids[,2]), " /// " )
tmp <- mapply( cbind, ids[,1], a )  #这是什么意思
ID2gene <- do.call(rbind.data.frame,tmp)
colnames( ID2gene ) <- c( "id", "gene" )
head(ID2gene)  ## 得到最终的探针id和geneID的配对信息
"HOX1" %in% ID2gene$gene
"HOX1I" %in% ID2gene$gene
"RUSAT1" %in% ID2gene$gene
target %in% ID2gene$gene
## 这一步根据ID2gene去除没有注释的探针
exprSet <- exprSet[ exprSet[,1] %in% ID2gene[ , 1 ], ]
ID2gene <- ID2gene[ match(exprSet[,1], ID2gene[ , 1 ] ), ]
dim( exprSet )
dim( ID2gene )
tail( sort( table( ID2gene[ , 2 ] ) ), n = 12L )

## 修正表达矩阵，添加行名
row.names <- exprSet[,1]
exprSet <- exprSet[,-1]
rownames(exprSet) <- row.names
exprSet[1:4,1:4]

##相同基因的表达数据取最大值，五万多个探针，这一步相对会运行较长时间
MAX <- by( exprSet, ID2gene[ , 2 ],function(x) rownames(x)[ which.max( rowMeans(x) ) ] )
MAX <- as.character(MAX)
exprSet <- exprSet[ rownames(exprSet) %in% MAX , ]
rownames( exprSet ) <- ID2gene[ match( rownames( exprSet ), ID2gene[ , 1 ] ), 2 ]





## 自动判断表达谱数据是否需要log2处理
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{ print("log2 transform not needed") }

## 根据"GSE2658" 找到互作的lncRNA-target，提取表达值，计算pearson相关系数以及p值
lncRNA.target.RID <- Data.accession[which(Data.accession$Data.accession == "GSE2658"),1]
regulator <- LncTarD$Regulator[which(LncTarD$RID == lncRNA.target.RID)]
target <- LncTarD$Target[which(LncTarD$RID == lncRNA.target.RID)]
c(regulator,target) %in% rownames(exprSet)


