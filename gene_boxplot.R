#生成mysql boxplot表格,三列数据如下:
# gene :  ENSGID  
# genename : TP53
# cancer: {"y":[0.6,0.7,0.3,0.6,0.5,0.7,0.9,0.5,0.8,0.7,0.2,0.6,0.7,0.3,0.6,0.5,0.7,0.9,0.5,0.8],"x":["BRCA","BRCA","BRCA","BRCA","BRCA","COAD","COAD","COAD","COAD","COAD","READ","READ","READ","READ","READ","UCEC","UCEC","UCEC","UCEC","UCEC"],"name":"cancer","boxpoints":"all","pointpos":0,"jitter":0.3,"marker":{"color":"#FF4136"},"type":"box"} 
# normal: {"y":[0.1,0.3,0.1,0.9,0.6,0.6,0.9,1,0.3,0.6,0.8,0.5],"x":["BRCA","BRCA","BRCA","COAD","COAD","COAD","READ","READ","READ","UCEC","UCEC","UCEC"],"name":"normal","boxpoints":"all","pointpos":0,"jitter":0.3,"marker":{"color":"#FF851B"},"type":"box"}
# yrange: [0,]
ExpressionMatrix <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_7_19/lncpcgMir.TCGA33cancerExpMatrix.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
SampleAnnotation <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_7_19/lncpcgMir.TCGA33cancerBarcode.txt",sep="\t",header=F,fill=T,quote=NULL,stringsAsFactors=F)
DiffExpression <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/processingData2019_7_19/lncpcgMir.TCGA24diff.table.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
head(SampleAnnotation)
SampleAnnotation$sampletype <- ifelse(grepl("normal",SampleAnnotation$V2),"normal",SampleAnnotation$V2)
#SampleAnnotation$sampletype1 <- ifelse(SampleAnnotation$V2 == SampleAnnotation$sampletype,SampleAnnotation$V2,paste0(SampleAnnotation$V2,SampleAnnotation$sampletype))
#colnames(SampleAnnotation) <- c("sampleid","cancertype","sampletype","sampletype1")
colnames(SampleAnnotation) <- c("sampleid","cancertype","sampletype")
library(stringr)
SampleAnnotation$cancertype <- str_split(SampleAnnotation$cancertype,"_",simplify = T)[,1]
class(ExpressionMatrix) #data.frame
tempMatrix <- ExpressionMatrix[,c(-1,-2)]
pos <- match(gsub("\\.","-",colnames(ExpressionMatrix[,c(-1,-2)])),SampleAnnotation$sampleid)
SampleAnnotation <- SampleAnnotation[pos,] #调整sampleannotation信息顺序
normal.x <- SampleAnnotation$cancertype[which(SampleAnnotation$sampletype == "normal")]
cancer.x <- SampleAnnotation$cancertype[which(SampleAnnotation$sampletype != "normal")]
unique(normal.x) #只有24个癌症有正常样本
boxplot.json <- data.frame(gene=character(0),genename=character(0),cancer=character(0),normal=character(0),yrange=character(0),stringsAsFactors = F)
library(dplyr)
library(rjson)  
for(i in 1:nrow(tempMatrix)){
  #normal.y <- as.numeric(tempMatrix[i,which(SampleAnnotation$sampletype == "normal")])
  #cancer.y <- as.numeric(tempMatrix[i,which(SampleAnnotation$sampletype != "normal")])
  normal.y <- log2(as.numeric(tempMatrix[i,which(SampleAnnotation$sampletype == "normal")])+1)
  cancer.y <- log2(as.numeric(tempMatrix[i,which(SampleAnnotation$sampletype != "normal")])+1)
  #var cancer = {"y":[],x:[],"name":"cancer","boxpoints":'suspectedoutliers',"pointpos":0,"jitter":0.6,"stream":{"maxpoints":50},"boxmean": "sd","marker":{"color":"rgb(7,40,89)","size":3},"type":"box"};
  #yrange <- range(boxplot.stats(c(normal.y,cancer.y))$stats)
  normal.list <- split(normal.y,normal.x) 
  normal.range.list <- lapply(normal.list,function(x){
    range(boxplot.stats(x)$stats)
  })
  #normal.range <- range(unlist((normal.range.list))) 
  cancer.list <- split(cancer.y,cancer.x) 
  cancer.range.list <- lapply(cancer.list,function(x){
    range(boxplot.stats(x)$stats)
  })
  #cancer.range <- range(unlist((cancer.range.list))) 
  yrange <- range(c(unlist(normal.range.list),unlist(cancer.range.list)))
  #????????????֢??*
  pos <- which(DiffExpression$GeneName==ExpressionMatrix[i,2])
  colnam <- colnames(DiffExpression[pos,c(-1,-2)])
  judgeVector <- as.numeric(DiffExpression[pos,c(-1,-2)])
  canpos <- seq(1,70,by=3)
  fcpos <- seq(2,71,by=3)
  fdrpos <- seq(3,72,by=3)
  flag <- ifelse(abs(judgeVector[fcpos]) >= 1 & judgeVector[fdrpos] <=0.05,T,F)
  sigCan <- colnam[canpos][flag]
  cancer.x.temp <- cancer.x
  cancer.x.temp[cancer.x.temp %in% sigCan] <- paste(cancer.x.temp[cancer.x.temp %in% sigCan],"*",sep = "")
  normal.x.temp <- normal.x
  normal.x.temp[normal.x.temp %in% sigCan] <- paste(normal.x.temp[normal.x.temp %in% sigCan],"*",sep = "")
  cancer <- list( y= cancer.y, x = cancer.x.temp, 
                  name = 'cancer', 
                  boxpoints = 'all',
                  pointpos = 0,
                  jitter = 0.6,
                  boxmean = TRUE,
                  marker = list(color='rgb(7,40,89)',size=2),
                  type = 'box')
  normal <- list( y= normal.y, x = normal.x.temp, 
                  name = 'normal', 
                  boxpoints = 'all',
                  pointpos = 0,
                  jitter = 0.6,
                  boxmean = TRUE,
                  marker = list(color='#3D9970',size=2),
                  type = 'box')
  cancer.json <- toJSON(cancer) #x是向量或者列??? 
  normal.json <- toJSON(normal)
  yrange <- toJSON(yrange)
  boxplot.json <- add_row(boxplot.json,gene=ExpressionMatrix[i,1],genename=ExpressionMatrix[i,2],cancer=cancer.json,normal=normal.json,yrange=yrange)
}
write.table(boxplot.json,file="D:/Rsources/Project/DataBase/HDfregulatory/processingData2019_7_19/boxplotdata.txt",sep="\t",quote=F,row.names=F,col.names=F)


