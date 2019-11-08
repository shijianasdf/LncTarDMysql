#'------------------------------------
#'
#'   处理boxplot
#'   生成mysql boxplot表格,三列数据如下:
#'   id : 0
#'   rid : RID00001
#    gene :  ENSGID  
#    genename : TP53
#    cancer: {"y":[0.6,0.7,0.3,0.6,0.5,0.7,0.9,0.5,0.8,0.7,0.2,0.6,0.7,0.3,0.6,0.5,0.7,0.9,0.5,0.8],"x":["BRCA","BRCA","BRCA","BRCA","BRCA","COAD","COAD","COAD","COAD","COAD","READ","READ","READ","READ","READ","UCEC","UCEC","UCEC","UCEC","UCEC"],"name":"cancer","boxpoints":"all","pointpos":0,"jitter":0.3,"marker":{"color":"#FF4136"},"type":"box"} 
#    normal: {"y":[0.1,0.3,0.1,0.9,0.6,0.6,0.9,1,0.3,0.6,0.8,0.5],"x":["BRCA","BRCA","BRCA","COAD","COAD","COAD","READ","READ","READ","UCEC","UCEC","UCEC"],"name":"normal","boxpoints":"all","pointpos":0,"jitter":0.3,"marker":{"color":"#FF851B"},"type":"box"}
#    yrange: [0,]
#'------------------------------------

cancer <- list( y= cancer.y, x = cancer.x.temp, 
                name = 'cancer', 
                boxpoints = 'all',
                pointpos = 0,
                jitter = 0.6,
                boxmean = TRUE,
                marker = list(color='rgb(7,40,89)',size=2),
                type = 'box' )
library(rjson)
library(dplyr)
## 读入boxplot数据
diffExp_boxplot <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_9_26/diffExp_boxplot.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
head(diffExp_boxplot)
colnames(diffExp_boxplot)

geoboxplot.json <- data.frame(id=character(0),rid=character(0),type=character(0),
           genename=character(0),cancer=character(0),normal=character(0),yrange=character(0))
for(i in 1:length(diffExp_boxplot$RID)){
  id <- i
  rid <- diffExp_boxplot[i,]$RID
  type <- diffExp_boxplot[i,]$type
  genename <- diffExp_boxplot[i,]$geneName
  pos <- which(unlist(strsplit(diffExp_boxplot[i,]$caseExp,";")) != "NA")
  cancer.y <- as.numeric( unlist(strsplit(diffExp_boxplot[i,]$caseExp,";"))[pos] )
  cancer.x.temp <- as.character(unlist(strsplit(diffExp_boxplot[i,]$caseDataset,";"))[pos])
  cancer <- list( y= cancer.y, x = cancer.x.temp, 
                  name = 'cancer', 
                  boxpoints = 'all',
                  pointpos = 0,
                  jitter = 0.6,
                  boxmean = TRUE,
                  marker = list(color='rgb(7,40,89)',size=2),
                  type = 'box' )
  pos1 <- which( unlist(strsplit(diffExp_boxplot[i,]$controlExp,";")) != "NA")
  normal.y <- as.numeric( unlist(strsplit(diffExp_boxplot[i,]$controlExp,";"))[pos1] )
  normal.x.temp <- as.character( unlist(strsplit(diffExp_boxplot[i,]$controlData,";"))[pos1] )
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
  yrange <- toJSON( range(c(range(boxplot.stats(cancer.y)$stats), range(boxplot.stats(normal.y)$stats) )))
  #yrange <- toJSON(yrange)
  geoboxplot.json <- add_row(geoboxplot.json,id=id,rid=rid,type=type,genename=genename,cancer=cancer.json,normal=normal.json,yrange=yrange)
}


write.table(geoboxplot.json,file="D:/Rsources/Project/DataBase/HDfregulatory/proceessingData2019_9_26/geoboxplotdata.txt",sep="\t",quote=F,row.names=F,col.names=F)



pos <- which( unlist(strsplit(diffExp_boxplot[47,]$caseExp,";")) != "NA" )
cancer.y <- as.numeric( unlist(strsplit(diffExp_boxplot[47,]$caseExp,";"))[pos] )
cancer.x.temp <- as.character( unlist(strsplit(diffExp_boxplot[47,]$caseDataset,";"))[pos] )
cancer <- list( y= cancer.y, x = cancer.x.temp, 
                name = 'cancer', 
                boxpoints = 'all',
                pointpos = 0,
                jitter = 0.6,
                boxmean = TRUE,
                marker = list(color='rgb(7,40,89)',size=2),
                type = 'box' )
pos1 <- which( unlist(strsplit(diffExp_boxplot[47,]$controlExp,";")) != "NA")
normal.y <- as.numeric( unlist(strsplit(diffExp_boxplot[47,]$controlExp,";"))[pos1] )
normal.x.temp <- as.character( unlist(strsplit(diffExp_boxplot[47,]$controlData,";"))[pos1] )
normal <- list( y= normal.y, x = normal.x.temp, 
                name = 'normal', 
                boxpoints = 'all',
                pointpos = 0,
                jitter = 0.6,
                boxmean = TRUE,
                marker = list(color='#3D9970',size=2),
                type = 'box')
















