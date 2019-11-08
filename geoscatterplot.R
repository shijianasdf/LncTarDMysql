#'------------------------------------
#'
#'   处理scatter plot
#'   生成mysql boxplot表格,三列数据如下:
#'   
#'--------------------------------------

#var trace1 = {
  #   x: [1, 2, 3, 4, 5], #regulator在GBM中疾病样本的表达值
  #   y: [1, 6, 3, 6, 1], #target在GBM中疾病样本的表达值
  #   mode: 'markers+text',
  #   type: 'scatter',
  #   name: 'GBM',
  #   text: ['TCGA SAMPLE1', 'TCGA SAMPLE2', 'TCGA SAMPLE3', 'TCGA SAMPLE4', 'TCGA SAMPLE5'],
  #   marker: { size: 12 }
  # };

## 读入散点图数据以及热图数据
targetPairs_scatterplot <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_9_26/targetPairs.PCC.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
head(targetPairs_scatterplot)

## 生成散点图数据
# id
# rid
# regulator
# target
# gseid
# scatterdata
geoscatterplot.json <- data.frame(id=character(0),rid=character(0),regulator=character(0),
                                  target=character(0),gseid=character(0),scatterdata=character(0))
library(rjson)
library(dplyr)
for(i in 1:length(targetPairs_scatterplot$RID)){
  id <- i
  RID <- targetPairs_scatterplot$RID[i]
  regulator <- targetPairs_scatterplot$Regulator[i]
  target <- targetPairs_scatterplot$Target[i]
  gseid <- targetPairs_scatterplot$GSEID[i]
  RTpos <- which(unlist(strsplit(targetPairs_scatterplot[i,]$expRegulator,";")) != "NA" & unlist(strsplit(targetPairs_scatterplot[i,]$expTarget,";")) != "NA")
  regulatorExpression <- unlist(strsplit(targetPairs_scatterplot[i,]$expRegulator,";"))[RTpos]
  targetExpression <- unlist(strsplit(targetPairs_scatterplot[i,]$expTarget,";"))[RTpos]
  samples <- unlist(strsplit(targetPairs_scatterplot[i,]$GSMID,";"))[RTpos]
  cancertype <- list( x = regulatorExpression ,
                      y = targetExpression,  
                      text = samples,   #
                      name = gseid,   # GSE00001
                      mode = 'markers',
                      marker =list(size = 10),
                      type = 'scatter')
  cancertype <- toJSON(cancertype)
  geoscatterplot.json <- add_row(geoscatterplot.json,id=id,rid=RID,regulator=regulator,target=target,gseid=gseid,scatterdata=cancertype)
}
write.table(geoscatterplot.json,file="D:/Rsources/Project/DataBase/HDfregulatory/proceessingData2019_9_26/geoScatterData.txt",sep="\t",quote=F,row.names=F,col.names=T)


## 生成热图矩阵
head(targetPairs_scatterplot)











