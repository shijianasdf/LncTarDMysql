#regulator和target散点图
#生成mysql scatter表格,几列数据如下:
# rid
# gene.regulator: ENSG  regulator Ensemble id
# genename.regulator : regulator gene symbol 
# gene.target: ENSG  target Ensemble id
# genename.target : target gene symbol 
# GBM:
# LUAD:
# ...
# var trace1 = {
#   x: [1, 2, 3, 4, 5], #regulator在GBM中疾病样本的表达值
#   y: [1, 6, 3, 6, 1], #target在GBM中疾病样本的表达值
#   mode: 'markers+text',
#   type: 'scatter',
#   name: 'GBM',
#   text: ['TCGA SAMPLE1', 'TCGA SAMPLE2', 'TCGA SAMPLE3', 'TCGA SAMPLE4', 'TCGA SAMPLE5'],
#   marker: { size: 12 }
# };
ExpressionMatrix <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_7_19/lncpcgMir.TCGA33cancerExpMatrix.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
SampleAnnotation <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_7_19/lncpcgMir.TCGA33cancerBarcode.txt",sep="\t",header=F,fill=T,quote=NULL,stringsAsFactors=F)
#colnames(SampleAnnotation) <- c("sampleid","cancertype","sampletype","sampletype1")
SampleAnnotation$sampletype <- ifelse(grepl("normal",SampleAnnotation$V2),"normal",SampleAnnotation$V2)
colnames(SampleAnnotation) <- c("sampleid","cancertype","sampletype")
correlationHeatmap <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_7_19/PCC.lncRNAtarget.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)

list( x = ,
      y= ,  
      mode = 'marker', 
      name = ,   #
      text = ,   #
      marker =list(size = 12),
      type = 'scatter')

pos <- match(gsub("\\.","-",colnames(ExpressionMatrix[,c(-1,-2)])),SampleAnnotation$sampleid)
SampleAnnotation <- SampleAnnotation[pos,] #调整sampleannotation信息顺序
rightExpressionMatrix <- ExpressionMatrix[,c(-1,-2)][,which(SampleAnnotation$sampletype != "normal")] #删除正常样本
SampleAnnotation <- SampleAnnotation[SampleAnnotation$sampletype != "normal",] #删除正常样本

#依据样本癌症类型，对样本分组
cancertype.list <- split(colnames(rightExpressionMatrix),SampleAnnotation$cancertype)
        #对表达矩阵进行拆分，分成33种TCGA癌症的表达矩阵
ExpressionMatrixWithoutNormal <- cbind.data.frame(ExpressionMatrix[,1:2],rightExpressionMatrix)
library(rjson)
library(dplyr)
#rid,genename.regulator,genename.target
scatterplot.data <- data.frame(ACC=character(0),
                                BLCA=character(0),
                                BRCA=character(0),
                                CESC=character(0),
                                CHOL=character(0),
                                COAD=character(0),
                                DLBC=character(0),
                                ESCA=character(0),
                                GBM=character(0),
                                HNSC=character(0),
                                KICH=character(0),
                                KIRC=character(0),
                                KIRP=character(0),
                                LAML=character(0),
                                LGG=character(0),
                                LIHC=character(0),
                                LUAD=character(0),
                                LUSC=character(0),
                                MESO=character(0),
                                OV=character(0),
                                PAAD=character(0),
                                PCPG=character(0),
                                PRAD=character(0),
                                READ=character(0),
                                SARC=character(0),
                                SKCM=character(0),
                                STAD=character(0),
                                TGCT=character(0),
                                THCA=character(0),
                                THYM=character(0),
                                UCEC=character(0),
                                UCS=character(0),
                                UVM=character(0),stringsAsFactors = F)
for(i in 1:nrow(correlationHeatmap)){
  
  tempScatterPlot.data  <- list()
  for(j in 4:ncol(correlationHeatmap)){
    if(is.na(correlationHeatmap[i,j])){
      tempScatterPlot.data[[(j-3)]] <- NA
      print(paste0("j =========== ",j))
    }else{
      regulatorPos <- match(correlationHeatmap[i,2],ExpressionMatrixWithoutNormal$GeneName) #regulator pos
      targetPos <- match(correlationHeatmap[i,3],ExpressionMatrixWithoutNormal$GeneName) # target pos
      pos <- match(colnames(correlationHeatmap)[j],names(cancertype.list))
      regulatorExpression <- as.numeric(ExpressionMatrixWithoutNormal[regulatorPos,match(cancertype.list[[pos]],colnames(ExpressionMatrixWithoutNormal))])
      targetExpression <- as.numeric(ExpressionMatrixWithoutNormal[targetPos,match(cancertype.list[[pos]],colnames(ExpressionMatrixWithoutNormal))])
      cancertype <- list( x = regulatorExpression ,
                          y= targetExpression,  
                          text = cancertype.list[[pos]],   #
                          name = colnames(correlationHeatmap)[j],   #
                          mode = 'markers',
                          marker =list(size = 10),
                          type = 'scatter')
      cancertype <- toJSON(cancertype)
      tempScatterPlot.data[[(j-3)]] <- cancertype
      print(paste0("j =========== ",j))
    }
  }
  scatterplot.data <- add_row(scatterplot.data,ACC=tempScatterPlot.data[[1]],
                              BLCA=tempScatterPlot.data[[2]],
                              BRCA=tempScatterPlot.data[[3]],
                              CESC=tempScatterPlot.data[[4]],
                              CHOL=tempScatterPlot.data[[5]],
                              COAD=tempScatterPlot.data[[6]],
                              DLBC=tempScatterPlot.data[[7]],
                              ESCA=tempScatterPlot.data[[8]],
                              GBM=tempScatterPlot.data[[9]],
                              HNSC=tempScatterPlot.data[[10]],
                              KICH=tempScatterPlot.data[[11]],
                              KIRC=tempScatterPlot.data[[12]],
                              KIRP=tempScatterPlot.data[[13]],
                              LAML=tempScatterPlot.data[[14]],
                              LGG=tempScatterPlot.data[[15]],
                              LIHC=tempScatterPlot.data[[16]],
                              LUAD=tempScatterPlot.data[[17]],
                              LUSC=tempScatterPlot.data[[18]],
                              MESO=tempScatterPlot.data[[19]],
                              OV=tempScatterPlot.data[[20]],
                              PAAD=tempScatterPlot.data[[21]],
                              PCPG=tempScatterPlot.data[[22]],
                              PRAD=tempScatterPlot.data[[23]],
                              READ=tempScatterPlot.data[[24]],
                              SARC=tempScatterPlot.data[[25]],
                              SKCM=tempScatterPlot.data[[26]],
                              STAD=tempScatterPlot.data[[27]],
                              TGCT=tempScatterPlot.data[[28]],
                              THCA=tempScatterPlot.data[[29]],
                              THYM=tempScatterPlot.data[[30]],
                              UCEC=tempScatterPlot.data[[31]],
                              UCS=tempScatterPlot.data[[32]],
                              UVM=tempScatterPlot.data[[33]])
  print(paste0("i =========== ",i))
}
# regulatorPos <- match(correlationHeatmap[1,2],ExpressionMatrixWithoutNormal$GeneName) #regulator pos
# targetPos <- match(correlationHeatmap[1,3],ExpressionMatrixWithoutNormal$GeneName) # target pos
# pos <- match(colnames(correlationHeatmap)[6],names(cancertype.list))
# regulatorExpression <- as.numeric(ExpressionMatrixWithoutNormal[regulatorPos,match(cancertype.list[[pos]],colnames(ExpressionMatrixWithoutNormal))])
# targetExpression <- as.numeric(ExpressionMatrixWithoutNormal[targetPos,match(cancertype.list[[pos]],colnames(ExpressionMatrixWithoutNormal))])
# cancertype <- list( x = regulatorExpression ,
#                     y= targetExpression,  
#                     mode = 'markers', 
#                     name = colnames(correlationHeatmap)[6],   #
#                     text = cancertype.list[[pos]],   #
#                     marker =list(size = 10),
#                     type = 'scatter')
# cancertype <- toJSON(cancertype)
# writeLines(cancertype, "D:/Rsources/testScatterPlot.json")

test <- bind_cols(correlationHeatmap[,1:3],scatterplot.data)
write.table(test,file="D:/Rsources/Project/DataBase/HDfregulatory/processingData2019_7_19/ScatterData.txt",sep="\t",quote=F,row.names=F,col.names=T)

