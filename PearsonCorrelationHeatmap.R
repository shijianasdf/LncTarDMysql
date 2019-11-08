#regulator和target表达相关性热???
# rid 
# heatdata
correlationHeatmap <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_7_19/PCC.lncRNAtarget.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)
correlationHeatmap.p <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_7_19/PCCpvalue.lncRNAtarget.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F)

# var data = [
#     {
#      z: [[1, 0.3, 30, 0.5, 1,67], [20, 1, 60, 80, 30,34], [30, 60, 1, -10, 20,-100]], //一个矩???
#      x: ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday','Sunday'], //x轴题???
#      y: ['Morning', 'Afternoon', 'Evening'], //y轴题???
#      type: 'heatmap',
#      colorscale: 'YlGnBu', 
#      reversescale: true, 
#     }
#   ];
library(dplyr)
library(rjson)
HeatmapData <- data.frame(rid=character(0),heatdata=character(0),stringsAsFactors = F)
for(i in 1:nrow(correlationHeatmap)){
  colNames <- colnames(correlationHeatmap[i,c(-1,-2,-3)])
  pos <- !is.na(correlationHeatmap[i,c(-1,-2,-3)])
  x <-  colNames[pos]
  if(length(x) == 1){
    x <- list(x)
  }
  y <- c("correlation p value",paste(paste(correlationHeatmap[i,2:3],collapse = "_"),"PCC"))
  if(length(round(na.omit(as.numeric(correlationHeatmap.p[i,c(-1,-2,-3)])),digits = 3)) == 1 & length(round(na.omit(as.numeric(correlationHeatmap[i,c(-1,-2,-3)])),digits = 3)) == 1){
    z <- list(list(round(na.omit(as.numeric(correlationHeatmap.p[i,c(-1,-2,-3)])),digits = 3)),list(round(na.omit(as.numeric(correlationHeatmap[i,c(-1,-2,-3)])),digits = 3)))
  }
  else{
    z <- list(round(na.omit(as.numeric(correlationHeatmap.p[i,c(-1,-2,-3)])),digits = 3),round(na.omit(as.numeric(correlationHeatmap[i,c(-1,-2,-3)])),digits = 3))
  }
  testHeatmap <- list(z=z,
                      x=x, #TCGA cancer type
                      y=y, #regulator_target
                      type='heatmap',
                      zmin = -1,
                      zmid = 0,
                      zmax = 1, 
                      colorscale= list(c(-1, '#c4463a'),c(0,'#fffbbc'),c(0.9, '#c4463a'),c(1, '#3060cf')) #'YlGnBu'
                      #reversescale=TRUE
                      )
  jsonHeatmap <- toJSON(testHeatmap)
  HeatmapData <- add_row(HeatmapData,rid=correlationHeatmap[i,1],heatdata=jsonHeatmap)
}
write.table(HeatmapData,file="D:/Rsources/Project/DataBase/HDfregulatory/processingData2019_7_19/HeatmapData.txt",sep="\t",quote=F,row.names=F,col.names=T)

correlationHeatmap[which(correlationHeatmap$Regulator == "miR-21"),]
i <- 71
colNames <- colnames(correlationHeatmap[i,c(-1,-2,-3)])
pos <- !is.na(correlationHeatmap[i,c(-1,-2,-3)])
x <-  colNames[pos]
if(length(x) == 1){
  x <- list(x)
}
y <- c("correlation p value",paste(paste(correlationHeatmap[i,2:3],collapse = "_"),"PCC"))
if(length(round(na.omit(as.numeric(correlationHeatmap.p[i,c(-1,-2,-3)])),digits = 3)) == 1 & length(round(na.omit(as.numeric(correlationHeatmap[i,c(-1,-2,-3)])),digits = 3)) == 1){
  z <- list(list(round(na.omit(as.numeric(correlationHeatmap.p[i,c(-1,-2,-3)])),digits = 3)),list(round(na.omit(as.numeric(correlationHeatmap[i,c(-1,-2,-3)])),digits = 3)))
}else{
  z <- list(round(na.omit(as.numeric(correlationHeatmap.p[i,c(-1,-2,-3)])),digits = 3),round(na.omit(as.numeric(correlationHeatmap[i,c(-1,-2,-3)])),digits = 3))
}
testHeatmap <- list(z=z,
                    x=x, #TCGA cancer type
                    y=y, #regulator_target
                    type='heatmap',
                    zmin = -1,
                    zmid = 0,
                    zmax = 1, 
                    colorscale= list(c(-1, '#c4463a'),c(0,'#fffbbc'),c(0.9, '#c4463a'),c(1, '#3060cf')) #'YlGnBu'
                    #reversescale=TRUE
)
jsonHeatmap <- toJSON(testHeatmap)
HeatmapData <- add_row(HeatmapData,rid=correlationHeatmap[i,1],heatdata=jsonHeatmap)
             