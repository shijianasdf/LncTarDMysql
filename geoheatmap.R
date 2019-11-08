#'-------------------------------------
#'
#'   GEO heatmap data
#'
#'-------------------------------------

# {"z":[[0,0.702,0.893,0.504,0.516,0.901,0.628,0.902,0.583,0.891,0.818,0.559,0.39,0.804,0.838,0.486,0.597,0.026,0.119,0.568,0.897],[0.104,-0.022,0.006,-0.099,0.048,0.01,0.021,-0.005,0.032,0.01,-0.01,-0.03,-0.037,-0.011,-0.01,0.031,0.033,0.109,-0.125,-0.025,0.006]],"x":["BRCA","CESC","COAD","DLBC","ESCA","GBM","HNSC","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PRAD","SARC","STAD","TGCT","THCA","UCEC"],"y":["correlation p value","LINC00313_miR-4429 PCC"],"type":"heatmap","zmin":-1,"zmid":0,"zmax":1,"colorscale":[["-1","#c4463a"],["0","#fffbbc"],["0.9","#c4463a"],["1","#3060cf"]]}

##  读入heatmap data
targetPairs_scatterplot <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/rawData2019_9_26/targetPairs.PCC.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) #记录转录因子结合位点以及TFClass
head(targetPairs_scatterplot)

##  对heatmap data初步处理

library(dplyr)
library(rjson)
HeatmapData <- data.frame(rid=character(0),heatdata=character(0),stringsAsFactors = F)
for(i in 1:length(unique(targetPairs_scatterplot$RID))){
  tt <- targetPairs_scatterplot[which(targetPairs_scatterplot$RID == unique(targetPairs_scatterplot$RID)[i]),]
  if(dim(tt)[1] == 1){
    z <- list(list(round(tt$p,digits = 3)),list(round(tt$PCC,digits=3)))
  }else{
    z <- list(round(tt$p,digits = 3),round(tt$PCC,digits=3))
  }
  x <- tt$GSEID
  if(length(tt$GSEID) == 1){
    x <- list(tt$GSEID)
  }
  y <- c("correlation p value",paste0(tt$Regulator[1],"_",tt$Target[1]," pcc")) 
  testHeatmap <- list(z=z, #[[0,0.702,0.893,0.504,0.516,0.901,0.628,0.902,0.583,0.891,0.818,0.559,0.39,0.804,0.838,0.486,0.597,0.026,0.119,0.568,0.897],[0.104,-0.022,0.006,-0.099,0.048,0.01,0.021,-0.005,0.032,0.01,-0.01,-0.03,-0.037,-0.011,-0.01,0.031,0.033,0.109,-0.125,-0.025,0.006]]
                      x=x, #TCGA cancer type ["BRCA","CESC","COAD","DLBC","ESCA","GBM","HNSC","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PRAD","SARC","STAD","TGCT","THCA","UCEC"]
                      y=y, #regulator_target ["correlation p value","LINC00313_miR-4429 PCC"]
                      type='heatmap',
                      zmin = -1,
                      zmid = 0,
                      zmax = 1, 
                      colorscale= list(c(-1, '#c4463a'),c(0,'#fffbbc'),c(0.9, '#c4463a'),c(1, '#3060cf')) #'YlGnBu'
                      #reversescale=TRUE
                     )
  ttHeatmap <- toJSON(testHeatmap)
  HeatmapData <- add_row(HeatmapData,rid=unique(tt$RID),heatdata=ttHeatmap)
}

write.table(HeatmapData,file="D:/Rsources/Project/DataBase/HDfregulatory/proceessingData2019_9_26/geoheatmap.txt",sep="\t",quote=F,row.names=F,col.names=F)



