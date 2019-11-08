#'--------------------------------------------------
#'
#' @author : shijian
#' @details :批量下载GEO数据              
#'
#'--------------------------------------------------
## 读入GSE以及TCGA大表
Data.accession <- read.table("D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/data.accession.new.txt",sep="\t",header=T,fill=T,quote=NULL,stringsAsFactors=F) 
head(Data.accession)
Data.accession <- Data.accession[grep("GSE|TCGA",Data.accession$Data.accession),]
write.table(Data.accession,file="D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/data.accession.new1.txt",sep="\t",quote=F,row.names=F,col.names=T)

library("GEOquery") #library(help="GEOquery")
## 获取GSE.soft.gz,得到对于数据的描述信息gse    head(Meta(gse))
DownloadGSE <- function(GSEid,RID,filepath){  #,Data.accession
  #' @param  filepath 头地址  e.g.  "D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSETCGA/"
  #' @param  GSEid GSEid号 建目录
  #' @param  RID  互作的id号，用于标记文件名
  #' @param  Data.accession 整理好的GSE大表
  #' @return 整合路径存储地址的大表 规律: 头地址+ GSEid目录 + 文件名(rid-GSEid) 匹配到Data.accession大表上面
  library("GEOquery")
  if(grepl("GSE",GSEid)){   #如果是GSE开头
    gse <- getGEO(GSEid,GSEMatrix = F) #获取soft信息
    if(Meta(gse)$type == "Expression profiling by array" || Meta(gse)$type == "Non-coding RNA profiling by array" || grepl("array",Meta(gse)$type)) # 判断是芯片数据还是测序数据   
    { 
      # 如果是芯片数据并且存在
      cat(paste(GSEid," is microarray"))
      dir.name <- paste0(filepath,GSEid)
      if( !file.exists(dir.name) ){
        dir.create( dir.name,recursive = T )
        gset <- getGEO( GSEid,getGPL = T,destdir="D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSEmicroarray") #下载数据
        options( stringsAsFactors = F )
        gset = gset[[1]]
        exprSet = exprs( gset ) ## “GEOquery”包中的exprs函数用来取出表达矩阵
        pdata = pData( gset ) ## “GEOquery”包中的pData函数用来取出样本信息
        GPL = gset@featureData@data  ##探针注释信息
        temp <- c()
        if( dim(exprSet)[1] > 0 ){
          expression.filename <- paste0( RID, "-", GSEid, "-expression.MATRIX.txt" )
          write.table(exprSet,file=paste0( filepath, GSEid,"/" , expression.filename), sep="\t", quote=F, row.names=T, col.names=T)
          temp <- append(temp,paste0( filepath, GSEid,"/" , expression.filename))
        }
        if( dim(pdata)[1] > 0 ){
          pData.filename <- paste0( RID, "-", GSEid, "-pData.txt" )
          write.table(pdata,file=paste0( filepath, GSEid,"/" , pData.filename), sep="\t", quote=F, row.names=T, col.names=T)
          temp <- append(temp,paste0( filepath, GSEid,"/" , pData.filename))
        }
        if( dim(GPL)[1] > 0 ){
          GPL.filename <- paste0( RID, "-", GSEid, "-ProbeAnnotation.txt" )
          write.table(GPL,file=paste0( filepath, GSEid,"/" , GPL.filename ), sep="\t", quote=F, row.names=T, col.names=T)
          temp <- append(temp,paste0( filepath, GSEid,"/" , GPL.filename))
        }
        #dataSite <- paste(temp,collapse=";")
      }
    }else if( grepl("sequencing",Meta(gse)$type) & grepl("Expression profiling",Meta(gse)$type) ){ ##是否是测序数据
      ## 首先判断是否提供处理好的数据
      gset <- getGEO( GSEid,getGPL = T,destdir="D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSEmicroarray") 
      options( stringsAsFactors = F )
      gset = gset[[1]]
      exprSet = exprs( gset ) ## “GEOquery”包中的exprs函数用来取出series matrix矩阵
      if( dim(exprSet)[1] > 0 ){
        expression.filename <- paste0( RID, "-", GSEid, "-expression.MATRIX.txt" )
        if( !file.exists(dir.name) ){
          dir.name <- paste0(filepath,GSEid)
          dir.create(dir.name,recursive = T)
        }
        write.table(exprSet,file=paste0( filepath, GSEid,"/" , expression.filename), sep="\t", quote=F, row.names=T, col.names=T)
        #dataSite <- paste0( filepath, GSEid,"/" , expression.filename)
        cat( paste0(GSEid," processed RNA-seq raw data") )
      }else{
        cat( paste0(GSEid," only RNA-seq raw data") ) ## RNA-seq原始数据就不下载了
        #dataSite <- paste0(GSEid," only RNA-seq raw data") 
      }
    }
  }else{ ## TCGA 数据
      cat(" TCGA,手动注释")
      #dataSite <- "TCGA"
  }
  #return( dataSite )
}

DownloadGSE <- function(GSEid,filepath){  #,Data.accession
  #' @param  filepath 头地址  e.g.  "D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSETCGA/"
  #' @param  GSEid GSEid号 建目录
  #' @return 整合路径存储地址的大表 规律: 头地址+ GSEid目录 + 文件名(rid-GSEid) 匹配到Data.accession大表上面
  library("GEOquery")
  if(grepl("GSE",GSEid)){   #如果是GSE开头
    gse <- getGEO(GSEid,GSEMatrix = F) #获取soft信息
    if(Meta(gse)$type == "Expression profiling by array" || Meta(gse)$type == "Non-coding RNA profiling by array" || grepl("array",Meta(gse)$type)) # 判断是芯片数据还是测序数据   
    { 
      # 如果是芯片数据并且存在
      cat(paste(GSEid," is microarray"))
      dir.name <- paste0(filepath,GSEid)
      if( !file.exists(dir.name) ){
        dir.create( dir.name,recursive = T )
        gset <- getGEO( GSEid,getGPL = T,destdir="D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSEmicroarray") #下载数据
        options( stringsAsFactors = F )
        gset = gset[[1]]
        exprSet = exprs( gset ) ## “GEOquery”包中的exprs函数用来取出表达矩阵
        pdata = pData( gset ) ## “GEOquery”包中的pData函数用来取出样本信息
        GPL = gset@featureData@data  ##探针注释信息
        temp <- c()
        if( dim(exprSet)[1] > 0 ){
          expression.filename <- paste0( GSEid, "-expression.MATRIX.txt" )
          write.table(exprSet,file=paste0( filepath, GSEid,"/" , expression.filename), sep="\t", quote=F, row.names=T, col.names=T)
          temp <- append(temp,paste0( filepath, GSEid,"/" , expression.filename))
        }
        if( dim(pdata)[1] > 0 ){
          pData.filename <- paste0( GSEid, "-pData.txt" )
          write.table(pdata,file=paste0( filepath, GSEid,"/" , pData.filename), sep="\t", quote=F, row.names=T, col.names=T)
          temp <- append(temp,paste0( filepath, GSEid,"/" , pData.filename))
        }
        if( dim(GPL)[1] > 0 ){
          GPL.filename <- paste0( GSEid, "-ProbeAnnotation.txt" )
          write.table(GPL,file=paste0( filepath, GSEid,"/" , GPL.filename ), sep="\t", quote=F, row.names=T, col.names=T)
          temp <- append(temp,paste0( filepath, GSEid,"/" , GPL.filename))
        }
        #dataSite <- paste(temp,collapse=";")
      }
    }else if( grepl("sequencing",Meta(gse)$type) & grepl("Expression profiling",Meta(gse)$type) ){ ##是否是测序数据
      ## 首先判断是否提供处理好的数据
      gset <- getGEO( GSEid,getGPL = T,destdir="D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSEmicroarray") 
      options( stringsAsFactors = F )
      gset = gset[[1]]
      exprSet = exprs( gset ) ## “GEOquery”包中的exprs函数用来取出series matrix矩阵
      if( dim(exprSet)[1] > 0 ){
        expression.filename <- paste0( GSEid, "-expression.MATRIX.txt" )
        if( !file.exists(dir.name) ){
          dir.name <- paste0(filepath,GSEid)
          dir.create(dir.name,recursive = T)
        }
        write.table(exprSet,file=paste0( filepath, GSEid,"/" , expression.filename), sep="\t", quote=F, row.names=T, col.names=T)
        #dataSite <- paste0( filepath, GSEid,"/" , expression.filename)
        cat( paste0(GSEid," processed RNA-seq raw data") )
      }else{
        cat( paste0(GSEid," only RNA-seq raw data") ) ## RNA-seq原始数据就不下载了
        #dataSite <- paste0(GSEid," only RNA-seq raw data") 
      }
    }
  }else{ ## TCGA 数据
    cat(" TCGA,手动注释")
    #dataSite <- "TCGA"
  }
  #return( dataSite )
}

sink(file="D:/Rsources/Project/DataBase/HDfregulatory/log.txt")
for(i in 1:length(Data.accession$RID)){
  print(i)
  DownloadGSE( Data.accession$Data.accession[i], Data.accession$RID[i], "D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSETCGA/microarray/" )
  Data.accession$DataSite[i] <- dataSite
}
sink()

sink(file="D:/Rsources/Project/DataBase/HDfregulatory/log1.txt")
for(i in 123:length(Data.accession$RID)){
  print(i)
  DownloadGSE( Data.accession$Data.accession[i], Data.accession$RID[i], "D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSETCGA/microarray/" )
  #Data.accession$DataSite[i] <- dataSite
}
sink()

sink(file="D:/Rsources/Project/DataBase/HDfregulatory/log2.txt")
for(i in 192:length(Data.accession$RID)){
  print(i)
  DownloadGSE( Data.accession$Data.accession[i], Data.accession$RID[i], "D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSETCGA/microarray/" )
  #Data.accession$DataSite[i] <- dataSite
}
sink()

head(Data.accession)

## 所有待下载的GSE
allfiles <- unique( Data.accession$Data.accession[grep("GSE",Data.accession$Data.accession)] )
## 已经下载的GSE
files <- list.files( "D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSETCGA/microarray" )
## 需要补充下载的GSE
LeftGSE <- allfiles[!(allfiles %in% files)]
sink( file="D:/Rsources/Project/DataBase/HDfregulatory/log3.txt" )
for( i in 1:length(LeftGSE) ){
  print(i)
  DownloadGSE( LeftGSE[i], "D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSETCGA/microarray/" )
  #Data.accession$DataSite[i] <- dataSite
}
sink()







