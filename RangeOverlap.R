#'------------------------------------------------------------------------
#'
#'    @description:基于lncRNA的区间信息，计算与已有注释系统的overlap信息
#'  
#'------------------------------------------------------------------------

## 导入计算overlap函数
source("/pub6/Temp/XAM/Rcode/General/BioInformatics/SequencingBasis/Ranges/BasicOperation/rangeOverlap.R")

## test.data
lncRNA.bed <- data.frame(chr=character(0),start=numeric(0),end=numeric(0),strand=character(0),name=character(0),stringsAsFactors = F)
annotation.bed <- data.frame(chr=character(0),start=numeric(0),end=numeric(0),strand=character(0),name=character(0),stringsAsFactors = F)
chr <- c("chr1","chr2","chr3","chr4","chr1","chr2","chr3","chr4")
start <- c(12345,23456,34567,45678,12340,23450,34560,45670)
end <- c(34567,45678,56789,67890,34560,45670,56780,67880)
strand <- c("*","*","*","*","*","*","*","*")
name <- c("HOTAIR","MEG3","MALAT1","GAS5","CERNA2","ZXF2","HGTY","FGSE")
lncRNA.bed <- cbind.data.frame(chr,start,end,strand,name)
annotation.bed <- lncRNA.bed
{
  write.table(lncRNA.bed,file="/pub6/Temp/sj/lncRNA.bed",sep="\t",quote=F,row.names=F,col.names=F)
  write.table(annotation.bed,file="/pub6/Temp/sj/annotation.bed",sep="\t",quote=F,row.names=F,col.names=F)
  system("/pub5/xiaoyun/BioSoftware/bedtools2/bin/sortBed  -i /pub6/Temp/sj/lncRNA.bed > /pub6/Temp/sj/lncRNA.sort.bed") 
  system("/pub5/xiaoyun/BioSoftware/bedtools2/bin/sortBed  -i /pub6/Temp/sj/annotation.bed > /pub6/Temp/sj/annotation.sort.bed") 
  system("/pub5/xiaoyun/BioSoftware/bedtools2/bin/intersectBed -a /pub6/Temp/sj/lncRNA.sort.bed -b /pub6/Temp/sj/annotation.sort.bed -wa -wb  > /pub6/Temp/sj/result.bed")
}

overlap.result <- RegionOverlapping.gr(lncRNA.bed,annotation.bed)
head(overlap.result)
## 删除自己和自己overlap的情况
pos <- which(overlap.result$Qindex == overlap.result$Sindex)
overlap.result <-  overlap.result[-pos,]
## 筛选overlap > 90% 的情况，如果lncRNA A的90%区间都与lncRNA B重合，那么lncRNA A就应该去掉，如果两个都大于90%，那么两个随便取一个lncRNA
pos_all <- which(overlap.result$OLpercQ  >= 90 & overlap.result$OLpercS >= 90)
anyName <- overlap.result[pos_all,c(6,20)] ##哪个名字都行
pos_Q <- which(overlap.result$OLpercQ  >= 90 & overlap.result$OLpercS < 90)
SName <- overlap.result[pos_Q,c(6,20)] ##使用第20列基因名字
pos_S <- which(overlap.result$OLpercQ  < 90 & overlap.result$OLpercS >= 90)
QName <- overlap.result[pos_S,c(20,6)] ##使用第6列基因名字
result <- list(any=anyname,SName=SName,QName=QName)
return(result)

OverLap.fun <- function(lncRNA.bed,annotation.bed){
  #' @param: lncRNA.bed  chr start end strand name
  #' @param: annotation.bed chr start end strand name
  #' @return : 返回一个列表，列表里有两个数据框，每个数据框都只有两列基因
  overlap.result <- RegionOverlapping.gr(lncRNA.bed,annotation.bed)
  pos <- which(overlap.result$Qindex == overlap.result$Sindex)
  overlap.result <-  overlap.result[-pos,]
  ## 筛选overlap > 90% 的情况，如果lncRNA A的90%区间都与lncRNA B重合，那么lncRNA A就应该去掉，如果两个都大于90%，那么两个随便取一个lncRNA
  pos_all <- which(overlap.result$OLpercQ  >= 90 & overlap.result$OLpercS >= 90)
  anyName <- overlap.result[pos_all,c(6,20)] ##哪个名字都行
  pos_Q <- which(overlap.result$OLpercQ  >= 90 & overlap.result$OLpercS < 90)
  SName <- overlap.result[pos_Q,c(6,20)] ##使用第20列基因名字
  pos_S <- which(overlap.result$OLpercQ  < 90 & overlap.result$OLpercS >= 90)
  QName <- overlap.result[pos_S,c(20,6)] ##使用第6列基因名字
  oneName <- rbind.data.frame(SName,QName)
  result <- list(any=anyname,oneName=oneName)
  return(result)
}




