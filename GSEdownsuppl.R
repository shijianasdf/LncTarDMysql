#'------------------------------------
#'   批量下载GSE的查缺补漏
#'
#'------------------------------------
"GSE65801"
"GSE61270"
library(GEOquery)
?getGEO
getGEO(filename ="D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSEmicroarray/GSE65801_series_matrix.txt.gz",getGPL=T)
gset <- getGEO( "GSE61270",getGPL = T,destdir="D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSEmicroarray") #下载数据

list.files("D:/Rsources/Project/DataBase/HDfregulatory/GSEdownload/GSETCGA/microarray/GSE4290")


