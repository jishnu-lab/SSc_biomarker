library(Seurat)
library(dplyr)



## LOAD IN DATA

load("harmonskinV6_full.rdata")

response<-read.csv("SkinScore_LafyatisAndBaselineTofa_MRSS.csv", fileEncoding = 'UTF-8-BOM')

ID<-response$orig.ident




## FILTER CLUSTERS
clusters<-new_levels[-c(19:35)]

harmonskin.top20 <- harmonskin.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)






## GET TOP 20 DEGs PER CLUSTER



# make empty matrix
topDEGs<-matrix(nrow = 24, ncol = 360)

# iterate through patients
for (i in 1:24){
  vec<-c()
  
  # iterate through clusters
  for (j in 1:18){
    # get a list of top 20 DEGs
    tmpList<-harmonskin.top20$gene[harmonskin.top20$cluster==new_levels[j]]
    # get seurat object from list
    tmp<-subset(harmonskin, idents = new_levels[j], features=tmpList, cells = WhichCells(subset(x = harmonskin, subset = orig.ident == ID[i])))
    # average across cluster
    cluster.averages <- AverageExpression(tmp, return.seurat = TRUE)
    #make df
    tmp_df<-as.data.frame(cluster.averages@assays[["RNA"]]@data)
    # add to vector
    vec<-c(vec,tmp_df[,1])
  }
  topDEGs[i,]<-vec
}

DEG.df<-data.frame(topDEGs)



## NAMING GENES
names<-c()
clstrs<-substr(harmonskin.top20$cluster, 1,2)
names<-paste(clstrs,harmonskin.top20$gene,sep=".")
names<-names[-c(361:700)]

names<-gsub(" ", "", names)
names<-gsub("-", "", names)
names<-gsub("/", "", names)
names<- paste("c", names, sep=".")
colnames(DEG.df)<-names







## GET TOP 20 HIGH VARIANCE GENES PER CLUSTER


# define function for row variance
RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}



varList<-c()
dfVar<-NULL
for (i in 1:18){
  matS<-matrix(nrow=33538,ncol=24)
  for (j in 1:1){
    vec<-c()
    #first isolate patient group
    tmp<-subset(harmonskin, idents = clusters[i] ,subset = orig.ident == ID[j])
    
    cluster.averages.temp <- AverageExpression(tmp, return.seurat = TRUE)
    tmp_df<-as.data.frame(cluster.averages.temp@assays[["RNA"]]@data)
    
    matS[,j]<-tmp_df[,1]
  }
  dfBIG<-as.data.frame(matS)
  # repeat gene :/
  dfBIG<-dfBIG[-22664,]
  
  names<-c()
  clstrs<-substr(clusters[i], 1,2)
  names<-paste(clstrs, rownames(tmp_df),sep=".")
  names<-gsub(" ", "", names)
  names<-gsub("-", "", names)
  names<-gsub("/", "", names)
  names<- paste("c", names, sep=".")
  names<-names[-22664]
  rownames(dfBIG)<-names
  
  dfBIG$rowVar<- RowVar(dfBIG)
  
  dfBIG<-dfBIG[order(dfBIG$rowVar,decreasing=T),]
  
  varOrder<-head(dfBIG,20)
  if (i==1){
    dfVar<-varOrder
  } else {
    dfVar<-rbind(dfVar, varOrder)
  }
  varList<-c(varList, row.names(varOrder))
  
}

Var.df<-as.data.frame(t(dfVar[,-25]))

write.csv(dfVar2, "Var.df.csv")







## SUBSET TOP 20 GENES PER CLUSTER USING SICK VS HEALTHY CONTROL LOG FOLD CHANGE

SickvsHC<-matrix(nrow = 24, ncol = 360)
for (i in 1:24){
  vec<-c()
  names<-c()
  for (j in 0:17){
    clusterN<-paste("cluster", j,".csv", sep = "")
    tmp_df<-read.csv(clusterN)
    tmp_df<-tmp_df[order(tmp_df$avg_log2FC,decreasing=T),]
    tmpList<-head(tmp_df$X,20)
    tmp<-subset(harmonskin, idents = new_levels[j+1], features=tmpList, cells = WhichCells(subset(x = harmonskin, subset = orig.ident == ID[i])))
    cluster.averages <- AverageExpression(tmp, return.seurat = TRUE)
    tmp_df<-as.data.frame(cluster.averages@assays[["RNA"]]@data)
    vec<-c(vec,tmp_df[,1])
    tmpList<-paste(new_levels[j+1],tmpList, sep=".")
    tmpList<-gsub(" ", "", tmpList)
    names<-c(names, tmpList)
  }
  SickvsHC[i,]<-vec
}

C.vs.Scl<-data.frame(SickvsHC)
colnames(C.vs.Scl)<-names
names2<-gsub("-", ".", names)
colnames(C.vs.Scl)<-names2

write.csv(C.vs.Scl, "C.vs.Scl.csv")










## TOP 300 GENES PER CELL GROUP

keratinocytes<-new_levels[c(1,2,12)]
fibro<-new_levels[c(4,7)]
mye<-new_levels[c(8,10)]

RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

grouping<-mye
varList<-c()
dfVar<-NULL
for (i in grouping){
  matS<-matrix(nrow=33538,ncol=24)
  names<-c()
  for (j in 1:24){
    vec<-c()
    tmp<-subset(harmonskin, idents = i ,subset = orig.ident == ID[j])
    
    cluster.averages.temp <- AverageExpression(tmp, return.seurat = TRUE)
    tmp_df<-as.data.frame(cluster.averages.temp@assays[["RNA"]]@data)
    
    matS[,j]<-tmp_df[,1]
  }
  dfBIG<-as.data.frame(matS)
  dfBIG<-dfBIG[-22664,]
  
  names<-paste(i, rownames(tmp_df),sep=".")
  names<-gsub(" ", "", names)
  names<-gsub("-", "", names)
  names<-gsub("/", "", names)
  names<- paste("c", names, sep=".")
  names<-names[-22664]
  rownames(dfBIG)<-names
  
  dfBIG$rowVar<- RowVar(dfBIG)
  
  dfBIG<-dfBIG[order(dfBIG$rowVar,decreasing=T),]
  
  varOrder<-head(dfBIG,300/length(grouping))
  if (i==1){
    dfVar<-varOrder
  } else {
    dfVar<-rbind(dfVar, varOrder)
  }
  varList<-c(varList, row.names(varOrder))
  
}

m_DF<-as.data.frame(t(dfVar[,-25]))


write.csv(m_DF,"mye.csv")




















## SUBSET PREVIOUSLY IDENTIFIED GENES FROM BULK


bulk<-c("CD14","SERPINE1","IL13RA1","CTGF","OSMR","THBS1", "COMP", "SIGLEC1", "IFI44", "CCL2", "IGFBP3", "ADAM12","THY1", "WIF1", "MS4A4A", "CD163")




bulkM<-matrix(nrow = 24, ncol = length(bulk)*18)
for (i in 1:24){
  vec<-c()
  for (j in 1:18){
    tmp<-subset(harmonskin, idents = new_levels[j], features=bulk, cells = WhichCells(subset(x = harmonskin, subset = orig.ident == ID[i])))
    cluster.averages <- AverageExpression(tmp, return.seurat = TRUE)
    tmp_df<-as.data.frame(cluster.averages@assays[["RNA"]]@data)
    vec<-c(vec,tmp_df[,1])
  }
  bulkM[i,]<-vec
}

Bulk.df<-data.frame(bulkM)



bulk2<-rep(bulk,18)





names<-c()
for (i in 1:18){
  names<-c(names, paste(new_levels[i], bulk, sep = "."))
  names<-gsub(" ", "", names)
  names<-gsub("-", "", names)
  names<-gsub("/", "", names)
}
names<- paste("c", names, sep=".")

colnames(Bulk.df)<-names

write.csv(Bulk.df,"bulk.csv")










## GENERATE DOT PLOTS

harmonskin.top2 <- harmonskin.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
harmonskin.top2<-harmonskin.top2[harmonskin.top2$cluster %in% clusters[-15],]
harmonskin.top2 <- harmonskin.top2[harmonskin.top2$gene %in% unique(harmonskin.top2$gene),]


DotPlot(harmonskin, features = unique(harmonskin.top2$gene),idents = clusters[-15]) + RotatedAxis()




varF<-c("IFI27", "SAA1", "CXCL3", "MT1X", "MT2A", "KRT17", "S100A9", "IGFBP5")
DotPlot(harmonskin, features = varF,idents = clusters[-15]) + RotatedAxis()


DEGF<-c("KRT6A", "KRT5", "CCL19", "RGS13", "HLA-DPB1", "TINAGL1", "PLN")
DotPlot(harmonskin, features = DEGF,idents = clusters[-15]) + RotatedAxis()








