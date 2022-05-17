#early version code used to generate bulk heatmaps !!!OUT OF DATE!!!


library(ggplot2)
library(dplyr)
library(reshape)
library(forcats)
library(pheatmap)
library(RColorBrewer)
library(viridis)


#bulk<-select(Var.df, c.9MonocyteDCNeutrophilTHBS1, c.6AdventitialFibroblastCOMP, c.7MacrophageCCL2, c.12MastCellCCL2, c.13SmoothMuscleCellCCL2, c.8DermalSheathandPapillaeIGFBP3)



bulk<-read.csv("bulk.csv")[,-1]


bulk<- bulk %>% select(-contains("c.14"))  
df<-f_DF


mat<-matrix(nrow = length(bulk), ncol = length(variables))

for (i in 1:length(bulk)){
  for (j in 1:length(variables)){
    mat[i,j]<-cor(bulk[,i],df[,variables[j]])
  }
}


rownames(mat)<-colnames(bulk)
colnames(mat)<-variables
mat2<-mat[apply(mat,1,max)>.5,]


bulk_mrss<-c()
for (i in 1:length(bulk)){
  bulk_mrss<-c(bulk_mrss,cor(bulk[,i],Y))
}








paletteLength <- 50
myColor <- colorRampPalette(c("#B600D9","#C73EE2","white","white","#2ECE76", "#00c458"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths




dfPlot<-as.data.frame(mat2)


dfPlot<-dfPlot[rowSums(is.na(dfPlot))==0,]

colnames(dfPlot)<-variables
rownames(dfPlot)<-colnames(bulk)



corList<-matrix(nrow = length(variables), ncol=2)
for (i in 1:length(variables)){
  corList[i,1]<-variables[i]
  corList[i,2]<-cor(df[,colnames(df)%in% variables[i]],Y,method = "spearman")
}
corList<-as.data.frame(corList)
ord<-corList$V1[order(corList$V2, decreasing = T)]
dfPlot <- subset(dfPlot, select=ord)



colnames(dfPlot)
colnames(dfPlot)<-gsub("0TerminalKeratinocyte","1.",colnames(dfPlot))
colnames(dfPlot)<-gsub("1BasalKeratinocyte","2.",colnames(dfPlot))
colnames(dfPlot)<-gsub("2Endothelial","3.",colnames(dfPlot))
colnames(dfPlot)<-gsub("3SFRP2Fibroblast","4.",colnames(dfPlot))
colnames(dfPlot)<-gsub("4RERGLPericyte","5.",colnames(dfPlot))
colnames(dfPlot)<-gsub("5STEAP4Pericyte","6.",colnames(dfPlot))
colnames(dfPlot)<-gsub("6AdventitialFibroblast","7.",colnames(dfPlot))
colnames(dfPlot)<-gsub("7Macrophage","8.",colnames(dfPlot))
colnames(dfPlot)<-gsub("8DermalSheathandPapillae","9.",colnames(dfPlot))
colnames(dfPlot)<-gsub("9MonocyteDCNeutrophil","10.",colnames(dfPlot))
colnames(dfPlot)<-gsub("10THelperCell","11.",colnames(dfPlot))
colnames(dfPlot)<-gsub("11ProliferatingKeratinocyte","12.",colnames(dfPlot))
colnames(dfPlot)<-gsub("12MastCell","13.",colnames(dfPlot))
colnames(dfPlot)<-gsub("13SmoothMuscleCell","14.",colnames(dfPlot))
colnames(dfPlot)<-gsub("15CA6AQP5SecretoryEpithelial","15.",colnames(dfPlot))
colnames(dfPlot)<-gsub("16MerkelCell","16.",colnames(dfPlot))
colnames(dfPlot)<-gsub("17CD8TCell","17.",colnames(dfPlot))
colnames(dfPlot)<-gsub("c.","",colnames(dfPlot))
dfPlotT<-t(dfPlot)



rownames(dfPlot)
rownames(dfPlot)<-gsub("0TerminalKeratinocyte","1.",rownames(dfPlot))
rownames(dfPlot)<-gsub("1BasalKeratinocyte","2.",rownames(dfPlot))
rownames(dfPlot)<-gsub("2Endothelial","3.",rownames(dfPlot))
rownames(dfPlot)<-gsub("3SFRP2Fibroblast","4.",rownames(dfPlot))
rownames(dfPlot)<-gsub("4RERGLPericyte","5.",rownames(dfPlot))
rownames(dfPlot)<-gsub("5STEAP4Pericyte","6.",rownames(dfPlot))
rownames(dfPlot)<-gsub("6AdventitialFibroblast","7.",rownames(dfPlot))
rownames(dfPlot)<-gsub("7Macrophage","8.",rownames(dfPlot))
rownames(dfPlot)<-gsub("8DermalSheathandPapillae","9.",rownames(dfPlot))
rownames(dfPlot)<-gsub("9MonocyteDCNeutrophil","10.",rownames(dfPlot))
rownames(dfPlot)<-gsub("10THelperCell","11.",rownames(dfPlot))
rownames(dfPlot)<-gsub("11ProliferatingKeratinocyte","12.",rownames(dfPlot))
rownames(dfPlot)<-gsub("12MastCell","13.",rownames(dfPlot))
rownames(dfPlot)<-gsub("13SmoothMuscleCell","14.",rownames(dfPlot))
rownames(dfPlot)<-gsub("15CA6AQP5SecretoryEpithelial","15.",rownames(dfPlot))
rownames(dfPlot)<-gsub("16MerkelCell","16.",rownames(dfPlot))
rownames(dfPlot)<-gsub("17CD8TCell","17.",rownames(dfPlot))
rownames(dfPlot)<-gsub("c.","",rownames(dfPlot))
rownames(dfPlot)<-gsub("\\..",".",rownames(dfPlot))
dfPlotT<-t(dfPlot)





myBreaks <- c(seq(-.7, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(0, .7, length.out=floor(paletteLength/2)+2)[-1])

p<-pheatmap(dfPlotT,cluster_rows=F,cluster_cols = F, color = myColor,breaks=myBreaks,cellwidth = 150/length(bulk),cellheight = 200/length(variables),fontsize_row = 12, fontsize_col = 12)

p<-pheatmap(dfPlotT,cluster_rows=F,cluster_cols = F, color = myColor,breaks=myBreaks)


save_pheatmap_pdf <- function(x, filename, width=15, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p, "bulkH_F.pdf")





