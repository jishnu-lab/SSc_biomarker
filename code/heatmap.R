library(ggplot2)
library(dplyr)
library(reshape)
library(forcats)
library(pheatmap)
library(RColorBrewer)
library(viridis)

df<-k_DF %>% select(-contains("c.14"))

resp<-read.csv("SkinScore_LafyatisAndBaselineTofa_MRSS.csv",header=T)
Y <- resp[,2]
ID<-resp[,1]
rownames(df)<-ID

dfPlot<-df[,colnames(df) %in% variables]

dfPlot$ID<-ID
dfPlot$ID <- factor(dfPlot$ID, levels = dfPlot$ID[order(Y)])
dfPlot<-dfPlot[order(dfPlot$ID),]



corList<-matrix(nrow = length(variables), ncol=2)
for (i in 1:length(variables)){
  corList[i,1]<-variables[i]
  corList[i,2]<-cor(df[,colnames(df)%in% variables[i]],Y,method = "spearman")
}
corList<-as.data.frame(corList)
ord<-corList$V1[order(corList$V2, decreasing = T)]
dfPlot <- subset(dfPlot, select=ord)




paletteLength <- 50
myColor <- colorRampPalette(c("#52adeb","white", "#fc2771"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths




dfPlot<-scale(dfPlot)


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
colnames(dfPlot)<-gsub("3","4",colnames(dfPlot))
colnames(dfPlot)<-gsub("6","7",colnames(dfPlot))
#colnames(dfPlot)<-gsub("\\..","\\.",colnames(dfPlot))
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
rownames(dfPlot)<-gsub("3","4",rownames(dfPlot))
rownames(dfPlot)<-gsub("6","7",rownames(dfPlot))
#rownames(dfPlot)<-gsub("\\..","\\.",rownames(dfPlot))
dfPlotT<-t(dfPlot)






myBreaks <- c(seq(min(dfPlot), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(dfPlot)/paletteLength, max(dfPlot), length.out=floor(paletteLength/2)))

p<-pheatmap(dfPlotT,cluster_rows=F,cluster_cols = F, color = myColor,breaks=myBreaks)








save_pheatmap_pdf <- function(x, filename, width=15/2, height=7/2) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p, "HK.pdf")









save_pheatmap_png <- function(x, filename, width=1500, height=700, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(p, "HDEG.png")



















dfM<-as.matrix(scale(df))

# Transform the matrix in long format
dfPlot <- melt(dfM)
colnames(dfPlot) <- c("ID", "y", "value")
dfPlot$ID<-ID


dfPlot<-dfPlot[dfPlot$y %in% variables,]



corList<-matrix(nrow = length(variables), ncol=2)
for (i in 1:length(variables)){
  corList[i,1]<-variables[i]
  corList[i,2]<-cor(df[,colnames(df)%in% variables[i]],Y)
}
corList<-as.data.frame(corList)
ORDER<-corList$V1[order(corList$V2, decreasing = F)]


dfPlot$ID <- factor(dfPlot$ID, levels = dfPlot$ID[order(Y)])

dfPlot$y <- factor(dfPlot$y, levels = ORDER)
dfPlotm<-as.matrix(dfPlot)

pheatmap(dfPlot)

ggplot(dfPlot,aes(x = ID, y = y, fill = value))+geom_tile()
