library(dplyr)
library(psych)


df<-f_DF %>% select(-contains("c.14"))
#df<-df[,-which( colnames(df)=="Y" )]


cormat<-matrix(nrow = length(variables)*length(df), ncol=6)
c<-1

for (i in 1:length(variables)){
  for (j in 1:length(df)){
    cormat[c,1]<-variables[i]
    cormat[c,2]<-colnames(df)[j]
    cormat[c,3]<-abs(cor(df[,variables[i]],df[,j], method='spearman'))
    cormat[c,4]<-ifelse(cor(df[,variables[i]],df[,j], method='spearman')<0,"Negative","Positive")
    ct<-corr.test(df[,variables[i]],df[,j], method='spearman',adjust = "fdr")
    cormat[c,5]<-ifelse(ct$p.adj<.05,1,0)
    cormat[c,6]<-ifelse(cor(df[,j], Y, method = "spearman")>0,1,0)
    c<-c+1
  }
}

corDF<-as.data.frame(cormat)
corDF$V3<-as.numeric(corDF$V3)


corDF2<-corDF[corDF$V3 >.5 & corDF$V3 !=1 & corDF$V5!=0,]

corDF2<-corDF2[!duplicated(corDF2$V3),]

l<-corDF2 %>% 
  group_by(V2) %>%
  summarise(n = n())
l<-l[l$n>1,]

df<-corDF2
df<-merge(corDF2, l, by = 'V2')





df$V1<-gsub("0TerminalKeratinocyte","1.",df$V1)
df$V1<-gsub("1BasalKeratinocyte","2.",df$V1)
df$V1<-gsub("2Endothelial","3.",df$V1)
df$V1<-gsub("3SFRP2Fibroblast","4.",df$V1)
df$V1<-gsub("4RERGLPericyte","5.",df$V1)
df$V1<-gsub("5STEAP4Pericyte","6.",df$V1)
df$V1<-gsub("6AdventitialFibroblast","7.",df$V1)
df$V1<-gsub("7Macrophage","8.",df$V1)
df$V1<-gsub("8DermalSheathandPapillae","9.",df$V1)
df$V1<-gsub("9MonocyteDCNeutrophil","10.",df$V1)
df$V1<-gsub("10THelperCell","11.",df$V1)
df$V1<-gsub("11ProliferatingKeratinocyte","12.",df$V1)
df$V1<-gsub("12MastCell","13.",df$V1)
df$V1<-gsub("13SmoothMuscleCell","14.",df$V1)
df$V1<-gsub("15CA6AQP5SecretoryEpithelial","15.",df$V1)
df$V1<-gsub("16MerkelCell","16.",df$V1)
df$V1<-gsub("17CD8TCell","17.",df$V1)
df$V1<-gsub("c.","",df$V1)

df$V2<-gsub("0TerminalKeratinocyte","1.",df$V2)
df$V2<-gsub("1BasalKeratinocyte","2.",df$V2)
df$V2<-gsub("2Endothelial","3.",df$V2)
df$V2<-gsub("3SFRP2Fibroblast","4.",df$V2)
df$V2<-gsub("4RERGLPericyte","5.",df$V2)
df$V2<-gsub("5STEAP4Pericyte","6.",df$V2)
df$V2<-gsub("6AdventitialFibroblast","7.",df$V2)
df$V2<-gsub("7Macrophage","8.",df$V2)
df$V2<-gsub("8DermalSheathandPapillae","9.",df$V2)
df$V2<-gsub("9MonocyteDCNeutrophil","10.",df$V2)
df$V2<-gsub("10THelperCell","11.",df$V2)
df$V2<-gsub("11ProliferatingKeratinocyte","12.",df$V2)
df$V2<-gsub("12MastCell","13.",df$V2)
df$V2<-gsub("13SmoothMuscleCell","14.",df$V2)
df$V2<-gsub("15CA6AQP5SecretoryEpithelial","15.",df$V2)
df$V2<-gsub("16MerkelCell","16.",df$V2)
df$V2<-gsub("17CD8TCell","17.",df$V2)
df$V2<-gsub("c.","",df$V2)





write.csv(df, "corF4.csv", row.names = F, col.names = F)
