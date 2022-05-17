library(glmnet)
library(pROC)
library(dplyr)
library(randomForest)
library(caret)
library(pls)
library(ggplot2)
library(e1071)
library(foreach)
library(doParallel)

registerDoParallel()
getDoParWorkers()

system.time( foreach(i=1:10000) %dopar% sum(tanh(1:i)) )



df<-k_DF %>% select(-contains("c.14"))


#df<-df[,-which( colnames(df)=="Y" )]

resp<-read.csv("SkinScore_LafyatisAndBaselineTofa_MRSS.csv",header=T)
ID<-resp$orig.ident

Y <- resp[,2]
df<-as.matrix(df)
variables<-c()

set.seed(405)#101 405 for M,V,K
lasso.cheat.cv <- cv.glmnet(df, Y, alpha=1, nfolds = 10)

lam<-lasso.cheat.cv$lambda.min
lam<- lam #can be changed

lasso.cheat<-glmnet(df, Y, alpha=1, lambda = lam)

coefs<-coef(lasso.cheat)

notNull<-which(coefs!=0)

  #remove intercept
variables <- tail(row.names(coefs)[notNull], -1)


variables<-variables[abs(cor(df[,variables],Y, method = "spearman"))>.1]


print(variables)

#variables<-c(variables,newVariables)


#freq = sort(table(variables),decreasing=TRUE)
#print(freq)

#stable_Var = names(which((freq)>=0.5))

vars <- paste(variables, collapse="+")
 
form_cheat<-as.formula(paste("Y ~ ",vars,sep = ""))







#PLS high low
df<-as.data.frame(df)
scl.DEG.var<-as.data.frame(scale(df))
MRSS_highlow<-ifelse(Y>median(Y),1,0)


formPLS<-as.formula(paste("MRSS_highlow ~ ",vars,sep = ""))
#formPLS<-as.formula(paste("MRSS_highlow ~ ",stable_Var,sep = ""))

modPLS<-plsr(formPLS, data=df, scale=T, validation="CV")


plsDF<-data.frame(modPLS$scores[,1],modPLS$score[,2],ifelse(as.factor(MRSS_highlow)==0,"Low","High"))

colnames(plsDF)<-c("Component1", "Component2","MRSS")

pls<-ggplot(plsDF, aes(Component1, Component2))+geom_point(aes(colour=MRSS))+xlab("Scores on LV1")+ylab("Scores on LV 2")+
  scale_color_manual(values = c("#00D962","#B600D9"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("plsDEG.pdf",device = "pdf")















modPLS<-plsr(formPLS, data=Var.df, scale=TRUE, valdiation="CV")

scoreplot(modPLS, pch=19, col=pal(Y), comps = 1:2) 
scoreplot(modPLS, pch=19, col=(c("darkred","darkblue")), comps = 1:2) 
loadingplot(modPLS, pch=19, col=(c("darkred","darkblue")), comps = 1:2)
summary(modPLS)



#PLS Full


scl.DEG.var<-as.data.frame(scale(DEG.Var.df))
modPLSFull<-plsr(form_cheat, data=scl.DEG.var, scale=TRUE, valdiation="CV")

pal = colorRampPalette(c("white","darkblue"))
pal = scale_color_viridis(discrete = TRUE, option = "D")

scoreplot(modPLS, pch=19, col=pal(Y), comps = 1:2)
scoreplot(modPLSFull, pch=19, col=(c("lightblue","darkblue")), comps = 1:2)

hist(Y,col = c('pink','violet','purple'))























df.cheat<-as.data.frame(df[,variables])
df.cheat$Y<-Y

corList<-matrix(nrow = length(df.cheat)-1, ncol=2)
for (i in 1:length(df.cheat)-1){
  corList[i,1]<-colnames(df.cheat)[i]
  corList[i,2]<-cor(df.cheat[,i],df.cheat$Y,method = "spearman")
}

tmp_df<-as.data.frame(corList)
colnames(tmp_df)<-c("variablesSaved", "cor")

#tmp_df2<-merge(dftab, tmp_df, by= "variablesSaved", sort=F)

write.csv(tmp_df, "freq.csv")
