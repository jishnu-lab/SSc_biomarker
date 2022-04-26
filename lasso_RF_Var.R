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
library(RColorBrewer)

registerDoParallel()
getDoParWorkers()

system.time( foreach(i=1:10000) %dopar% sum(tanh(1:i)) )

Var.df<-read.csv("varDF.csv")[,-1]
df<-Var.df %>% select(-contains("c.14"))
#df<-Var.df %>% select(contains("IGFBP5")|contains("c.14"))

resp<-read.csv("SkinScore_LafyatisAndBaselineTofa_MRSS.csv",header=T)
Y <- resp[,2]
df$Y<-Y
set.seed(101)


#df$test<-df$c.0TerminalKeratinocyteLOR+6
#remove redundant

cormat<-cor(df[,1:length(df)-1])

highlyCorrelated <- findCorrelation(cormat, cutoff=0.9)

df<-df[,-highlyCorrelated]

#FOLDS
k <- 24



#CREATE VECTORS FOR LOOP
variables<-c()
variablesSaved<-c()
trueY<-c()
pred<-c()
trueY<-c()
RFpred<-c()
rf.MSE.Var<-c()
rf.R.Var<-c()
svm.R<-c()
trueYsave<-c()
RFpredsave.Var<-c()
svm.pred<-c()
RFpredperm<-c()
rf.R.Var.perm<-c()
RFpredsavePerm.Var<-c()
svm.predsave<-c()
#rfFMat<-matrix(nrow=20, ncol=10)


#OUTER LOOP
for (i in 1:10){
  #MAKE K FOLDS
  folds <- createFolds(y = df$Y, k=k, list = FALSE, returnTrain = FALSE)
  #INNER LOOP
  for (NoF in 1:k){  
    fold = which(folds == NoF)
    
    #SPLIT TRAIN TEST
    train <- df[-fold, ]
    test <- df[fold, ]
    #set.seed(405)
    
    
    testY <- test[,ncol(df)]
    
    #NOTE TRUE Y
    trueY<-c(trueY,testY)
    
    
    #train<-train[,abs(cor(train[,colnames(train)],train$Y,method = "spearman"))>.1]
    #stable_Var<-stable_Var[cor(train[,stable_Var],trainY, method = "spearman")>.1]
    
    
    {
      #INNER VALIDATION
      variables<-c()
      reps<-20
      #foldsI <- createFolds(y = train$Y, k=reps, list = FALSE, returnTrain = FALSE)
      
      
      for (j in 1:reps){
        
        #foldI = which(foldsI == j)
        trainX<-as.matrix(train[,-ncol(train)])
        trainY <- train$Y
        
        
        lasso.cv <- cv.glmnet(trainX, trainY, alpha=1, nfolds = 10)
        
        lam<-lasso.cv$lambda.min
        lam<- lam* .8 #can be changed
        
        lasso.fit<-glmnet(trainX, trainY, alpha=1, lambda = lam)
        
        coefs <- coef(lasso.fit)
        
        notNull<-which(coefs!=0)
        
        newVariables <- tail(row.names(coefs)[notNull], -1)
        
        #print(tmp_variables)
        
        variables<-c(variables,newVariables)
      }
      
      
      freq = sort(table(variables),decreasing=TRUE)/(reps)
      #print(freq)
      
      if (length(which(freq>=.8))>10){
        stable_Var = names(which((head(freq,10))>=.8))
      }
      else {
        stable_Var = names(which((freq)>=0.8))
      }
      variablesSaved<-c(variablesSaved, stable_Var)
      
      print(i)
      print(NoF)
      #print(freq)
      
      
      stable_Var<-stable_Var[abs(cor(train[,stable_Var],train$Y, method = "spearman"))>.1]
      #train<-train[,abs(cor(train[,colnames(train)],train$Y,method = "spearman"))>.1]
      
      
      #RANDOM FOREST
      
      
      vars <- paste(stable_Var, collapse="+")
      
      form<-as.formula(paste("Y ~ ",vars,sep = ""))
      
      formPerm<-as.formula(paste("Yperm ~ ",vars,sep = ""))
      Yperm<-sample(train$Y)
    } #test
    
    
    rf.fit<-randomForest(form, train, ntree=500, importance=T)
    
    
    
    rf.fit.Perm<-randomForest(formPerm, train, importance=F)
    #importance(rf.test)
    
    test<-test[colnames(train)]
    RFpred<-c(RFpred,predict(rf.fit, test))
    
    RFpredperm<-c(RFpredperm,predict(rf.fit.Perm,test))
    
    
    
  }
  
  rf.R.Var <- c(rf.R.Var, cor(RFpred, trueY, method="spearman"))
  rf.R.Var.perm <- c(rf.R.Var.perm, cor(RFpredperm, trueY, method = "spearman"))
  
  rf.MSE.Var<-c(rf.MSE.Var, mean((trueY - RFpred) ^ 2))
  
  trueYsave<-c(trueYsave,trueY)
  RFpredsave.Var<-c(RFpredsave.Var, RFpred)
  RFpredsavePerm.Var<-c(RFpredsavePerm.Var, RFpredperm)
  
  print(rf.R.Var)
  
  trueY<-c()
  RFpred<-c()
  RFpredperm<-c()
  
}



freqsave = sort(table(variablesSaved),decreasing=TRUE)/(k*10)
print(freqsave)

dftab<-as.data.frame(freqsave)

corList<-matrix(nrow = length(df)-1, ncol=2)
for (i in 1:length(df)-1){
  corList[i,1]<-colnames(df)[i]
  corList[i,2]<-cor(df[,i],df$Y)
}

tmp_df<-as.data.frame(corList)
colnames(tmp_df)<-c("variablesSaved", "cor")

tmp_df2<-merge(dftab, tmp_df, by= "variablesSaved", sort=F)

write.csv(tmp_df2, "Varfreq.csv")





boxplot(rf.MSE.Var, rf.MSE., ylab="Spearman Correlation" ,col="steelblue", names = c("MSE Variance Model", "MSE Combined Model"))


boxplot(rf.R.Var, rf.R.Var.perm, ylab="Spearman Correlation" ,col="steelblue", names = c("Random Forest", "Permuted Response"))+
  expand_limits(y=c(-1,1))

corr<-c(rf.R.Var, rf.R.Var.perm)
var<-as.data.frame(corr)
var$group<-1
var$group[1:10]<-"Actual"
var$group[11:20]<-"Permuted"
var$group <- factor(var$group, levels = c("Actual", "Permuted"))

pal = colorRampPalette(c("#52adeb", "#fc2771"))

varBP<-ggplot(var, aes(group,corr,color=group))+geom_boxplot()+expand_limits(y=c(-.7,.7))+xlab("")+ylab("Spearman Correlation")+
  scale_color_manual(values=c("#D91100", "#949494"))+ scale_y_continuous(breaks = seq(-.6, .6, by = .4))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggsave('varBP.pdf',varBP, device = "pdf",dpi=500)






RFpredsave.Var

ggDF<-data.frame(RFpredsave.Var, trueYsave)

ggplot(ggDF, aes(RFpredsave.Var, trueYsave))+geom_point(color= 'blue')+
  geom_abline(slope=1, intercept = 0, color='red')+xlab("Random Forest Predicted Skin Score")+ylab("Observed Skin Score")


