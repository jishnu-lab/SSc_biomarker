#pipeline for fibroblasts

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

set.seed(405)
### IMPORT MANIPULATED DATA


k_DF<-read.csv("keratinocytes.csv")


df<-k_DF

df$Y<-resp$MRSS



cormat<-cor(df[,1:length(df)-1])

highlyCorrelated <- findCorrelation(cormat, cutoff=0.95)

df<-df[,-highlyCorrelated]

#FOLDS
k <- 24



#CREATE VECTORS FOR LOOP
Variables<-c()
VariablesSaved<-c()
trueY<-c()
pred<-c()
trueY<-c()
RFpred<-c()
rf.MSE.Ker<-c()
rf.R.Ker<-c()
svm.R<-c()
trueYsave<-c()
RFpredsave.Ker<-c()
svm.pred<-c()
RFpredperm<-c()
rf.R.Ker.perm<-c()
RFpredsavePerm.Ker<-c()
svm.predsave<-c()
#rfFMat<-matrix(nrow=20, ncol=10)

set.seed(101)
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
    
    #trainX <- as.matrix(train[, -ncol(df)]) #drop Y
    #testX <- as.matrix(test[, -ncol(df)])
    
    #trainY <- train[,ncol(df)]
    testY <- test[,ncol(df)]
    
    #NOTE TRUE Y
    trueY<-c(trueY,testY)
    
    
    {
      #INNER VALIDATION
      Variables<-c()
      reps<-15
      #foldsI <- createFolds(y = train$Y, k=reps, list = FALSE, returnTrain = FALSE)
      
      for (j in 1:reps){
        
        #foldI = which(foldsI == j)
        trainX<-as.matrix(train[,-ncol(df)])
        trainY <- as.matrix(train[,ncol(df)])
        
        
        lasso.cv <- cv.glmnet(trainX, trainY, alpha=1, nfolds = 10)
        
        lam<-lasso.cv$lambda.min
        lam<- lam* .7 #can be changed
        
        lasso.fit<-glmnet(trainX, trainY, alpha=1, lambda = lam)
        
        coefs <- coef(lasso.fit)
        
        notNull<-which(coefs!=0)
        
        newVariables <- tail(row.names(coefs)[notNull], -1)
        
        #print(tmp_Variables)
        
        Variables<-c(Variables,newVariables)
      }
      
      
      freq = sort(table(Variables),decreasing=TRUE)/(reps)
      #print(freq)
      
      if (length(which(freq>=.8))>10){
        stable_Var = names(which((head(freq,10))>=.8))
      }
      else {
        stable_Var = names(which((freq)>=0.8))
      }
      stable_Var<-stable_Var[abs(cor(train[,stable_Var],train$Y, method = "spearman"))>.1]
      VariablesSaved<-c(VariablesSaved, stable_Var)
      
      print(i)
      print(NoF)
      #print(freq)
      
      #RANDOM FOREST
      
      
      Vars <- paste(stable_Var, collapse="+")
      
      form<-as.formula(paste("Y ~ ",Vars,sep = ""))
      
      formPerm<-as.formula(paste("Yperm ~ ",Vars,sep = ""))
      Yperm<-sample(train$Y)
    } #test
    
    
    rf.fit<-randomForest(form, train, ntree=500, importance=T)
    
    
    
    rf.fit.Perm<-randomForest(formPerm, train, importance=F)
    #importance(rf.test)
    
    RFpred<-c(RFpred,predict(rf.fit, test))
    
    RFpredperm<-c(RFpredperm,predict(rf.fit.Perm,test))
    
    
    
  }
  
  rf.R.Ker <- c(rf.R.Ker, cor(RFpred, trueY, method="spearman"))
  rf.R.Ker.perm <- c(rf.R.Ker.perm, cor(RFpredperm, trueY, method = "spearman"))
  
  rf.MSE.Ker<-c(rf.MSE.Ker, mean((trueY - RFpred) ^ 2))
  
  trueYsave<-c(trueYsave,trueY)
  RFpredsave.Ker<-c(RFpredsave.Ker, RFpred)
  RFpredsavePerm.Ker<-c(RFpredsavePerm.Ker, RFpredperm)
  
  print(rf.R.Ker)
  
  trueY<-c()
  RFpred<-c()
  RFpredperm<-c()
  
}



corr<-c(rf.R.Ker, rf.R.Ker.perm)
ker<-as.data.frame(corr)
ker$group<-1
ker$group[1:10]<-"Actual"
ker$group[11:20]<-"Permuted"
ker$group <- factor(ker$group, levels = c("Actual", "Permuted"))


kerBP<-ggplot(ker, aes(group,corr,color=group))+geom_boxplot()+expand_limits(y=c(-.7,.7))+xlab("")+ylab("Spearman Correlation")+
  scale_color_manual(values=c("#D91100", "#949494"))+ scale_y_continuous(breaks = seq(-.6, .6, by = .4))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave('kerBP.pdf',kerBP,device="pdf",dpi=500)


