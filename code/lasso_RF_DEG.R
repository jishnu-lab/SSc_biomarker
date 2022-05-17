#RUN pipeline on DEGs


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
## IMPORT MANIPULATED DATA


resp<-read.csv("SkinScore_LafyatisAndBaselineTofa_MRSS.csv",header=T)
ID<-resp$orig.ident

DEG.df<-read.csv("DEGdf.csv")



Y <- resp[,2]
DEG.df$Y<-Y

DEG.Var.df$Y<-Y


#remove highly correlated data
cormatDEG<-cor(DEG.df[,1:length(DEG.df)-1])

highlyCorrelated <- findCorrelation(cormatDEG, cutoff=0.85)

#manually remove neural
DEG.df<-DEG.df %>% select(-contains("c.14"))


df<-DEG.df[,-highlyCorrelated]
#FOLDS
k <- 24



#CREATE VECTORS FOR LOOP
variables<-c()
variablesSaved<-c()
trueY<-c()
pred<-c()
trueY<-c()
RFpred<-c()
rf.MSE.DEG<-c()
rf.R.DEG<-c()
svm.R<-c()
trueYsave<-c()
RFpredsave<-c()
svm.pred<-c()
RFpredperm<-c()
rf.R.DEG.perm<-c()
RFpredsavePerm<-c()
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
          
          #trainX <- as.matrix(train[, -ncol(df)]) #drop Y
          #testX <- as.matrix(test[, -ncol(df)])
          
          #trainY <- train[,ncol(df)]
          testY <- test[,ncol(df)]
          
          #NOTE TRUE Y
          trueY<-c(trueY,testY)
          
          
          {
          #INNER VALIDATION
          variables<-c()
          reps<-30
          #foldsI <- createFolds(y = train$Y, k=reps, list = FALSE, returnTrain = FALSE)
          
          for (j in 1:reps){
          
            #foldI = which(foldsI == j)
            trainX<-as.matrix(train[,-ncol(df)])
            trainY <- as.matrix(train[,ncol(df)])
      
      
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
          
          if (length(which(freq>=.6))>8){
            stable_Var = names(which((head(freq,8))>=.6))
          }
          else {
            stable_Var = names(which((freq)>=0.6))
          }
          
          stable_Var<-stable_Var[abs(cor(train[,stable_Var],train$Y, method = "spearman"))>.1]
          
          variablesSaved<-c(variablesSaved, stable_Var)
          
          
          
          print(i)
          print(NoF)
          #print(freq)
          
           #RANDOM FOREST
          
          
          vars <- paste(stable_Var, collapse="+")
          
          form<-as.formula(paste("Y ~ ",vars,sep = ""))
          
          formPerm<-as.formula(paste("Yperm ~ ",vars,sep = ""))
          Yperm<-sample(train$Y)
          } #test
          
          
          rf.fit<-randomForest(form, train, ntree=500, importance=T)
      
          
          
         rf.fit.Perm<-randomForest(formPerm, train, importance=F)
          #importance(rf.test)
          
          RFpred<-c(RFpred,predict(rf.fit, test))
          
          RFpredperm<-c(RFpredperm,predict(rf.fit.Perm,test))
          
          
          
        }
        
        rf.R.DEG <- c(rf.R.DEG, cor(RFpred, trueY, method="spearman"))
        rf.R.DEG.perm <- c(rf.R.DEG.perm, cor(RFpredperm, trueY, method = "spearman"))
        rf.MSE.DEG<-c(rf.MSE.DEG, mean((trueY - RFpred) ^ 2))
      
        
        trueYsave<-c(trueYsave,trueY)
        RFpredsave<-c(RFpredsave, RFpred)
        RFpredsavePerm<-c(RFpredsavePerm, RFpredperm)
        
        print(rf.R.DEG)
        print(rf.R.DEG.perm)
        
        trueY<-c()
        RFpred<-c()
        RFpredperm<-c()
        
      }

  
  




rf.MSE<-mean((trueYsave - RFpredsave) ^ 2)
rf.cor<-cor(trueYsave, RFpredsave)



freqsave = sort(table(variablesSaved),decreasing=TRUE)/(k*10)
print(freqsave)

TotalStableVar = names(which(freq>0.29))  #### the 0.6 can be changed ######



#print(stable_Var)





dftab<-as.data.frame(freqsave)

corList<-matrix(nrow = length(df)-1, ncol=2)
for (i in 1:length(df)-1){
  corList[i,1]<-colnames(df)[i]
  corList[i,2]<-cor(df[,i],df$Y)
}

tmp_df<-as.data.frame(corList)
colnames(tmp_df)<-c("variablesSaved", "cor")

tmp_df2<-merge(dftab, tmp_df, by= "variablesSaved", sort=F)

write.csv(tmp_df, "freq.csv")





boxplot(rf.R.DEG, rf.R.DEG.perm, ylab="Spearman Correlation" ,col="steelblue", names = c("Random Forest", "Permuted Response"))

RFpredsave

ggDF<-data.frame(RFpredsave, trueYsave)

ggplot(ggDF, aes(RFpredsave, trueYsave))+geom_point(color= 'blue')+
  geom_abline(slope=1, intercept = 0, color='red')+xlab("Random Forest Predicted Skin Score")+ylab("Observed Skin Score")



ggplot(ggDF, aes(RFpredsavePerm, trueYsave))+geom_point(color= 'blue')+
  geom_abline(slope=1, intercept = 0, color='red')+xlab("Random Forest Predicted Skin Score Using Permutated Data")+ylab("Observed Skin Score")

cor(RFpred,trueY)




corr<-c(rf.R.DEG, rf.R.DEG.perm)
DEG<-as.data.frame(corr)
DEG$group<-1
DEG$group[1:10]<-"Actual"
DEG$group[11:20]<-"Permuted"
DEG$group <- factor(DEG$group, levels = c("Actual", "Permuted"))


degBP<-ggplot(DEG, aes(group,corr,color=group))+geom_boxplot()+expand_limits(y=c(-.7,.7))+xlab("")+ylab("Spearman Correlation")+
  scale_color_manual(values=c("#D91100", "#949494"))+ scale_y_continuous(breaks = seq(-.6, .6, by = .4))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave('DEGBP.pdf',degBP,device="pdf",dpi=500)





