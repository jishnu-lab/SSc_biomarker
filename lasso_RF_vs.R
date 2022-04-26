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

set.seed(405)
### IMPORT MANIPULATED DATA


resp<-read.csv("Scl/SkinScore_LafyatisAndBaselineTofa_MRSS.csv",header=T)
ID<-resp$orig.ident

C.vs.Scl<-read.csv("C.vs.Scl.csv")



Y <- resp[,2]
C.vs.Scl$Y<-Y



cormat.vs<-cor(C.vs.Scl[,1:length(C.vs.Scl)-1])

highlyCorrelated <- findCorrelation(cormat.vs, cutoff=0.99)

df<-C.vs.Scl[,-highlyCorrelated]

#df<-C.vs.Scl
#FOLDS
k <- 24



#CREATE VECTORS FOR LOOP
variables<-c()
variablesSaved<-c()
trueY<-c()
pred<-c()
trueY<-c()
RFpred<-c()
rf.MSE.vs<-c()
rf.R.vs<-c()
svm.R<-c()
trueYsave<-c()
RFpredsave<-c()
svm.pred<-c()
RFpredperm<-c()
rf.R.vs.perm<-c()
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
      reps<-20
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
      
      if (length(which(freq>=.6))>10){
        stable_Var = names(which((head(freq,10))>=.6))
      }
      else {
        stable_Var = names(which((freq)>=0.6))
      }
      variablesSaved<-c(variablesSaved, stable_Var)
      
      print(i)
      print(NoF)
      print(freq)
      
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
  
  rf.R.vs <- c(rf.R.vs, cor(RFpred, trueY, method="spearman"))
  rf.R.vs.perm <- c(rf.R.vs.perm, cor(RFpredperm, trueY, method = "spearman"))
  rf.MSE.vs<-c(rf.MSE.vs, mean((trueY - RFpred) ^ 2))
  
  
  trueYsave<-c(trueYsave,trueY)
  RFpredsave<-c(RFpredsave, RFpred)
  RFpredsavePerm<-c(RFpredsavePerm, RFpredperm)
  
  print(rf.R.vs)
  print(rf.R.vs.perm)
  
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

write.csv(tmp_df2, "freq.csv")





boxplot(rf.R.vs, rf.R.vs.perm, ylab="Spearman Correlation" ,col="steelblue", names = c("Random Forest", "Permuted Response"))





RFpredsave

ggDF<-data.frame(RFpredsave, trueYsave)

ggplot(ggDF, aes(RFpredsave, trueYsave))+geom_point(color= 'blue')+
  geom_abline(slope=1, intercept = 0, color='red')+xlab("Random Forest Predicted Skin Score")+ylab("Observed Skin Score")



ggplot(ggDF, aes(RFpredsavePerm, trueYsave))+geom_point(color= 'blue')+
  geom_abline(slope=1, intercept = 0, color='red')+xlab("Random Forest Predicted Skin Score Using Permutated Data")+ylab("Observed Skin Score")

cor(RFpred,trueY)








#EACH MODEL BOX PLOT
#Change reps to 20
#Include permute

#one slide with boxplots comp
# one dot plot
# one selected features
# shuffle y with respect to x







#plsr

set.seed(405)

freq = sort(table(variablesSaved),decreasing=TRUE)/(k*10*10)
print(freq)

stable_Var = names(which(freq>0.4))


Svars <- paste(stable_Var, collapse="+")
formPLS<-as.formula(paste("Y ~ ",Svars,sep = ""))



modPLS<-plsr(formPLS, data=df, scale=TRUE, valdiation="CV")


#modPLS<-plsr(Y~., data=df, scale=TRUE, valdiation="CV")

summary(modPLS)

validationplot(modPLS,val.type="MSEP")


plot(modPLS, plottype = "scores")

plsDF<-data.frame(modPLS)

PLSdf<-data.frame(modPLS$scores)


pal = colorRampPalette(c("blue", "red"))

highlow<-ifelse(Y>median(Y),1,0)

plot(modPLS$scores, pch=19, col=pal(Y)) 


image(1, modPLS$scores, t(seq_along(modPLS$scores)), col=pal, axes=FALSE)
axis(4)

ggplot(modPLS, aes(scores))

modPLS$ncomp

pls.MSE<-mean((test$price - predict(modPLS, test,12)) ^ 2)

summary(modPLS)


































