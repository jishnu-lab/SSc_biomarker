## NOT OPERATIONAL


m_DF<-read.csv("mye.csv")[-1]
df<-m_DF
df$Y<-resp$MRSS

set.seed(405)
k<-24

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







boxplot(rf.R.Mye, rf.R.Mye.perm, ylab="Spearman Correlation" ,col="steelblue", names = c("Random Forest", "Permuted Response"))

