#look at features w/ whole dataset; most important

rm(list = ls())
cat("\014")

#3 groups?

# 0 <- healthy
# .5 <- survived
# 1 <- died

library(glmnetUtils) #cv.glmnet()
library(pROC) #auc()
library(caret) #createDataPartition(), createFolds()
library(e1071) #svm()
library(randomForest)
library(kernlab) #tuning svm train

setwd("/bgfs/djishnu/Trirupa/ABomics.Prj/")

dataSet=read.csv("KidneyTransplant_DataProbes_RCA_SNAunimputed.csv", header=TRUE, sep=",")

# dataSet <- read.csv("fcProfile_nonCanonical_Healthy_Covid_Updated.csv")
# dataSet <- read.csv("fcProfile_Healthy_Covid_C_NC_Updated.csv")
# dataSet <- read.csv("fcProfile_Healthy_Covid_C_NC_without_FCRs_Glycans.csv")
# dataSet <- read.csv("fcProfile_Healthy_Covid_3_Flu_Updated.csv")
#dataSet <- read.csv("fcProfile_Healthy_Covid_2_Flu_Updated.csv")

#healthy <- 2  survived <- 0    died <- 1
#healthy <- 0  survived <- 0.5  died <- 1

# all 0's to 0.5    then all 2's to 0
##dataSet["Y"][dataSet["Y"] == 0] <- 0.5
##dataSet["Y"][dataSet["Y"] == 2] <- 0

# data <- dataSet
# 
# SampleID <- data[, 1]
# Y <- data[, 2]
# threshold <- 50
# data <- data[, !sapply(data, function(x) mean(x)) < threshold]
# data <- cbind.data.frame(SampleID, Y, data)

Y <- dataSet[, 1]
# Normalize using z score
scaled <- scale(dataSet[,-c(1)])
dataSet <- cbind(Y, scaled)


aucSVM <- c()
aucRF <- c()
replicateCV <- 10
for(cv in 1:replicateCV) 
{
  
  numk <- 10
  
  folds <- createFolds(dataSet[,"Y"], k=numk, list=FALSE, returnTrain = FALSE) #n folds with row indexes to be included in fold
  trainFolds <- cbind(folds, dataSet)
  
  repl.Y <- c()
  repl.svm.Y <- c()
  repl.rf.Y <- c()
  numVarsUsed <- c()
  
  for(n in 1:numk)
  {
    tempTest <- trainFolds[which(trainFolds$folds==n), ] #nth fold #! for names & - for col indexes
    tempTrain <- trainFolds[which(trainFolds$folds!=n),] #all folds except nth fold #need matrix format
    
    cv.lasso <- cv.glmnet(x=as.matrix(tempTrain[,!colnames(tempTrain) %in% c("folds","SampleID","Y")]), y=tempTrain[,"Y"],alpha=1, nfolds=10)
    lambda <- cv.lasso$lambda.min   #lambda <- cv.lasso$lambda.1se 
    lasso <- glmnet(x=as.matrix(tempTrain[,!colnames(tempTrain) %in% c("folds","SampleID","Y")]), y=tempTrain[,"Y"], alpha=1, lambda = lambda*0.1)
    #6/7/10 features lambda 0.7  
    
    
    allVars <- data.frame(data.matrix(lasso[["beta"]]))
    selectedVars <- row.names(subset(allVars, s0!=0.000000e+00))
    vars <- c()
  
    if (length(selectedVars)==1) ##If no variables
      {
      data <- tempTrain[,!colnames(tempTrain) %in% c("folds","SampleID","Y")]
      randomSelect <- sample(ncol(data), 3) ##
      vars <- sapply(randomSelect, function(x){vars <- append(colnames(data)[x], vars) }) ##
      selectedVars <- vars
     }
    
    numVarsUsed <- append(numVarsUsed, length(selectedVars))
    
    tempVarsTrain <- data.frame(cbind(tempTrain[,c("Y")],tempTrain[,colnames(tempTrain) %in% selectedVars]))
    tempVarsTest <- data.frame(cbind(tempTest[,c("Y")] ,tempTest[,colnames(tempTest) %in% selectedVars]))
    
    colnames(tempVarsTrain)[1] <- "Y"
    colnames(tempVarsTest)[1] <- "Y"
    
    repl.Y <- append(repl.Y, tempVarsTest$Y)
    
    #svm
  
    #beforehand tune; cost=10/100; kernel: radial/linear; make scale true if not normalized
    
    svm.model <- svm(Y ~., data=tempVarsTrain, kernal = "linear", cost=10, scale=FALSE)
    svm.predY <- predict(svm.model, newdata = tempVarsTest[,-1])
    repl.svm.Y <- append(repl.svm.Y, svm.predY)
    
    #rf
    rf.model <- randomForest(Y~.,data=tempVarsTrain, ntree=100, importance = TRUE)
    rf.predY <- predict(rf.model, newdata = tempVarsTest[,-1])
    repl.rf.Y <- append(repl.rf.Y, rf.predY)
  
  }

#aucSVM[cv] <- auc(repl.Y, repl.svm.Y)
aucSVM[cv] <- (multiclass.roc(repl.Y, repl.svm.Y, levels=c(0, .5, 1)))$auc
aucRF[cv] <- (multiclass.roc(repl.Y, repl.rf.Y, levels=c(0, .5, 1)))$auc
#aucRF[cv] <- auc(repl.Y, repl.rf.Y)

}

# numVarsUsed
aucSVM
aucRF

median(aucSVM)
median(aucRF)

# write.csv(aucRF, "AUC_Lasso_RF_ETP_fcProfile_3levels.csv")
# write.csv(aucSVM, "AUC_Lasso_SVM_ETP_fcProfile_3levels.csv")

# write.csv(aucRF, "AUC_Lasso_RF_fcProfile_3levels_without_FCRs_Glycans.csv")
# write.csv(aucSVM, "AUC_Lasso_SVM_ETP_fcProfile_3levels_without_FCRs_Glycans.csv")

# write.csv(aucRF, "AUC_Lasso_RF_fcProfile_3levels_3Flu.csv")
# write.csv(aucSVM, "AUC_Lasso_SVM_fcProfile_3levels_3Flu.csv")

# write.csv(aucRF, "AUC_Lasso_RF_fcProfile_3levels_2Flu.csv")
# write.csv(aucSVM, "AUC_Lasso_SVM_fcProfile_3levels_2Flu.csv")


#actual_aucInfo <- data.frame(result)
#permuted_aucInfo <- data.frame(result)

#actual_aucs <- as.vector(unlist(actual_aucInfo$auc.RF))
#actual_aucs <- as.vector(rbind(unlist(actual_aucInfo$auc.SVM),unlist(actual_aucInfo$auc.RF)))

#permuted_aucs <- as.vector(unlist(permuted_aucInfo$auc.RF))
#permuted_aucs <- as.vector(rbind(unlist(permuted_aucInfo$auc.SVM), unlist(permuted_aucInfo$auc.RF)))


#pdf("boxplot with svm.pdf")
#boxplot(actual_aucs, permuted_aucs, names = c("actual", "permuted"))
#dev.off()


#lasso whole dataset

#Train & Test
#splitIndex <- createDataPartition(dataSet$Y, times = 1, p = 0.67, list = FALSE, group = 2) #splitIndex stores vector of values to be included in train based on p=.67
#train <- dataSet[splitIndex,] # selects rows @ all row indexes specified in vector
#test <- dataSet[-splitIndex,] # [-c(),] = drops/excludes rows (by row indexes) specified in vector

cv.lasso <- cv.glmnet(x=as.matrix(dataSet[,-c(1,2)]), y=dataSet[,"Y"], alpha=1, nfolds=10)
lambda <- cv.lasso$lambda.min   #lambda <- cv.lasso$lambda.1se 
# lasso <- glmnet(x=as.matrix(dataSet[,-c(1,2)]), y=dataSet[,"Y"], alpha=1, lambda = lambda*.95)
lasso <- glmnet(x=as.matrix(dataSet[,-c(1,2)]), y=dataSet[,"Y"], alpha=1, lambda = lambda)
everyVar <- data.frame(data.matrix(lasso[["beta"]]))
impVars <- row.names(subset(allVars, s0!=0.000000e+00))
impVars

cv.lasso <- cv.glmnet(x=as.matrix(train[,-c(1,2)]), y=train$Y, alpha=1, nfolds=10)
lambda <- cv.lasso$lambda.min   #lambda <- cv.lasso$lambda.1se
lasso <- glmnet(x=as.matrix(train[,-c(1,2)]), y=train$Y, alpha=1, lambda = lambda*.75)
everyVar <- data.frame(data.matrix(lasso[["beta"]]))
impVars <- row.names(subset(allVars, s0!=0.000000e+00))

#[1] "IgG2.Nsp13" "IgM.Nsp13"  "FcR1.Orf3a" "FcR2b.Nsp3" "SNACy3.M"


#tuning svm with whole dataSet

#svm.model <- svm(as.factor(Y) ~., data=dataSet, kernal = "radial", type = "C-classification", cost=.25, sigma=0.01783144, scale=FALSE) 
#tune.svm = tune(method="svm", Y~., data = dataSet, kernal = "radial", type = "C-classification", ranges = list(cost = c(10), gamma = c(1e-8), epsilon = c(1e-04)))
#trnControl <- trainControl(method = "repeatedcv", number = 10, repeats = 100)
#tuneModel.svm <- train(as.factor(Y)~., data=dataSet, method = "svmRadialSigma", trControl = trnControl)
#tune$best.parameters
#names(getModelInfo())

#tuning rf with whole dataSet

#rf.model <- randomForest(Y~.,data=dataSet, ntree=10, mtry=50, importance = TRUE)
#tune.rf = tune(method="rf", Y~., data = tempVarsTrain, ntree = 1)
#trainCntrl <- trainControl(method = "repeatedcv", number = 10, repeats = 30)
#tuneModel.rf <- train(as.factor(Y)~., data=dataSet, method = "rf", trControl = trainCntrl)
#showing mtry=50
