rm(list = ls())
cat("\014")

# setwd("/Users/sar210/Box/MSD_data_ABMR/")
setwd("/bgfs/djishnu/Trirupa/ABomics.Prj/")

library(e1071)
library(caret)
library(glmnet)
library(randomForest)
library(tree)
library(gbm)
library(matrixStats)
library(readxl)
library(cvAUC)
library(pROC)

data=read.csv("KidneyTransplant_DataProbes_RCA_SNAunimputed.csv", header=TRUE, sep=",")


# stable vs TCMR
# data <- read.csv("MSD_data_raw_13.csv")
# 0.62
# stable vs ABMR
# data <- read.csv("MSD_data_raw_14.csv")
# 0.64

# change this (slide - 15)
# data <- read.csv("MSD_data_raw_3_EL.csv")
# data <- data[,-c(1,2,3,4)]
# 0.78
#  (slide - 14)
data <- read.csv("MSD_data_raw_4_EL.csv")
data <- data[,-c(1,2,3,4)]
# 0.64

# data <- read.csv("MSD_data_raw_34_early.csv")
# 0.96

# (check slide - 13)
# data <- read.csv("MSD_data_raw_34_late.csv")
# 0.68

# data <- read.csv("MSD_data_raw_134.csv")
# data <- read.csv("MSD_data_raw_034.csv")
# data <- read.csv("MSD_data_raw_0134.csv")

# data <- data[,-c(1,2,4,5)]

Y <- data[,1]
# threshold <- 100
# data <- data[, !sapply(data, function(x) mean(x)) < threshold] 
# data <- cbind.data.frame(Y, data)

scaledData <- scale(data[,-1])

# Y <- data[,1]
data <- data[, -1]
data <- cbind.data.frame(Y, scaledData)


k <- 10
# k <- nrow(data) - 1
ACC <- c() ##can remove this
aucRF <- c()
auctree <- c() ##can remove this
aucSVM <- c()

dF <- c(1:nrow(data))

# replicateCV <- 1
for (replicateCV in 1:10)
{
  
    folds <- createFolds(y = data$Y, k=k, list = FALSE, returnTrain = FALSE)
    myData <- cbind(data, folds)
    # FOLDS[replicateCV, ] <- folds
    ##print(dim(myData))
    ##print(colnames(myData))
    
    Y <- c()
    append.RF <- c()
    append.SVM <- c()
    append.tree <- c()
    append.lasso <- c()
    selectedVars <- c()
    
    # NoF <- 1
    for (NoF in 1:k)
    {
      
      fold = which(folds == NoF)
      ##print(fold)
      ##print(folds)
      
      train <- myData[-fold, ]
      test <- myData[fold, ]
      
      train <- train[, -ncol(train)]
      test <- test[, -ncol(test)]
      ##print(train)
      ##print(test)
      
      X = as.matrix(train[, -1])
      y = train$Y
  
      glmnet1 <- cv.glmnet(X, y=train$Y, alpha=1, nfolds = k, family="multinomial")
      
      lambda <- glmnet1$lambda.min
      # lambda <- glmnet1$lambda.1se
      lambda <- lambda*0.75
      
      glmnet2 <- glmnet(X, y=train$Y, alpha=1, lambda = lambda, family="multinomial")
      c <- coef(glmnet2)
      
      inds<-which(c!=0)
      variables<-row.names(c)[inds]
      len <- length(variables)
      
      if (len == 1)
      {
        randomSelect <- sample(ncol(data), 3)
        variables <- row.names(c)[randomSelect]
      } else
      {
        variables<-row.names(c)[inds]
        variables <- variables[2:len]
      }
      selectedVars <- append(selectedVars, len)
        
      tempTr <- train[, (names(train) %in% variables)]
      tempTrain <- cbind.data.frame(train$Y, tempTr)
      tempTr <- test[, (names(test) %in% variables)]
      tempTest <- cbind.data.frame(test$Y, tempTr)
      
      colnames(tempTrain)[1] <- "Y"
      colnames(tempTest)[1] <- "Y"
          
      Y <- append(Y, tempTest$Y)
      
      # svm
      svmfit = svm(Y~ ., data = tempTrain , kernel="linear", cost=10, scale=FALSE)
      yhat.SVM = predict(svmfit, newdata = tempTest[, -1])
      append.SVM <- append(append.SVM, yhat.SVM)

      # RF
      RFfit <- randomForest(Y~ ., data = tempTrain, importance=TRUE, ntree = 100)
      yhat.RF = predict(RFfit, newdata = tempTest[, -1])
      append.RF <- append(append.RF, yhat.RF)
      
      
    }
    
    # tempdF <- cbind.data.frame(Y, append.SVM)
    # ACC[replicateCV] <- length(which(Y == append.SVM))/nrow(data)
    aucRF[replicateCV] <- auc(Y, append.RF)
    aucSVM[replicateCV] <- auc(Y, append.SVM)
    # auctree[replicateCV] <- auc(Y, append.tree)
    # aucLasso[replicateCV] <- auc(Y, append.lasso)
    
    # dF <- cbind.data.frame(dF, tempdF)

}
selectedVars

aucRF
median(aucRF)
aucSVM
median(aucSVM)

# boxplot(aucSVM)
# boxplot(aucRF)
sdf=head(data)


# write.csv(aucRF, "AUC_MSD_data_raw_34_early_RF.csv")
# write.csv(aucSVM, "AUC_MSD_data_raw_34_early_SVM.csv")

# write.csv(aucRF, "AUC_MSD_data_raw_34_late_RF.csv")
# write.csv(aucSVM, "AUC_MSD_data_raw_34_late_SVM.csv")

# write.csv(aucRF, "AUC_MSD_data_raw_3_EL_RF.csv")
# write.csv(aucSVM, "AUC_MSD_data_raw_3_EL_SVM.csv")

# write.csv(aucRF, "AUC_MSD_data_raw_4_EL_RF.csv")
# write.csv(aucSVM, "AUC_MSD_data_raw_4_EL_SVM.csv")

# write.csv(aucRF, "AUC_MSD_data_std_14_RF.csv")
# write.csv(aucSVM, "AUC_MSD_data_std_14_SVM.csv")

# write.csv(aucRF, "AUC_MSD_data_std_3_EL_RF.csv")
# write.csv(aucSVM, "AUC_MSD_data_std_3_EL_SVM.csv")

# write.csv(aucRF, "AUC_MSD_data_std_4_EL_RF.csv")
# write.csv(aucSVM, "AUC_MSD_data_std_4_EL_SVM.csv")

# write.csv(aucRF, "AUC_MSD_data_std_34_early_RF.csv")
# write.csv(aucSVM, "AUC_MSD_data_std_34_early_SVM.csv")

# write.csv(aucRF, "AUC_MSD_data_std_34_late_RF.csv")
# write.csv(aucSVM, "AUC_MSD_data_std_34_late_SVM.csv")


# write.csv(aucRF, "AUC_MSD_data_std_134_late_RF.csv")
# write.csv(aucSVM, "AUC_MSD_data_std_134_late_SVM.csv")

# write.csv(aucRF, "AUC_MSD_data_std_034_late_RF.csv")
# write.csv(aucSVM, "AUC_MSD_data_std_034_late_SVM.csv")

# write.csv(aucRF, "AUC_MSD_data_std_0134_late_RF.csv")
# write.csv(aucSVM, "AUC_MSD_data_std_0134_late_SVM.csv")
