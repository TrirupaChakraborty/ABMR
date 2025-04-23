rm(list = ls())
cat("\014")

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

#data <- read.csv("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/miRNA_data/miRNA_features_labels.csv")
data= read.csv("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/miRNA_data/forLASSO/XY_maxEXP_sub3pt6_var25_cov1_filtered_miRNA_17042023.csv")
Y <- data[,1]
scaledData <- scale(data[,-1])
data <- cbind.data.frame(Y, scaledData)


dF <- c(1:nrow(data))

train <- data


X = as.matrix(train[, -1]) ### selecting features from the whole data
y = Y

## running LASSO with tuned lambda_factor to get features 50 times
selectedVars <- c()
#lambda=0.04199923 ##for tuning factor of 0.5
for (num_fs in 1:50){
  glmnet1 <- cv.glmnet(X, y=Y, alpha=1, nfolds = 10)
  
  lambda <- glmnet1$lambda.min
  lambda <- lambda*0.5
  
  glmnet2 <- glmnet(X, y=Y, alpha=1, lambda = lambda)
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
  selectedVars <- append(selectedVars, variables)
}
print(paste0("lambda=",lambda))
##selecting the final set of stable features
feat_freq=as.data.frame(sort(table(selectedVars),decreasing=TRUE))

feat_select = feat_freq[c(which(((feat_freq$Freq)/50)>=0.25)),]
feat_names = as.character(feat_select$selectedVars)
print("selected vars: ")
print(feat_names)

