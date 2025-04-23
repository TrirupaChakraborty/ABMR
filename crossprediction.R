###this code does cross-prediction between C1 cohort and validation cohort (VC). 

rm(list = ls())
cat("\014")

#setwd("/Users/trirupachakraborty/Desktop/Ab profiling project/inputData/miRNA_data/")
#setwd("/Users/trirupachakraborty/Desktop/Ab profiling project/inputData/")
setwd("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/miRNA_data/")

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

##loading datasets
Dataset_C1=read.csv("./miRNA_features_labels.csv")
Dataset_VC=read.csv("./ValidationCohort1/miRNA_data/input_data/unfiltmiR_XY.csv")

##features selected in C1
#feat=c("miR.129.1.3p", "miR.326","miR.3162.5p", "miR.431.5p", "miR.28.5p", 
       #"miR.4737", "miR.486.5p", "miR.495.3p", "miR.210.3p", 
       #"miR.543", "miR.192.5p",
       #"miR.6079", "miR.181a.3p", "miR.6865.5p", "miR.3692.3p", "miR.154.3p", "miR.134.5p")

#feat_new=c("miR.129.1.3p", "miR.495.3p" , "miR.431.5p" , "miR.3162.5p" , "miR.28.5p" , "miR.4737" , "miR.486.5p", "miR.210.3p" , "miR.192.5p" , "miR.154.3p" , "miR.6865.5p" , "miR.6782.5p" , "miR.134.5p")

feat=c("miR.129.1.3p", "miR.192.5p", "miR.3162.5p", "miR.4737", "miR.486.5p", "miR.495.3p", "miR.6865.5p", "miR.93.3p")
###first performing sanity check to see that auc_X1 is close to 1 with the feats selected from X1
# subsetting C1 based on feats selected
X1_feat= which(colnames(Dataset_C1)%in% feat )
data_X1=Dataset_C1[,X1_feat] ### X1 with lasso selected features
Y1=Dataset_C1[,1]


#calculating auc_X1 (should be close to 1)
#for svm
model_svmc1 <- svm(x=as.matrix(data_X1),y=Y1, kernel="radial",cost=20, scale=FALSE) ## training svm on C1
predY1_svm=predict(model_svmc1, data_X1) ##predicting on C1 
aucY1_svm=auc(Y1,predY1_svm)
aucY1_svm

#for RF
model_RFc1 <- randomForest(x=as.matrix(data_X1),y=Y1, importance=TRUE, ntree = 10)
predY1_rf=predict(model_RFc1,data_X1)
aucY1_rf=auc(Y1,predY1_rf)
aucY1_rf

###Now recording performance on VC
# subsetting VC based on features selected in Cohort C1
data_X2=Dataset_VC[,X1_feat]
Y2=Dataset_VC[,1]

# use model trained on C1 to check performance on (auc_X2) 
#for SVM
predY2_svm=predict(model_svmc1, data_X2) ##predicting on C1 
aucY2_svm=auc(Y2,predY2_svm)
aucY2_svm

#for RF
predY2_rf=predict(model_RFc1,data_X2)
aucY2_rf=auc(Y2,predY2_rf)
aucY2_rf

###plotting the ROC curve
roc_score = roc(Y2,predY2_svm)  # AUC score
plot(roc_score, main="ROC curve -- Logistic Regression ")

roc_score = roc(Y2,predY2_rf)  # AUC score
plot(roc_score, main="ROC curve -- Logistic Regression ")


