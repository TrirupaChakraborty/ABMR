
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



#################Importing z-score data ######################


data <- read.csv("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/miRNA_data/miRNA_features_labels.csv")
Y_data <- data[,1]
colnames(data)[1] <- "Y" 
#data <- data[-c(1)]
#X_data=ori_myData[,-1]
X_data=data[,-1]
################Permutation Test Outside##############
# data <- read.csv("/Users/hanxix/Desktop/Research/Latent_Model/Skin3PrimeV2_NormalSSC_Lafyatis/seurat_new_DEG/Lasso/lasso_fibro.csv", header = FALSE)
# #data <- read.csv("/Users/hanxix/Desktop/Research/Latent_Model/Skin3PrimeV2_NormalSSC_Lafyatis/seurat_new_DEG/Lasso/lasso_myeloids.csv", header = FALSE)
# Y <- data[,1]
# colnames(data)[1] <- "Y"
# permutedY <- sample(Y)
# 
# data$Y <- permutedY

################ Set K for K-fold CV  #######################
k <- 10

ACC <- c()
aucRF <- c()
aucRF_permut <- c()
auctree <- c()
aucSVM <- c()
aucSVM_permut <- c()


dF <- c(1:nrow(data))
permutedY <- sample(Y_data) ##sampling Y from the main dataset
#replicateCV=1
for (replicateCV in 1:3){
    print("in Replication")
    print(replicateCV)
    ################### Scramble The Response Vector ###############
    permutedY <- sample(Y_data) 
    
    ################## Get Folds for the Data ######################
    folds <- createFolds(y = data$Y, k=k, list = FALSE, returnTrain = FALSE)
    # at each fold, we make sure to start with the true y label
    ori_myData <- cbind(data, folds)
    
    Y <- c() ##not in use
    trueY <- c()
    permutY <- c()
    
    append.RF <- c()
    append.SVM <- c()
    
    append.RF_permut <- c()
    append.SVM_permut <- c()
    
    append.tree <- c() ##not in use
    append.lasso <- c() ##not in use
    
    selectedVars <- c() 
    selectedVars_permut <- c() ##not in use

    #for each fold of the data
    #NoF=1
    for (NoF in 1:k){
      fold = which(folds == NoF)
      #print(dim(ori_myData))
      myData <- ori_myData
      
      
      #1 is regular, 2 is permute 
      for (m in 1:2){
        # we set y as regular
        if (m == 1){
          #print('Regular Lasso Running...')
          train <- myData[-fold, ]
          test <- myData[fold, ]
          
          # delete the fold column
          train <- train[, -ncol(train)]
          test <- test[, -ncol(test)]
          
          X = as.matrix(train[, -1]) # X stays the same for permut and non-permut
          y_train = train$Y
        }
        # the second loop, we use scrambled y
        if (m == 2){
          #print('Permutation Test Lasso Running...')
          #print(dim(myData))
          myData$Y <- permutedY
          train <- myData[-fold, ]
          test <- myData[fold, ]
          
          # delete the fold column
          train <- train[, -ncol(train)]
          test <- test[, -ncol(test)]
          
          X = as.matrix(train[, -1]) # X stays the same for permute and non-permute
          y_train = train$Y
        }
        
        
        #feature selection
        variables <- c()
        #print("selecting features")
        ########################### Loop for Feature Selection ###############
        #num_FS=1
        for (num_FS in 1:5){
          #print(num_FS)
          ######################### Run CV.Glmnet #############################
          glmnet1 <- cv.glmnet(X, y=train$Y, alpha=1, nfolds = k)
          
          # tuning lambda
          lambda <- glmnet1$lambda.min
          lambda <- lambda*0.7
          #print(lambda)
          
          
          ###################### Using Selected Lambda to Run Lasso ############
          glmnet2 <- glmnet(X, y=train$Y, alpha=1, lambda = lambda)
          c <- coef(glmnet2)
          
          ##################### Saving Selected Features From Lasso ############
          inds<-which(c!=0)
          
          #remove intercept from the list
          tmp_variables <- tail(row.names(c)[inds], -1)
          #print(tmp_variables)
          
          variables<-c(variables,tmp_variables)
        }
        ################# Selecting Features Based On Frequency ############
        freq = sort(table(variables),decreasing=TRUE)/num_FS
        #print(freq)
        
        stable_Var = names(which(freq>0.2))  #### the 0.6 can be changed ######
        #print(stable_Var)
        
        #selectedVars <- append(selectedVars, length(stable_Var))
        selectedVars <- append(selectedVars, stable_Var)
        ################## Use Stable Features on RF & SVM ##################### 
        tempTr <- train[, (names(train) %in% stable_Var)]
        
       
        tempTrain <- cbind.data.frame(train$Y, tempTr)
        
        
        tempTr <- test[, (names(test) %in% stable_Var)]
        
        
        tempTest <- cbind.data.frame(test$Y, tempTr)
        
        
        colnames(tempTrain)[1] <- "Y"
        colnames(tempTest)[1] <- "Y"
        
        ######################## Run SVM and RF #########################
        # append all test Y for AUC calculation   
        #Y <- append(Y, tempTest$Y)
        #For true Y and permut Y, we put them in separate lists but name both Y at the end for readability
        if (m == 1){
          trueY <- append(trueY, tempTest$Y)
          Y_train <- train$Y ##this Y will be used to train SVM and RF
        }
        if (m == 2){
          permutY <- append(permutY, tempTest$Y)
          Y_train <- train$Y ##this Y will be used to train SVM and RF
        }
        
        #svm
        colnames(tempTrain)[1]<- "Y_train"
        svmfit = svm(Y_train~ ., data = tempTrain , kernel="linear", cost=19, scale=FALSE,tolerance = 0.05)
        yhat.SVM = predict(svmfit, newdata = tempTest[, -1])
        
        #rf
        RFfit <- randomForest(Y_train~ ., data = tempTrain, importance=TRUE, ntree = 10)
        yhat.RF = predict(RFfit, newdata = tempTest[, -1])
        
        #################### Store yhat in different array depend on the mode #############
        if (m == 1){
          append.SVM <- append(append.SVM, yhat.SVM)
          append.RF <- append(append.RF, yhat.RF)
        }
        if (m == 2){
          append.SVM_permut <- append(append.SVM_permut, yhat.SVM)
          append.RF_permut <- append(append.RF_permut, yhat.RF)
        }
      }
    }

    # get auc score for each replication
    #print("printing AUCs")
    aucRF[replicateCV] <- auc(trueY, append.RF)
    #print(aucRF[replicateCV])
    aucSVM[replicateCV] <- auc(trueY, append.SVM)
    #print(aucSVM[replicateCV])
    
    aucRF_permut[replicateCV] <- auc(permutY, append.RF_permut)
    #print(aucRF_permut[replicateCV])
    aucSVM_permut[replicateCV] <- auc(permutY, append.SVM_permut)
    #print(aucSVM_permut[replicateCV])
}



print("selected vars: ")
print(selectedVars)
#check_df=as.data.frame(X[,c(selectedVars)])

#aucRF
#median(aucRF)
#print(median(aucRF))

#aucSVM
#median(aucSVM)
print(paste0("median aucSVM: ", median(aucSVM)))

#aucSVM_permut
#median(aucSVM_permut)
print(paste0("median aucSVM_permute: ",median(aucSVM_permut)))

print(paste0("median aucRF: ",median(aucRF)))
print(paste0("median aucRF_permute: ", median(aucRF_permut)))
#aucRF_permut
#median(aucRF_permut)
#print(median(aucRF_permut))

f=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/Lasso/miRNA/feats/feats_l0.6.csv")
rf_df=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/Lasso/miRNA/RF_auc_0.6.csv")
#write.csv(aucRF, "/ix/djishnu/Trirupa/ABomics.Prj/Lasso/miRNA/RF_auc_0.7.csv")
#write.csv(aucSVM, "/ix/djishnu/Trirupa/ABomics.Prj/Lasso/miRNA/SVM_auc_0.7.csv")

#write.csv(aucRF_permut, "/ix/djishnu/Trirupa/ABomics.Prj/Lasso/miRNA/RF_aucPermute_0.7.csv")
#write.csv(aucSVM_permut, "/ix/djishnu/Trirupa/ABomics.Prj/Lasso/miRNA/SVM_aucPermute_0.7.csv")

##preparing for plotting boxplots
#library(tidyr)
##loading the datasets
#RF_csv=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/Lasso/miRNA/RF_auc_lambda0.3.csv")
#SVM_csv=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/Lasso/miRNA/SVM_auc_lambda0.3.csv")

##changing column names
#RF_df=as.data.frame(RF_csv[-21,-1]) ###removing the first column that has repCV number
#SVM_df=as.data.frame(SVM_csv[-21,-1]) ## also removing the 21st element since it was the column header while concatenating actual +permuted files

##now we need to label the Actual and Permuted data points as a second column
#Type_actual <-rep("Actual",times=20)
#Type_permuted =rep("Permuted",times=20)
#RF_df$Type=c(Type_actual,Type_permuted)
#SVM_df$Type=c(Type_actual,Type_permuted)
#colnames(RF_df)[1]<- "AUC"
#colnames(SVM_df)[1]<- "AUC"

### if in any run no features are selected then select any 3 at random so that loop doesn't break####
if (len == 1) 
{
  randomSelect <- sample(ncol(data), 3)
  tmp_variables <- row.names(c)[randomSelect]
} else
{
  variables<-row.names(c)[inds]
  tmp_variables <- variables[2:len]
}