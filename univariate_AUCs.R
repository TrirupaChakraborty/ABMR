#### This code gets AUCs of features picked up by ER selected LFs in both discovery and validation cohorts ### 

library(dplyr)
library(pROC)
library(ROCR)
library(matrixStats)
library(ggplot2)
library(glmnet)


#### for EARLY ABOMICS features for VC#####
#VC_Early_x = read.csv("/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/early_x.csv",row.names=1)
#VC_Early_y =  read.csv("/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/early_y.csv",row.names=1)
Early_x= read.csv("/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/x.csv", row.names=1)
Early_y= read.csv("/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/y.csv", row.names=1)

## loading the features from the sig. LFs ###
featlist_dir="/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/SLIDE/"
Early_sigLF_feats= readRDS(paste0(featlist_dir,"plotSigGenes_data.RDS"))

#Zfeats= read.csv(paste0(featlist_dir,"feature_list_Z1.txt"))
univar_AUC.fx= function(x,y,sigLF_feats){
  result_df=data.frame(name=character(),pval=numeric(), corXiY=numeric(), auc_canada=numeric(), auc_pitt=numeric())
  
  for (feat in sigLF_feats$names){
    
    mannU_pval=wilcox.test(x[,feat]~as.matrix(y))$p.value
    auc_VC=(pROC:::roc(as.factor(y$Y),x[,feat], direction="<"))$auc
    auc_DC=sigLF_feats$AUCs[sigLF_feats$names == feat]
    result_df= rbind(result_df, data.frame(name=feat, pval=mannU_pval,
                                           corXiY= cor(x[,feat],y),auc_canada=auc_VC,auc_pitt=auc_DC))
  }
  print(result_df)
  #sorted_result_df <- result_df[order(result_df$auc_pitt,decreasing=TRUE),]
  return(result_df)
}

#univar_AUC.fx <- function(x, y, sigLF_feats) {
  # Initialize an empty data frame with the correct column types
  result_df <- data.frame(
    name = character(),
    pval = numeric(),
    corXiY = numeric(),
    auc_canada = numeric(),
    auc_pitt = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (feat in sigLF_feats$names){
    mannU_pval <- wilcox.test(x[[feat]] ~ as.matrix(y))$p.value
    auc_VC <- pROC::roc(as.factor(y$Y), x[[feat]], direction = "<")$auc
    auc_DC <- sigLF_feats$AUCs[sigLF_feats$names == feat]
    corrXiY <- cor(x[[feat]], y$Y)
    
    result_df <- rbind(result_df, data.frame(
      name = feat,
      pval = mannU_pval,
      corXiY = corrXiY,
      auc_canada = auc_VC,
      auc_pitt = auc_DC
    ))
  }
  
  sorted_result_df <- result_df[order(result_df$auc_pitt, decreasing = TRUE), ]
  #saveRDS(sorted_result_df, "/ix/djishnu/Trirupa/ABomics.Prj/Validation/Early_univariateAUCs.RDS")
  return(result_df)
#}

Early_uniAUC= univar_AUC.fx(x=Early_x, y= Early_y, sigLF_feats= Early_sigLF_feats)


