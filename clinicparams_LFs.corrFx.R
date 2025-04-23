############## fx that will be used in correlation plots of Zs with clinical parameters ##############
library(EssReg)
library(ggplot2)
library(yaml)
library(dplyr)
library(pROC)
library(ROCR)
library(matrixStats)
library(ggpubr)
library(reshape2)
library(SLIDE)

###### significant Z extraction #####
predZ <- function(x, er_res) {
  A_hat <- er_res$A
  C_hat <- er_res$C
  Gamma_hat <- er_res$Gamma
  Gamma_hat <- ifelse(Gamma_hat == 0, 1e-10, Gamma_hat)
  Gamma_hat_inv <- diag(Gamma_hat ** (-1))
  G_hat <- crossprod(A_hat, Gamma_hat_inv) %*% A_hat + solve(C_hat)
  Z_hat <- x %*% Gamma_hat_inv %*% A_hat %*% MASS::ginv(G_hat)
  return (Z_hat)
}

getLFdf<- function(SLIDE_res, z_matrix, interactions = TRUE){
  # marginal data
  sigK <- SLIDE_res$marginal_vals
  sigMargData <- z_matrix[,sigK]
  
  if(!interactions){
    return(as.data.frame(sigMargData))
  }
  
  # interaction
  sigIn <- SLIDE_res$SLIDE_res$interaction_vars
  if(!is.null(sigIn)){
    IntData <- pairwiseInteractions(sigK,z_matrix)
    sigIntData <- IntData$interaction[ ,sigIn]
    final_mtx <- cbind(as.data.frame(sigMargData), as.data.frame(sigIntData))
    return(final_mtx)
  }
  return(as.data.frame(sigMargData))
}

load_valSigZ <- function(valX_path,valY_path,yaml_path){
  val_x=as.matrix(read.csv(valX_path, row.names=1))
  val_y=as.matrix(read.csv(valY_path, row.names=1))
  input <- read_yaml(yaml_path)
  
  train_x <- as.matrix(read.csv(input$x_path, row.names=1))
  train_y <- as.matrix(read.csv(input$y_path, row.names=1))
  
  # if validation dataset has aditional features, remove them
  val_not_in_train <- colnames(val_x)[which(!colnames(val_x) %in% colnames(train_x))]
  val_x <- scale(as.matrix(val_x))
  train_x <- scale(as.matrix(train_x))
  
  
  er_results = readRDS(list.files(input$out_path, full.names = T,  pattern = "final_"))
  slide_res = readRDS(list.files(input$out_path,recursive = T, full.names = T, pattern = "slide_res")[1])
  
  #### predicting the LFs in the Validation cohort using the group structure (in er_result) from discover cohort ###
  val_z <- predZ(val_x, er_results)
  colnames(val_z) <- paste0("Z", c(1:ncol(val_z)))
  ##### getting the significant LFs (val_sigZ) #####  
  val_sigZ = getLFdf(slide_res, val_z, interactions=TRUE)
  
  return(val_sigZ)
}

#### plotting function ######
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr) 

create.LF_clin.df= function(meta_df, sigZ_df, clin_param,rename_param){
  
  sigZ_patid <- sigZ_df 
  sigZ_patid['pat_id']=rownames(sigZ_patid) ### since rownames of sigZdf are the patient ids
  clinparam_df = meta_df[meta_df$`patient id`%in% sigZ_patid$pat_id ,
                         c('patient id', "current study group", clin_param)]
  names(clinparam_df)[names(clinparam_df)== 'patient id'] = "pat_id"
  clinparam_sigZdf= base::merge(sigZ_patid, clinparam_df, by ="pat_id")
  names(clinparam_sigZdf)[names(clinparam_sigZdf)== clin_param] = rename_param
  Z_param.col_lis= c(colnames(sigZ_df), rename_param)
  sub.clinparam_sigZdf = clinparam_sigZdf[ , colnames(clinparam_sigZdf) %in% Z_param.col_lis]
  
  return(sub.clinparam_sigZdf)
}

cor_scatterplt= function(Z_param.df,Z_name,clinc_param, plt_title, x_title, y_title,xaxis_len){
  ggplot( Z_param.df, aes_string( x=Z_name, y=clinc_param))+ 
    geom_point(color = "#809c13", size = 3, alpha = 1.0)+ 
    geom_smooth(method = "lm", se = FALSE, color = "#ab0000") +
    stat_cor(method = "spearman", label.x = xaxis_len, label.y = 30) +
    labs(title = plt_title,
         x = x_title,
         y = y_title) + theme_classic() +
    theme(plot.title = element_text(hjust=0.5))}