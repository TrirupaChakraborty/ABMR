rm(list = ls())
cat("\014")

library(EssReg)
library(doParallel)
library(dplyr)
library(pROC)
library(ROCR)
library(SLIDE)
registerDoParallel(detectCores())

y=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/miRNA_data/Ylabels_maxExpression_miRNA.csv", row.names = 1)
z=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/EssReg-main/ER_kidneyTransplant/miRNA/Filtered_miRNA/maxExpression_filtered/fdr1.0/del0.02_lam1.0/zmatrix_maxEXP.miRNA.delta0.02_lam1.0.corr.csv", row.names=1)

elbow_vars<- function(z, y, niter = 1000, fdr = 0.1, parallel = TRUE) {
  ## run second order knockoffs
  results <- secondKO(z = z,
                      y = y,
                      fdr = fdr,
                      niter = niter,
                      parallel = parallel)
  
  ## find most frequently selected variable
  #ii <- which.max(results$tab_data)
  selected_vars <- names(results$tab_data)
  
  return (results)
}

elbow_vars_res=elbow_vars(z,y)
elbow_vars_res
