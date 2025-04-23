setrm(list = ls())
cat("\014")

library(EssReg)
library(doParallel)
library(dplyr)
library(pROC)
library(ROCR)
library(SLIDE)

registerDoParallel(detectCores())

x=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/miRNA_data/X_maxEXP.sub3.6.var25.filtered.miRNA.csv", row.names=1)
y=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/miRNA_data/Ylabels_maxExpression_miRNA.csv", row.names = 1)
x_std <- scale(x, T, T)

##loading er result
er_result=readRDS("final_delta_0.01_lambda_0.5.rds")
##getting zs
zs <- predZ(x_std, er_result)
colnames(zs)=paste0("Z",seq(1,ncol(zs)))
#setwd("/ix/djishnu/Trirupa/ABomics.Prj/EssReg-main/ER_kidneyTransplant/miRNA_ABMR/maxExpression_filtered/fdr1.0/del0.065_lam1.0/")
write.csv(zs,file="zmatrix_maxEXP.sub3.6.var25.filtered.miRNA.d0.01_lam0.5.auc.csv",row.names=TRUE)
#zs=read.csv("zmatrix_maxEXP.miRNA_ABMR.delta0.065_lam1.0.corr.csv",row.names=1)

## running SLIDE (without interaction terms) ##significant latent factors
SLIDE_marginal_res =SLIDE(zs,y,niter=5000,do_interacts = F, spec=0.2)
saveRDS(SLIDE_marginal_res, file="hierER_maxEXPmirna_ABMR.del0.065_lam1.0.rds")
#### getting marginals at different loading thresholds
A_abs <-  apply(er_result$A,c(1,2),abs)
A_k <- as.data.frame(A_abs) %>% select(31)
var_name <- rownames(subset(A_k,A_k>0.1)) 
var_name
