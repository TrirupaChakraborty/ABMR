### to calculate the univariate pvals with Wilcox test
rm(list = ls())
cat("\014")
library(plyr)

Data=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/data/class2_SEASchistoAg_withLabel.csv")
x=Data[,-1]
y=Data[,1]
#y=Data[,1]
#y=mapvalues(Data[,"Y"],from=c(2,1), to=c(0,1))

#y=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/miRNA_data/Ylabels_maxExpression_miRNA.csv", row.names = 1)
var_name=colnames(x)

result_df=data.frame(name=NULL,pval=NULL)
for (i in var_name){res=wilcox.test(x[,i]~as.matrix(y)); result_df<-rbind(result_df,data.frame(name=i,pval=res$p.value))}
ii <- order(as.numeric(result_df$pval),decreasing = F)
result_df_sorted<- result_df[ii,]

## finding the univariate pvals of high loaded ER features
slide_result=readRDS(file="hierER.class2.SchsitoAg.MFI50.d0.1_lam1.spec0.1.auc.rds")
er_result=readRDS("final_delta_0.1_lambda_1.rds")
A_abs <-  apply(er_result$A,c(1,2),abs)
A_k <- as.data.frame(A_abs) %>% select(10)
marginal_vars <- rownames(subset(A_k,A_k > 0.08)) 

vars_check=lapply(X = marginal_vars, FUN = function(t) gsub(pattern = ".", replacement = "-", x = t, fixed = TRUE))
vars_check
max_expX=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/miRNA_data/ValidationCohort1/miRNA_data/ER_input/X_maxExp_miRNA.csv", row.names=1)
t=which(colnames(max_expX)=="miR.941")
pval_var=wilcox.test(max_expX[,t]~as.matrix(y))
pval_var

VCresult_df=data.frame(name=NULL,pval=NULL)
for (i in vars_check){res=wilcox.test(x[,i]~as.matrix(y)); VCresult_df<-rbind(VCresult_df,data.frame(name=i,pval=res$p.value))}
ii <- order(as.numeric(result_df$pval),decreasing = F)
result_df_sorted<- result_df[ii,]