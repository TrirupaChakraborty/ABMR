rm(list = ls())
cat("\014")

setwd("/ix/djishnu/Trirupa/ABomics.Prj/")

library(EssReg)
library(doParallel)
registerDoParallel(detectCores())

##checking datasets 
x_features=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/X_prescale_features_Kidney_0vs2.csv", header=TRUE, sep=",")
y_labels=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/EssReg/Y_labels_KidneyTransplant.csv", header=TRUE)

#####running ER
yaml_path="/ix/djishnu/Trirupa/ABomics.Prj/EssReg-main/ER_kidneyTransplant/class0vs2/fdr0.05/pipeline1_kidney.yaml"
pipelineER1(yaml_path,"all")

yaml_path2="/ix/djishnu/Trirupa/ABomics.Prj/EssReg-main/ER_kidneyTransplant/class0vs2/fdr0.05/pipeline2_kidney.yaml"
pipelineER2(yaml_path2)

yaml_path3="/ix/djishnu/Trirupa/ABomics.Prj/EssReg-main/ER_kidneyTransplant/class0vs2/fdr0.05/pipeline3_kidney.yaml"
pipelineER3(yaml_path3)

plainER_yaml=""
plainER_out=parseRun(yaml_path)

saveRDS(plainER_out, file = "/ix/djishnu/Trirupa/ABomics.Prj/EssReg-main/BenchMark_plainER/plainER_output_repcv200fdr0.1/er_output_0.1.rds")


##steps to extract the Zs from the betas
x=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/X_prescale_features_KidneyTransplant.csv", row.names=1) ##loading x
#standardising x
x_std <- scale(x, T, T) #mean=0, std dev=1
#to get the top 10 betas
beta_top10=sigBetas(final_delta_0.1_lambda_0.01$beta, 0.1)
#to get all Zs
zs <- predZ(x_std, final_delta_0.1_lambda_0.01) ##er_res=output of plain ER, in this case the output of final delta and lambda
#to get Zs for top betas
top_betas=c(beta_top10$pos_sig,beta_top10$neg_sig)
top_zs <- zs[, top_betas]
#to rename columns
colnames(top_zs) <- paste0("Z", top_betas)
#to get cluster memberships: 
clusters <- readER(final_delta_0.1_lambda_0.01)
#to see a given Z (ex: Z31): 
Z5 <- clusters$clusters[[5]]
Z5_colind=c(Z5$pos, Z5$neg) #to go from list to vector of column indices from Z31
##to go from indices to column names: 
Z5_names <- indName(Z5_colind, colnames(x), F)

#to get just pure variables in Z31:
a_mat <- final_delta_0.1_lambda_0.01$A
pure_vars <- which(abs(a_mat[, 5]) >= 1) ## these are the indices. abs = absolute values
pure_vars_names <- indName(pure_vars, colnames(x), F)


#######
