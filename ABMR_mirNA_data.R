rm(list = ls())
cat("\014")

setwd("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/")
data_path="KidneyTransplant_DataProbes_v3_after_imputation_class2vs1_withPatID.csv"
data=read.csv(data_path,header=TRUE,sep=",")
Y <- data[,1]
pat_id <- data[,2]
scaledData <- scale(data[,3:ncol(data)])
final_data <- cbind.data.frame(Y,pat_id,scaledData)

write.csv(final_data,"scaled_KidneyTransplant_DataProbes_v3_after_imputation_class2vs1_withPatid.csv")

miRNA_x=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/miRNA_data/X_prescale_miRNA_class1vs2.csv", row.names=1)
scaled_miRNA_x=scale(miRNA_x, T, T)
