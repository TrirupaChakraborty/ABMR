metadf=read.csv("PittCohort_kidneytransplant_ABMR_patientID_clinicalData.csv")
mydf=metadf[,c(1,2,5,9)]
mydf
hist(mydf[,3], xlab="time to tx (mo)")
DSAdf= mydf[(mydf$GRP == "DSA"), ]
ABMRdf=mydf[(mydf$GRP == "ABMR"), ]

hist(DSAdf[,3], xlab="time to tx (mo)", ylab="no. of patients",main="DSA only")
hist(ABMRdf[,3], xlab="time to tx (mo)", ylab="no. of patients",main="ABMR only")
summary(DSAdf)
summary(ABMRdf)

DSAdf_longTx=DSAdf[DSAdf["Time.to.Tx..mo."]>=median((mydf$Time.to.Tx..mo.)),]
DSAdf$NewGrp=ifelse(DSAdf$Time.to.Tx..mo.< 1.65, "DSA_short", "DSA_long")

#DSAdf_shortTx=DSAdf[DSAdf["Time.to.Tx..mo."]<median((mydf$Time.to.Tx..mo.)),]

ABMRdf_longTx=ABMRdf[ABMRdf["Time.to.Tx..mo."]>=median((mydf$Time.to.Tx..mo.)),]
ABMRdf_shortTx=ABMRdf[ABMRdf["Time.to.Tx..mo."]<median((mydf$Time.to.Tx..mo.)),]
ABMRdf$NewGrp=ifelse(ABMRdf$Time.to.Tx..mo.< 1.65, "ABMR_short", "ABMR_long")

fulldf=rbind(DSAdf,ABMRdf)
value_map <- c("DSA_long" = 0, "DSA_short" = 1, "ABMR_long" = 2, "ABMR_short" = 3)
fulldf$Y <- value_map[match(fulldf$NewGrp, names(value_map))]

abmir_x=read.csv("/ix/djishnu/Marisa/ABomics/AB_miRNA/input/AB_miRNA_x.csv",row.names=1)
abID=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/AB_patID_tempdata_1sep2023.csv")
abID=abID[,-c(3,4)]

### first merging abMIR data with the temp data ###
abmir_patID_df <- merge(abID,abmir_x, by = "IgG_ClassI.1.6.", all = FALSE)
### next merging abMIR_patID data with the new grps in fulldf dataset ###
fulldf_subset=fulldf[,c(1,3,4,5,6)]
colnames(fulldf_subset)[colnames(fulldf_subset) == "CUSTOMER_ID"] <- "Customer.ID"
abmir_patIDNewGrp_df = merge(fulldf_subset,abmir_patID_df, by="Customer.ID", all=FALSE)
### the groups are based on median Tx on teh whole data (long >1.65, short<1.65) ###
write.csv(abmir_patIDNewGrp_df,"/ix/djishnu/Trirupa/ABomics.Prj/input_data/abmir_patID_TxGrping_1sep2023.csv")
#checkdf=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/input_data/abmir_patID_TxGrping_1sep2023.csv")


