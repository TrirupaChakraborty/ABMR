#### this code gets an odds ration between the enrichment of memory vs de novo DSAs #########
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
source("/ix/djishnu/Trirupa/ABomics.Prj/scripts/DSA_hyperGeom_Fx.R")
source("/ix/djishnu/Trirupa/ABomics.Prj/scripts/Ag_hyperg_Fx.R")
source("/ix/djishnu/Marisa/ABomics/pitt_cohort/antigens/Plotting_helper_fxns.R")

setwd("/ix/djishnu/Trirupa/ABomics.Prj")

hist <- read_excel("/ix/djishnu/Marisa/ABomics/pitt_cohort/antigens/5_28/Historical_ABMR_DSA_CLEAN.xlsx",sheet = 2)


##################################################################
# Historical information filtering based on Txn date
##################################################################

# WITHIN 3 MONTH TIMEFRAME of Tx  --------------------memory DSA
hist_tx_3mo <- hist[hist$DAYS_TX_TO_RESULT <= 90,]

pt_tx_3mo = hist_tx_3mo[,c("PAT_ID", "Status","SEROTYPE")]
pt_tx_3mo= na.omit(pt_tx_3mo)
                          
                          
# POST 3 MONTH TIMEFRAME of Tx -----------------------de-novo DSA
hist_tx_post <- hist[hist$DAYS_TX_TO_RESULT > 90,]

pt_tx_post = hist_tx_post[,c("PAT_ID", "Status","SEROTYPE")]
pt_tx_post=na.omit(pt_tx_post)


######## dividing the patient DSAs into de-novo and memory for ONLY REJECTORS ###
denovoDSA = pt_tx_3mo[pt_tx_3mo$Status == "ABMR", -1]
denovoDSA= separate_rows(denovoDSA, SEROTYPE, sep = ",")
denovoDSA$SEROTYPE = str_trim(denovoDSA$SEROTYPE)


memDSA= pt_tx_post[pt_tx_post$Status == "ABMR", -1] 
memDSA= separate_rows(memDSA, SEROTYPE, sep = ",")
memDSA$SEROTYPE = str_trim(memDSA$SEROTYPE)

##### 
#commonDSA= unique(intersect(denovoDSA$SEROTYPE,memDSA$SEROTYPE))
onlyin_memDSA = anti_join(memDSA,denovoDSA, by= "SEROTYPE")
onlyin_denovoDSA = anti_join(denovoDSA,memDSA, by= "SEROTYPE")
#print(commonDSA)

print(paste0("# of denovo=",length(unique(denovoDSA$SEROTYPE))))
print(paste0("# of mem=",length(unique(memDSA$SEROTYPE))))
print(paste0("#common=", length(unique(commonDSA))))
print(paste0("only in mem DSA= ",length(unique(onlyin_memDSA$SEROTYPE))))
print(paste0("only in denovo DSA= ",length(unique(onlyin_denovoDSA$SEROTYPE))))


###### assay probe-ag mapping 
ag_mappingMT= readRDS("./MARISA/ag_mapping_matrix.rds")
ag_mappingMT <- ag_mappingMT %>%
  select(-matches("Mica"))
rows_to_keep <- !grepl("MICA", rownames(ag_mappingMT))
ag_mappingMT= ag_mappingMT[rows_to_keep,]
ag_mappingMT$Ag_class= rownames(ag_mappingMT) 


####### compiling all the sig LFs and the features in them ######
## for Early Ab-omics model###
E_Zlist=c("Z1","Z12") 
E_Zfeat_path= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/SLIDE/"
E_allZ_df= load_allZ(Z_list= E_Zlist, FeatList_path=E_Zfeat_path)
colnames(E_allZ_df)[colnames(E_allZ_df)== "Ag_pool"] <- "Ag_class"



##### first performing analysis for late model 
####### compiling all the sig LFs and the features in them ######
## for Late Ab-omics model###
L_Zlist=c("Z1","Z2","Z7","Z10") 
L_Zfeat_path= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/SLIDE/"
L_allZ_df= load_allZ(Z_list= L_Zlist, FeatList_path= L_Zfeat_path)
colnames(L_allZ_df)[colnames(L_allZ_df)== "Ag_pool"] <- "Ag_class"





######### for the ODDS Ratio ###
## Pmem = # of mem DSA in LF / # of mem DSA in pool
## Pdnovo = # of de novo DSA in LF / # of de novo DSA in pool

common_DSA = unique(intersect(mem_DSAdf$SEROTYPE,dnovo_df$SEROTYPE))
common.ABMR.DSA= unique(intersect(mem.ABMR_agdf$SEROTYPE,dnovo.ABMR_agdf$SEROTYPE))

                        