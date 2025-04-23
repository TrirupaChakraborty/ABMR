########## this script performs hypergeometric test for DSAs found in the VC #####
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(purrr)
library(xlsx)
setwd("/ix/djishnu/Trirupa/ABomics.Prj/")
source("/ix/djishnu/Trirupa/ABomics.Prj/scripts/DSA_hyperGeom_Fx.R")
source("/ix/djishnu/Trirupa/ABomics.Prj/scripts/Ag_hyperg_Fx.R")


sero.grpdf= read_xlsx("./input_data/Toronto-DSA_editedAg_V1.xlsx", sheet=3)
DSA_agdf <- separate_rows(sero.grpdf, DSA, sep = ",")
DSA_agdf$DSA <- str_trim(DSA_agdf$DSA)
DSA_agdf$DSA= gsub("Cw1", "C1",DSA_agdf$DSA )
DSA_agdf$DSA= gsub("Cw2", "C2",DSA_agdf$DSA )

ABMR_agdf= DSA_agdf[DSA_agdf$Group == "ABMR",]
NR_agdf= DSA_agdf[DSA_agdf$Group == "NR",]


####### compiling all the sig LFs and the features in them ######
## for Early Ab-omics model###
E_Zlist=c("Z1","Z12") 
E_Zfeat_path= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/SLIDE/"
E_allZ_df= load_allZ(Z_list= E_Zlist, FeatList_path=E_Zfeat_path)
colnames(E_allZ_df)[colnames(E_allZ_df)== "Ag_pool"] <- "Ag_class"

## for Late Ab-omics model###
L_Zlist=c("Z1","Z2","Z7","Z10") 
L_Zfeat_path= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/SLIDE/"
L_allZ_df= load_allZ(Z_list= L_Zlist, FeatList_path= L_Zfeat_path)
colnames(L_allZ_df)[colnames(L_allZ_df)== "Ag_pool"] <- "Ag_class"

##### performing hypergeometric test ###################################

#### for EARLY ####
E.hyperg_NR= Ag_hyperg.fx(allZ_df= E_allZ_df, grp_agdf= NR_agdf)
E.probe_agsumdf_NR= unique(E.hyperg_NR %>% arrange(desc(numQueryAg_probe)))
E.probe_agsumdf_NR= E.probe_agsumdf_NR[E.probe_agsumdf_NR$pval <= 0.05,]
write.csv(E.probe_agsumdf_NR,"./Figures/Abomics_Early/VC.Ag_EarlyNR.hypergeomV2.csv", row.names=F, quote=F)


E.hyperg_ABMR= Ag_hyperg.fx(allZ_df= E_allZ_df, grp_agdf= ABMR_agdf)
E.probe_agsumdf_ABMR= unique(E.hyperg_ABMR %>% arrange(desc(numQueryAg_probe)))
E.probe_agsumdf_ABMR= E.probe_agsumdf_ABMR[E.probe_agsumdf_ABMR$pval <= 0.05,]
write.csv(E.probe_agsumdf_ABMR,"./Figures/Abomics_Early/VC.Ag_EarlyABMR.hypergeomV2.csv", row.names=F, quote=F)


#### for LATE ####
L.Ag_hyperg_NR= Ag_hyperg.fx(allZ_df= L_allZ_df, grp_agdf= NR_agdf)
L.probe_agsumdf_NR= unique(L.Ag_hyperg_NR %>% arrange(desc(numQueryAg_probe)))
L.probe_agsumdf_NR= L.probe_agsumdf_NR[L.probe_agsumdf_NR$pval <= 0.05,]
write.csv(L.probe_agsumdf_NR,"./Figures/Abomics_Late/VC.Ag_LateNR.hypergeomV2.csv", row.names=F, quote=F)


L.Ag_hyperg_ABMR= Ag_hyperg.fx(allZ_df= L_allZ_df, grp_agdf= ABMR_agdf)
L.probe_agsumdf_ABMR= unique(L.Ag_hyperg_ABMR %>% arrange(desc(numQueryAg_probe)))
L.probe_agsumdf_ABMR= L.probe_agsumdf_ABMR[L.probe_agsumdf_ABMR$pval <= 0.05,]
write.csv(L.probe_agsumdf_ABMR,"./Figures/Abomics_Late/VC.Ag_LateABMR.hypergeomV2.csv", row.names=F, quote=F)


