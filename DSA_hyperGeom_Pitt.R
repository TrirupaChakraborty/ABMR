### DSA specificity analysis for DC (Pitt) cohort ####

## loading libraries and fxs
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
source("/ix/djishnu/Trirupa/ABomics.Prj/scripts/DSA_hyperGeom_Fx.R")

setwd("/ix/djishnu/Trirupa/ABomics.Prj")
####### Loading datasets and dfs #####
patid.serotype_df = read_excel("./input_data/PITT_Historical_DSA_Full_5_20.xlsx", 
                               sheet = "Patid_Serotype")
new_cols= c("PAT_ID","all", "DSA")
colnames(patid.serotype_df)= new_cols
patid.serotype_df= patid.serotype_df %>% select(-"all")

meta_df= read.csv("/ix/djishnu/Marisa/ABomics/pitt_cohort/meta_36.csv")
patid.status_df= meta_df[,c("ID","GRP")]
new_cols2= c("PAT_ID", "Group")
colnames(patid.status_df) = new_cols2

## combining pat_id-status(group)-serotype 
pat.grp.sero_df= merge(patid.status_df, patid.serotype_df, by="PAT_ID")
#colnames(pat.grp.sero_df)[colnames(pat.grp.sero_df) == "Status"] = "Group"
pat.grp.sero_df$Group[pat.grp.sero_df$Group == "DSA"] = "NR"
pat.grp.sero_df=pat.grp.sero_df %>% select(-"PAT_ID")
#write.csv(pat.grp.sero_df, file="./input_data/Pitt_serotype.status.csv",row.names=FALSE)

DSA_agdf <- separate_rows(pat.grp.sero_df, DSA, sep = ",")
DSA_agdf$DSA <- str_trim(DSA_agdf$DSA)
ABMR_agdf= DSA_agdf[DSA_agdf$Group == "ABMR",]
NR_agdf= DSA_agdf[DSA_agdf$Group == "NR",]
#VC_HLAdf= VC_HLAdf[!is.na(VC_HLAdf$`HTG #`),c("DSA","Group")]

##### Loading the Antigen set that was used in the assay #######
BeadAg_path <- "./input_data/onelambda_beadantigen_regionsEDITED.xlsx"
c1c2.Ag_df= loadc1c2df(BeadAg_path)

################### Now looking at the enrichment of AbMR or NR associated DSA-Ags in all sig LFs #######

####### compiling all the sig LFs and the features in them ######
## for Early Ab-omics model###
E_Zlist=c("Z1","Z12") 
E_Zfeat_path= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/SLIDE/"
E_allZ_df= load_allZ(Z_list= E_Zlist, FeatList_path=E_Zfeat_path)

## for Late Ab-omics model###
L_Zlist=c("Z1","Z2","Z7","Z10") 
L_Zfeat_path= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/SLIDE/"
L_allZ_df= load_allZ(Z_list= L_Zlist, FeatList_path= L_Zfeat_path)

##### performing hypergeometric test ###################################
## parameters ##
# no. of DSA Ags = m
# total no.of Ags tested = x
# no. of non-DSA Ags = n
# no. of Ags in LF = k
# no. of DSA Ags OF the no. of Ags in LF = q

### for Early
E_ABMR_hypergdf = run_hypergeom1(ABMR_agdf, c1c2.Ag_df, E_allZ_df)
E_NR_hypergdf = run_hypergeom1(NR_agdf, c1c2.Ag_df, E_allZ_df)

### for Late 
L_ABMR_hypergdf = run_hypergeom1(ABMR_agdf, c1c2.Ag_df, L_allZ_df)
L_NR_hypergdf = run_hypergeom1(NR_agdf, c1c2.Ag_df, L_allZ_df)

check_df= read.csv("/ix/djishnu/Marisa/ABomics/pitt_cohort/abmr_36/inputs/abmr_x.csv")
Echeck_abdf= read.csv("/ix/djishnu/Marisa/ABomics/pitt_cohort/ab_36_cp/early_y.csv")
Lcheck_abdf= read.csv("/ix/djishnu/Marisa/ABomics/pitt_cohort/ab_36_cp/late_y.csv")