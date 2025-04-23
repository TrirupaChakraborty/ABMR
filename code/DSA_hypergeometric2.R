### DSA specificity analysis for VC cohort ####

## loading libraries and fxs
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
source("/ix/djishnu/Trirupa/ABomics.Prj/scripts/DSA_hyperGeom_Fx.R")

####### Loading datasets and dfs #####
VC_HLApath = "/ix/djishnu/Trirupa/ABomics.Prj/input_data/Toronto-DSA_editedAg_V1.xlsx"
VC_HLAdf = read_excel(VC_HLApath, sheet=4)
DSA_agdf <- separate_rows(VC_HLAdf, DSA, sep = ", ")
ABMR_agdf= DSA_agdf[DSA_agdf$Group == "ABMR",]
NR_agdf= DSA_agdf[DSA_agdf$Group == "NR",]
#VC_HLAdf= VC_HLAdf[!is.na(VC_HLAdf$`HTG #`),c("DSA","Group")]

##### Loading the Antigen set that was used in the assay #######
BeadAg_path <- "/ix/djishnu/Trirupa/ABomics.Prj/input_data/onelambda_beadantigen_regionsEDITED.xlsx"
c1c2.Ag_df= loadc1c2df(BeadAg_path)

## NOTE: this analysis ignores DQA1 and DQA5 (5 entries in VC_HLAdf) since c1c2.Ag_df don't have DQA1, 5.

################### Now looking at the enrichment of AbMR or NR associated DSA-Ags in all sig LFs #######

####### compiling all the sig LFs and the features in them ######
## for Late Ab-omics model###
L_Zlist=c("Z1","Z2","Z7","Z10") 
L_Zfeat_path= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/SLIDE/"
L_allZ_df= load_allZ(Z_list= L_Zlist, FeatList_path= L_Zfeat_path)

## for Early Ab-omics model###
E_Zlist=c("Z1","Z12") 
E_Zfeat_path= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/SLIDE/"
E_allZ_df= load_allZ(Z_list= E_Zlist, FeatList_path=E_Zfeat_path)

##### performing hypergeometric test ###################################
## parameters ##
# no. of DSA Ags = m
# total no.of Ags tested = x
# no. of non-DSA Ags = n
# no. of Ags in LF = k
# no. of DSA Ags OF the no. of Ags in LF = q

### for Late 
L_ABMR_hypergdf = run_hypergeom1(ABMR_agdf, c1c2.Ag_df, L_allZ_df)
L_NR_hypergdf = run_hypergeom1(NR_agdf, c1c2.Ag_df, L_allZ_df)


### for Early
E_ABMR_hypergdf = run_hypergeom1(ABMR_agdf, c1c2.Ag_df, E_allZ_df)
E_NR_hypergdf = run_hypergeom1(NR_agdf, c1c2.Ag_df, E_allZ_df)


############# Calculating ENRICHMENT ###################







