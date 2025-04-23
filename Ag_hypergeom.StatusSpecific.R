### this script uses patient DSAs divided by early and late and for each status (early/ late) finds the 
## enrichment for non-DSA HLAs. 

library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(purrr)
setwd("/ix/djishnu/Trirupa/ABomics.Prj/")
#source("/ix/djishnu/Trirupa/ABomics.Prj/scripts/DSA_hyperGeom_Fx.R")
source("/ix/djishnu/Trirupa/ABomics.Prj/scripts/Ag_hyperg_Fx.R")

####### compiling all the sig LFs and the features in them ######
## for Early Ab-omics model###
E_Zlist=c("Z1","Z12") 
E_Zfeat_path= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/SLIDE/"
E_allZ_df= load_allZ(Z_list= E_Zlist, FeatList_path=E_Zfeat_path)
colnames(E_allZ_df)[colnames(E_allZ_df)== "Ag_pool"] <- "Ag_class"

#### for EARLY ####
early_df= create_agdf("./input_data/Pitt_EARLY_serotype.csv")

Ag_hyperg_NR= Ag_hyperg.fx(allZ_df= E_allZ_df, grp_agdf= early_df$NR)
probe_agsumdf_NR= unique(Ag_hyperg_NR %>% arrange(desc(numQueryAg_probe)))
probe_agsumdf_NR= probe_agsumdf_NR[probe_agsumdf_NR$pval <= 0.05,]
write.csv(probe_agsumdf_NR,"./Figures/Abomics_Early/Ag_EarlyNR.hypergeomV5_earlypatONLY.csv", row.names=F, quote=F)

Ag_hyperg_ABMR= Ag_hyperg.fx(allZ_df= E_allZ_df, grp_agdf= early_df$ABMR)
probe_agsumdf_ABMR= unique(Ag_hyperg_ABMR %>% arrange(desc(numQueryAg_probe)))
probe_agsumdf_ABMR= probe_agsumdf_ABMR[probe_agsumdf_ABMR$pval <= 0.05,]
write.csv(probe_agsumdf_ABMR,"./Figures/Abomics_Early/Ag_EarlyABMR.hypergeomV5_earlypatONLY.csv", row.names=F, quote=F)

Ag_hyperg_nonDSA= Ag_hyperg.fx(allZ_df= E_allZ_df, grp_agdf= early_df$nonDSA)
probe_agsumdf_nonDSA= unique(Ag_hyperg_nonDSA %>% arrange(desc(numQueryAg_probe)))
probe_agsumdf_nonDSA= probe_agsumdf_nonDSA[probe_agsumdf_nonDSA$pval <= 0.05,]

probe_agsumdf_nonDSA.final=probe_agsumdf_nonDSA[probe_agsumdf_nonDSA$numQueryAg_probe >=2,   ]
write.csv(probe_agsumdf_nonDSA.final,"./Figures/Abomics_Early/Ag_EarlyNonDSA.hypergeomV5_earlypatONLY.csv", row.names=F, quote=F)



####### compiling all the sig LFs and the features in them ######
## for Early Ab-omics model###
L_Zlist=c("Z1","Z2","Z7","Z10") 
L_Zfeat_path= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/SLIDE/"
L_allZ_df= load_allZ(Z_list= L_Zlist, FeatList_path= L_Zfeat_path)
colnames(L_allZ_df)[colnames(L_allZ_df)== "Ag_pool"] <- "Ag_class"

#### for LATE ####
late_df= create_agdf("./input_data/Pitt_LATE_serotype.csv")
Ag_hyperg_NR= Ag_hyperg.fx(allZ_df= L_allZ_df, grp_agdf= late_df$NR)
probe_agsumdf_NR= unique(Ag_hyperg_NR %>% arrange(desc(numQueryAg_probe)))
probe_agsumdf_NR= probe_agsumdf_NR[probe_agsumdf_NR$pval <= 0.05,]
probe_agsumdf_NR.final=probe_agsumdf_NR[probe_agsumdf_NR$numQueryAg_probe >2,   ]
write.csv(probe_agsumdf_NR.final,"./Figures/Abomics_Late/Ag_LateNR.hypergeomV5_latepatONLY.csv", row.names=F, quote=F)

Ag_hyperg_ABMR= Ag_hyperg.fx(allZ_df= L_allZ_df, grp_agdf= late_df$ABMR)
probe_agsumdf_ABMR= unique(Ag_hyperg_ABMR %>% arrange(desc(numQueryAg_probe)))
probe_agsumdf_ABMR= probe_agsumdf_ABMR[probe_agsumdf_ABMR$pval <= 0.05,]
probe_agsumdf_ABMR.final=probe_agsumdf_ABMR[probe_agsumdf_ABMR$numQueryAg_probe >2,   ]
write.csv(probe_agsumdf_ABMR.final,"./Figures/Abomics_Late/Ag_LateABMR.hypergeomV5_latepatONLY.csv", row.names=F, quote=F)

Ag_hyperg_nonDSA= Ag_hyperg.fx(allZ_df= L_allZ_df, grp_agdf= late_df$nonDSA)
probe_agsumdf_nonDSA= unique(Ag_hyperg_nonDSA %>% arrange(desc(numQueryAg_probe)))
probe_agsumdf_nonDSA= probe_agsumdf_nonDSA[probe_agsumdf_nonDSA$pval <= 0.05,]
probe_agsumdf_nonDSA.final=probe_agsumdf_nonDSA[probe_agsumdf_nonDSA$numQueryAg_probe >2,   ]
write.csv(probe_agsumdf_nonDSA.final,"./Figures/Abomics_Late/Ag_LateNonDSA.hypergeomV5_latepatONLY.csv", row.names=F, quote=F)
#----------------##### for non-DSA HLAs ######---------------------------#


EARLY_hyperg_nonDSA= Ag_hyperg.fx(allZ_df= E_allZ_df, grp_agdf= nonDSA_df)
probe_agsumdf_nonDSA= unique(EARLY_hyperg_nonDSA %>% arrange(desc(numQueryAg_probe)))
probe_agsumdf_nonDSA= probe_agsumdf_nonDSA[probe_agsumdf_nonDSA$pval <= 0.05,]


library(ggplot2)
library(dplyr)


df = probe_agsumdf_ABMR
# Define color categories
df <- df %>%
  mutate(Category = case_when(
    pval <= 0.02 ~ "Rejector",
    pval > 0.02 & pval <= 0.05 ~ "Non-Rejector",
    TRUE ~ "Absent"
  ))

# Define colors
fill_colors <- c("Rejector" = "#C95A76", "Non-Rejector" = "#E7C46C", "Absent" = "white")

# Create the plot
ggplot(df, aes(x = QueryAg, y = Probe, fill = Category)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = fill_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Patient donor-specific antigens",
    y = "",
    fill = "Enrichment"
  )


### Two things to be done:
# 1. ceheck why B51 is coming up for three groups in early. Drop the ones that show up in all three, and keep only specific ones
# 2. Put a threshold on teh numQuery in LF (i.e the numerator) >1. 