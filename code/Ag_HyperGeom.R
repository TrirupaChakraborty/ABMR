library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
source("/ix/djishnu/Trirupa/ABomics.Prj/scripts/DSA_hyperGeom_Fx.R")

setwd("/ix/djishnu/Trirupa/ABomics.Prj")

sero.grpdf= read.csv("./input_data/Pitt_serotype.status.csv")
DSA_agdf <- separate_rows(sero.grpdf, DSA, sep = ",")
DSA_agdf$DSA <- str_trim(DSA_agdf$DSA)
ABMR_agdf= DSA_agdf[DSA_agdf$Group == "ABMR",]
NR_agdf= DSA_agdf[DSA_agdf$Group == "NR",]

##### Loading the Antigen set that was used in the assay #######
BeadAg_path <- "./input_data/onelambda_beadantigen_regionsEDITED.xlsx"
c1c2.Ag_df= loadc1c2df(BeadAg_path)

####### compiling all the sig LFs and the features in them ######
## for Early Ab-omics model###
E_Zlist=c("Z1","Z12") 
E_Zfeat_path= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/SLIDE/"
E_allZ_df= load_allZ(Z_list= E_Zlist, FeatList_path=E_Zfeat_path)
colnames(E_allZ_df)[colnames(E_allZ_df)== "Ag_pool"] <- "Ag_class"

####### compiling all the sig LFs and the features in them ######
## for Late Ab-omics model###
L_Zlist=c("Z1","Z2","Z7","Z10") 
L_Zfeat_path= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/SLIDE/"
L_allZ_df= load_allZ(Z_list= L_Zlist, FeatList_path= L_Zfeat_path)
colnames(L_allZ_df)[colnames(L_allZ_df)== "Ag_pool"] <- "Ag_class"


##### performing hypergeometric test ###################################

Ag_hyperg.fx= function(allZ_df, grp_agdf, c1c2.Ag_df){
  ## parameters ##
  # no. of Ag1 = m
  # total no. of Ags tested= x
  # no. of Ags that are not Ag1 = n
  #
  # no. of Ags against a given probe in  LF = k
  # no. of Ag1s against a given probe in LF = q
  
  #NOTE!!! for total number of Ags tested we should not do a unique 
  probe_agsumdf= data.frame(Probe=character(), QueryAg=character(), QueryAg_inpool= numeric(),
                            numQueryAg_probe= numeric(), numAg_probe= numeric(), Total_Ags_inpool= numeric(),
                            num_minip= numeric(), pval= numeric())
  
  #allZ_df= E_allZ_df
  allpatAGs=unique(grp_agdf$DSA)
  tot_allpatAgs= length(allpatAGs)
  tot_Agstested = length(unlist(lapply(c1c2.Ag_df[-1], function(x) x[x != "."]))) ## x
  probe_list= unique(allZ_df$Probe)
  
  for (queryAg in allpatAGs){
    
    numof_queryAg= sum(c1c2.Ag_df == queryAg) ## m
    
    num_nonqueryAg= tot_Agstested - numof_queryAg ## n
    
    for (probe in probe_list){
      
      probe_df= allZ_df[allZ_df$Probe == probe,c("Probe", "Ag_class")]
      
      class_listinLF= unique(probe_df$Ag_class)
      
      Ags_inclassdf= c1c2.Ag_df[c1c2.Ag_df$Ag_class %in% class_listinLF, ]
      Ags_againstPrb = length(unlist(lapply(Ags_inclassdf[-1], function(x) x[x != "."]))) ## k
      queryAg_agianstPrb = sum(Ags_inclassdf == queryAg, na.rm = TRUE )  ## q
      
      if(queryAg_agianstPrb != 0){
        ########### running hypergeometric test #########
        # phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
        p_value <- phyper(queryAg_agianstPrb, numof_queryAg, num_nonqueryAg, Ags_againstPrb, lower.tail = FALSE)
        
        probe_agsumdf= rbind(probe_agsumdf, data.frame(Probe= probe, QueryAg= queryAg,
                                                       QueryAg_inpool= numof_queryAg,
                                                       numQueryAg_probe= queryAg_agianstPrb,
                                                       numAg_probe= Ags_againstPrb,
                                                       Total_Ags_inpool= tot_Agstested,
                                                       num_minip= length(class_listinLF),
                                                       pval= round(p_value,2)))
      }
      
    }
  } 
  
  return(probe_agsumdf)
}

  
#### for EARLY ####
Ag_hyperg_NR= Ag_hyperg.fx(allZ_df= E_allZ_df, grp_agdf= NR_agdf, c1c2.Ag_df= c1c2.Ag_df)
probe_agsumdf_NR= unique(Ag_hyperg_NR %>% arrange(desc(numQueryAg_probe)))
probe_agsumdf_NR= probe_agsumdf_NR[probe_agsumdf_NR$pval <= 0.05,]

write.csv(probe_agsumdf_NR,"./Figures/Abomics_Early/Ag_EarlyNR.hypergeom.csv", sep="\t", row.names=F, quote=F)

#### for LATE ####
Ag_hyperg_ABMR= Ag_hyperg.fx(allZ_df= L_allZ_df, grp_agdf= ABMR_agdf, c1c2.Ag_df= c1c2.Ag_df)
probe_agsumdf_ABMR= unique(Ag_hyperg_ABMR %>% arrange(desc(numQueryAg_probe)))
probe_agsumdf_ABMR= probe_agsumdf_ABMR[probe_agsumdf_ABMR$pval <= 0.05,]

write.csv(probe_agsumdf_ABMR,"./Figures/Abomics_Early/Ag_EarlyABMR.hypergeom.csv", sep="\t", row.names=F, quote=F)

