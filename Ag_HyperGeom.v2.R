library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(purrr)
setwd("/ix/djishnu/Trirupa/ABomics.Prj/")
source("/ix/djishnu/Trirupa/ABomics.Prj/scripts/DSA_hyperGeom_Fx.R")


ag_mappingMT= readRDS("./MARISA/ag_mapping_matrix.rds")
ag_mappingMT <- ag_mappingMT %>%
  select(-matches("Mica"))
rows_to_keep <- !grepl("MICA", rownames(ag_mappingMT))
ag_mappingMT= ag_mappingMT[rows_to_keep,]
ag_mappingMT$Ag_class= rownames(ag_mappingMT) 

#unique_elements <- lapply(ag_mappingMT, unique) ## the dataframe has only 0,1 

#####creating a df that has teh total counts of each Ag in the assay
#Ag_counts_df= data.frame(unlist(lapply(ag_mappingMT, sum)))
#Ag_counts_df$Ag= rownames(Ag_counts_df)
#colnames(Ag_counts_df)[1] = "Ag_count"
#Ag_counts_df=Ag_counts_df %>% select(Ag,Ag_count)



sero.grpdf= read.csv("./input_data/Pitt_serotype.status.csv")
DSA_agdf <- separate_rows(sero.grpdf, DSA, sep = ",")
DSA_agdf$DSA <- str_trim(DSA_agdf$DSA)
ABMR_agdf= DSA_agdf[DSA_agdf$Group == "ABMR",]
NR_agdf= DSA_agdf[DSA_agdf$Group == "NR",]
#NR_uniqAg= setdiff(NR_agdf$DSA, ABMR_agdf$DSA)
#ABMR_uniqAg= setdiff(ABMR_agdf$DSA,NR_agdf$DSA)


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

Ag_hyperg.fx= function(allZ_df, grp_agdf){
  ## parameters ##
  # no. of Ag1 = m
  # total no. of Ags tested= x
  # no. of Ags that are not Ag1 = n
  #
  # no. of Ags against a given probe in  LF = k
  # no. of Ag1s against a given probe in LF = q
  probe_agsumdf= data.frame(Probe=character(), QueryAg=character(), QueryAg_inpool= numeric(),
                          numQueryAg_probe= numeric(), numAg_probe= numeric(), Total_Ags_inpool= numeric(),
                          num_mpools= numeric(), pval= numeric())
  allpatAGs=unique(grp_agdf$DSA)
  tot_allpatAgs= length(allpatAGs)
  
  ##### MAKING A DF THAT HAS THE TOTAL COUNTS OF ALL THE AGS TESTED IN THE ASSAY
  Agscounts <- t((ag_mappingMT %>% summarise_all(~ sum(. == 1))))
  Agscounts = data.frame(Agscounts)
  Agscounts$Ag= rownames(Agscounts)
  Agscounts=Agscounts[,c("Ag","Agscounts")]
  colnames(Agscounts)= c("Ag","counts")
  tot_Agstested = sum(Agscounts$counts) ###### x
  
  probe_list= unique(allZ_df$Probe)
  
  for (queryAg in allpatAGs){
    if(queryAg %in% colnames(ag_mappingMT)){
      #print(queryAg)
      numof_queryAg= Agscounts[Agscounts$Ag == queryAg,]$counts ## m
      
      num_nonqueryAg= tot_Agstested - numof_queryAg ## n
      
      for (probe in probe_list){
        probe_df= allZ_df[allZ_df$Probe == probe,c("Probe", "Ag_class")]
        class_listinLF= unique(probe_df$Ag_class)
        inLFclass_ag_map= ag_mappingMT[ag_mappingMT$Ag_class %in% class_listinLF ,]
        inLFclass_ag_map <- inLFclass_ag_map %>% select(-Ag_class)
        Ags_inclassdf= t((inLFclass_ag_map %>% 
                            summarise_all(~ sum(. == 1)))) 
        Ags_inclassdf= data.frame(Ags_inclassdf)
        Ags_inclassdf$Ag= rownames(Ags_inclassdf) 
        colnames(Ags_inclassdf)= c("counts","Ag")
        Ags_againstPrb = sum(Ags_inclassdf$counts) ##k
        
        queryAg_agianstPrb = Ags_inclassdf[Ags_inclassdf$Ag == queryAg, ]$counts ##q
        queryAg_agianstPrb=as.numeric(queryAg_agianstPrb)
        print(queryAg_agianstPrb)
        if(queryAg_agianstPrb != 0 |length(queryAg_agianstPrb) <1){
          ########### running hypergeometric test #########
          # phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
          p_value <- phyper(queryAg_agianstPrb, numof_queryAg, num_nonqueryAg, Ags_againstPrb, lower.tail = FALSE)
          
          probe_agsumdf= rbind(probe_agsumdf, data.frame(Probe= probe, QueryAg= queryAg,
                                                         QueryAg_inpool= numof_queryAg,
                                                         numQueryAg_probe= queryAg_agianstPrb,
                                                         numAg_probe= Ags_againstPrb,
                                                         Total_Ags_inpool= tot_Agstested,
                                                         num_mpools=length(class_listinLF),
                                                         pval= round(p_value,2)))
        }
        }
      }
    }
    return(probe_agsumdf)
  }
  
    
#### for EARLY ####
Ag_hyperg_NR= Ag_hyperg.fx(allZ_df= E_allZ_df, grp_agdf= NR_agdf)
probe_agsumdf_NR= unique(Ag_hyperg_NR %>% arrange(desc(numQueryAg_probe)))
probe_agsumdf_NR= probe_agsumdf_NR[probe_agsumdf_NR$pval <= 0.05,]
write.csv(probe_agsumdf_NR,"./Figures/Abomics_Early/Ag_EarlyNR.hypergeomV2.csv", row.names=F, quote=F)


Ag_hyperg_ABMR= Ag_hyperg.fx(allZ_df= E_allZ_df, grp_agdf= ABMR_agdf)
probe_agsumdf_ABMR= unique(Ag_hyperg_ABMR %>% arrange(desc(numQueryAg_probe)))
probe_agsumdf_ABMR= probe_agsumdf_ABMR[probe_agsumdf_ABMR$pval <= 0.05,]
write.csv(probe_agsumdf_ABMR,"./Figures/Abomics_Early/Ag_EarlyABMR.hypergeomV2.csv", row.names=F, quote=F)


#### for LATE ####
L.Ag_hyperg_NR= Ag_hyperg.fx(allZ_df= L_allZ_df, grp_agdf= NR_agdf, ag_mappingMT= ag_mappingMT)
L.probe_agsumdf_NR= unique(L.Ag_hyperg_NR %>% arrange(desc(numQueryAg_probe)))
L.probe_agsumdf_NR= L.probe_agsumdf_NR[L.probe_agsumdf_NR$pval <= 0.05,]
write.csv(L.probe_agsumdf_NR,"./Figures/Abomics_Late/Ag_LateNR.hypergeomV2.csv", row.names=F, quote=F)


L.Ag_hyperg_ABMR= Ag_hyperg.fx(allZ_df= L_allZ_df, grp_agdf= ABMR_agdf, ag_mappingMT= ag_mappingMT)
L.probe_agsumdf_ABMR= unique(L.Ag_hyperg_ABMR %>% arrange(desc(numQueryAg_probe)))
L.probe_agsumdf_ABMR= L.probe_agsumdf_ABMR[L.probe_agsumdf_ABMR$pval <= 0.05,]
write.csv(L.probe_agsumdf_ABMR,"./Figures/Abomics_Late/Ag_LateABMR.hypergeomV2.csv", row.names=F, quote=F)
