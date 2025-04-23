library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(purrr)
library(readr)
setwd("/ix/djishnu/Trirupa/ABomics.Prj/")
source("/ix/djishnu/Trirupa/ABomics.Prj/scripts/DSA_hyperGeom_Fx.R")

ag_mappingMT= readRDS("./MARISA/ag_mapping_matrix.rds")
ag_mappingMT <- ag_mappingMT %>%
  select(-matches("Mica"))
rows_to_keep <- !grepl("MICA", rownames(ag_mappingMT))
ag_mappingMT= ag_mappingMT[rows_to_keep,]
ag_mappingMT$Ag_class= rownames(ag_mappingMT) 
#write.csv(ag_mappingMT, "./input_data/ag_mapping_mt.TC.csv", quote=F)

create_agdf= function(sero_path){
  ag_mappingMT=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/input_data/ag_mapping_mt.TC.csv", row.names=1)
  sero.grpdf= read_csv(sero_path)
  colnames(sero.grpdf)[colnames(sero.grpdf)=="serotype"] = "DSA"
  DSA_agdf <- separate_rows(sero.grpdf, DSA, sep = ",")
  DSA_agdf$DSA <- str_trim(DSA_agdf$DSA)
  DSA_agdf=DSA_agdf[,-1]
  DSA_agdf=unique(DSA_agdf)
  ABMR_agdf= DSA_agdf[DSA_agdf$Group == "ABMR",]
  NR_agdf= DSA_agdf[DSA_agdf$Group == "NR",]
  
  ## if an antigen occurs in both NR and ABMR then keep its instance only in ABMR
  common_ag= intersect(ABMR_agdf$DSA,NR_agdf$DSA)
  NR_agdf=NR_agdf[-(NR_agdf$DSA %in% common_ag),]
  
  nonDSA_df= data.frame(colnames(ag_mappingMT)[(!(colnames(ag_mappingMT) %in% (DSA_agdf$DSA)))])
  colnames(nonDSA_df)="DSA"
  nonDSA_df$grp= "nonDSA"
  
  ### keeping only the Ags unique to each status (i.e ABMR, NR and nonDSA)
  unique_ABMR_dsa <- setdiff(ABMR_agdf$DSA, union(NR_agdf$DSA, nonDSA_df$DSA))
  ABMR_agdf=ABMR_agdf[ABMR_agdf$DSA %in% unique_ABMR_dsa,]
  
  unique_NR_dsa <- setdiff(NR_agdf$DSA, union(ABMR_agdf$DSA, nonDSA_df$DSA))
  NR_agdf <- NR_agdf[NR_agdf$DSA %in% unique_NR_dsa,]
  
  unique_nonDSA_dsa <- setdiff(nonDSA_df$DSA, union(ABMR_agdf$DSA, NR_agdf$DSA))
  nonDSA_df <- nonDSA_df[nonDSA_df$DSA %in% unique_nonDSA_dsa,]
  
  return(list(ABMR = ABMR_agdf, NR = NR_agdf, nonDSA = nonDSA_df))
  
}

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
        if (length(queryAg_agianstPrb) > 0 && !is.na(queryAg_agianstPrb) && queryAg_agianstPrb > 0){
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

### change ag_mappingMT to long version and then re-write this code.
get.mpool_ag <- function(grp_df, allZ_df) {
  mpool_ag_list <- list() 
  for (ag in unique(grp_df$QueryAg)){
    ab_list <- grp_df$Probe[grp_df$QueryAg == ag]
    print(paste0("Antigen: ", ag))
    print(paste0("Probe: ",paste(ab_list, collapse=", ")))
    
    for (ab in ab_list){
      mpool_lis <- unique(allZ_df$Ag_class[allZ_df$Probe %in% ab])
      print(paste0("mpool_list= ",paste(mpool_lis, collapse= ", ")))
      # Checking which minipools contain the antigen
      filtered_rows.1 <- ag_mappingMT[ag_mappingMT$Ag_class %in% mpool_lis,]
      print(paste("filt_rows.1: ", paste(rownames(filtered_rows.1), collapse= ", ")))
      filtered_rows.2 = filtered_rows.1[filtered_rows.1[[ag]] == 1,]
      mpool_with_ag <- rownames(filtered_rows.2)
      print(paste("mpool_with_ag: ",paste(mpool_with_ag, collapse= ", ")))
      
      if (length(mpool_with_ag) > 0) {
        new_rows <- data.frame(
          ab_mpool = paste0(ab, "_", mpool_with_ag),
          QueryAg = ag,
          stringsAsFactors = FALSE
        )
        mpool_ag_list[[length(mpool_ag_list) + 1]] <- new_rows
      }
    }
    mpool_ag.df <- do.call(rbind, mpool_ag_list)      
  }
  return(mpool_ag.df)
}

## fx for plotting enirchment
plot_enrich= function(enrich_df, save_path){
  full_data <- enrich_df %>% # Ensure all Probe-QueryAg pairs exist, even if missing
    complete(ab_mpool, QueryAg, fill = list(status = NA))  # Adds missing combinations
  
  enrich.plt= ggplot(full_data, aes(x = QueryAg, y = ab_mpool, fill = status)) +
    geom_tile(color = "grey", size = 0.5) +  # Black grid lines everywhere
    scale_fill_manual(values = c("ABMR" = "#ca4f73","NR" = "#eeca59", "nonDSA" = "#6b98c0"), na.value = "white") +  # White for missing values
    theme_minimal(base_size = 14) +  
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  
      panel.grid.major = element_blank(),  # Remove ggplot default grid
      panel.grid.minor = element_blank()  
    ) +
    labs( x = "Query Antigens", y = "Probes")

  ggsave(save_path,enrich.plt)
  return(enrich.plt)
}