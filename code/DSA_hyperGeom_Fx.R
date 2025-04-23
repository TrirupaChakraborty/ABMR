library(dplyr)
library(tidyr)
library(stringr)

#### this fx combines class I and class II Ags that were tested in the assay as a combined df ######
loadc1c2df <- function(BeadAg_path){
  c1.Ag_df <- read_excel(BeadAg_path, sheet = "ClassI", col_names=F)
  c1.Ag_df= c1.Ag_df[-1,]
  c1_colnames <- c("Ag_class","A", "A", "B", "B", "B", "B", "C", "C")
  colnames(c1.Ag_df) = c1_colnames
  for (i in 2:length(c1.Ag_df)){
    class_name=colnames(c1.Ag_df[i])
    c1.Ag_df[i] = lapply(c1.Ag_df[i], function(x) ifelse(x !=".",paste0(class_name,x), x))
  }
  
  c2.Ag_df <- read_excel(BeadAg_path, sheet = "ClassII", col_names=F)
  c2.Ag_df= c2.Ag_df[-1,]
  c2_colnames <- c("Ag_class","DR", "DR", "DR", "DR", "DQ", "DQ", "DP", "DP")
  colnames(c2.Ag_df) = c2_colnames
  for (i in 2:length(c2.Ag_df)){
    class_name=colnames(c2.Ag_df[i])
    c2.Ag_df[i] = lapply(c2.Ag_df[i], function(x) ifelse(x !=".",paste0(class_name,x), x))
  }
  
  #### merging infor from both c1 and c2 antigens ####
  colnames(c2.Ag_df) <- colnames(c1.Ag_df) ## just to combine the two dfs making the column names same
  c1c2.Ag_df = rbind(c1.Ag_df,c2.Ag_df)
  c1c2.Ag_df$Ag_class= gsub("-",".",c1c2.Ag_df$Ag_class)
  c1c2.Ag_df$Ag_class= gsub("\\s+","",c1c2.Ag_df$Ag_class)
  
  return(c1c2.Ag_df)
}

###### this fx loads all the sig LFs as a combined df #####
load_allZ <- function(Z_list, FeatList_path, A_thresh){
  allZ_df= data.frame()
  for (Znum in Z_list){
    Z_df=read.table(paste0(FeatList_path,"feature_list_",Znum,".txt"), header=T)
    colnames(Z_df)=c("Profile","A_loading", "AUCs","corrs","color")
    # Split the "names" column into Ab and Ag 
    Z_newdf <- separate(Z_df, "Profile", into = c("Probe", "Ag_pool"), sep = "_")
    allZ_df=rbind(allZ_df,Z_newdf)
  }
  
  allZ_df=na.omit(allZ_df)
  filt_allZ_df= allZ_df[allZ_df$A_loading >= A_thresh | allZ_df$AUCs >= 0.6 | allZ_df$AUCs <= 0.4 , ] # previously we took A_thresh=0.2
  #### dropping any characters after the first instance of . in minipools eg. ClassI.3.15.
  clean_mpool_name <- function(string){
    modified_str <- sub("^([^\\.]+\\.[^\\.]+)\\..*", "\\1", string)
    return(modified_str)
  }
  
  filt_allZ_df$Ag_pool=sapply(filt_allZ_df$Ag_pool,clean_mpool_name)
  
  return(filt_allZ_df)
}

##### performing hypergeometric test ###################################

#### this fx runs hypergeometric to check enrichment of patient specific DSAs in the LFs by probe.
## parameters ##
# no. of DSA Ags = m
# total no.of Ags tested = x
# no. of non-DSA Ags = n
# no. of Ags in LF = k
# no. of DSA Ags OF the no. of Ags in LF = q

run_hypergeom1 <- function(Ag_df, HLAclass_df,filt_allZ_df){
  total_DSAag = length(unique(Ag_df$DSA)) ### m 
  
  ##### total no.of Ags tested for = x ########
  total_Ags= length(unique(unlist(lapply(HLAclass_df[-1], function(x) x[x != "."])))) ### X
  
  ##### no. of non-DSA Ags = n ##########
  nonDSA_Ag= total_Ags - total_DSAag ### n 
  
  ###### computing values "k" and "q" for each probe ######
  ## for a given probe ##
  probe_list=unique(filt_allZ_df$Probe)
  probe_DSAenrich <- data.frame(probe = character(), p_value = numeric())
  param_df=data.frame(Probe=character(),Total_DSA.Ags=numeric(),NonDSA_Ags=numeric(),
                      c1c2Ags_inLF=numeric(),numDSAags_inLF=numeric(), p_value = numeric())
  
  class_list= unique(HLAclass_df$Ag_class)
  
  for (probe in probe_list) {
    subset_probe <- filt_allZ_df[filt_allZ_df$Probe == probe,] 
    LFprobe_list <- unique(subset_probe$Ag_pool[subset_probe$Ag_pool %in% class_list])  # Filter LFprobe_list
    
    Ags_inLF= c()
    for (j in LFprobe_list) {
      class_df <- HLAclass_df[HLAclass_df$Ag_class == j,]
      uniq_ags <- unique(unlist(lapply(class_df[-1], function(x) x[x != "."]))) ### class_df needs to be defined
      Ags_inLF <- c(Ags_inLF, uniq_ags)
    }
    
    Ags_inLF <- Ags_inLF[!is.na(Ags_inLF)]
    numAgs_inLF <- length(unique(Ags_inLF))  # k
    numDSA_inLF <- length(intersect(Ags_inLF, unique(Ag_df$DSA)))  # q 
    
    # calculating p-value from hypergeometric test to check enrichment of DSAs in all LFs combined
    p_value <- phyper(numDSA_inLF, total_DSAag, nonDSA_Ag, numAgs_inLF, lower.tail = FALSE)
    
    # Append the results to the dataframe
    probe_DSAenrich <- rbind(probe_DSAenrich, data.frame(probe = probe, p_value = p_value))
    param_df <- rbind(param_df, data.frame(Probe = probe, Total_DSA.Ags = total_DSAag,
                                           NonDSA_Ags = nonDSA_Ag, c1c2Ags_inLF = numAgs_inLF,
                                           numDSAags_inLF = numDSA_inLF, 
                                           p_value = round(p_value,digits=3)))
  }
  
  return(param_df)
}

