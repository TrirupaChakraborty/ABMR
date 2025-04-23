library(ggplot2)
library(tidyr)
library(dplyr)

setwd("/ix/djishnu/Trirupa/ABomics.Prj/")
source("/ix/djishnu/Trirupa/ABomics.Prj/scripts/Ag_hyperg_Fx.R")

### FIRST, GETTING all_Z_df (get Ab associated with QueryAg) #######################
#L_Zlist=c("Z1","Z2","Z7","Z10") 
L_Zlist= "Z1"
L_Zfeat_path= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/SLIDE/"
L_allZ_df= load_allZ(Z_list= L_Zlist, FeatList_path= L_Zfeat_path)
colnames(L_allZ_df)[colnames(L_allZ_df)== "Ag_pool"] <- "Ag_class"
L_allZ_df$probe_mpool= paste0(L_allZ_df$Probe,"_", L_allZ_df$Ag_class)

######### SECOND, GETTING PROBE-AG ENRICHMENT MAP ################# 
select_cols= c("Probe","QueryAg")

L.NR_df = read.csv("./Figures/Abomics_Late/Ag_LateNR.hypergeomV4_latepatONLY.csv")
filt.L.NR_df= L.NR_df[,colnames(L.NR_df) %in% select_cols]
filt.L.NR_df$status = "NR"

L.ABMR_df = read.csv("./Figures/Abomics_Late/Ag_LateABMR.hypergeomV4_latepatONLY.csv")
filt.L.ABMR_df= L.ABMR_df[,colnames(L.ABMR_df) %in% select_cols]
filt.L.ABMR_df$status = "ABMR"

L.nonDSA_df = read.csv("./Figures/Abomics_Late/Ag_LateNonDSA.hypergeomV4_latepatONLY.csv")
filt.L.nonDSA_df= L.nonDSA_df[,colnames(L.nonDSA_df) %in% select_cols]
filt.L.nonDSA_df$status = "nonDSA"

L_enrich.df= rbind(filt.L.NR_df,filt.L.ABMR_df,filt.L.nonDSA_df)

L.mpool_ag.df= get.mpool_ag(L_enrich.df, L_allZ_df) #### resolve issue with get.mpool_ag. check why IgG3_Class1.3 showing up- 11th March'25
L.mpool_ag.df= L.mpool_ag.df %>% separate(ab_mpool, into = c("Probe", "mpool"), sep = "_", remove =FALSE)
L.mpool_ag.df= L.mpool_ag.df[L_allZ_df$probe_mpool %in% L.mpool_ag.df$ab_mpool,]
L_enrich.df= L_enrich.df %>% left_join(L.mpool_ag.df, by = c("Probe", "QueryAg"))
L_enrich.df <- L_enrich.df %>% select(-mpool)

L_enrich.df_C1= L_enrich.df[grepl("^[ABC]", L_enrich.df$QueryAg), ]
L_enrich.df_C1$QueryAg <- factor(L_enrich.df_C1$QueryAg, levels = unique(L_enrich.df_C1$QueryAg[order(L_enrich.df_C1$status)]))

L_enrich.df_C2= L_enrich.df[grepl("^D", L_enrich.df$QueryAg), ]
L_enrich.df_C2$QueryAg <- factor(L_enrich.df_C2$QueryAg, levels = unique(L_enrich.df_C2$QueryAg[order(L_enrich.df_C2$status)]))

L.C1_path= "./Figures/Abomics_Late/L.Ag_enrich.ClassI.v2.pdf"
L.C1_plt= plot_enrich(L_enrich.df_C1, L.C1_path)

L.C2_path= "./Figures/Abomics_Late/L.Ag_enrich.ClassII.pdf"
L.C2_plt= plot_enrich(L_enrich.df_C2, L.C2_path)


#############################################
#################### FOR EARLY ##############
#############################################

E.NR_df = read.csv("./Figures/Abomics_Early/Ag_EarlyNR.hypergeomV4_earlypatONLY.csv")
filt.E.NR_df= E.NR_df[,colnames(E.NR_df) %in% select_cols]
filt.E.NR_df$status = "NR"

E.ABMR_df = read.csv("./Figures/Abomics_Early/Ag_EarlyABMR.hypergeomV4_earlypatONLY.csv")
filt.E.ABMR_df= E.ABMR_df[,colnames(E.ABMR_df) %in% select_cols]
filt.E.ABMR_df$status = "ABMR"

E.nonDSA_df = read.csv("./Figures/Abomics_Early/Ag_EarlyNonDSA.hypergeomV4_earlypatONLY.csv")
filt.E.nonDSA_df= E.nonDSA_df[,colnames(E.nonDSA_df) %in% select_cols]
filt.E.nonDSA_df$status = "nonDSA"

E_enrich.df= rbind(filt.E.NR_df,filt.E.ABMR_df,filt.E.nonDSA_df)


E_Zlist=c("Z1","Z12") 
E_Zfeat_path= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/SLIDE/"
E_allZ_df= load_allZ(Z_list= E_Zlist, FeatList_path=E_Zfeat_path)
colnames(E_allZ_df)[colnames(E_allZ_df)== "Ag_pool"] <- "Ag_class"


E.mpool_ag.df= get.mpool_ag(E_enrich.df, E_allZ_df)
E.mpool_ag.df= E.mpool_ag.df %>% separate(ab_mpool, into = c("Probe", "mpool"), sep = "_", remove =FALSE)
E_enrich.df= E_enrich.df %>% left_join(E.mpool_ag.df, by = c("Probe", "QueryAg"))
E_enrich.df <- E_enrich.df %>% select(-mpool)
## grouping by status ###
E_enrich.df$QueryAg <- factor(E_enrich.df$QueryAg, levels = unique(E_enrich.df$QueryAg[order(E_enrich.df$status)]))
E.plt_path= "./Figures/Abomics_Early/E.Ag_enrich.pdf"
E.plt= plot_enrich(E_enrich.df, E.plt_path)






#################

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)  # For reordering factors

# Define the custom antigen order
antigen_order <- c("A", "B", "C", "DP", "DQ", "DR")

# Function to extract the first letter for ordering
get_antigen_group <- function(ag) {
  ag_letter <- substr(ag, 1, 1)  # Get the first character (A, B, C, DP, DQ, DR)
  if (ag_letter %in% c("A", "B", "C")) {
    return(ag_letter)
  } else if (grepl("^DP", ag)) {
    return("DP")
  } else if (grepl("^DQ", ag)) {
    return("DQ")
  } else if (grepl("^DR", ag)) {
    return("DR")
  }
  return(ag)  # Default fallback
}

# Add antigen groups for sorting
L_enrich.df <- L_enrich.df %>%
  mutate(
    antigen_group = sapply(QueryAg, get_antigen_group),  # Assign groups (A, B, C, DP, DQ, DR)
    antigen_group = factor(antigen_group, levels = antigen_order),  # Convert to factor with order
    QueryAg = fct_reorder(QueryAg, as.numeric(antigen_group))  # Reorder within groups
  )

# Ensure all Probe-QueryAg pairs exist
full_data <- L_enrich.df %>%
  complete(Probe, QueryAg, fill = list(status = NA))  # Adds missing combinations

# Plot heatmap with Probe-based grouping for ab_mpool
late_enrich.plt <- ggplot(full_data, aes(x = QueryAg, y = ab_mpool, fill = status)) +
  geom_tile(color = "grey", size = 0.5) +  # Black grid lines
  scale_fill_manual(values = c("ABMR" = "#ca4f73", "NR" = "#eeca59", "nonDSA" = "#6b98c0"), na.value = "white") +  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Probe-QueryAg Heatmap", x = "Query Antigens", y = "ab_mpool") +
  facet_grid(rows = vars(Probe), scales = "free_y", space = "free_y")  # Group by Probe

# Show plot
late_enrich.plt

