library(readxl)
library(ggplot2)
library(dplyr)

metadf=read.csv("PittCohort_kidneytransplant_ABMR_patientID_clinicalData.csv")
subset_metadf=metadf[,c(1,2,5)]
subset_metadf$CUSTOMER_ID= gsub(".*-0*([1-9]\\d+)-.*", "\\1", subset_metadf$CUSTOMER_ID)
Pitt_patid_status <- read_excel("Pitt_patid_status_serotype.xlsx")

## first filtering out those patients who weren't included in the final study
subset_metadf= subset_metadf[subset_metadf$CUSTOMER_ID %in% Pitt_patid_status$ID,]

median_t= median(subset_metadf$Time.to.Tx..mo.)
subset_metadf$Time_to_Tx_weeks= subset_metadf$Time.to.Tx..mo.*4.345 #(converting months to week)
median_t_weeks= median(subset_metadf$Time_to_Tx_weeks)


#ABMR_df= subset_metadf[subset_metadf$GRP == "ABMR", ]
#DSA_df= subset_metadf[subset_metadf$GRP == "DSA", ]
#median_t_abmr= median(ABMR_df$Time.to.Tx..mo.)
#median_t_dsa= median(DSA_df$Time.to.Tx..mo.)

latedf= subset_metadf[subset_metadf$Time.to.Tx..mo. >= 1.65, ]
latedf= subset_metadf[subset_metadf$Time_to_Tx_weeks >= 5.214, ]
earlydf= subset_metadf[subset_metadf$Time_to_Tx_weeks <= 5.214, ]

################################################################
c4meta= read_xlsx("/ix/djishnu/Trirupa/ABomics.Prj/input_data/C4dscores_newPitt_patients2020-2024.xlsx", sheet=1)
c4meta=c4meta[,c("STUDY ID","TX_DATE","BIOPSY_DATE","GROUP","C4D")]
c4meta$TX_DATE = as.Date(c4meta$TX_DATE)
c4meta$BIOPSY_DATE = as.Date(c4meta$BIOPSY_DATE)
c4meta$Time_to_Tx.weeks= as.numeric(difftime(c4meta$BIOPSY_DATE, c4meta$TX_DATE, units = "weeks")) ## 1 month has avg of 30.44 days

c4_late= c4meta[c4meta$Time_to_Tx.weeks >= 6, ]
c4_late= na.omit(c4_late)
## counting how many rejectors versus nonrejectors in the late category
grp.counts <- c4_late %>% group_by(GROUP) %>% summarize(Count=n())

# Count elements in each group where C4d <= 2
# Count how many in each category have C4d ≤ 2
c4d.counts <- c4_late %>%
  group_by(GROUP) %>%
  summarize(count_high_C4d = sum(C4D >= 2))

c4_late$highlight = ifelse(c4_late$C4D >= 2, ">=2", "<2")

# Perform the Wilcoxon test (Mann-Whitney U test)
wilcox_test <- wilcox.test(C4D ~ GROUP, data = c4_late)

# Extract the p-value from the test result
p_value <- wilcox_test$p.value

# Create the violin plot with ggplot
c4d_plot= ggplot(c4_late, aes(x = GROUP, y = C4D, fill = GROUP)) +
  geom_violin(alpha = 0.55) +  # Violin plot with transparency and no border
  geom_jitter(aes(color = highlight), width = 0.2, height=0.03, size = 3, alpha = 0.55) +  # Dots
  scale_color_manual(values = c(">=2" = "#DC381F", "<2" = "#ADADC9")) +  # Dot colors
  scale_fill_manual(values = c("ABMR+DSA+" = "#ca4f73", "ABMR-DSA+" = "#eeca59")) +  # Violin colors
  labs(title = "C4d Distribution by Group", x = "GROUP", y = "C4d Score", color = "C4d Value") +
  annotate("text", x = 1.5, y = max(c4_late$C4D) + 0.5, 
           label = paste("p =", format(p_value, digits = 2)), 
           size = 3, color = "black") +  # Add the p-value annotation
  theme_minimal() +
  theme(legend.position = "bottom", 
        axis.line = element_line(color = "black"))  # Make the axis lines black

ggsave("/ix/djishnu/Trirupa/ABomics.Prj/Figures/Abomics_Late/c4d_plot_median6wks.pdf",c4d_plot ,width=6, height=5 )



