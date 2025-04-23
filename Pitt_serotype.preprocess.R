library(readxl)
library(dplyr)
### loading the serotypes for full pitt cohort ###
full_cohort= read_excel("/ix/djishnu/Trirupa/ABomics.Prj/input_data/Pitt_patid_status_serotype.xlsx", sheet=1)

## separating the early and the late patients ###
x_df.E=read.csv("/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/x.csv")### early ##
x_df.L=read.csv("/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/x.csv")### early ##
### 

early_sero= full_cohort[full_cohort$ID %in% x_df.E$X,]
write.csv(early_sero, "/ix/djishnu/Trirupa/ABomics.Prj/input_data/Pitt_EARLY_serotype.csv", quote= F, row.names= F)

late_sero= full_cohort[full_cohort$ID %in% x_df.L$X,]
write.csv(late_sero, "/ix/djishnu/Trirupa/ABomics.Prj/input_data/Pitt_LATE_serotype.csv", quote= F, row.names= F)
