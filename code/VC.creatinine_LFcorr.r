############## this code calls fxs from clinicparams_LFs.corrFx.R and plots corr plots between clinical params and LFs ############# 
source("/ix/djishnu/Trirupa/ABomics.Prj/scripts/clinicparams_LFs.corrFx.R")
setwd("/ix/djishnu/Trirupa/ABomics.Prj/Figures/")

###### loading data #####
canada_meta=read_excel("/ix/djishnu/Trirupa/ABomics.Prj/input_data/Dec2023_ClinicalandHistologyMetaData_Canada.xlsx")

############# loading the significant Zs for EARLY Abomics in Canada cohort #########
early_cp_yaml = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/input.yaml"
VC_early_x = "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/early_x.csv"
VC_early_y =  "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/early_y.csv"
earlyVC_sigZ <- load_valSigZ(valX_path=VC_early_x, valY_path=VC_early_y, yaml_path=early_cp_yaml)
earlyVC_sigZ_lis <- colnames(earlyVC_sigZ)

################ calculating correlation between different clinical params with the sig Zs ##########

########## Z1 #########
earlyZ_creatdf= create.LF_clin.df(meta_df= canada_meta, sigZ_df = earlyVC_sigZ,
                             clin_param = 'creatinine at the time of biopsy',
                             rename_param = "creat.atBx")

earlyZ_creatdf

#### Z1 ######
Z1_creatdf= earlyZ_creatdf[, c("Z1","creat.atBx" )]
Z1_creatFiltdf = Z1_creatdf[-nrow(Z1_creatdf),] ##### dropping the outlier 
Z1_creat.plot= cor_scatterplt(Z1_creatFiltdf, Z_name="Z1", clinc_param = "creat.atBx",
                              plt_title = "Correlation between Z1 and Creatinine at Bx",
                              x_title = "Z1 values", y_title = "Creatinine at Bx", xaxis_len = -1.5)
Z1_creat.plot
file_name = "./Abomics_Early/Z1_creatinine.corrplot.pdf"
ggsave(filename=file_name,plot=Z1_creat.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)

######### Z12 #########
Z12_creatdf= earlyZ_creatdf[, c("Z12","creat.atBx" )]
#Z_creatFiltdf = Z1_creatdf[-nrow(Z1_creatdf),] ##### dropping the outlier 
Z12_creat.plot= cor_scatterplt(earlyZ_creatdf, Z_name="Z12", clinc_param = "creat.atBx",
                              plt_title = "Correlation between Z12 and Creatinine at Bx",
                              x_title = "Z12 values", y_title = "Creatinine at Bx", xaxis_len = -4.0)
Z12_creat.plot
file_name = "./Abomics_Early/Z12_creatinine.corrplot.pdf"
ggsave(filename=file_name,plot=Z12_creat.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)
#########  XXXXXXXXXXXXXXXXXXXXXXXX ############### XXXXXXXXXXXXXXXXXXXXX ###################################

############# loading the significant Zs for LATE Abomics in Canada cohort ############
late_cp_yaml = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/input.yaml"
VC_late_x = "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/late_x.csv"
VC_late_y =  "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/late_y.csv"

###### first getting the significant LFs from the VC (from cross-prediction) #####
VC_late.sigZ <- load_valSigZ(valX_path=VC_late_x, valY_path=VC_late_y, yaml_path=late_cp_yaml)
VC_late.sigZ_lis <- colnames(VC_late.sigZ) 

####### calculating correlation between different clinical params with the sig Zs ##########
lateZ_creatdf= create.LF_clin.df(meta_df= canada_meta, sigZ_df = VC_late.sigZ,
                                  clin_param = 'creatinine at the time of biopsy',
                                  rename_param = "creat.atBx")

lateZ_creatdf
lateZ_creatFiltdf = lateZ_creatdf[!(lateZ_creatdf$creat.atBx == 2290),] ## dropping the outlier

#### Z1 ######
lateZ1_creatdf= lateZ_creatFiltdf[, c("Z1","creat.atBx" )]
lateZ1_creat.plot= cor_scatterplt(lateZ1_creatdf, Z_name="Z1", clinc_param = "creat.atBx",
                              plt_title = "Correlation between Z1 and Creatinine at Bx",
                              x_title = "Z1 values", y_title = "Creatinine at Bx", xaxis_len = -1.0)
lateZ1_creat.plot
file_name = "./Abomics_Late/lateZ1_creatinine.corrplot.pdf"
ggsave(filename=file_name,plot=lateZ1_creat.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)

##### Z2 ####
lateZ2_creatdf= lateZ_creatFiltdf[, c("Z2","creat.atBx" )]
lateZ2_creat.plot= cor_scatterplt(lateZ2_creatdf, Z_name="Z2", clinc_param = "creat.atBx",
                                  plt_title = "Correlation between Z2 and Creatinine at Bx",
                                  x_title = "Z2 values", y_title = "Creatinine at Bx", xaxis_len = -1.0)
lateZ2_creat.plot
file_name = "./Abomics_Late/lateZ2_creatinine.corrplot.pdf"
ggsave(filename=file_name,plot=lateZ2_creat.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)

##### Z7 ####
lateZ7_creatdf= lateZ_creatFiltdf[, c("Z7","creat.atBx" )]
lateZ7_creatFiltdf = lateZ7_creatdf[!(lateZ7_creatdf$Z7 == max(lateZ7_creatdf$Z7)), ]
lateZ7_creat.plot= cor_scatterplt(lateZ7_creatFiltdf, Z_name="Z7", clinc_param = "creat.atBx",
                                  plt_title = "Correlation between Z7 and Creatinine at Bx",
                                  x_title = "Z7 values", y_title = "Creatinine at Bx", xaxis_len = -1.0)
lateZ7_creat.plot
file_name = "./Abomics_Late/lateZ7_creatinine.corrplot.pdf"
ggsave(filename=file_name,plot=lateZ7_creat.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)

##### Z10 ####
lateZ10_creatdf= lateZ_creatFiltdf[, c("Z10","creat.atBx" )]
#lateZ10_creatFiltdf = lateZ10_creatdf[!(lateZ10_creatdf$Z7 == max(lateZ7_creatdf$Z7)), ]
lateZ10_creat.plot= cor_scatterplt(lateZ10_creatdf, Z_name="Z10", clinc_param = "creat.atBx",
                                  plt_title = "Correlation between Z10 and Creatinine at Bx",
                                  x_title = "Z10 values", y_title = "Creatinine at Bx", xaxis_len = -1.0)
lateZ10_creat.plot
file_name = "./Abomics_Late/lateZ10_creatinine.corrplot.pdf"
ggsave(filename=file_name,plot=lateZ10_creat.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)
