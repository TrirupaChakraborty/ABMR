source("/ix/djishnu/Trirupa/ABomics.Prj/scripts/clinicparams_LFs.corrFx.R")
setwd("/ix/djishnu/Trirupa/ABomics.Prj/Figures/")

canada_meta=read_excel("/ix/djishnu/Trirupa/ABomics.Prj/input_data/Dec2023_ClinicalandHistologyMetaData_Canada.xlsx")

####### interpreting Banff scores #######
## Activity = g+ptc+v+c4d
## Chronicity= (cg*2)+ci+ct+cv

canada_meta['Activity'] = (canada_meta$g) + (canada_meta$ptc) + (canada_meta$v) + (canada_meta$c4d)
canada_meta['Chronicity'] = (2*canada_meta$cg) + (canada_meta$ci) +(canada_meta$ct) + (canada_meta$cv)

############# loading the significant Zs for EARLY Abomics in Canada cohort #########
early_cp_yaml = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/input.yaml"
VC_early_x = "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/early_x.csv"
VC_early_y =  "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/early_y.csv"
earlyVC_sigZ <- load_valSigZ(valX_path=VC_early_x, valY_path=VC_early_y, yaml_path=early_cp_yaml)
earlyVC_sigZ_lis <- colnames(earlyVC_sigZ)

################ calculating correlation between different clinical params with the sig Zs ##########
##### with Activity scores ########
earlyZ_Actdf= create.LF_clin.df(meta_df= canada_meta, sigZ_df = earlyVC_sigZ,
                                  clin_param = 'Activity',
                                  rename_param = "Activity")

earlyZ_Actdf

#### Z1 ######  
Z1_Actdf= earlyZ_Actdf[, c("Z1","Activity" )]
Z1_ActFiltdf = Z1_Actdf[!(Z1_Actdf$Z1 == max(Z1_Actdf$Z1)),] ##### dropping outlier 
Z1_Act.plot= cor_scatterplt(Z1_ActFiltdf, Z_name="Z1", clinc_param = "Activity",
                              plt_title = "Correlation between Z1 and Activity score",
                              x_title = "Z1 values", y_title = "Activity", xaxis_len = -1.2)
Z1_Act.plot
file_name = "./Abomics_Early/Z1_activity.corrplot.pdf"
ggsave(filename=file_name,plot=Z1_Act.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)

#### Z12 ######
Z12_Actdf= earlyZ_Actdf[, c("Z12","Activity" )]
#Z12_ActFiltdf = Z12_Actdf[!(Z12_Actdf$Z1 == max(Z12_Actdf$Z1)),] ##### dropping outlier 
Z12_Act.plot= cor_scatterplt(Z12_Actdf, Z_name="Z12", clinc_param = "Activity",
                            plt_title = "Correlation between Z12 and Activity score",
                            x_title = "Z12 values", y_title = "Activity", xaxis_len = -1.5)
Z12_Act.plot
file_name = "./Abomics_Early/Z12_activity.corrplot.pdf"
ggsave(filename=file_name,plot=Z12_Act.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)


#### with Chronicity scores ########
earlyZ_Chrondf= create.LF_clin.df(meta_df= canada_meta, sigZ_df = earlyVC_sigZ,
                                clin_param = 'Chronicity',
                                rename_param = "Chronicity")

earlyZ_Chrondf

##### Z1 ####
Z1_Chrondf= earlyZ_Chrondf[, c("Z1","Chronicity" )]
#Z1_ChronFiltdf = Z1_Chrondf[!(Z1_Chrondf$Z1 == max(Z1_Chrondf$Z1)),] ##### dropping outlier 
Z1_Chron.plot= cor_scatterplt(Z1_ChronFiltdf, Z_name="Z1", clinc_param = "Chronicity",
                            plt_title = "Correlation between Z1 and Chronicity score",
                            x_title = "Z1 values", y_title = "Chronicity", xaxis_len = -1.2)
Z1_Chron.plot
file_name = "./Abomics_Early/Z1_chronicity.corrplot.pdf"
ggsave(filename=file_name,plot=Z1_Chron.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)

#### Z12 ######
Z12_Chrondf= earlyZ_Chrondf[, c("Z12","Chronicity" )]
#Z12_ActFiltdf = Z12_Actdf[!(Z12_Actdf$Z1 == max(Z12_Actdf$Z1)),] ##### dropping outlier 
Z12_Chron.plot= cor_scatterplt(Z12_Chrondf, Z_name="Z12", clinc_param = "Chronicity",
                             plt_title = "Correlation between Z12 and Chronicity score",
                             x_title = "Z12 values", y_title = "Chronicity", xaxis_len = -1.5)
Z12_Chron.plot
file_name = "./Abomics_Early/Z12_chronicity.corrplot.pdf"
ggsave(filename=file_name,plot=Z12_Chron.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)


############# loading the significant Zs for EARLY Abomics in Canada cohort #########
late_cp_yaml = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/input.yaml"
VC_late_x = "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/late_x.csv"
VC_late_y =  "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/late_y.csv"
VC_late.sigZ <- load_valSigZ(valX_path=VC_late_x, valY_path=VC_late_y, yaml_path=late_cp_yaml)
VC_late.sigZ_lis <- colnames(VC_late.sigZ) 

################ calculating correlation between different clinical params with the sig Zs ##########
##### with Activity scores ########
lateZ_Actdf= create.LF_clin.df(meta_df= canada_meta, sigZ_df = VC_late.sigZ,
                                clin_param = 'Activity',
                                rename_param = "Activity")

lateZ_Actdf

#### Z1 ######  
lateZ1_Actdf= lateZ_Actdf[, c("Z1","Activity" )]
#lateZ1_ActFiltdf = lateZ1_Actdf[!(lateZ1_Actdf$Z1 == max(lateZ1_Actdf$Z1)),] ##### dropping outlier 
lateZ1_Act.plot= cor_scatterplt(lateZ1_Actdf, Z_name="Z1", clinc_param = "Activity",
                            plt_title = "Correlation between Z1 and Activity score",
                            x_title = "Z1 values", y_title = "Activity", xaxis_len = -1.0)
lateZ1_Act.plot
file_name = "./Abomics_Late/lateZ1_activity.corrplot.pdf"
ggsave(filename=file_name,plot=lateZ1_Act.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)

#### Z2 ######  
lateZ2_Actdf= lateZ_Actdf[, c("Z2","Activity" )]
#lateZ1_ActFiltdf = lateZ1_Actdf[!(lateZ1_Actdf$Z1 == max(lateZ1_Actdf$Z1)),] ##### dropping outlier 
lateZ2_Act.plot= cor_scatterplt(lateZ2_Actdf, Z_name="Z2", clinc_param = "Activity",
                                plt_title = "Correlation between Z2 and Activity score",
                                x_title = "Z2 values", y_title = "Activity", xaxis_len = -2.0)
lateZ2_Act.plot
file_name = "./Abomics_Late/lateZ2_activity.corrplot.pdf"
ggsave(filename=file_name,plot=lateZ2_Act.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)

#### Z7 ######  
lateZ7_Actdf= lateZ_Actdf[, c("Z7","Activity" )]
#lateZ1_ActFiltdf = lateZ1_Actdf[!(lateZ1_Actdf$Z1 == max(lateZ1_Actdf$Z1)),] ##### dropping outlier 
lateZ7_Act.plot= cor_scatterplt(lateZ7_Actdf, Z_name="Z7", clinc_param = "Activity",
                                plt_title = "Correlation between Z7 and Activity score",
                                x_title = "Z7 values", y_title = "Activity", xaxis_len = -1.0)
lateZ7_Act.plot
file_name = "./Abomics_Late/lateZ7_activity.corrplot.pdf"
ggsave(filename=file_name,plot=lateZ7_Act.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)

#### Z10 ######  
lateZ10_Actdf= lateZ_Actdf[, c("Z10","Activity" )]
#lateZ1_ActFiltdf = lateZ1_Actdf[!(lateZ1_Actdf$Z1 == max(lateZ1_Actdf$Z1)),] ##### dropping outlier 
lateZ10_Act.plot= cor_scatterplt(lateZ10_Actdf, Z_name="Z10", clinc_param = "Activity",
                                plt_title = "Correlation between Z10 and Activity score",
                                x_title = "Z10 values", y_title = "Activity", xaxis_len = -1.0)
lateZ10_Act.plot
file_name = "./Abomics_Late/lateZ10_activity.corrplot.pdf"
ggsave(filename=file_name,plot=lateZ10_Act.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)


##### with Chronicity scores ########
lateZ_Chrondf= create.LF_clin.df(meta_df= canada_meta, sigZ_df = VC_late.sigZ,
                               clin_param = 'Chronicity',
                               rename_param = "Chronicity")

lateZ_Chrondf

#### Z1 ######  
lateZ1_Chrondf= lateZ_Chrondf[, c("Z1","Chronicity" )]
#lateZ1_ActFiltdf = lateZ1_Actdf[!(lateZ1_Actdf$Z1 == max(lateZ1_Actdf$Z1)),] ##### dropping outlier 
lateZ1_Chron.plot= cor_scatterplt(lateZ1_Chrondf, Z_name="Z1", clinc_param = "Chronicity",
                                plt_title = "Correlation between Z1 and Chronicity score",
                                x_title = "Z1 values", y_title = "Chronicity", xaxis_len = -1.0)
lateZ1_Chron.plot
file_name = "./Abomics_Late/lateZ1_chronicity.corrplot.pdf"
ggsave(filename=file_name,plot=lateZ1_Chron.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)

#### Z2 ######  
lateZ2_Chrondf= lateZ_Chrondf[, c("Z2","Chronicity" )]
#lateZ1_ActFiltdf = lateZ1_Actdf[!(lateZ1_Actdf$Z1 == max(lateZ1_Actdf$Z1)),] ##### dropping outlier 
lateZ2_Chron.plot= cor_scatterplt(lateZ2_Chrondf, Z_name="Z2", clinc_param = "Chronicity",
                                plt_title = "Correlation between Z2 and Chronicity score",
                                x_title = "Z2 values", y_title = "Chronicity", xaxis_len = -2.0)
lateZ2_Chron.plot
file_name = "./Abomics_Late/lateZ2_chronicity.corrplot.pdf"
ggsave(filename=file_name,plot=lateZ2_Chron.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)

#### Z7 ######  
lateZ7_Chrondf= lateZ_Chrondf[, c("Z7","Chronicity" )]
#lateZ1_ActFiltdf = lateZ1_Actdf[!(lateZ1_Actdf$Z1 == max(lateZ1_Actdf$Z1)),] ##### dropping outlier 
lateZ7_Chron.plot= cor_scatterplt(lateZ7_Chrondf, Z_name="Z7", clinc_param = "Chronicity",
                                plt_title = "Correlation between Z7 and Chronicity score",
                                x_title = "Z7 values", y_title = "Chronicity", xaxis_len = -1.0)
lateZ7_Chron.plot
file_name = "./Abomics_Late/lateZ7_chronicity.corrplot.pdf"
ggsave(filename=file_name,plot=lateZ7_Chron.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)

#### Z10 ######  
lateZ10_Chrondf= lateZ_Chrondf[, c("Z10","Chronicity" )]
#lateZ1_ActFiltdf = lateZ1_Actdf[!(lateZ1_Actdf$Z1 == max(lateZ1_Actdf$Z1)),] ##### dropping outlier 
lateZ10_Chron.plot= cor_scatterplt(lateZ10_Chrondf, Z_name="Z10", clinc_param = "Chronicity",
                                 plt_title = "Correlation between Z10 and Chronicity score",
                                 x_title = "Z10 values", y_title = "Chronicity", xaxis_len = -1.0)
lateZ10_Chron.plot
file_name = "./Abomics_Late/lateZ10_chronicity.corrplot.pdf"
ggsave(filename=file_name,plot=lateZ10_Chron.plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)
