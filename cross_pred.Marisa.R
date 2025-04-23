
source("/ix/djishnu/Marisa/ABomics/cross_preds/cp_SLIDE.R")
# cross prediction function 


# slide yaml file
early_yaml = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/input.yaml"
#new_outpath = "/ix/djishnu/Trirupa/ABomics.Prj/cross_pred/DC_Early.VC_Late.out/"

#validation dataset x and y 
val_x_path = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/x.csv"
val_y_path =  "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/y.csv"


# scale needs to be set to FALSE
slide_crossPred_Early<- load_data_cp(yaml_path = early_yaml, 
                        val_x = val_x_path,
                        val_y = val_y_path,
                        interactions = TRUE, 
                        scale = TRUE, 
                        new = NULL)

late_yaml = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/input.yaml"
#new_outpath = "/ix/djishnu/Trirupa/ABomics.Prj/cross_pred/DC_Early.VC_Late.out/"

#validation dataset x and y 
val_x_path = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/x.csv"
val_y_path =  "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/y.csv"


# scale needs to be set to FALSE
slide_crossPred_Early<- load_data_cp(yaml_path = late_yaml, 
                                     val_x = val_x_path,
                                     val_y = val_y_path,
                                     interactions = TRUE, 
                                     scale = TRUE, 
                                     new = NULL)


##### Early R vs Late R ######
rejector_yaml= ""

rejector_x_path= "/ix/djishnu/Marisa/ABomics/canada_cohort/SLIDE_inputs/rejector_x.csv"
rejector_y_path= "/ix/djishnu/Marisa/ABomics/canada_cohort/SLIDE_inputs/rejector_y.csv"