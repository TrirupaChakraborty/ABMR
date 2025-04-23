source("/ix/djishnu/Trirupa/ABomics.Prj/scripts/Heatmap_fx.R")

############################# EARLY AB-omics ##############################

X_file="/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/x_id.csv"
Y_file= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/y.csv"
Y0_name ="DSA+ AbMR-"
Y1_name ="DSA+ AbMR+"
colors_list=list(Y=c("DSA+ AbMR-"="#EECA59","DSA+ AbMR+"="#CA4F73"))
savePath="/ix/djishnu/Trirupa/ABomics.Prj/Figures/Abomics_Early/EarlyAbomics_heatmap_PatId.pdf"
RDS_file= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/SLIDE/plotSigGenes_data.RDS"

EarlyAb_plot=plot_heatmap(X_file,Y_file, Y0_name=Y0_name, Y1_name=Y1_name, colors_list,
                        savePath=savePath, plotSigGenesRDS=RDS_file)

############################# LATE AB-omics ##############################

X_file="/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/x_id.csv"
Y_file= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/y.csv"
Y0_name ="DSA+ AbMR-"
Y1_name ="DSA+ AbMR+"
colors_list=list(Y=c("DSA+ AbMR-"="#EECA59","DSA+ AbMR+"="#CA4F73"))
savePath="/ix/djishnu/Trirupa/ABomics.Prj/Figures/Abomics_Late/LateAbomics_heatmap_PatId.pdf"
RDS_file= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/SLIDE/plotSigGenes_data.RDS"

LateAb_plot=plot_heatmap(X_file,Y_file, Y0_name=Y0_name, Y1_name=Y1_name, colors_list,
                        savePath=savePath, plotSigGenesRDS=RDS_file)


