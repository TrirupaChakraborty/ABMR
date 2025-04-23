library(ggplot2)
library(ggpubr)

############### loading the plotting function ####################
dot_plt = function(df_name,x_name,y_name){
  
  dot_plot <- ggplot(data = df_name, aes_string(x = x_name, y = y_name, color = x_name)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    scale_color_manual(values = colour_list) +
    ggtitle(colnames(df_name)[2]) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.spacing = unit(0.5, "lines"),
          panel.border = element_rect(fill = NA, color = "black", linetype = "dashed"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    stat_summary(fun = median, geom = "errorbar", aes(ymin = ..y.., ymax = ..y.., group = GRP),
                 width = 0.2,  # Controls the width of the median lines
                 color = "black") +
    stat_compare_means(symnum.args = signif_labels,
                       comparison = list(c("DSA", "ABMR")),
                       aes(label = paste0("p = ", after_stat(p.format))),
                       label.x.npc = "middle", label.y.npc = "top")
  
  return(dot_plot)
}

################# loading data ####################
pitt_meta=read.csv("/ix/djishnu/Marisa/ABomics/pitt_cohort/meta_36.csv")

colour_list = c("DSA"="#EECA59","ABMR"="#CA4F73")
signif_labels <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                      symbols = c("****", "***", "**", "*", "ns"))

################# univariate dot plots ####################

################## eGFR.at.Bx ###################
egfrBx_df=pitt_meta[,c("GRP","eGFR.at.Bx")]
n_60= sum(grepl(">60",egfrBx_df$eGFR.at.Bx))
## >60 value was set to 80 in eGFR_atBx. and eGFR_atSample
egfrBx_df$eGFR.at.Bx = gsub(">60",80, egfrBx_df$eGFR.at.Bx) 
egfrBx_df$eGFR.at.Bx=as.numeric(egfrBx_df$eGFR.at.Bx)

egfr.atBxdot_plot= dot_plt(df_name=egfrBx_df, x_name = "GRP", y_name="eGFR.at.Bx")
print(egfr.atBxdot_plot)

file_name = "/ix/djishnu/Trirupa/ABomics.Prj/Figures/eGFR_atBx_dotplot.pdf"
ggsave(filename=file_name,plot=egfr.atBxdot_plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)


################## eGFR.at.sample ###################
egfrSampl_df=pitt_meta[,c("GRP","eGFR.at.sample")]
n_60= sum(grepl(">60",egfrSampl_df$eGFR.at.sample)) ### 5 AbMR patients had >60 value
## >60 value was set to 80 in eGFR_atBx. and eGFR_atSample
egfrSampl_df$eGFR.at.sample = gsub(">60",80, egfrSampl_df$eGFR.at.sample) 
egfrSampl_df = egfrSampl_df[!egfrSampl_df$eGFR.at.sample == "cyclosporine",] #### dropping the DSA only patient that had "cyclosporine" in place of egFR value
egfrSampl_df$eGFR.at.sample = gsub("<2",0.5, egfrSampl_df$eGFR.at.sample) ### setting <2 to 0.5 value
egfrSampl_df$eGFR.at.sample=as.numeric(egfrSampl_df$eGFR.at.sample)

egfr.atSdot_plot= dot_plt(df_name=egfrSampl_df, x_name = "GRP", y_name="eGFR.at.sample")
print(egfr.atSdot_plot)

file_name = "/ix/djishnu/Trirupa/ABomics.Prj/Figures/eGFR_atSample_dotplot.pdf"
ggsave(filename=file_name,plot=egfr.atSdot_plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)


