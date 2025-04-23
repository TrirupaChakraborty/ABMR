rm(list = ls())
cat("\014")

library(ggplot2)
library(dplyr)
library(ggpubr)
bp_df=readRDS("pipeline_step5.rds")
bp_ER=bp_df %>% filter(!(method %in% c("lasso","lasso_y")))  
Actual=as.vector(bp_ER[which(bp_ER$method == "plainER"),2]) ## these are unused in the code. The code automatically does it in wilcox.test()
Permuted=as.vector(bp_ER[which(bp_ER$method == "plainER_y"),2])
signif_labels=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
lambda_boxplot=ggplot2::ggplot(data = bp_ER,
                                 ggplot2::aes(x = method,
                                              y = auc,fill=method)) +
  ggplot2::geom_boxplot(width=0.20) + scale_x_discrete(labels=c('Actual', 'Permuted')) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                              panel.background = element_blank(), axis.line = element_line(colour = "black"))+ stat_compare_means(symnum.args = signif_labels,aes(label = paste0("p = ", after_stat(p.format))),label.x.npc="middle",label.y.npc="top")+
  #ggplot2::labs(fill = "Method") +
  ggtitle("del=0.05 lam=1.0")+theme(plot.title = element_text(hjust = 0.5, face="bold"))
  ggplot2::scale_alpha(guide = 'none') 

lambda_boxplot 
                                                                  

#lambda_boxplot+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      #panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_x_discrete(labels=c('Actual', 'Permuted'))
#stat_compare_means(method = wilcox.test, comparisons = (“plainER”, “plainER_y”)
ggsave(file="EvsNE.del0.05_lam1.0_pipe3Plot.pdf", lambda_boxplot, width = 4, height = 4, units = "in")
