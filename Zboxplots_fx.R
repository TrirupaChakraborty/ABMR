#### This script has all functions used for generating Z-level boxplots ###

library(EssReg)
library(ggplot2)
library(yaml)
library(dplyr)
library(pROC)
library(ROCR)
library(matrixStats)
library(ggpubr)
library(reshape2)
library(SLIDE)
############ the following functions will be used for plotting Validation cohort Z-boxplots ##############

## predZ ################################################################################################################
predZ <- function(x, er_res) {
  A_hat <- er_res$A
  C_hat <- er_res$C
  Gamma_hat <- er_res$Gamma
  Gamma_hat <- ifelse(Gamma_hat == 0, 1e-10, Gamma_hat)
  Gamma_hat_inv <- diag(Gamma_hat ** (-1))
  G_hat <- crossprod(A_hat, Gamma_hat_inv) %*% A_hat + solve(C_hat)
  Z_hat <- x %*% Gamma_hat_inv %*% A_hat %*% MASS::ginv(G_hat)
  return (Z_hat)
}

### Pairwise interaxtions #########################################################
pairwiseInteractions <- function(index_list, mat) {
  num_cols <- ncol(mat)
  index_combinations <- expand.grid(seq_len(num_cols),index_list,stringsAsFactors=F)
  temp <- index_combinations$Var1
  index_combinations$Var1 <- index_combinations$Var2
  index_combinations$Var2 <- temp
  
  col_names <- paste0(colnames(mat)[as.numeric(index_combinations[, 1])], ".", colnames(mat)[as.numeric(index_combinations[, 2])])
  interaction_mat <- mat[, as.numeric(index_combinations[, 1])] * mat[, as.numeric(index_combinations[, 2])]
  if(is.null(dim(interaction_mat))){
    interaction_mat <- matrix(interaction_mat,nrow=1)}
  
  colnames(interaction_mat) <- col_names
  return(list(interaction=as.data.frame(interaction_mat)))
}

######## fucntion to get sig LF df #########################################################
getLFdf<- function(SLIDE_res, z_matrix, interactions = TRUE){
  # marginal data
  sigK <- SLIDE_res$marginal_vals
  sigMargData <- z_matrix[,sigK]
  
  if(!interactions){
    return(as.data.frame(sigMargData))
  }
  
  # interaction
  sigIn <- SLIDE_res$SLIDE_res$interaction_vars
  if(!is.null(sigIn)){
    IntData <- pairwiseInteractions(sigK,z_matrix)
    sigIntData <- IntData$interaction[ ,sigIn]
    final_mtx <- cbind(as.data.frame(sigMargData), as.data.frame(sigIntData))
    return(final_mtx)
  }
  return(as.data.frame(sigMargData))
}

####### function to get significant Zs (LFs) from the validation cohort ########
load_valSigZ <- function(valX_path,valY_path,yaml_path){
  val_x=as.matrix(read.csv(valX_path, row.names=1))
  val_y=as.matrix(read.csv(valY_path, row.names=1))
  input <- read_yaml(yaml_path)
  
  train_x <- as.matrix(read.csv(input$x_path, row.names=1))
  train_y <- as.matrix(read.csv(input$y_path, row.names=1))
  
  # if validation dataset has aditional features, remove them
  val_not_in_train <- colnames(val_x)[which(!colnames(val_x) %in% colnames(train_x))]
  val_x <- scale(as.matrix(val_x))
  train_x <- scale(as.matrix(train_x))
  
  
  er_results = readRDS(list.files(input$out_path, full.names = T,  pattern = "final_"))
  slide_res = readRDS(list.files(input$out_path,recursive = T, full.names = T, pattern = "slide_res")[1])
  
  #### predicting the LFs in the Validation cohort using the group structure (in er_result) from discover cohort ###
  val_z <- predZ(val_x, er_results)
  colnames(val_z) <- paste0("Z", c(1:ncol(val_z)))
  ##### getting the significant LFs (val_sigZ) #####  
  val_sigZ = getLFdf(slide_res, val_z, interactions=TRUE)
  
  return(val_sigZ)
}

################# this is the general function for plotting Zs ######
Zbox_plot <- function(Znum,colour_list,dfsigZ,dfY_path,upper_limit,lower_limit){
  #dfZ=data.frame(scale(dfZ_input))
  dfY=read.csv(dfY_path, row.names=1)
  ##renaming Ys
  dfZ_Y=cbind(dfY,dfsigZ)
  dfZ_Y$Y[dfZ_Y$Y == 1] <- "DSA+AbMR+"
  dfZ_Y$Y[dfZ_Y$Y == 0] <- "DSA+AbMR-"
  
  sigZ_lis= colnames(dfsigZ)
  
  value <- dfsigZ[ ,Znum]
  df1 <- as.data.frame(value)
  df1["Status"] <- dfZ_Y[, 1]
  
  signif_labels=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
  zplot=ggplot2::ggplot(data=df1,aes(x=Status,y=value))+
    ggplot2::geom_boxplot(width=0.40, aes(fill=as.factor(Status))) + 
    scale_fill_manual(values=colour_list)+
    ggtitle(Znum) +theme(plot.title = element_text(hjust = 0.5),
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), 
                         panel.spacing = unit(0.5, "lines"),
                         panel.border = element_rect(fill = NA, color = "black", linetype = "dashed"),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black")) +
    stat_compare_means(method = "wilcox.test", label = "p.signif",
                       symnum.args = signif_labels,
                       comparison = list(c("DSA+AbMR+", "DSA+AbMR-")),
                       aes(label = paste0("p = ", after_stat(p.format))),
                       label.x.npc = "middle", label.y.npc = "top")
  
  if (!is.null(upper_limit)) {
    zplot <- zplot + ylim(min(df1$value), upper_limit)
  }
  if (!is.null(lower_limit)){
    zplot <- zplot + ylim(lower_limit,max(df1$value))
  }
  return(zplot)
}

Zdot_plot <- function(Znum, colour_list, dfsigZ, dfY_path, upper_limit, lower_limit) {
  # Read Y dataframe
  dfY <- read.csv(dfY_path, row.names = 1)
  
  # Rename Ys
  dfZ_Y <- cbind(dfY, dfsigZ)
  dfZ_Y$Y[dfZ_Y$Y == 1] <- "DSA+AbMR+"
  dfZ_Y$Y[dfZ_Y$Y == 0] <- "DSA+AbMR-"
  
  # Extract Z values
  value <- dfsigZ[, Znum]
  df1 <- as.data.frame(value)
  df1$Status <- dfZ_Y[, 1]
  
  # Define significance labels
  signif_labels <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
  
  # Create dot plot
  zplot <- ggplot(data = df1, aes(x = Status, y = value, color = Status)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    scale_color_manual(values = colour_list) +
    ggtitle(Znum) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.spacing = unit(0.5, "lines"),
          panel.border = element_rect(fill = NA, color = "black", linetype = "dashed"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    stat_summary(fun = median, geom = "errorbar", aes(ymin = ..y.., ymax = ..y.., group = Status),
                 width = 0.2,  # Controls the width of the median lines
                 color = "black")+
    stat_compare_means(method = "wilcox.test", label = "p.signif",
                       symnum.args = signif_labels,
                       comparison = list(c("DSA+AbMR+", "DSA+AbMR-")),
                       aes(label = paste0("p = ", after_stat(p.format))),
                       label.x.npc = "middle", label.y.npc = "top") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)),
                       breaks = seq(floor(min(df1$value)), ceiling(max(df1$value)), by = 0.5))  # Adjust y-tick intervals
  
  # Adjust ylim if upper_limit or lower_limit is specified
  if (!is.null(upper_limit)) {
    zplot <- zplot + ylim(min(df1$value), upper_limit)
  }
  if (!is.null(lower_limit)) {
    zplot <- zplot + ylim(lower_limit, max(df1$value))
  }
  
  return(zplot)
}


