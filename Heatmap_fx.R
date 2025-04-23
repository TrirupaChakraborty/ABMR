library(ggplot2) 
library(dplyr)
library(plyr)
library(tidyr) 
library(stringr) 
library(RColorBrewer)
library(sjPlot)
library(pheatmap)
library(genefilter)



plot_heatmap <- function(X_id_df, Y_df, Y0_name,Y1_name,annot_colors, savePath, plotSigGenesRDS){
  X_df=read.csv(X_id_df,row.names=1)
  Y_df=read.csv(Y_df, row.names=1)
  XY_df=cbind(Y_df,X_df) #### note: the X and Y dfs rows should corespond to teh same PatID
  
  #### extracting the set of features in all the significant Zs
  sigFeats_df=readRDS(plotSigGenesRDS)
  sigFeats_df=na.omit(sigFeats_df)
  sigFeats=unique(sigFeats_df$names)
  
  Y0_df=XY_df[XY_df$Y ==0,]
  Y1_df=XY_df[XY_df$Y ==1,]
  
  ### calculating the fold change for each feature between the 2 groups of Y
  ### This order will determine the order of rows in the heatmap, with highest fold-change (FC) at the top
  
  FC_df=data.frame(Feature=as.character(), Fold_change=as.numeric())
  for (feat in sigFeats){
    mean_Y0 <- mean(Y0_df[[feat]], na.rm = TRUE)
    mean_Y1 <- mean(Y1_df[[feat]], na.rm = TRUE)
    FC <- (abs(mean_Y1-mean_Y0))/mean_Y0
    new_row=data.frame(Feature=feat, Fold_change=FC)
    FC_df=rbind(FC_df,new_row)
  }
  
  FC_df_ordered <- FC_df[order(FC_df$Fold_change, decreasing = TRUE), ]
  sigfeat_XY_df=cbind(Y_df,XY_df[,FC_df_ordered$Feature]) ### this is the master df we will be subsetting from
  
  ##################### the following code is for plotting the heatmap ##############
  
  ###creating features X patId matrix
  #sigfeat_XY_df$PatID <- 1:(nrow(XY_df))
  Pat_id = row.names(X_df)
  sigfeat_XY_df$PatID= as.numeric(Pat_id)
  sigfeat_XY_df <- sigfeat_XY_df[, c("PatID", names(sigfeat_XY_df)[names(sigfeat_XY_df) != "PatID"])]
  patID_X <- subset(sigfeat_XY_df, select = -Y)
  all_feats=colnames(subset(patID_X, select = -PatID))
  X_feat_patID=as.data.frame(t(patID_X))
  tpmMat=X_feat_patID[-1,]
  
  ##creating patID X group labels matrix. This is colAnnot
  ## sort the features in the feature matrix based on heatmap_featMFIsort.py script and then use the input ##
  patID_grp <- subset(sigfeat_XY_df, select = c(Y))
  rownames(patID_grp) <- sigfeat_XY_df$PatID
  feat_list=data.frame(target_id=all_feats,feat_name=all_feats)
  results <- list()
  
  # Z-transform
  tpmMat <- as.matrix(tpmMat)
  colnames(tpmMat)= rownames(patID_grp)
  row_means <- rowMeans(tpmMat)
  row_stds <- rowSds(tpmMat)
  zMatrix <- (tpmMat - row_means) / row_stds
  
  colnames(zMatrix)=paste0("id_",colnames(zMatrix))
  rownames(patID_grp)=paste0("id_",rownames(patID_grp))
  
  # subset z-matrix for selected feats
  subZMat <- data.frame(zMatrix[rownames(zMatrix) %in% feat_list$target_id, ])
  subZMat_ordered <- cbind(subZMat[feat_list$target_id, ], "feat_name"= feat_list$feat_name)
  
  rownames(subZMat_ordered) <- subZMat_ordered$feat_name
  subZMat2 <- subZMat_ordered[,-ncol(subZMat_ordered)] ##
  results[["zMat"]] <- subZMat2
  
  # Define color scale
  margins <- c(min(subZMat2), max(subZMat2)) 
  breaksList=NA
  if(sum(abs(margins) > 1.5) > 0){
    breaksList <- seq(-1.5, 1.5, by = 0.03)
  } 
  
  feat_list=data.frame(target_id=rownames(subZMat2),feat_name=rownames(subZMat2))
  
  patID_grp[patID_grp$Y == 1, "Y"] <- Y1_name
  patID_grp[patID_grp$Y == 0, "Y"] <- Y0_name
  patID_grp$Y=as.factor(patID_grp$Y)
  colAnnot=patID_grp["Y"]
  
  ##function for saving heatmap
  save.heatmap <- function(x, filename,width=6,height=2.5) {
    pdf(filename,width=width,height=height) 
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  ### heatmap
  hMap <- pheatmap(subZMat2, color=colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100), 
                   breaks=breaksList, show_rownames=TRUE, show_colnames=TRUE,
                   cluster_rows = FALSE, rownames=rownames(subZMat2), cluster_cols = FALSE, 
                   annotation_col = colAnnot, fontsize_row = 4, fontsize_col = 4, annotation_colors=annot_colors,
                   annotation_legend=TRUE, legend=TRUE, border_color = "#f0f0f0", fontsize = 6)
  
  save.heatmap(x = hMap, filename = savePath)
  results[["hMap"]] <- hMap
  
  return(hMap)
}







