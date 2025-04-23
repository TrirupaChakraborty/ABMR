############ DEG analysis between Bcells in low risk vs high risk
library(Matrix)
library(spacexr)
library(spatstat.explore)
library(Seurat)
library(ggplot2)
library(dplyr)

low_risk= "14808"
high_risk= "15965"
base_dir= "/ix/djishnu/Trirupa/CHERUKURI/output_ref4apr25/"

high_risk.obj= readRDS(paste0(base_dir, high_risk,"/step12_",high_risk,"_16um.rds"))
low_risk.obj= readRDS(paste0(base_dir, low_risk,"/step12_",low_risk,"_16um.rds"))

##### subsetting B cells ##
Bcell.high_risk.obj= subset(high_risk.obj, idents ="B Cell")
Bcell.low_risk.obj= subset(low_risk.obj, idents ="B Cell")

seurat_high_risk=Bcell.high_risk.obj
seurat_low_risk=Bcell.low_risk.obj

# Preprocess high-risk dataset
seurat_high_risk <- NormalizeData(seurat_high_risk)
seurat_high_risk <- FindVariableFeatures(seurat_high_risk, selection.method = "vst", nfeatures = 18085)
seurat_high_risk <- ScaleData(seurat_high_risk)

# Preprocess low-risk dataset
seurat_low_risk <- NormalizeData(seurat_low_risk)
seurat_low_risk <- FindVariableFeatures(seurat_low_risk, selection.method = "vst", nfeatures = 18085) ## 18085 genes
seurat_low_risk <- ScaleData(seurat_low_risk)

seurat_high_risk$sample_origin <- "high_risk"
seurat_low_risk$sample_origin <- "low_risk"

# Create a list of Seurat objects
seurat_list <- list(high_risk = seurat_high_risk, low_risk = seurat_low_risk)

# Select features for integration
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 18085) 

# Find integration anchors
integration_anchors <- FindIntegrationAnchors( 
  object.list = seurat_list,
  anchor.features = features,
  dims = 1:10
)
############### breaking here since there are <50 cells in low-risk samples
# Integrate the datasets
combined.obj <- IntegrateData(
  anchorset = integration_anchors,
  dims = 1:10
)

# Set default assay and preprocess
DefaultAssay(combined.obj) <- "integrated"
combined.obj <- ScaleData(combined.obj)
combined.obj <- RunPCA(combined.obj, npcs = 30)
combined.obj <- RunUMAP(combined.obj, reduction = "pca", dims = 1:30)

# Plot UMAP with sample origin labels
DimPlot(combined.obj, reduction = "umap", group.by = "sample_origin", pt.size = 0.5) +
  ggtitle("UMAP of Integrated Data by Sample Origin")


# Set the identity to sample origin
Idents(combined.obj) <- "sample_origin"

# Perform DEG analysis
# Find DEGs between high_risk and low_risk samples
deg_results <- FindMarkers(
  combined.obj,
  ident.1 = "high_risk",
  ident.2 = "low_risk",
  group.by = "sample_origin",
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# View the results
head(deg_results)

top10_genes <- rownames(deg_results) %>% head(30)

# Plot heatmap
htmap= DoHeatmap(combined.obj, features = top10_genes, group.by = "sample_origin")
ggsave(htmap,file= paste0(base_dir,"15965vs14808.Trm.DEGheatmap.18kfeats.pdf"))
saveRDS(combined.obj,paste0(base_dir,"15965vs14808.Trm_combined.seurat.18kfeats.rds"))
saveRDS(deg_results,paste0(base_dir,"15965vs14808.Trm_DEGresults.18kfeats.rds"))

##################################################################################################
##################################################################################################