# Install necessary packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("ReactomePA")
#BiocManager::install("KEGGREST")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(readxl)
#library(KEGGREST)
setwd("/ix/djishnu/Trirupa/ABomics.Prj/SUHANA/")

extract_TFpathway = function(tf_path){
  tf_list <- tf_path$Symbol
  # Convert gene symbols to Entrez IDs
  gene_ids <- bitr(tf_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Perform KEGG pathway enrichment analysis
  kegg_enrich <- enrichKEGG(gene = gene_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
  kegg_enrich_readable <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  kegg_allres= as.data.frame(kegg_enrich_readable@result)
  
  # Perform Reactome pathway enrichment analysis
  reactome_enrich <- enrichPathway(gene = gene_ids$ENTREZID, organism = 'human', pvalueCutoff = 0.05)
  reactome_enrich_readable <- setReadable(reactome_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  reactome_allres = as.data.frame(reactome_enrich_readable@result)
  
  return(list(KEGG = kegg_allres, Reactome = reactome_allres))
  
}


##### for the 4/7 miRs that were common across all LFs in DSA vs AbMR analysis. The target genes were not thresholded 
# for the number of miRs that target that gene. i.e. these genes had atleast 1 miR targetting them.
commonmiR.inLFs_path= read_excel("/ix/djishnu/Trirupa/ABomics.Prj/SUHANA/LFcommon_miR_TFs.AnimalTFDB.xls")
commonmiR.inLFs_allres= extract_TFpathway(commonmiR.inLFs_path)

kegg_allres= commonmiR.inLFs_path$KEGG
reactome_allres= commonmiR.inLFs_path$Reactome

write.csv(kegg_allres, "kegg_commonmiR_TF.allpathways.csv", quote=F, row.names = F)
write.csv(reactome_allres, "reactome_commonmiR_TF.allpathways.csv", quote=F, row.names = F)

###### for all the miRs present in the LFs. the genes were thresholded at >=5 miRs targeting that gene. 
miRs_19.inLFs_path= read_excel("/ix/djishnu/Trirupa/ABomics.Prj/SUHANA/19miR.target_TFs.AnimalTFDB.xls")
miRs_19.inLFs_allres = extract_TFpathway(miRs_19.inLFs_path)

miRs_19.kegg_allres= miRs_19.inLFs_allres$KEGG
miRs_19.reactome_allres= miRs_19.inLFs_allres$Reactome

miRs_19.kegg_allres[miRs_19.kegg_allres$ID == "R-HSA-2559583",]
