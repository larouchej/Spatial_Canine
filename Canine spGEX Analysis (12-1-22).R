# Canine spGEX Analysis.R
# Jacqueline Larouche
# Fall 2022

# Contains the code used for spGEX data analysis of canine VML tissues

# Workflow: 
# 1. Initialize variables
# 2. Load ST Data and Create List of Seurat Objects, Perform Seurat Label Transfer
# 3. Generate seurat object with all ST datasets
# 4. Plot (UMAPs, Overlays)
# 5. GSEA on Seurat Clusters

#### 0. Initialize Environment ####
library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(Matrix)
library(data.table)
library(RColorBrewer)
library(tibble)
library(cowplot)
library(reshape)
library(amap)
library(hdf5r)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(AnnotationDbi)
library(org.Cf.eg.db)
library(dittoSeq) 
library(scales)
library(clusterProfiler)

setwd("~/Documents/UMichigan/Aguilar_Lab/ResearchProjects/TA_VML/SpatialTranscriptomics/Visium/Data/spGEX")
#### 1. Initialize Variables ####
sample_dict <- list()
sample_dict[["Sample_5513-JL-S2-A_CGCGCACT-AGAATACA"]] = "D7_Mid_1"
sample_dict[["Sample_5513-JL-S2-B_CCTGTCAG-GTTACGGG"]] = "D7_Edge_1"
sample_dict[["Sample_5513-JL-S2-C_GTCCTTCG-CTGTGCAT"]] = "D14_Mid_1"
sample_dict[["Sample_5513-JL-S2-D_AATGTATC-TAAGCTCA"]] = "D14_Edge_1"
sample_dict[["Sample_6539-JL-S1-A_GTGGATCAAA-CAGGGTTGGC"]] = "D0_Intact_R1"
sample_dict[["Sample_6833-JL-S1-A_TGTCCCAA-TGGACATC"]] = "D0_Intact_R2"
sample_dict[["Sample_6834-JL-S1-C-GEX_TCCCAAGG-AAAGGTAG"]] = "D7_Mid_2"
sample_dict[["Sample_6834-JL-S1-D-GEX_GCCTTCGG-AAATCGTT"]] = "D14_Mid_2"
samples <- names(sample_dict)

st_list <- readRDS("Data/Canis_ST_list.RDS")
st_all <- readRDS(file = "Data/canis_merged_D0_D7_D14.RDS")

#### 2. Load ST Data and Create List of Seurat Objects ####
st_list_1 <- sapply(samples[1:4], USE.NAMES = TRUE, function(geo) {
  print(sample_dict[[geo]])
  st_se <- Load10X_Spatial(paste0(getwd(), sprintf("/10xOutputData/Canine/%s", geo)), 
                           assay = "Spatial", filter.matrix = TRUE, 
                           slice = sprintf("%s", sample_dict[[geo]]), 
                           filename = "filtered_feature_bc_matrix.h5")
  # remove spots with <1 gene detected
  st_se <- subset(st_se, subset = nFeature_Spatial > 0)
  st_se$orig.ident <- sample_dict[[geo]]
  st_se$timepoint <- strsplit(sample_dict[[geo]], '_')[[1]][1]
  st_se$location <- strsplit(sample_dict[[geo]], '_')[[1]][2]
  st_se$animal <- strsplit(sample_dict[[geo]], '_')[[1]][3]
  
  st_se <- SCTransform(st_se, assay = "Spatial", verbose = FALSE)
  st_se <- st_se %>% RunPCA(assay = "SCT", verbose = FALSE) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% 
    FindClusters(verbose = FALSE) %>% 
    RunUMAP(reduction = "pca", dims = 1:30)
  return(list(st_se))
})

st_list_2 <- sapply(samples[5:8], USE.NAMES = TRUE, function(geo) {
  print(sample_dict[[geo]])
  st_se <- Load10X_Spatial(paste0(getwd(), sprintf("/10xOutputData/Canine/%s", geo)), 
                           assay = "Spatial", filter.matrix = TRUE, 
                           slice = sprintf("%s", sample_dict[[geo]]), 
                           filename = "filtered_feature_bc_matrix.h5")
  
  # Convert coordinates of the spatial image to integers from characters.
  st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["tissue"]] <- as.integer(st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["tissue"]])
  st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["row"]] <- as.integer(st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["row"]])
  st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["col"]] <- as.integer(st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["col"]])
  st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["imagerow"]] <- as.integer(st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["imagerow"]])
  st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["imagecol"]] <- as.integer(st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["imagecol"]])

  st_se <- subset(st_se, subset = nFeature_Spatial > 0)
  
  st_se$orig.ident <- sample_dict[[geo]]
  st_se$timepoint <- strsplit(sample_dict[[geo]], '_')[[1]][1]
  st_se$location <- strsplit(sample_dict[[geo]], '_')[[1]][2]
  st_se$animal <- strsplit(sample_dict[[geo]], '_')[[1]][3]
  
  st_se <- SCTransform(st_se, assay = "Spatial", verbose = FALSE)
  st_se <- st_se %>% RunPCA(assay = "SCT", verbose = FALSE) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% 
    FindClusters(verbose = FALSE) %>% 
    RunUMAP(reduction = "pca", dims = 1:30)
  return(list(st_se))
})

st_list <- c(st_list_1, st_list_2)
saveRDS(st_list, file = "Data/Canis_ST_list.RDS")

#### 3. Generate seurat object with all ST datasets ####
st_1 <- st_list[[samples[1]]]; st_2 <- st_list[[samples[2]]] 
st_3 <- st_list[[samples[3]]]; st_4 <- st_list[[samples[4]]]
st_5 <- st_list[[samples[5]]]; st_6 <- st_list[[samples[6]]];
st_7 <- st_list[[samples[7]]]; st_8 <- st_list[[samples[8]]];

st_all <- merge(st_1, c(st_2, st_3, st_4, st_5, st_6, st_7, st_8))
table(st_all$timepoint, st_all$location)
table(st_all$timepoint, st_all$animal)

st_all <- SCTransform(st_all, assay = "Spatial", verbose = FALSE)
st_all <- st_all %>% RunPCA(assay = "SCT", verbose = FALSE) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(verbose = FALSE, res = 0.2) %>% 
  RunUMAP(reduction = "pca", dims = 1:30)

saveRDS(st_all, file = "Data/canis_merged_D0_D7_D14.RDS")

#### 4. Plot (UMAPs, Overlays, Heatmap) ####
pdf(file="Plots/Finalized/canis_all_umap_dimplot.pdf", onefile = TRUE)
DimPlot(st_all, reduction = "umap", group.by = c('timepoint'), label = TRUE, pt.size = 3) + NoLegend()
DimPlot(st_all, reduction = "umap", group.by = c('orig.ident'), label = TRUE, pt.size = 3) + NoLegend()
DimPlot(st_all, reduction = "umap", group.by = c('seurat_clusters'), label = TRUE, pt.size = 3) + NoLegend()
dev.off()

st_all <- FindClusters(st_all, verbose = FALSE, res = 0.3)
pdf("Plots/Finalized/canis_all_umap_clusters.pdf", width = 3, height = 3)
dittoDimPlot(st_all, "seurat_clusters", color.panel = c(hue_pal()(14)), do.label = TRUE) + NoLegend()
dev.off()

pdf("Plots/Finalized/canis_all_umap_timepoint.pdf", width = 3, height = 3)
dittoDimPlot(st_all, "timepoint", color.panel = c(hue_pal()(3)), do.label = TRUE) + NoLegend()
dev.off()

pdf("Plots/Finalized/Canine_clusters_highres_noLabel.pdf", width = 12, height = 5)
SpatialDimPlot(st_all, group.by = 'seurat_clusters', images = c('D0_Intact_R2', 'D7_Mid_1', 'D14_Mid_1'), pt.size.factor = 1, crop = FALSE, label = FALSE)
dev.off()

pdf("Plots/Finalized/Canine_clusters_highres_Label.pdf", width = 12, height = 5)
SpatialDimPlot(st_all, group.by = 'seurat_clusters', images = c('D0_Intact_R2', 'D7_Mid_1', 'D14_Mid_1'), pt.size.factor = 1, crop = FALSE, label = TRUE)
dev.off()

pdf(file="Plots/Finalized/canis_all_umap_qc_filtered.pdf", width = 10, height = 4)
p1 <- FeaturePlot(st_all, reduction = "umap", features = 'nCount_Spatial', pt.size = 1)
p2 <- FeaturePlot(st_all, reduction = "umap", features = 'nFeature_Spatial', pt.size = 1)
p1 + p2 + plot_layout(ncol = 2)
dev.off()

pdf("Plots/Finalized/canine_adgre1_cd68_overlay.pdf", onefile = TRUE, width = 12, height = 5)
SpatialFeaturePlot(st_all, features = c("ADGRE1", "CD68"), 
                   images = c('D0_Intact_R2', 'D7_Mid_1', 'D14_Mid_1'),
                   min.cutoff = 0, max.cutoff = 3, alpha = c(0,1),
                   pt.size.factor = 2, crop = FALSE)
dev.off()

pdf("Plots/canine_fibrotic_mac_markers_overlay.pdf", onefile = TRUE, width = 12, height = 6)
SpatialFeaturePlot(st_all, features = c('TREM2', 'SPP1'),  
                   images = c('D0_Intact_R2', 'D7_Mid_1', 'D14_Mid_1'),
                   min.cutoff = 0, max.cutoff = 4, alpha = c(0,1),
                   pt.size.factor = 2, crop = FALSE)
dev.off()

pdf("Plots/Finalized/canine_tgfb_overlays.pdf", onefile = TRUE, width = 12, height = 5)
SpatialFeaturePlot(st_all, features = c("TGFB1"),  
                   images = c('D0_Intact_R2', 'D7_Mid_1', 'D14_Mid_1'),
                   pt.size.factor = 1, crop = FALSE, min.cutoff = 0, max.cutoff = 2)
SpatialFeaturePlot(st_all, features = c("TGFBR1"),  
                   images = c('D0_Intact_R2', 'D7_Mid_1', 'D14_Mid_1'),
                   pt.size.factor = 1, crop = FALSE, min.cutoff = 0, max.cutoff = 2)
SpatialFeaturePlot(st_all, features = c("TGFBR2"),  
                   images = c('D0_Intact_R2', 'D7_Mid_1', 'D14_Mid_1'),
                   pt.size.factor = 1, crop = FALSE, min.cutoff = 0, max.cutoff = 2)
dev.off()

DefaultAssay(st_all) <- "SCT"
Idents(st_all) <- 'orig.ident'
degs <- FindAllMarkers(st_all, only.pos = TRUE)
top10 <- degs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf("Plots/Finalized/canine_deg_heatmap_tissue.pdf", width = 8, height = 6)
DoHeatmap(st_all, features = top10$gene)
dev.off()

#### 5. GSEA on Seurat Clusters ####
st_all <- FindClusters(st_all, res = 0.3)
DimPlot(st_all, group.by = c('seurat_clusters', 'orig.ident'))

# PATHWAY ANALYSIS ON CLUSTERS
Idents(st_all) <- 'seurat_clusters'
genes <- !grepl(pattern = "^RP[L|S]|MT", x = rownames(st_all))
genes <- !grepl(pattern = "^RP[L|S]|MT|ENSCAFG", x = rownames(st_all))
genes <- rownames(st_all)[genes]
de_markers_cluster <- FindAllMarkers(st_all, only.pos = T, features = genes)

top5 <- de_markers_cluster %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

pdf(file="Plots/Canis_Heatmap_Cluster_DEGs.pdf", width = 6, height = 8)
DoHeatmap(st_all, features = unique(top5$gene), slot = 'scale.data') +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(text = element_text(color = 'black'))
dev.off()

de_markers_cluster$entrez <- mapIds(org.Cf.eg.db, rownames(de_markers_cluster),'ENTREZID', 'SYMBOL')

## GOTerm analysis
go_cluster_comp_bp <- compareCluster(entrez ~ cluster, data = de_markers_cluster, fun = 'enrichGO', 
                                  OrgDb = "org.Cf.eg.db", ont = c("BP"), maxGSSize = 5000, 
                                  qvalueCutoff = 1)
go_cluster_comp_cc <- compareCluster(entrez ~ cluster, data = de_markers_cluster, fun = 'enrichGO', 
                                     OrgDb = "org.Cf.eg.db", ont = c("CC"), maxGSSize = 5000, 
                                     qvalueCutoff = 1)
go_cluster_comp <- compareCluster(entrez ~ cluster, data = de_markers_cluster, fun = 'enrichGO', 
                                     OrgDb = "org.Cf.eg.db", ont = c("ALL"), maxGSSize = 5000, 
                                     qvalueCutoff = 1)
write.table(go_cluster_comp@compareClusterResult, file = 'Canine_GOTerm_Cluster_Enrichment_Results.csv', sep = ',', append = FALSE)

pdf(file="Plots/Finalized/Canine_Cluster_enrichGO.pdf", width = 12, height = 4)
dotplot(go_cluster_comp_bp, showCategory = c('oxidative phosphorylation',
                                             'inflammatory response', 'immune response',
                                             'extracellular matrix organization',
                                             'macrophage derived foam cell differentiation',
                                             'coagulation')) +
dotplot(go_cluster_comp_cc, showCategory = c('myofibril', 'sarcomere',
                                             'phagocytic vesicle membrane', 
                                             'endocytic vesicle membrane',
                                             'collagen-containing extracellular matrix',
                                             'high-density lipoprotein particle',
                                             'endopeptidase complex', 'proteosome complex'
                                             ))
dev.off()

Idents(st_all) <- 'timepoint'
my_levels <- c('D0', 'D7', 'D14')
Idents(st_all) <- factor(Idents(st_all), levels= my_levels)

pdf('Plots/Canine_Violin_Myogenic_genes.pdf', width = 5, height = 5)
VlnPlot(st_all, features = c('MYOG', 'MYOD1', 'MYMK', 'MYH3', 'DES'))
dev.off()
