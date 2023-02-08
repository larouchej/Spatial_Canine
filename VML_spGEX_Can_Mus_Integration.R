# spGEX Canine Mouse Integration.R
# Jacqueline Larouche
# Fall 2022

# Integrate spGEX datasets across species (canine and mouse), analyze, and plot.

#### 0. Set up environment ####
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
library(corrplot)
library(amap)
library(hdf5r)
library(patchwork)
library(clusterProfiler)
library(tidyverse)
library(readxl)
library(ggpubr)
library(AnnotationDbi)
library(org.Cf.eg.db)
library(dittoSeq)
library(scales)

setwd("~/Documents/UMichigan/Aguilar_Lab/ResearchProjects/TA_VML/SpatialTranscriptomics/Visium/Data/spGEX")

#### 1. Initialize Variables ####
sample_dict <- list()
sample_dict[["PBS_599L"]] = "3765-JL"
sample_dict[["ITD1_599R"]] = "3765-JL"
sample_dict[["ITD1_600L"]] = "3765-JL"
sample_dict[["PBS_600R"]] = "3765-JL"
sample_dict[["PBS_1197L"]] = "5347-JL"
sample_dict[["PBS_1203L"]] = "5347-JL"
sample_dict[["ITD1_1203R"]] = "5347-JL"
sample_dict[["ITD1_1202L"]] = "5347-JL"
sample_dict[["Sample_5513-JL-S2-A_CGCGCACT-AGAATACA"]] = "D7_Mid_1"
sample_dict[["Sample_5513-JL-S2-B_CCTGTCAG-GTTACGGG"]] = "D7_Edge_1"
sample_dict[["Sample_5513-JL-S2-C_GTCCTTCG-CTGTGCAT"]] = "D14_Mid_1"
sample_dict[["Sample_5513-JL-S2-D_AATGTATC-TAAGCTCA"]] = "D14_Edge_1"
sample_dict[["Sample_6539-JL-S1-A_GTGGATCAAA-CAGGGTTGGC"]] = "D0_Intact_1"
sample_dict[["Sample_6833-JL-S1-A_TGTCCCAA-TGGACATC"]] = "D0_Intact_2"
sample_dict[["Sample_6834-JL-S1-C-GEX_TCCCAAGG-AAAGGTAG"]] = "D7_Mid_2"
sample_dict[["Sample_6834-JL-S1-D-GEX_GCCTTCGG-AAATCGTT"]] = "D14_Mid_2"
sample_dict[["Sample_6834-JL-S1-A-GEX_GCGGGTAA-CTTAGTGC"]] = "D14_Mouse_1"
sample_dict[["Sample_6834-JL-S1-B-GEX_CCTATCCT-GTTAGTAT"]] = "D14_Mouse_2"
samples <- names(sample_dict)

st_list <- readRDS("Data/Integrated_ST_list.RDS")
st_all <- readRDS("Data/Integrated_Merged_spGEX_SeuratObj_v3.RDS")

#### 2. Re-make canine list and convert to mouse orthologs, then merge with mouse list ####
gene.map <- read.table("Data/mus_can_mart_export.txt", sep = "\t", header = TRUE)

st_list_1 <- sapply(samples[13:16], USE.NAMES = TRUE, function(geo) {
  print(sample_dict[[geo]])
  st_se <- Load10X_Spatial(paste0(getwd(), sprintf("/10xOutputData/Canine/%s", geo)), 
                           assay = "Spatial", filter.matrix = TRUE, 
                           slice = sprintf("%s", sample_dict[[geo]]), 
                           filename = "filtered_feature_bc_matrix.h5")
  # convert to gene names
  symbol <- mapIds(org.Cf.eg.db, keys = st_se@assays$Spatial@counts@Dimnames[[1]],keytype = 'SYMBOL', 'ENSEMBL', multiVals = "first")
  # remove genes that don't map to symbols
  to.keep <- st_se@assays$Spatial@counts@Dimnames[[1]][!is.na(symbol)]
  st_se <- subset(st_se, features = to.keep)
  
  # convert to mouse gene names
  features.list <- st_se@assays$Spatial@counts@Dimnames[[1]]
  to.keep <- features.list %in% c(gene.map$Gene.stable.ID, gene.map$Gene.name)
  st_se <- subset(st_se, features = features.list[to.keep])
  features.list.2 <- st_se@assays$Spatial@counts@Dimnames[[1]]
  mapped.features <- features.list.2
  for (i in 1:length(features.list.2)){
    feature <- features.list.2[i]
    # if the feature is a gene name, map canine gene name to mouse gene name.
    if (feature %in% gene.map$Gene.name) {
      mapped.features[i] <- gene.map$Mouse.gene.name[gene.map$Gene.name == feature]
    }
    else {
      mapped.features[i] <- gene.map$Mouse.gene.name[gene.map$Gene.stable.ID == feature]
      }
    }

  st_se@assays$Spatial@counts@Dimnames[[1]] <- mapped.features
  st_se@assays$Spatial@data@Dimnames[[1]] <- mapped.features
  to.remove <- mapped.features == ''
  st_se <- subset(st_se, features = mapped.features[!to.remove])
  
  # Convert coordinates of the spatial image to integers from characters.
  st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["tissue"]] <- as.integer(st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["tissue"]])
  st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["row"]] <- as.integer(st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["row"]])
  st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["col"]] <- as.integer(st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["col"]])
  st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["imagerow"]] <- as.integer(st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["imagerow"]])
  st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["imagecol"]] <- as.integer(st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["imagecol"]])
  
  st_se$orig.ident <- sample_dict[[geo]]
  st_se$timepoint <- strsplit(sample_dict[[geo]], '_')[[1]][1]
  st_se$location <- strsplit(sample_dict[[geo]], '_')[[1]][2]
  st_se$animal <- strsplit(sample_dict[[geo]], '_')[[1]][3]
  st_se$species <- 'Canis'
  st_se$comp <- st_se$timepoint
  
  st_se <- SCTransform(st_se, assay = "Spatial", verbose = FALSE)
  return(list(st_se))
})

st_list_2 <- sapply(samples[9:12], USE.NAMES = TRUE, function(geo) {
  print(sample_dict[[geo]])
  st_se <- Load10X_Spatial(paste0(getwd(), sprintf("/10xOutputData/Canine/%s", geo)), 
                           assay = "Spatial", filter.matrix = TRUE, 
                           slice = sprintf("%s", sample_dict[[geo]]), 
                           filename = "filtered_feature_bc_matrix.h5")
  
  # convert to mouse gene names
  features.list <- st_se@assays$Spatial@counts@Dimnames[[1]]
  to.keep <- features.list %in% c(gene.map$Gene.stable.ID, gene.map$Gene.name)
  st_se <- subset(st_se, features = features.list[to.keep])
  features.list.2 <- st_se@assays$Spatial@counts@Dimnames[[1]]
  mapped.features <- features.list.2
  for (i in 1:length(features.list.2)){
    feature <- features.list.2[i]
    # if the feature is a gene name, map canine gene name to mouse gene name.
    if (feature %in% gene.map$Gene.name) {
      mapped.features[i] <- gene.map$Mouse.gene.name[gene.map$Gene.name == feature]
    }
    else {
      mapped.features[i] <- gene.map$Mouse.gene.name[gene.map$Gene.stable.ID == feature]
    }
  }
  
  st_se@assays$Spatial@counts@Dimnames[[1]] <- mapped.features
  st_se@assays$Spatial@data@Dimnames[[1]] <- mapped.features
  
  to.keep <- mapped.features != ''
  st_se <- subset(st_se, features = mapped.features[to.keep])
  
  st_se$orig.ident <- sample_dict[[geo]]
  st_se$timepoint <- strsplit(sample_dict[[geo]], '_')[[1]][1]
  st_se$location <- strsplit(sample_dict[[geo]], '_')[[1]][2]
  st_se$animal <- strsplit(sample_dict[[geo]], '_')[[1]][3]
  st_se$species <- 'Canis'
  st_se$comp <- st_se$timepoint
  
  st_se <- SCTransform(st_se, assay = "Spatial", verbose = FALSE)
  return(list(st_se))
})

st_list_mus <- readRDS("Data/Mus_ST_list_label_transfer_1201.RDS")

st_list <- c(st_list_1, st_list_2, st_list_mus)
saveRDS(st_list, "Data/Integrated_ST_list.RDS")

#### 3. Re-make merged seurat object for joint analysis ####
# Use only the PBS treated tissues
st_1 <- st_list[[samples[1]]]; st_4 <- st_list[[samples[4]]]; 
st_5 <- st_list[[samples[5]]]; st_6 <- st_list[[samples[6]]];
st_9 <- st_list[[samples[9]]]; st_10 <- st_list[[samples[10]]];
st_11 <- st_list[[samples[11]]]; st_12 <- st_list[[samples[12]]];
st_13 <- st_list[[samples[13]]]; st_14 <- st_list[[samples[14]]];
st_15 <- st_list[[samples[15]]]; st_16 <- st_list[[samples[16]]];
st_17 <- st_list[[samples[17]]]; st_18 <- st_list[[samples[18]]];

st_mus <- merge(st_1, c(st_4, st_5, st_6, st_17, st_18)) %>% SCTransform(assay = 'Spatial')
names(st_mus@images)
st_mus$location <- st_mus$Zone
# D0_R2 removed because low quality data
st_canis <- merge(st_9, c(st_10, st_11, st_12, st_13, st_14, st_15, st_16)) %>% SCTransform(assay = 'Spatial')
names(st_canis@images)
st_all <- merge(st_mus, st_canis)
names(st_all@images)
table(st_all$location)
# Formatting
Idents(st_all) <- 'location'
current.cluster.ids <- c('', 'Defect', 'Intact', 'IntactMuscle', 'Transition', "Edge", "Mid")
corrected.cluster.ids <- c('', 'Defect', 'Intact', 'Intact', 'Transition', 'Edge', 'Mid')
st_all$location <- plyr::mapvalues(x = st_all$location, from = current.cluster.ids, to = corrected.cluster.ids)
st_all <- subset(st_all, idents = '', invert = TRUE)

#### 4. Integrate datasets by species using Seurat ####
# The below is a workaround for the duplicate images that are created during integration
species.list <- SplitObject(st_all, split.by = "species")
species.list[[1]] <- SCTransform(species.list[[1]], assay = 'Spatial')
species.list[[2]] <- SCTransform(species.list[[2]], assay = 'Spatial')
features <- SelectIntegrationFeatures(object.list = species.list, nfeatures = 3000)
species.list <- PrepSCTIntegration(object.list = species.list, anchor.features = features)
species.anchors <- FindIntegrationAnchors(object.list = species.list, 
                                         normalization.method = "SCT",
                                         anchor.features = features)
species.combined.sct <- IntegrateData(anchorset = species.anchors, 
                                      normalization.method = "SCT")
species.combined.sct <- RunPCA(species.combined.sct, verbose = FALSE) %>%
 FindNeighbors(dims = 1:30) %>%
 FindClusters(resolution = 0.1) %>%
 RunUMAP(reduction = "pca", dims = 1:30) 
DimPlot(species.combined.sct, group.by = c('species', 'comp', 'seurat_clusters'))
names(species.combined.sct@images)
# species.combined.sct@images <- species.combined.sct@images[-c(5:14)]

# Transfer integrated assay and metadata to st_all
st_all@assays$integrated <- GetAssay(species.combined.sct, assay = 'integrated')
DefaultAssay(st_all) <- 'integrated'
st_all <- RunPCA(st_all)
st_all <- FindNeighbors(st_all, dims = 1:30)
st_all <- FindClusters(st_all, resolution = 0.1)
st_all <- RunUMAP(st_all, reduction = "pca", dims = 1:30) 

saveRDS(st_all, "Data/Integrated_Merged_spGEX_SeuratObj_v3.RDS")

#### 5. Plot ####
pdf("Plots/Integrated/UMAP_species_comp_clusters.pdf", onefile = TRUE, width = 20, height = 5)
DimPlot(st_all, group.by = c('species', 'comp', 'seurat_clusters'))
dev.off()

pdf("Plots/Integrated/UMAP_clusters.pdf", onefile = TRUE, width = 2.5, height = 2.5)
dittoDimPlot(st_all, "seurat_clusters", color.panel = c(hue_pal()(7)), do.label = TRUE) + NoLegend()
dev.off()

pdf("Plots/Integrated/UMAP_comp_origident_splitby_species.pdf", onefile = TRUE, width = 5.5, height = 2.5)
dittoDimPlot(st_all, "comp", split.by = c("species"), color.panel = c(hue_pal()(9)), do.label = FALSE)
dittoDimPlot(st_all, "orig.ident", split.by = c("species"), color.panel = c(hue_pal()(14)))
dev.off()

pdf("Plots/Integrated/StackedBar_timepoint_location_clusters.pdf", width = 6, height = 3)
p1 <- dittoBarPlot(st_all, group.by = "timepoint", "seurat_clusters", split.by = 'species',
             color.panel = c(hue_pal()(7)), x.reorder = c(1, 3, 2), 
             cells.use = meta("species", st_all) == "Canis")
p2 <- dittoBarPlot(st_all, group.by = "comp", "seurat_clusters", split.by = 'species',
             color.panel = c(hue_pal()(7)), x.reorder = c(4,6,5,1,3,2),
             cells.use = meta("species", st_all) == "Mus")
p1 + p2 + plot_layout(ncol = 2)
dev.off()

st_all <- FindClusters(st_all, resolution = 0.1)
DimPlot(st_all, group.by = c('seurat_clusters', 'orig.ident'))

cluster.degs <- FindAllMarkers(st_all, only.pos = TRUE)
top10 <- cluster.degs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

Idents(st_canis) <- 'timepoint'
timepoint.degs <- FindAllMarkers(st_canis, only.pos = TRUE)
top10 <- timepoint.degs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
Idents(st_all) <- 'comp'
my_levels <- c('D7_IntactMuscle', 'D14_IntactMuscle', 'D0', 'D7_Transition', 'D14_Transition', 'D7', 'D7_Defect', 'D14_Defect', 'D14')
st_all@active.ident <- factor(x = st_all@active.ident, levels = my_levels)
col_pal <- c(hue_pal()(9))

pdf("Plots/Integrated/Heatmap_comp_degs.pdf", width = 7, height = 7)
DoHeatmap(st_all, features = top10$gene[c(21:30,1:20)],
          group.colors = col_pal[c(8,4,1,9,5,6,7,3,2)]) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(text = element_text(color = 'black'))
dev.off()

DefaultAssay(st_all) <- 'SCT'

pdf("Plots/Integrated/SpatialFeature_Myog_Aspn_Myh1_Adgre1.pdf", onefile = TRUE, width = 2.5, height = 12)
SpatialFeaturePlot(st_all, features = c('Myh1', 'Myog', 'Aspn', 'Adgre1'), images = c('PBS_1197L'), image.alpha = 1, ncol = 1, pt.size.factor = 2.5)
SpatialFeaturePlot(st_all, features = c('Myh1', 'Myog', 'Aspn', 'Adgre1'), images = c('D14_Mouse_1'), image.alpha = 1, ncol = 1, pt.size.factor = 2.5)
SpatialFeaturePlot(st_all, features = c('Myh1', 'Myog', 'Aspn', 'Adgre1'), images = c('D0_Intact_2'), image.alpha = 1, ncol = 1)
SpatialFeaturePlot(st_all, features = c('Myh1', 'Myog', 'Aspn', 'Adgre1'), images = c('D7_Mid_1'), image.alpha = 1, ncol = 1)
SpatialFeaturePlot(st_all, features = c('Myh1', 'Myog', 'Aspn', 'Adgre1'), images = c('D14_Mid_1'), image.alpha = 1, ncol = 1)
dev.off()

pdf("Plots/Integrated/SpatialFeature_Myh4_Tnnt3_Pax7_Myod1_Acta2_Col1a1_Csf2rb_S100a8_d0.pdf", onefile = TRUE, width = 18, height = 3)
SpatialFeaturePlot(st_all, features = c('Tnni2', 'Tnnt3', 'Pax7', 'Myod1', 'Acta2', 'Col1a1', 'Csf2rb', 'S100a8'), images = c('D0_Intact_2'), image.alpha = 1, ncol = 8)
dev.off()
pdf("Plots/Integrated/SpatialFeature_Myh4_Tnnt3_Pax7_Myod1_Acta2_Col1a1_Csf2rb_S100a8_d7.pdf", onefile = TRUE, width = 18, height = 3)
SpatialFeaturePlot(st_all, features = c('Tnni2', 'Tnnt3', 'Pax7', 'Myod1', 'Acta2', 'Col1a1', 'Csf2rb', 'S100a8'), images = c('D7_Mid_1'), image.alpha = 1, ncol = 8)
dev.off()
pdf("Plots/Integrated/SpatialFeature_Myh4_Tnnt3_Pax7_Myod1_Acta2_Col1a1_Csf2rb_S100a8_d14.pdf", onefile = TRUE, width = 18, height = 3)
SpatialFeaturePlot(st_all, features = c('Tnni2', 'Tnnt3', 'Pax7', 'Myod1', 'Acta2', 'Col1a1', 'Csf2rb', 'S100a8'), images = c('D14_Mid_1'), image.alpha = 1, ncol = 8)
dev.off()

