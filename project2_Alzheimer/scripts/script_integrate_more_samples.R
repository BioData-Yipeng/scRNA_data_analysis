## This script is used to integrate more samples to avoid batch effects.
## After integration, different expressed genes were calculated among clusters and between HC and AD.
## Clusters were annotated based on DE genes. 
## AD associated genes were explored between HC and AD in all clusters or subclusters.


#Load library
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(presto)
library(tidyverse)
library(multtest)
library(metap)
library(clustree)

#set working dirtory
getwd()
setwd("C:/R_project/project2_Alzheimer")

#creat seurat object
# Sample IDs
ad_samples <- c("AD_A1", "AD_A3")
hc_samples <- c("HC_A2", "HC_A4")
all_samples <- c(ad_samples, hc_samples)
# Initialize list to store Seurat objects
seurat_list <- list()
# Loop through each sample
for (sample in all_samples) {
  data_path <- paste0("data/", sample, "/")
  data <- Read10X(data.dir = data_path)
  seurat_obj <- CreateSeuratObject(counts = data, project = tolower(sample), 
                                   min.cells = 5, min.features = 200)
  seurat_list[[tolower(sample)]] <- seurat_obj
}

# optional: Loop through each Seurat object in the list and qc plot save for further explore
# vlnplot and featureplot are used to QC.
vlnplot_list <-list()
featureplot_list <- list()

for (sample_name in names(seurat_list)) {
  # Get the Seurat object
  seurat_obj <- seurat_list[[sample_name]]
  
  # Calculate percentage of mitochondrial genes
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  
  vlnplot_list[[sample_name]] <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  featureplot_list[[sample_name]] <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  # Subset the Seurat object based on criteria
  seurat_list[[sample_name]] <- subset(
    seurat_obj, 
    subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5
  )
}

# Check the number of cells in each subset Seurat object
cell_counts <- sapply(seurat_list, ncol)

# Merge all Seurat objects in the list with add.cell.ids and project name to faciliate the following steps.
seurat <- merge(
  x = seurat_list[[1]], 
  y = seurat_list[-1], 
  add.cell.ids = names(seurat_list),  # Add sample names as cell ID prefixes
  project = "Merged"  # Define the merged project name
)

rm(data, seurat_list, seurat_obj) # once merged, remove used data to save RAM

# Add a column in meta.data to group HC and AD
seurat$condition <- ifelse(seurat$orig.ident %in% tolower(hc_samples), "HC", "AD")

#Standard Analysis
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))

#Base on elbowplot pick up number of pc used for followng steps
ElbowPlot(seurat,ndims=50)
#Check variances for each dimision, help to determine how many pc used for following steps
DimHeatmap(seurat, dims = 1:9, cells = 500, balanced = TRUE)

seurat <- FindNeighbors(seurat,dims = 1:15) #15dims contains most of the variations
seurat <- FindClusters(seurat, resolution = 0.2) #try various resolutions 0.1,0.2, 0.4, 0.8 
#run clustree to see distribution of clusters at difference resolution.
clustree(merged_seurat)

#Check the UMAP to determine whether need integration or not
seurat <- RunUMAP(seurat, dims = 1:15)

#This is checking batch effects from different samples
p1 <- DimPlot(seurat, reduction = "umap", group.by="orig.ident") + ggtitle("unintegrated_UMAP")

#Because of the batch effects shown by the umap, considering integration methods
seurat <- IntegrateLayers(
  object = seurat, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE
)

seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:15)
seurat <- FindClusters(seurat, resolution = 0.2, cluster.name = "harmony_clusters")
seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:15, reduction.name = "harmony")

#draw UMAP after integration
p2 <- DimPlot(
  seurat,
  reduction = "harmony",
  group.by = "orig.ident"
) + ggtitle("intergated_UMAP")

p1+p2 #compare umap before and after integration

#Need to join all layers to find different expression genes
seurat <- JoinLayers(seurat)

View(seurat@meta.data)
#Find different expression genes in each cluster vs all other cluster. cluster markers
de_cluster <- wilcoxauc(seurat, 'seurat_clusters') # by presto method
#Find different expression genes between health control and AD
de_condition <- wilcoxauc(seurat, 'condition')

#Filter top genes in each cluster
top_markers_cluster <- de_cluster %>%
  group_by(group) %>%
  filter(auc > 0.7, logFC > 1, pct_in > 50, pct_out < 35, padj < 0.05) %>%
  slice_max(order_by = auc, n = 40)

#based on the top marker information finder out annotation for each cluster through SingleR or other resources.
seurat <- RenameIdents(seurat,
                       '0'='oligodendrocytes',
                       '1'='excitatory neurons',
                       '2'='inhibitory neurons',
                       '3'='Microglia',
                       '4'='Astrocytes',
                       '5'='Fibroblasts',
                       '6'='striatal medium spiny neurons',
                       '7'='vascular endothelial cells',
                       '8'='hippocampal neurons',
                       '9'='cortical layer-specific neurons',
                       '10'='neuronal-glial hybrid',
                       '11'='cerebellar neurons',
                       '12'='perivascular macrophages')
Idents(seurat)
#Creat one column containing cluster annotation info
seurat$celltype <- Idents(seurat)


DimPlot(
  seurat,
  reduction = "harmony",
  group.by = "celltype"
) + ggtitle("cluster annotation")

DimPlot(
  seurat,
  reduction = "harmony",
  group.by = "celltype",
  split.by = "condition",
) + ggtitle("compare each cluster in HC and AD")

#example markers for clusters
oligodendro <- FeaturePlot(seurat,
                           reduction = "harmony",
                           features = c('PLP1','MOBP'),
                           split.by = 'condition',
                           min.cutoff = 'q10')

excit_neuron <- FeaturePlot(seurat,
                            reduction = "harmony",
                            features = c('ADGRV1','RYR3'),
                            split.by = 'condition',
                            min.cutoff = 'q10')
microglia <- FeaturePlot(seurat,
                            reduction = "harmony",
                            features = c('INPP5D','CSF1R'),
                            split.by = 'condition',
                            min.cutoff = 'q10')
astrocytes <- FeaturePlot(seurat,
                            reduction = "harmony",
                            features = c('VCAN','PTPRZ1'),
                            split.by = 'condition',
                            min.cutoff = 'q10')

#Filter top genes in HC and AD condition
top_markers_condition <- de_condition %>%
  group_by(group) %>%
  filter(auc > 0.6, logFC > 0.7, padj < 0.05) %>%
  slice_max(order_by = auc, n = 40)

#example markers for AD
DE1 <- FeaturePlot(seurat,
                          reduction = "harmony",
                          features = c('RBFOX1','SYT1'),
                          split.by = 'condition',
                          min.cutoff = 'q10')
DE2 <- FeaturePlot(seurat,
                   reduction = "harmony",
                   features = c('NEAT1','FKBP5'),
                   split.by = 'condition',
                   min.cutoff = 'q10')
DE1
DE2

#analysis different gene expression AD vs HC in one subcluster like astrocytes
# For example, subset cells labeled "Microglia"
microglia_seurat <- subset(
  seurat,
  idents = "Microglia"
)

#find DE genes in microglia
de_condition_migroglia <- wilcoxauc(microglia_seurat, 'condition')

DE1_microglia <- FeaturePlot(seurat,
                   reduction = "harmony",
                   features = c('XIST','MAP1B'),
                   split.by = 'condition',
                   min.cutoff = 'q10')
VlnPlot(seurat, group.by = "celltype", features = 'MAP1B', split.by = "condition")
VlnPlot(seurat, group.by = "celltype", features = 'RBFOX1', split.by = "condition")
VlnPlot(seurat, group.by = "celltype", features = 'NEAT1', split.by = "condition") # normalized data as y-axis
VlnPlot(seurat, group.by = "celltype", features = 'NEAT1', split.by = "condition", slot="count") #raw count as y-axis

#Check how many cells in each cluster
CellinEachCluster <- table(seurat@meta.data$celltype)
CellinEachCluster
#Check how many cells in each cluster under differnt condition
CellinEachClusterPerCondition <- table(Idents(seurat), seurat$condition)
CellinEachClusterPerCondition
#Check how many cells in each cluster from different samples
CellinEachClusterPerSample <- table(Idents(seurat), seurat$orig.ident)
CellinEachClusterPerSample

#Write cell count info into files
write.csv(top_markers_cluster,"output/top20markersforeachgroup.csv")
write.csv(CellinEachCluster,"output/CellinEachCluster.csv")
write.csv(CellinEachClusterPerCondition,"output/CellinEachClusterPerCondition.csv")
write.csv(CellinEachClusterPerSample,"output/CellinEachClusterPerSample.csv")
write.csv(cell_counts,"output/cellcountpersample.csv")
