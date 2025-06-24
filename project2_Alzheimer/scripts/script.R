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

# Merge all Seurat objects in the list with add.cell.ids and project name
merged_seurat <- merge(
  x = seurat_list[[1]], 
  y = seurat_list[-1], 
  add.cell.ids = names(seurat_list),  # Add sample names as cell ID prefixes
  project = "Merged"  # Define the merged project name
)

rm(data, seurat_list, seurat_obj) # once merged, remove used data to save RAM

# Add a column in meta.data to group HC and AD
merged_seurat$condition <- ifelse(merged_seurat$orig.ident %in% tolower(hc_samples), "HC", "AD")

#Standard Analysis
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(object = merged_seurat))

#Base on elbowplot pick up number of pc used for followng steps
ElbowPlot(merged_seurat,ndims=50)
#Check variances for each dimision, help to determine how many pc used for following steps
DimHeatmap(merged_seurat, dims = 1:15, cells = 500, balanced = TRUE)

merged_seurat <- FindNeighbors(merged_seurat,dims = 1:15) #15dims contains most of the variations
merged_seurat <- FindClusters(merged_seurat, resolution = 0.2) #try various resolutions 0.1,0.2, 0.4, 0.8 
#run clustree to see distribution of clusters at difference resolution.
clustree(merged_seurat)

merged_seurat <- RunUMAP(merged_seurat, dims = 1:15)

#This is checking batch effects from different samples
DimPlot(merged_seurat, reduction = "umap", group.by="orig.ident") + ggtitle("merged")

#Because of the batch effects shown by the umap, considering integration methods

integrated_seurat <- IntegrateLayers(
  object = merged_seurat, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE
)

integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "harmony", dims = 1:15)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.2, cluster.name = "harmony_clusters")
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "harmony", dims = 1:15, reduction.name = "harmony")
DimPlot(
  integrated_seurat,
  reduction = "harmony",
  group.by = "orig.ident"
) + ggtitle("intergated_harmony")

DimPlot(
  integrated_seurat,
  reduction = "harmony",
  group.by = "harmony_clusters",
  split.by = "condition",
  label = TRUE,
) + ggtitle("intergated_harmony")

DimPlot(
  integrated_seurat,
  reduction = "harmony",
  group.by = "harmony_clusters",
  label = TRUE,
) + ggtitle("intergated_harmony")

View(integrated_seurat@meta.data)
#Check how many cells in each cluster
CellinEachCluster <- table(integrated_seurat@meta.data$celltype)
CellinEachCluster
#Check how many cells in each cluster under differnt condition
CellinEachClusterPerCondition <- table(Idents(integrated_seurat), integrated_seurat$condition)
CellinEachClusterPerCondition
#Check how many cells in each cluster from different samples
CellinEachClusterPerSample <- table(Idents(integrated_seurat), integrated_seurat$orig.ident)
CellinEachClusterPerSample
#Need to join all layers to find different expression genes
integrated_seurat <- JoinLayers(integrated_seurat)

#Find different expression genes in each cluster vs all other cluster. cluster markers
de_cluster <- wilcoxauc(integrated_seurat, 'seurat_clusters') # by presto method
#Find different expression genes between health control and AD
de_condition <- wilcoxauc(integrated_seurat, 'condition')

#Filter top genes in each cluster
top_markers_cluster <- de_cluster %>%
  group_by(group) %>%
  filter(auc > 0.7, logFC > 1, pct_in > 50, pct_out < 35, padj < 0.05) %>%
  slice_max(order_by = auc, n = 100)

cell_counts


# Look into each cluster for gene markers and annotation
# There are several ways to identify and annotated clusters, Deepseek and ChatGPT is helpful.
Idents(integrated_seurat)
integrated_seurat <- RenameIdents(integrated_seurat,
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

merged_seurat$cellgroup <- paste0(Idents(merged_seurat),'_',merged_seurat$condition)

#Creat one column containing cluster annotation info
integrated_seurat$celltype <- Idents(integrated_seurat)
DimPlot(
  integrated_seurat,
  reduction = "harmony",
  group.by = "celltype"
) + ggtitle("integrated_harmony")

DimPlot(
  integrated_seurat,
  reduction = "harmony",
  group.by = "celltype",
  split.by = "condition",
  label = TRUE,
) + ggtitle("intergated_harmony")



#Plot each genes in UMAP plot, change features with interesting genes
FeaturePlot(integrated_seurat,
            features = c('ST18'),
            split.by = 'condition',
            min.cutoff = 'q10')

#Explore interested gene expression among clusters
target_genes <- readLines("data/NominatedTargets.txt")
seurat_subset <- subset(integrated_seurat, features = target_genes)

# Ensure your Seurat object has been processed
seurat_subset <- NormalizeData(seurat_subset)
seurat_subset <- FindVariableFeatures(seurat_subset)
seurat_subset <- ScaleData(seurat_subset)

# Now plot the heatmap
DoHeatmap(seurat_subset, features = target_genes) #912X21171 matrix
DoHeatmap(seurat_subset, features = target_genes, cells = micro_cells,  group.by = "condition")
DoHeatmap(seurat_subset, features = target_genes, cells = astro_cells,  group.by = "condition")
astro_cells <- WhichCells(seurat_subset, idents = "Astrocytes")
#Write cell count info into files
write.csv(top_markers_cluster,"output/top20markersforeachgroup.csv")
write.csv(CellinEachCluster,"output/CellinEachCluster.csv")
write.csv(CellinEachClusterPerCondition,"output/CellinEachClusterPerCondition.csv")
write.csv(CellinEachClusterPerSample,"output/CellinEachClusterPerSample.csv")
