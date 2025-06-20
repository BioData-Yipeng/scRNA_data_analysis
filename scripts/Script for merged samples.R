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
setwd("C:/R_project/1_Alzheimer")

#creat seurat object
# Sample IDs
ad_samples <- c("AD02", "AD04", "AD12", "AD16", "AD17", "AD30")
hc_samples <- c("HC03", "HC07", "HC14", "HC19", "HC35", "HC37")
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

# Loop through each Seurat object in the list
for (sample_name in names(seurat_list)) {
  # Get the Seurat object
  seurat_obj <- seurat_list[[sample_name]]
  
  # Calculate percentage of mitochondrial genes
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
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
merged_seurat <- FindNeighbors(merged_seurat,dims = 1:15) #15dims contains most of the variations
merged_seurat <- FindClusters(merged_seurat, resolution = 0.1) #try various resolutions 0.1,0.2, 0.4, 0.8 
#run clustree to see distribution of clusters at difference resolution.
clustree(merged_seurat)

merged_seurat <- RunUMAP(merged_seurat, dims = 1:15)

#This is checking batch effects from different condition
DimPlot(merged_seurat, reduction = "umap", group.by="condition") + ggtitle("merged")

#This is checking batch effects from different samples
DimPlot(merged_seurat, reduction = "umap", group.by="orig.ident") + ggtitle("merged")

#Split HC and AD to check each clusters
DimPlot(merged_seurat, reduction = "umap", group.by="seurat_clusters", split.by="condition", label = TRUE) + ggtitle("merged")
gc()

#Check how many cells in each cluster
CellinEachCluster <- table(merged_seurat@meta.data$seurat_clusters)

#Check how many cells in each cluster under differnt condition
CellinEachClusterPerCondition <- table(Idents(merged_seurat), merged_seurat$condition)

#Check how many cells in each cluster from different samples
CellinEachClusterPerSample <- table(Idents(merged_seurat), merged_seurat$orig.ident)

#Need to join all layers to find different expression genes
merged_seurat <- JoinLayers(merged_seurat)

#Find different expression genes in each cluster vs all other cluster. cluster markers
de_cluster <- wilcoxauc(merged_seurat, 'seurat_clusters') # by presto method
#Find different expression genes between health control and AD
de_condition <- wilcoxauc(merged_seurat, 'condition')

#Filter top genes in each cluster
top_markers_cluster <- de_cluster %>%
  group_by(group) %>%
  filter(auc > 0.7, logFC > 1, pct_in > 85, pct_out < 45, padj < 0.05) %>%
  slice_max(order_by = auc, n = 20)

# Look into each cluster for gene markers and annotation
Idents(merged_seurat)
merged_seurat <- RenameIdents(merged_seurat,
                              '0'='oligodendrocytes',
                              '1'='astrocytes',
                              '2'='neuron-I',
                              '3'='neuron-II',
                              '4'='OPC',
                              '5'='microglia',
                              '6'='sensory neurons',
                              '7'='neuron-II',
                              '8'='unknown',
                              '9'='Interneurons')

merged_seurat$cellgroup <- paste0(Idents(merged_seurat),'_',merged_seurat$condition)

#Creat one column containing cluster annotation info
merged_seurat$celltype <- Idents(merged_seurat)

#Plot each genes in UMAP plot, change features with interesting genes
FeaturePlot(merged_seurat,
            features = c('MEGF11','VCAN'),
            split.by = 'condition',
            min.cutoff = 'q10')

Idents(merged_seurat)

#This is checking batch effects from different condition
DimPlot(merged_seurat, reduction = "umap", group.by="cellgroup", split.by = "condition", label = T) + ggtitle("merged")
DimPlot(merged_seurat, reduction = "umap", group.by="celltype", split.by = "condition", label = T) + ggtitle("merged")


#Write cell count info into files
write.csv(top_markers,"output/top50markersforeachgroup.csv")
write.csv(CellinEachCluster,"output/CellinEachCluster.csv")
write.csv(CellinEachClusterPerCondition,"output/CellinEachClusterPerCondition.csv")
write.csv(CellinEachClusterPerSample,"output/CellinEachClusterPerSample.csv")
