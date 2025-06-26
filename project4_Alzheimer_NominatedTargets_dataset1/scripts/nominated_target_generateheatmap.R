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
library(pheatmap) 

#set working dirtory
getwd()
setwd("C:/R_project/project4_Alzheimer_NominatedTargets_dataset1")

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

rm(data, seurat_list, seurat_obj)

gc()
# Add a column in meta.data to group HC and AD
merged_seurat$condition <- ifelse(merged_seurat$orig.ident %in% tolower(hc_samples), "HC", "AD")

merged_seurat <- JoinLayers(merged_seurat)

norminated_genes <- readLines("data/NominatedTargets.txt")

target_seurat <- subset(merged_seurat, features=norminated_genes) # 912 over 948 target genes fetched from seurat dataset.

#Standard Analysis
target_seurat <- NormalizeData(target_seurat)
target_seurat <- FindVariableFeatures(target_seurat)

VariableFeaturePlot(target_seurat)

target_seurat <- ScaleData(target_seurat)
# or ComplexHeatmap for more advanced options

View(target_seurat@meta.data)


# Extract scaled data (assuming RNA assay)
scaled_data <- GetAssayData(target_seurat, slot = "scale.data")
dim(scaled_data)

avg_exp <- AverageExpression(target_seurat,
                             assays = "RNA",
                             slot = "scale.data",
                             group.by = "condition")$RNA
write.csv(avg_exp, file = "output/average_expression_by_group.csv")

#generally check the average expression between HC and AD
hm <- pheatmap(avg_exp,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,    # Hide row names for clarity
         show_colnames = TRUE,    # Hide column names
         scale = "none",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Average Scaled Expression by Group")
#extract clustering assignments
row_cluster <- cutree(hm$tree_row, k=4)
row_cluster
#retrieve genes in each cluster
genes_by_cluster <- split(names(row_cluster), row_cluster)

cluster1_genes <- genes_by_cluster[[1]]
cluster2_genes <- genes_by_cluster[[2]]
cluster3_genes <- genes_by_cluster[[3]]
cluster4_genes <- genes_by_cluster[[4]]

cluster2_seurat <- subset(target_seurat, features=cluster2_genes)
cluster3_seurat <- subset(target_seurat, features=cluster3_genes)
cluster4_seurat <- subset(target_seurat, features=cluster4_genes)

DoHeatmap(
  cluster2_seurat,              
  group.by = "condition",  
  slot = "scale.data",
  raster = TRUE,
  draw.lines = TRUE
) + scale_y_discrete(
  breaks = function(x) x[seq(1, length(x), by = 5)]
)


#an alternative way to draw heatmap of specific genes or cells
DoHeatmap(
  target_seurat,
  features = cluster4_genes,
  group.by = "condition",  
  slot = "scale.data",
  raster = TRUE,
  draw.lines = TRUE
)


write.csv(cluster1_genes,"output/cluster1_genes.csv")
write.csv(cluster2_genes,"output/cluster2_genes.csv")
write.csv(cluster3_genes,"output/cluster3_genes.csv")
write.csv(cluster4_genes,"output/cluster4_genes.csv")
