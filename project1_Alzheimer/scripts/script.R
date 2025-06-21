library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(presto)
library(tidyverse)
library(multtest)
library(metap)


install.packages("presto")
devtools::install_github("immunogenomics/presto")

getwd()
setwd("C:/R_project/1_Alzheimer")

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

merged_seurat$condition <- ifelse(merged_seurat$orig.ident %in% tolower(hc_samples), "HC", "AD")

#Standard Analysis
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(object = merged_seurat))

ElbowPlot(merged_seurat,ndims=50)

merged_seurat <- FindNeighbors(merged_seurat,dims = 1:10)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.1)
merged_seurat <- RunUMAP(merged_seurat, dims = 1:10)

DimPlot(merged_seurat, reduction = "umap", group.by="condition") + ggtitle("merged")
DimPlot(merged_seurat, reduction = "umap", group.by="orig.ident") + ggtitle("merged")
DimPlot(merged_seurat, reduction = "umap", group.by="seurat_clusters", split.by="condition", label = TRUE) + ggtitle("merged")
gc()

integrated_seurat <- IntegrateLayers(object = merged_seurat,
                                   method = HarmonyIntegration,
                                   orig.reduction = "pca",
                                   new.reduction = "harmony",
                                   verbose=TRUE)

integrated_seurat <- FindNeighbors(integrated_seurat,
                                 reduction = "harmony",
                                 dims = 1:10)
integrated_seurat <- FindClusters(integrated_seurat,
                                resultion=0.1,
                                cluster.name = "harmonycluster")
integrated_seurat <- RunUMAP(integrated_seurat,
                           reduction = "harmony",
                           dims = 1:10,
                           reduction.name = "harmonyUMAP")
DimPlot(integrated_seurat,
        reduction = "harmonyUMAP",
        group.by = "condition",
        ) + ggtitle("integrated-harmony")
DimPlot(integrated_seurat,
        reduction = "harmonyUMAP",
        group.by = "orig.ident",
) + ggtitle("integrated-harmony")
DimPlot(integrated_seurat,
        reduction = "harmonyUMAP",
        group.by = "harmonycluster",
        split.by = "condition",
        label = TRUE,
) + ggtitle("integrated-harmony")


table(merged_seurat@meta.data$seurat_clusters)
table(integrated_seurat@meta.data$seurat_clusters)

table(Idents(merged_seurat), merged_seurat$condition)

table(Idents(integrated_seurat), integrated_seurat$condition)

library(clustree)
clustree(merged_seurat)
install.packages("clustree")



# 2nd intergration method
integrated_seurat <- IntegrateLayers(object = merged_seurat,
                                     method = CCAIntegration,
                                     orig.reduction = "pca",
                                     new.reduction = "CCA",
                                     verbose=TRUE)

integrated_seurat <- FindNeighbors(integrated_seurat,
                                   reduction = "CCA",
                                   dims = 1:10)
integrated_seurat <- FindClusters(integrated_seurat,
                                  resultion=0.5,
                                  cluster.name = "CCAcluster")
integrated_seurat <- RunUMAP(integrated_seurat,
                             reduction = "CCA",
                             dims = 1:10,
                             reduction.name = "CCAUMAP")
DimPlot(integrated_seurat,
        reduction = "CCAUMAP",
        group.by = "condition",
) + ggtitle("integrated-CCA")

DimPlot(integrated_seurat,
        reduction = "CCAUMAP",
        group.by = "CCAcluster",
        split.by = "condition",
        label = TRUE,
) + ggtitle("integrated-CCA")

# 3rd intergration method
integrated_seurat <- IntegrateLayers(object = merged_seurat,
                                     method = RPCAIntegration,
                                     orig.reduction = "pca",
                                     new.reduction = "RPCA",
                                     verbose=TRUE)

integrated_seurat <- FindNeighbors(integrated_seurat,
                                   reduction = "RPCA",
                                   dims = 1:10)
integrated_seurat <- FindClusters(integrated_seurat,
                                  resultion=0.5,
                                  cluster.name = "RPCAcluster")
integrated_seurat <- RunUMAP(integrated_seurat,
                             reduction = "RPCA",
                             dims = 1:10,
                             reduction.name = "RPCAUMAP")
DimPlot(integrated_seurat,
        reduction = "RPCAUMAP",
        group.by = "condition",
) + ggtitle("integrated-RPCA")

DimPlot(integrated_seurat,
        reduction = "RPCAUMAP",
        group.by = "RPCAcluster",
        split.by = "condition",
        label = TRUE,
) + ggtitle("integrated-RPCA")

integrated_seurat <- JoinLayers(integrated_seurat)
merged_seurat <- JoinLayers(merged_seurat)
de <- wilcoxauc(integrated_seurat, 'harmonycluster') # by presto method
de_merged <- wilcoxauc(merged_seurat, 'seuratcluster') # by presto method


top_markers <- de %>%
  group_by(group) %>%
  filter(auc > 0.7, logFC > 1, pct_in > 85, pct_out < 45, padj < 0.05) %>%
  slice_max(order_by = auc, n = 50)
markers_cluster2 <- FindConservedMarkers(integrated_seurat,
                     ident.1 = 2,
                     grouping.var = 'condition')
head(markers_cluster2)
FeaturePlot(integrated_seurat,reduction = "harmonyUMAP", features = c("GFAP"),min.cutoff = 'q10')
integrated_seurat$cellgroup <- paste0(integrated_seurat$harmonycluster,'_',integrated_seurat$condition)
View(integrated_seurat@meta.data)

Idents(integrated_seurat) <- integrated_seurat$cellgroup
hcad_cluster2 <- FindMarkers(integrated_seurat,ident.1 = '2_AD',ident.2 = '2_HC')
head(hcad_cluster2)

Idents(integrated_seurat)
integrated_seurat <- RenameIdents(integrated_seurat, '1_HC'='Astrocyte_HC')


FeaturePlot(integrated_seurat, 
            reduction = "harmonyUMAP",
            features = c('GFAP','GRM3'),
            split.by = 'condition',
            min.cutoff = 'q10')


install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')
rm(de)
gc()
View(de)
write.csv(top_markers,"output/top50markersforeachgroup.csv")
#method = CCAIntegration, method = HarmonyIntegration, method = RPCAIntegration, method = FastMNNIntegration, method = scVIIntegration
