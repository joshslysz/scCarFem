#### Carotid Fem Pipe #### 
# Check for missing R packages, install if needed
list_of_pkgs <- c("hdf5r", "tidyverse", "ggpubr", "Seurat", "enrichR",
                  "cowplot", "patchwork", "VennDiagram", "devtools",
                  "RColorBrewer", "pheatmap", "msigdbr")
install.packages(list_of_pkgs[! list_of_pkgs %in% rownames(installed.packages())])
# Check for missing BioC packages, install if needed
list_of_pkgs <- c("limma", "org.Hs.eg.db", "SingleR", "batchelor")
bioc2install <- list_of_pkgs[! list_of_pkgs %in% rownames(installed.packages())]
if (length(bioc2install) > 0) {
  BiocManager::install(bioc2install)
}

# Install SeuratWrappers (optional)
devtools::install_github('satijalab/seurat-wrappers')

# Install monocle3 (optional)
devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")
devtools::install_github("mw201608/msigdb")
devtools::install_github("romanhaa/cerebroApp")
install.packages("UpSetR")
install.packages("SCINA")
install.packages("harmony")


## Load libraries ##
library(Seurat)
library(tidyverse)
library(ggpubr)
library(limma)
library(SingleR)
library(enrichR)
library(cowplot)
library(patchwork)
library(VennDiagram)
library(pheatmap)
library(RColorBrewer)
library(monocle3)
library(msigdbr)
library(scmap)
library(scater)
library(msigdb)
library(SCINA)
library(devtools)
library(cerebroApp)
library(shiny)
library(SeuratWrappers)
library(UpSetR)
library(harmony)

# Define color palettes
WtOrRd_pal <- colorRampPalette(c("#FFFFFF", brewer.pal(9, "OrRd")))(100)

# Cleanup the environment to free up memory
rm(list=ls()); gc() 

## Get the data ##
fem1.data <- Read10X(data.dir = "femoral_outs")
car1.data <- Read10X(data.dir = "carotid_outs")

# Initialize the Femoral Seurat object with the raw (non-normalized data).
srtfemoral<- CreateSeuratObject(counts = fem1.data, project = "femoral")
srtfemoral

#Quality control on Femoral Seurat Object

# Add in the Mitochondrial PCT% information
srtfemoral$percent.mt <- PercentageFeatureSet(srtfemoral, pattern = "^MT-")

# nCount_RNA is the number of UMI counts in a cell
hist(srtfemoral$nCount_RNA)

# nFeature_RNA is the number of different genes that had any reads
hist(srtfemoral$nFeature_RNA)

# percent.mt is the percent mitochondrial reads
hist(srtfemoral$percent.mt)

# Make a violin plot of the QC columns
plt <- VlnPlot(srtfemoral, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               ncol = 3)
ggsave(filename = "QC_Femoral.png", plot = plt, width = 7, height = 3.5)

# Filter out unwanted cells from srtfemoral
srtfemoral <- subset(srtfemoral, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mt < 10)


# Initialize the Carotid Seurat object with the raw (non-normalized data).
srtcarotid <- CreateSeuratObject(counts = car1.data, project = "carotid")
srtcarotid

#Quality control on Carotid Seurat Object
# Add in the Mitochondrial PCT% information
srtcarotid$percent.mt <- PercentageFeatureSet(srtcarotid, pattern = "^MT-")

# nCount_RNA is the number of UMI counts in a cell
hist(srtcarotid$nCount_RNA)

# nFeature_RNA is the number of different genes that had any reads
hist(srtcarotid$nFeature_RNA)

# percent.mt is the percent mitochondrial reads
hist(srtcarotid$percent.mt)

# Make a violin plot of the QC columns
plt <- VlnPlot(srtcarotid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               ncol = 3)
ggsave(filename = "QC_Carotid.png", plot = plt, width = 7, height = 3.5)

# Filter out unwanted cells from srtcarotid
srtcarotid <- subset(srtcarotid, subset = nFeature_RNA > 200 & nFeature_RNA < 25000 & percent.mt < 10)

### END QUALIT CONTROL ON SEPARATE SAMPLES ###

### COMBINE SAMPLES ###

# combine carotid and femoral seurat objects
Carotid_Femoral_srt <- merge(x= srtfemoral, y= srtcarotid)

## Normalize and Scale the data ##
# Log-transform the counts
Carotid_Femoral_srt <- NormalizeData(Carotid_Femoral_srt)
# Find Variable Features
Carotid_Femoral_srt <- FindVariableFeatures(Carotid_Femoral_srt)
# Scale the data
Carotid_Femoral_srt <- ScaleData(Carotid_Femoral_srt)
# Run PCA
Carotid_Femoral_srt <- RunPCA(Carotid_Femoral_srt)
# Plot the PCA
DimPlot(Carotid_Femoral_srt, reduction = "pca", group.by = "orig.ident")
# Plot the loadings
VizDimLoadings(Carotid_Femoral_srt, dims = 1:2, reduction = "pca")
# Choose the number of principle components to keep
ElbowPlot(Carotid_Femoral_srt,ndims = 50)
# Find nearest neighbors and construct the graph
Carotid_Femoral_srt <- FindNeighbors(Carotid_Femoral_srt, k.param = 50, dims = 1:40)
# Find the clusters
Carotid_Femoral_srt <- FindClusters(Carotid_Femoral_srt, resolution = 0.35)
# Get the UMAP embedding
Carotid_Femoral_srt <- RunUMAP(Carotid_Femoral_srt, dims = 1:40)
# Plot the UMAP with clustering
DimPlot(Carotid_Femoral_srt, reduction = "umap", label = TRUE)
# Dim Plots
DimPlot(Carotid_Femoral_srt, reduction = "umap", group.by = "orig.ident") 
DimPlot(Carotid_Femoral_srt, group.by = "orig.ident") + DimPlot(Carotid_Femoral_srt, group.by = "seurat_clusters") 

# Compare Batch and Cluster ID
compTable <- table(Carotid_Femoral_srt$orig.ident, Carotid_Femoral_srt$seurat_clusters)
compTable <- (compTable / rowSums(compTable)) * 100
pheatmap(compTable, color = WtOrRd_pal)

#####INTEGRATE WITH HARMONY##################
Carotid_Femoral_srt.int.har <- RunHarmony(Carotid_Femoral_srt, group.by.vars = "orig.ident")
Carotid_Femoral_srt.int.har <- RunUMAP(Carotid_Femoral_srt.int.har, reduction = "harmony", dims = 1:20)
DimPlot(Carotid_Femoral_srt.int.har, reduction = "umap", group.by = "orig.ident")

# Find nearest neighbors and construct the graph
Carotid_Femoral_srt.int.har <- FindNeighbors(Carotid_Femoral_srt.int.har, dims = 1:20)
# Find the clusters
Carotid_Femoral_srt.int.har <- FindClusters(Carotid_Femoral_srt.int.har, resolution = 0.35)


DimPlot(Carotid_Femoral_srt.int.har, reduction = "umap", split.by = "orig.ident", label = TRUE)



#####INTEGRATE WITH SEURAT#####################
# Split the seurat object and integrate with CCA
Carotid_Femoral_srtList <- SplitObject(Carotid_Femoral_srt, split.by = "orig.ident")

# Normalize and identify variable features for each dataset independently
Carotid_Femoral_srtList <- lapply(X = Carotid_Femoral_srtList, SCTransform) 
features <- SelectIntegrationFeatures(object.list = Carotid_Femoral_srtList, 
                                      nfeatures = 3000)
Carotid_Femoral_srtList <- PrepSCTIntegration(object.list = Carotid_Femoral_srtList, 
                                              anchor.features = features)
Carotid_Femoral_srt.anchors <- FindIntegrationAnchors(object.list = Carotid_Femoral_srtList,
                                                      normalization.method = "SCT",
                                                      anchor.features = features)
Carotid_Femoral_srt.int <- IntegrateData(anchorset = Carotid_Femoral_srt.anchors, 
                                         normalization.method = "SCT")
# Re-run the pipeline on the integrated data

# Scale the data
Carotid_Femoral_srt.int <- ScaleData(Carotid_Femoral_srt.int)
# Run PCA
Carotid_Femoral_srt.int <- RunPCA(Carotid_Femoral_srt.int)
# Choose the number of principle components to keep
ElbowPlot(Carotid_Femoral_srt.int,ndims = 50)
# Find nearest neighbors and construct the graph
Carotid_Femoral_srt.int <- FindNeighbors(Carotid_Femoral_srt.int, k.param = 50, dims = 1:40)
# Find the clusters
Carotid_Femoral_srt.int <- FindClusters(Carotid_Femoral_srt.int, resolution = 0.5)
# Get the UMAP embedding
Carotid_Femoral_srt.int <- RunUMAP(Carotid_Femoral_srt.int, dims = 1:40)
Carotid_Femoral_srt.int <- RunTSNE(Carotid_Femoral_srt.int, dims = 1:40)
# Plot the UMAP with clustering
DimPlot(Carotid_Femoral_srt.int, reduction = "umap", split.by = "orig.ident", label = TRUE)
DimPlot(Carotid_Femoral_srt.int, reduction = "tsne", split.by = "orig.ident", label = TRUE) + NoLegend()
# Dim Plots
DimPlot(Carotid_Femoral_srt.int, reduction = "umap", group.by = "orig.ident") + DimPlot(Carotid_Femoral_srt, reduction = "umap", group.by = "orig.ident")
DimPlot(Carotid_Femoral_srt.int, group.by = "orig.ident") + DimPlot(Carotid_Femoral_srt.int, group.by = "seurat_clusters", label = TRUE)


DimPlot(Carotid_Femoral_srt.int, split.by = "orig.ident", group.by = "seurat_clusters", label = TRUE) + NoLegend()
DimPlot(Carotid_Femoral_srt.int, cells = unlist(CellsByIdentities(Carotid_Femoral_srt.int, idents = c(2,5,6,7,9,12,15))), label = TRUE) + NoLegend()

FeaturePlot(Carotid_Femoral_srt.int, features = "nCount_RNA")
FeaturePlot(Carotid_Femoral_srt.int, features = "nFeature_RNA")


# Compare Batch and Cluster ID
compTable <- table(Carotid_Femoral_srt.int$orig.ident, Carotid_Femoral_srt.int$seurat_clusters)
compTable <- (compTable / rowSums(compTable)) * 100
pheatmap(compTable, color = WtOrRd_pal)

