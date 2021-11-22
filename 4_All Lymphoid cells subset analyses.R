#Ensure Idents are set to the cluster
Idents(Carotid_Femoral_srt.int) <- Carotid_Femoral_srt.int$seurat_clusters

###############Subset CD14+ Clusters##############
Carotid_Femoral_srt.Alllymphoid <- subset(Carotid_Femoral_srt.int, idents = c(0,2,3,7,8,9,11,12,13,16))
DimPlot(Carotid_Femoral_srt.Alllymphoid, group.by = "seurat_clusters", split.by = "orig.ident", label = TRUE)

###############Extract raw counts from all Myeloid clusters and Re-cluster##############
###Split into Carotid and Femoral Samples###
Carotid_Femoral_srt.Alllymphoid.List <- SplitObject(Carotid_Femoral_srt.Alllymphoid, split.by = "orig.ident")


#Get assay data from Carotid sample
Alllymphoid.raw.carotid.data <- as.matrix(GetAssayData(Carotid_Femoral_srt.Alllymphoid.List[["carotid"]], slot = "counts"))

#Make new seurat object with carotid CD4+ T cells
Alllymphoid.carotid_srt <- CreateSeuratObject(counts = Alllymphoid.raw.carotid.data , project = "carotid")

#Get assay data from Femoral sample
Alllymphoid.raw.femoral.data <- as.matrix(GetAssayData(Carotid_Femoral_srt.Alllymphoid.List[["femoral"]], slot = "counts"))

#Make new seurat object with femoral CD4+ T cells
Alllymphoid.femoral_srt <- CreateSeuratObject(counts = Alllymphoid.raw.femoral.data, project = "femoral")

#merge the new CD4+ seurat objects
Alllymphoid.Carotid_Femoral_srt <- merge(x= Alllymphoid.carotid_srt, y= Alllymphoid.femoral_srt)

## Normalize and Scale the data ##
# Log-transform the counts
Alllymphoid.Carotid_Femoral_srt <- NormalizeData(Alllymphoid.Carotid_Femoral_srt)
# Find Variable Features
Alllymphoid.Carotid_Femoral_srt  <- FindVariableFeatures(Alllymphoid.Carotid_Femoral_srt)
# Scale the data
Alllymphoid.Carotid_Femoral_srt <- ScaleData(Alllymphoid.Carotid_Femoral_srt)
# Run PCA
Alllymphoid.Carotid_Femoral_srt  <- RunPCA(Alllymphoid.Carotid_Femoral_srt)
#RUN UMAP
Alllymphoid.Carotid_Femoral_srt <- RunUMAP(Alllymphoid.Carotid_Femoral_srt, dims = 1:50)
#Show Dimplot
DimPlot(Alllymphoid.Carotid_Femoral_srt, reduction = "umap", group.by = "orig.ident")

##INTEGRATE###
#Run Harmony
Alllymphoid.Carotid_Femoral_srt.int <- RunHarmony(Alllymphoid.Carotid_Femoral_srt, group.by.vars = "orig.ident")
#Run UMAP
Alllymphoid.Carotid_Femoral_srt.int <- RunUMAP(Alllymphoid.Carotid_Femoral_srt.int, reduction = "harmony", dims = 1:50)
#Show Dimplot
DimPlot(Alllymphoid.Carotid_Femoral_srt.int, reduction = "umap", group.by = "orig.ident")
# Find nearest neighbors and construct the graph
Alllymphoid.Carotid_Femoral_srt.int <- FindNeighbors(Alllymphoid.Carotid_Femoral_srt.int, dims = 1:30)
# Find the clusters
Alllymphoid.Carotid_Femoral_srt.int <- FindClusters(Alllymphoid.Carotid_Femoral_srt.int, resolution = 0.4)
# Plot the UMAP with clustering
DimPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", reduction = "umap", label = TRUE)
DimPlot(Alllymphoid.Carotid_Femoral_srt.int, reduction = "umap", label = TRUE)


############################ Cell Typing############################################
#################General lymphoid Markers######################################
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD3D")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD3D")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD8A")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD8A")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD4")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD4")

###
###Find T Cells ###
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD3D")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD3D")

###Find plasmacytoid DCs ###
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CLEC4C")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CLEC4C")

###Find Plasma Cells ###
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD27")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD27")



alllymphoid.compTable.abs <- table(Alllymphoid.Carotid_Femoral_srt.int$orig.ident, Alllymphoid.Carotid_Femoral_srt.int$seurat_clusters)
alllymphoid.compTable <- (allmyeloid.compTable / rowSums(allmyeloid.compTable)) * 100



