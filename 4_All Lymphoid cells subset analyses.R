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

#################### Cluster proportions and log odds ratio between carotid and Femoral ########################
#Absolute and Proportion Tables
alllymphoid.compTable.abs <- table(Alllymphoid.Carotid_Femoral_srt.int$orig.ident, Alllymphoid.Carotid_Femoral_srt.int$seurat_clusters)
alllymphoid.compTable.rel <- (alllymphoid.compTable.abs / rowSums(alllymphoid.compTable.abs)) * 100


###################################################################
################ Cell Typing########################
##################################################################

#################Find Markers #####################################
#Find Markers
alllymphoid.All.markers <- FindAllMarkers(Alllymphoid.Carotid_Femoral_srt.int, 
                                         only.pos = TRUE)
write_csv(alllymphoid.All.markers, file = "alllymphoid.All.markers.csv")


alllymphoid.All.sig_markers <- alllymphoid.All.markers %>% 
  filter(p_val_adj < .05)

write_csv(alllymphoid.All.sig_markers, file = "alllymphoid.All.sig_markers.csv")

################# T cells ######################################
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD3D")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD3D") # T cells and NKs
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD2")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD2") #T cells and NKs
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD69")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD69") #Exhausted T cells
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD44")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD44") #Distinguishes memory and effector Ts from Naive

################# Naive T cells ######################################
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD28")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD28")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CCR7")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CCR7")


###################################################################
################# CD8+ T cells ######################################
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD8A")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD8A")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD8B")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD8B")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "GZMK")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "GZMK")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "GZMB")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "GZMB")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "IFNG")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "IFNG")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "TNF")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "TNF")



#######################################################################
################# CD4+ T cells ######################################
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD4")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD4")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD28")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD28")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CCR7")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CCR7") # central memory T cells


### CD4+ T reg CElls ### CLUSTER 9#### TEST TEST TEST
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "IL2RA")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "IL2RA")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "Fox3P")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "FOX3P")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "Fox3P")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "FOX3P")
### CD4+ Th1 Cells ### DOES NOT SEEM TO PRESENT IN DATASET
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "TBX21")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "TBX21")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CXCR3")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CXCR3")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CCR5")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CCR5")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "IFNG")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "IFNG")


### CD4+ Th2 Cells ###DOES NOT SEE TO BE PRESENT IN DATASET
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "IL4")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "IL4")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "GATA3")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "GATA3") # only really expressed in carotid
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "IL5")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "IL5")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "IL13")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "IL13")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "IL10")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "IL10")

### CD4+ Th9 Cells ###DOES NOT SEE TO BE PRESENT IN DATASET
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "IL9")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "IL9")

### CD4+ Th17 Cells ###DOES NOT SEE TO BE PRESENT IN DATASET
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "RORC")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "RORC")


### CD4+ Th22 Cells 
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "IL22")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "IL22")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "AHR")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "AHR")

### CD4+ ThFH Cells 
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "BCL6")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "BCL6")

###Find Gamma sigma T cells  ### POssible cluster 14
#cholesterol-related pathway genes: Cheng, Hsin-Yuan, et al. "Increased cholesterol content in gammadelta (γδ) T lymphocytes differentially regulates their activation." PloS one 8.5 (2013): e63746.
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "ACAT1") #cholesterol esterification
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, features = "ACAT1")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "ACAT2")#cholesterol esterification
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, features = "ACAT2")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "TRERF1")#cholesterol utilization
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, features = "TRERF1")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "IL6")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, features = "IL6")


###Find plasmacytoid DCs ### CLUSTER 10
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CLEC4C")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CLEC4C")

###Find Plasma Cells ### 
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD27")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD27")

###Proliferating T Cells ### CLUSTER 15 & CLuster 8
Alllymphoid.Carotid_Femoral_srt.int <- CellCycleScoring(Alllymphoid.Carotid_Femoral_srt.int, search = TRUE,
                                                       s.features = cc.genes.updated.2019$s.genes,
                                                       g2m.features = cc.genes.updated.2019$g2m.genes)
DimPlot(Alllymphoid.Carotid_Femoral_srt.int, group.by = "Phase", split.by = "orig.ident")

###Natural Killer  Cells ### CLUSTER 4 
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "KLRC1")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "KLRC1")

### T reg CElls ### CLUSTER 10
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "IL2RA")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "IL2RA")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "Fox3P")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "FOX3P")


####CLUSTER 12 ####
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CCR7")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CCR7")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CCL19")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CCL19")
FeaturePlot(Alllymphoid.Carotid_Femoral_srt.int, features = "CD83")
VlnPlot(Alllymphoid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD83")


