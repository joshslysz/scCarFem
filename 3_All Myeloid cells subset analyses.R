#Ensure Idents are set to the cluster
Idents(Carotid_Femoral_srt.int) <- Carotid_Femoral_srt.int$seurat_clusters

###############Subset CD14+ Clusters##############
Carotid_Femoral_srt.AllMyeloid <- subset(Carotid_Femoral_srt.int, idents = c(1,4,5,6,10,14,15,17))
DimPlot(Carotid_Femoral_srt.AllMyeloid, group.by = "seurat_clusters", split.by = "orig.ident", label = TRUE)

###############Extract raw counts from all Myeloid clusters and Re-cluster##############
###Split into Carotid and Femoral Samples###
Carotid_Femoral_srt.AllMyeloid.List <- SplitObject(Carotid_Femoral_srt.AllMyeloid, split.by = "orig.ident")

#Get assay data from Carotid sample
allmyeloid.raw.carotid.data <- as.matrix(GetAssayData(Carotid_Femoral_srt.AllMyeloid.List[["carotid"]], slot = "counts"))

#Make new seurat object with carotid CD4+ T cells
allmyeloid.carotid_srt <- CreateSeuratObject(counts = allmyeloid.raw.carotid.data , project = "carotid")

#Get assay data from Femoral sample
allmyeloid.raw.femoral.data <- as.matrix(GetAssayData(Carotid_Femoral_srt.AllMyeloid.List[["femoral"]], slot = "counts"))

#Make new seurat object with femoral CD4+ T cells
allmyeloid.femoral_srt <- CreateSeuratObject(counts = allmyeloid.raw.femoral.data, project = "femoral")

#merge the new CD4+ seurat objects
allmyeloid.Carotid_Femoral_srt <- merge(x= allmyeloid.carotid_srt, y= allmyeloid.femoral_srt)

## Normalize and Scale the data ##
# Log-transform the counts
allmyeloid.Carotid_Femoral_srt <- NormalizeData(allmyeloid.Carotid_Femoral_srt)
# Find Variable Features
allmyeloid.Carotid_Femoral_srt  <- FindVariableFeatures(allmyeloid.Carotid_Femoral_srt)
# Scale the data
allmyeloid.Carotid_Femoral_srt <- ScaleData(allmyeloid.Carotid_Femoral_srt)
# Run PCA
allmyeloid.Carotid_Femoral_srt  <- RunPCA(allmyeloid.Carotid_Femoral_srt)
#RUN UMAP
allmyeloid.Carotid_Femoral_srt <- RunUMAP(allmyeloid.Carotid_Femoral_srt, dims = 1:50)
#Show Dimplot
DimPlot(allmyeloid.Carotid_Femoral_srt, reduction = "umap", group.by = "orig.ident")

##INTEGRATE###
#Run Harmony
allmyeloid.Carotid_Femoral_srt.int <- RunHarmony(allmyeloid.Carotid_Femoral_srt, group.by.vars = "orig.ident")
#Run UMAP
allmyeloid.Carotid_Femoral_srt.int <- RunUMAP(allmyeloid.Carotid_Femoral_srt.int, reduction = "harmony", dims = 1:50)
#Show Dimplot
DimPlot(allmyeloid.Carotid_Femoral_srt.int, reduction = "umap", group.by = "orig.ident")
# Find nearest neighbors and construct the graph
allmyeloid.Carotid_Femoral_srt.int <- FindNeighbors(allmyeloid.Carotid_Femoral_srt.int, dims = 1:50)
# Find the clusters
allmyeloid.Carotid_Femoral_srt.int <- FindClusters(allmyeloid.Carotid_Femoral_srt.int, resolution = 0.25)
# Plot the UMAP with clustering
DimPlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", reduction = "umap", label = TRUE)
DimPlot(allmyeloid.Carotid_Femoral_srt.int, reduction = "umap", label = TRUE)



#################### Cluster proportions and log odds ratio between carotid and Femoral ########################
#Absolute and Proportion Tables
allmyeloid.compTable.abs <- table(allmyeloid.Carotid_Femoral_srt.int$orig.ident, allmyeloid.Carotid_Femoral_srt.int$seurat_clusters)
allmyeloid.compTable <- (allmyeloid.compTable.abs / rowSums(allmyeloid.compTable.abs)) * 100

#Log Odds ratio table
compTable.myleloid.odds <- table(allmyeloid.Carotid_Femoral_srt.int$seurat_clusters, allmyeloid.Carotid_Femoral_srt.int$orig.ident)
compTable.myleloid.odds 
compDFmyl <- as.data.frame.matrix(compTable.myleloid.odds )

compDFmyl <- cbind(compDFmyl, oddscar = compDFmyl$carotid/(sum(compDFmyl$carotid)-compDFmyl$carotid))
compDFmyl<- cbind(compDFmyl, oddsfem = compDFmyl$femoral/(sum(compDFmyl$femoral)-compDFmyl$femoral))
compDFmyl <- cbind(compDFmyl, oddscarratio = compDFmyl$oddscar/compDFmyl$oddsfem)
compDFmyl

compDFmyl <- cbind(compDFmyl, log_oddscarratio = log10(compDFmyl$oddscarratio))
compDFmyl <- cbind(compDFmyl, oddsfemratio = compDFmyl$oddsfem/compDFmyl$oddscar)
compDFmyl <- cbind(compDFmyl, log_oddsfemratio = log10(compDFmyl$oddsfemratio))
compDFmyl

library("epitools")
compTable.myeloid <- table(allmyeloid.Carotid_Femoral_srt.int$seurat_clusters, allmyeloid.Carotid_Femoral_srt.int$orig.ident)
compTable.myeloid
compTable.myeloid.Wald <- as.matrix(compTable.myeloid)
compTable.myeloid.Wald
Myeloid.odds <- oddsratio.wald(compTable.myeloid.Wald)
Myeloid.odds
oddsratio.wald(compTable.myeloid.Wald)

##odds ratios and Fishers exact test ####
#4000 carotid myeloid cells
#579 femoral myeloid cells 
install.packages("epitools")
library("epitools")
cluster0_ORtable <- matrix(c(1307,2693,177,402), nrow = 2, ncol = 2)
oddsratio.wald(cluster0_ORtable)

cluster1_ORtable <- matrix(c(911,3089,108,471), nrow = 2, ncol = 2)
oddsratio.wald(cluster1_ORtable)

cluster2_ORtable <- matrix(c(596,3404,162,417), nrow = 2, ncol = 2)
oddsratio.wald(cluster2_ORtable)

cluster3_ORtable <- matrix(c(560,3440,87,492), nrow = 2, ncol = 2)
oddsratio.wald(cluster3_ORtable)

cluster6_ORtable <- matrix(c(126,3874,2,577), nrow = 2, ncol = 2)
oddsratio.wald(cluster6_ORtable)

cluster8_ORtable <- matrix(c(53,3947,3,576), nrow = 2, ncol = 2)
oddsratio.wald(cluster8_ORtable)

cluster9_ORtable <- matrix(c(43,3957,7,572), nrow = 2, ncol = 2)
oddsratio.wald(cluster9_ORtable)

#log 10 transformation


####################################### Stacked Bar graph ################################
library(RColorBrewer)
df1 <- data.frame(samp=rep(c("Carotid", "femoral"), each=12),
                  Cluster=rep(c("0: Monocytes","1: Monocyte-derived DCs","2: Inflammatory Macs","3: Non-inlfammatory Macs","4: Misclass",
                                "5: Mast cells","6: Neutrophils","7: Proliferating DCs","8: cDCs type I","9: Interferon Macs","10: cDCs type II","11: TREM2+ Macs"),2),
                  ncells=c(1307, 911, 596, 560, 154, 129, 126, 121, 53, 43, 0, 0,
                           177, 108, 152 , 87, 0, 0, 2, 0, 3, 7, 23, 21)
                  )

ggplot(data = df1, aes(fill=Cluster, x=samp, y=ncells)) +
  geom_bar(position = "fill", stat= "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set3") +
  xlab("") +
  ylab("Percent") +
  theme_bw()


#Find Markers
allmyeloid.All.markers <- FindAllMarkers(allmyeloid.Carotid_Femoral_srt.int, 
                                   only.pos = TRUE)
allmyeloid.All.sig_markers <- allmyeloid.All.markers %>% 
  filter(p_val_adj < .05)

write_csv(allmyeloid.All.sig_markers, file = "allmyeloid.All.sig_markers.csv")


################################################################################
###############################################################################
##############################################################################

###############################################################################
#################General Myelopid Markers######################################
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CD14")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD14")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CD68")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD68")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CCR2")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "CCR2")

###############################################################################
##################Cluster 0: Monocytes##########################################
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CD14")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "CD14")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CD68")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "CD68")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "IL32")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "IL32")


VlnPlot(allmyeloid.Carotid_Femoral_srt.int, idents = "0", split.by = "orig.ident", features = "IL1B")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, idents = "0", split.by = "orig.ident", features = "IL32")

#### Differential State Analyses ####
#over-expressed in fem versus car
fem_vs_car_MyeloidCluster0 <- FindMarkers(allmyeloid.Carotid_Femoral_srt.int, ident.1 = "femoral", ident.2 = "carotid", only.pos = TRUE,
                                          group.by = "orig.ident", subset.ident = 0)

fem_vs_car_MyeloidCluster0.sig <- fem_vs_car_MyeloidCluster0 %>% 
  filter(p_val_adj < .05) %>%
  rownames_to_column("gene")

write_csv(fem_vs_car_MyeloidCluster0.sig, file = "fem_vs_car_MyeloidCluster0.sig.csv")

#over-expressed in car versus fem
car_vs_fem_MyeloidCluster0 <- FindMarkers(allmyeloid.Carotid_Femoral_srt.int, ident.1 = "carotid", ident.2 = "femoral", only.pos = TRUE,
                                          group.by = "orig.ident", subset.ident = 0)
car_vs_fem_MyeloidCluster0.sig <- car_vs_fem_MyeloidCluster0 %>% 
  filter(p_val_adj < .05) %>%
  rownames_to_column("gene") 

write_csv(car_vs_fem_MyeloidCluster0.sig, file = "car_vs_fem_MyeloidCluster0.sig.csv")
################################################################################

###############################################################################
######################Cluster 2: Pro-inflammatory Macrophages #################
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "IL1B")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "IL1B")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CASP4")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "CASP4")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CASP1")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "CASP1")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "TNF")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "TNF")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "IL6")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "IL6")
#Cochain, C. et al. Single-cell RNA-seq reveals the transcriptional landscape and heterogeneity of aortic macrophages in murine atherosclerosis. Circ. Res. 122, 1661–1674 (2018).
#Winkels, H. et al. Atlas of the immune cell repertoire in mouse atherosclerosis defined by single- cel RNA- sequencing and mass cytometry. Circ. Res.  122, 1675–1688 (2018).

#### Differential State Analyses ####
#over-expressed in fem versus car
fem_vs_car_MyeloidCluster2 <- FindMarkers(allmyeloid.Carotid_Femoral_srt.int, ident.1 = "femoral", ident.2 = "carotid", only.pos = TRUE,
                                          group.by = "orig.ident", subset.ident = 2)

fem_vs_car_MyeloidCluster2.sig <- fem_vs_car_MyeloidCluster2 %>% 
  filter(p_val_adj < .05) %>%
  rownames_to_column("gene")

write_csv(fem_vs_car_MyeloidCluster2.sig, file = "fem_vs_car_MyeloidCluster2.sig.csv")

#over-expressed in car versus fem
car_vs_fem_MyeloidCluster2 <- FindMarkers(allmyeloid.Carotid_Femoral_srt.int, ident.1 = "carotid", ident.2 = "femoral", only.pos = TRUE,
                                          group.by = "orig.ident", subset.ident = 2)
car_vs_fem_MyeloidCluster2.sig <- car_vs_fem_MyeloidCluster2 %>% 
  filter(p_val_adj < .05) %>%
  rownames_to_column("gene") 

write_csv(car_vs_fem_MyeloidCluster2.sig, file = "car_vs_fem_MyeloidCluster2.sig.csv")

VlnPlot(allmyeloid.Carotid_Femoral_srt.int, idents = "2", split.by = "orig.ident", features = "IL1B")


###############################################################################
######################Cluster 3: Non-inflammatory Macrophages #################
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "IL1B")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "IL1B")

#### Differential State Analyses ####
#over-expressed in fem versus car
fem_vs_car_MyeloidCluster3 <- FindMarkers(allmyeloid.Carotid_Femoral_srt.int, ident.1 = "femoral", ident.2 = "carotid", only.pos = TRUE,
                                          group.by = "orig.ident", subset.ident = 3)

fem_vs_car_MyeloidCluster3.sig <- fem_vs_car_MyeloidCluster3 %>% 
  filter(p_val_adj < .05) %>%
  rownames_to_column("gene")

write_csv(fem_vs_car_MyeloidCluster3.sig, file = "fem_vs_car_MyeloidCluster3.sig.csv")

#over-expressed in car versus fem
car_vs_fem_MyeloidCluster3 <- FindMarkers(allmyeloid.Carotid_Femoral_srt.int, ident.1 = "carotid", ident.2 = "femoral", only.pos = TRUE,
                                          group.by = "orig.ident", subset.ident = 3)
car_vs_fem_MyeloidCluster3.sig <- car_vs_fem_MyeloidCluster3 %>% 
  filter(p_val_adj < .05) %>%
  rownames_to_column("gene") 

write_csv(car_vs_fem_MyeloidCluster3.sig, file = "car_vs_fem_MyeloidCluster3.sig.csv")

VlnPlot(allmyeloid.Carotid_Femoral_srt.int, idents = "3", split.by = "orig.ident", features = "IL1B")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, idents = "3", split.by = "orig.ident", features = "CD44")

###############################################################################
######################Cluster 9: Interferon Macrophages #################
#Interferon Macrophages
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "IFIT3")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "IFIT3")

#### Differential State Analyses ####
#over-expressed in fem versus car
fem_vs_car_MyeloidCluster9 <- FindMarkers(allmyeloid.Carotid_Femoral_srt.int, ident.1 = "femoral", ident.2 = "carotid", only.pos = TRUE,
                                          group.by = "orig.ident", subset.ident = 9)

fem_vs_car_MyeloidCluster9.sig <- fem_vs_car_MyeloidCluster9 %>% 
  filter(p_val_adj < .05) %>%
  rownames_to_column("gene")

write_csv(fem_vs_car_MyeloidCluster9.sig, file = "fem_vs_car_MyeloidCluster9.sig.csv")

#over-expressed in car versus fem
car_vs_fem_MyeloidCluster9 <- FindMarkers(allmyeloid.Carotid_Femoral_srt.int, ident.1 = "carotid", ident.2 = "femoral", only.pos = TRUE,
                                          group.by = "orig.ident", subset.ident = 9)
car_vs_fem_MyeloidCluster9.sig <- car_vs_fem_MyeloidCluster9 %>% 
  filter(p_val_adj < .05) %>%
  rownames_to_column("gene") 

write_csv(car_vs_fem_MyeloidCluster9.sig, file = "car_vs_fem_MyeloidCluster9.sig.csv")

###############################################################################
######################Cluster 11: TREM2+ Macrophages ########################
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "TREM2")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "TREM2")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CD9")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "CD9")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "FABP4")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "FABP4")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CTSD")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "CTSD")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "ABCG1")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "ABCG1")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "IL10")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "IL10")

###############################################################################
######################Cluster 1: Monocyte-derived Monocytes ####################
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "THBD") #CD141
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "THBD")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "XCR1") #CD172
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "XCR1")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CADM1")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "CADM1")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CLEC10A")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "CLEC10A")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "FCGR1A")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "FCGR1A")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD74")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD74")

#### Differential State Analyses ####
#over-expressed in fem versus car
fem_vs_car_MyeloidCluster1 <- FindMarkers(allmyeloid.Carotid_Femoral_srt.int, ident.1 = "femoral", ident.2 = "carotid", only.pos = TRUE,
                                          group.by = "orig.ident", subset.ident = 1)

fem_vs_car_MyeloidCluster1.sig <- fem_vs_car_MyeloidCluster1 %>% 
  filter(p_val_adj < .05) %>%
  rownames_to_column("gene")

write_csv(fem_vs_car_MyeloidCluster1.sig, file = "fem_vs_car_MyeloidCluster1.sig.csv")

#over-expressed in car versus fem
car_vs_fem_MyeloidCluster1 <- FindMarkers(allmyeloid.Carotid_Femoral_srt.int, ident.1 = "carotid", ident.2 = "femoral", only.pos = TRUE,
                                          group.by = "orig.ident", subset.ident = 1)
car_vs_fem_MyeloidCluster1.sig <- car_vs_fem_MyeloidCluster1 %>% 
  filter(p_val_adj < .05) %>%
  rownames_to_column("gene") 

write_csv(car_vs_fem_MyeloidCluster1.sig, file = "car_vs_fem_MyeloidCluster1.sig.csv")

###############################################################################
######################Cluster 10: cDCc type 2 ####################
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "ITGAX") #CD1C+ 
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "ITGAX")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CD1C") #CD1C+ 
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "CD1C")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CLEC9A") #CLEC9A-
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CLEC9A")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "HLA-DRA") #HLA-DR
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "HLA-DRA")

###############################################################################
######################Cluster 8: cDCc type 1 ###############################
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CLEC9A")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "CLEC9A")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "FLT3")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "FLT3")

#### Differential State Analyses ####
#over-expressed in fem versus car
fem_vs_car_MyeloidCluster8 <- FindMarkers(allmyeloid.Carotid_Femoral_srt.int, ident.1 = "femoral", ident.2 = "carotid", only.pos = TRUE,
                                          group.by = "orig.ident", subset.ident = 8)

fem_vs_car_MyeloidCluster8.sig <- fem_vs_car_MyeloidCluster8 %>% 
  filter(p_val_adj < .05) %>%
  rownames_to_column("gene")

write_csv(fem_vs_car_MyeloidCluster8.sig, file = "fem_vs_car_MyeloidCluster8.sig.csv")

#over-expressed in car versus fem
car_vs_fem_MyeloidCluster8 <- FindMarkers(allmyeloid.Carotid_Femoral_srt.int, ident.1 = "carotid", ident.2 = "femoral", only.pos = TRUE,
                                          group.by = "orig.ident", subset.ident = 8)
car_vs_fem_MyeloidCluster8.sig <- car_vs_fem_MyeloidCluster8 %>% 
  filter(p_val_adj < .05) %>%
  rownames_to_column("gene") 

write_csv(car_vs_fem_MyeloidCluster8.sig, file = "car_vs_fem_MyeloidCluster8.sig.csv")


######################################################################################
###### Neutrophils ##############

FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CSF3R")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CSF3R")

######################################################################################
######## #Mast Cells ########
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "KIT")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "KIT")

########################################################################
#Cluster 3: Misclass cluster
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CD3E")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "CD3E")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CD8A")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, features = "CD8A")
FeaturePlot(allmyeloid.Carotid_Femoral_srt.int, features = "CD3D")
VlnPlot(allmyeloid.Carotid_Femoral_srt.int, split.by = "orig.ident", features = "CD8A")
 #########################################################################################


#Cell Cycle
allmyeloid.Carotid_Femoral_srt.int <- CellCycleScoring(allmyeloid.Carotid_Femoral_srt.int, search = TRUE,
                                                 s.features = cc.genes.updated.2019$s.genes,
                                                 g2m.features = cc.genes.updated.2019$g2m.genes)
# Plot in UMAP
DimPlot(allmyeloid.Carotid_Femoral_srt.int, group.by = "Phase", split.by = "orig.ident")



