###Ensure Idents are set to the cluster
Idents(Carotid_Femoral_srt.int) <- Carotid_Femoral_srt.int$seurat_clusters

### Find Markers in all cells ###
All_markers_srt.int <- Seurat::FindAllMarkers(Carotid_Femoral_srt.int, assay = "RNA", only.pos = TRUE)

###Find Top 5 for each cluster and visualize
# Top 5 marker genes in each cluster ###
top5 <- All_markers_srt.int %>% group_by(cluster) %>%
  dplyr::slice_max(get(grep("^avg_log", colnames(All_markers_srt.int), value = TRUE)),
                   n = 5)
# Create the dot plot
Seurat::DotPlot(Carotid_Femoral_srt.int, features = unique(top5$gene)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,
                                                     size = 8, hjust = 1)) +
  Seurat::NoLegend()

# Create the heatmap
Seurat::DoHeatmap(Carotid_Femoral_srt.int, features = unique(top5$gene)) +
  Seurat::NoLegend() +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5))


###Find significant Markers and write to Excel
All_sig_markers_srt.int <- All_markers_srt.int %>% 
  filter(p_val_adj < .05)

write_csv(All_sig_markers_srt.int, file = "All_sigmarkers_srt_int.csv")


#####FIND GENERAL CELL TYPES#####
###Find B Cells ###
FeaturePlot(Carotid_Femoral_srt.int, features = "CD79A")
VlnPlot(Carotid_Femoral_srt.int, features = "CD79A")

###Find CD14+ Myeloid Cells ###
FeaturePlot(Carotid_Femoral_srt.int, features = "CD14")
VlnPlot(Carotid_Femoral_srt.int, features = "CD14")
FeaturePlot(Carotid_Femoral_srt.int, features = "CD68")
VlnPlot(Carotid_Femoral_srt.int, features = "CD68")

#Neutrophils ##
FeaturePlot(Carotid_Femoral_srt.int, features = "CSF3R")
VlnPlot(Carotid_Femoral_srt.int, features = "CSF3R")

#Mast Cells ##
FeaturePlot(Carotid_Femoral_srt.int, features = "KIT")
VlnPlot(Carotid_Femoral_srt.int, features = "KIT")

###Find T Cells ###
FeaturePlot(Carotid_Femoral_srt.int, features = "CD3D")
VlnPlot(Carotid_Femoral_srt.int, features = "CD3D")

###Find plasmacytoid DCs ###
FeaturePlot(Carotid_Femoral_srt.int, features = "CLEC4C")
VlnPlot(Carotid_Femoral_srt.int, features = "CLEC4C")

###Find Plasma Cells ###
FeaturePlot(Carotid_Femoral_srt.int, features = "CD27")
VlnPlot(Carotid_Femoral_srt.int, features = "CD27")

###Find DCs ###
FeaturePlot(Carotid_Femoral_srt.int, features = "CD86")
VlnPlot(Carotid_Femoral_srt.int, features = "CD86")
FeaturePlot(Carotid_Femoral_srt.int, features = "HLA-DRA")
VlnPlot(Carotid_Femoral_srt.int, features = "HLA-DRA")
FeaturePlot(Carotid_Femoral_srt.int, features = "HLA-DQB1")
VlnPlot(Carotid_Femoral_srt.int, features = "HLA-DQB1")
FeaturePlot(Carotid_Femoral_srt.int, features = "HLA-DPB1")
VlnPlot(Carotid_Femoral_srt.int, features = "HLA-DPB1")
FeaturePlot(Carotid_Femoral_srt.int, features = "HLA-DQA1")
VlnPlot(Carotid_Femoral_srt.int, features = "HLA-DQA1")
FeaturePlot(Carotid_Femoral_srt.int, features = "IRF8")
VlnPlot(Carotid_Femoral_srt.int, features = "IRF8")




compTable.allcells <- table(Carotid_Femoral_srt.int$orig.ident, Carotid_Femoral_srt.int$seurat_clusters)
compTable.allcells <- (compTable.allcells / rowSums(compTable.allcells)) * 100
