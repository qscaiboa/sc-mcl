library(Seurat)
library(Seurat)
library(ggplot2)
library(readr)
library(readxl)
library(Matrix)
library(tidyr)
library(dplyr)
library(readxl)
library(cowplot)

Apt <- readRDS ("./Data/Object.PatientA&PBMC.human_only.RDS")
Apt <- subset(Apt,       subset = seurat_clusters %in% c(0:8,10,11,13,14,16) )
 Apt <- NormalizeData(object = Apt, normalization.method = "LogNormalize", scale.factor = 10000)  #, verbose = FALSE)
 Apt <- FindVariableFeatures(object = Apt, selection.method = "vst", nfeatures = 2000)
 Apt <- ScaleData(Apt,verbose = FALSE)
 Apt <- RunPCA(Apt, features = VariableFeatures(object = Apt))
 Apt <- FindNeighbors(Apt, dims = 1:10)
 Apt <- FindClusters(Apt, resolution = 0.5)
 Apt <- RunUMAP(Apt, dims = 1:30)
 Apt <- RunTSNE(Apt, dims = 1:30)
saveRDS (Apt, file= "./Data/Object.PatientA&PBMC.Bcell_only.RDS")

 png("ZZZ.1020_Apt.Bcell_only.tsne.png", width = 600, height = 600)
 DimPlot(Apt, reduction = "tsne", label = TRUE)
 dev.off()

 png("ZZZ.1020_Apt.Bcell_only.tsne.Finalname.png", width = 900, height = 600)
 DimPlot(Apt, reduction = "tsne", label = TRUE, group.by= "Finalname")
 dev.off()

 png("ZZZ.1020.Apt.Bcell_only.tsne.splitbyFinalname.groupbyPt.png", width = 1200, height = 800)
 DimPlot(Apt, reduction = "tsne", label = TRUE, split.by= "Finalname",group.by= "species", ncol=8)
 dev.off()
 
 png("ZZZ.1020.Apt.Bcell_only.tsne.splitbyFinalname.png", width = 1200, height = 800)
 DimPlot(Apt, reduction = "tsne", label = TRUE, split.by= "Finalname", ncol=8)
 dev.off()
 
  
 png("ZZZ.1020.Apt.Bcell_only.tsne.byPt.png", width = 600, height = 600)
 DimPlot(Apt, reduction = "tsne", label = TRUE, group.by= "Pt")
 dev.off()
 
 
 png("ZZZ.1020_Apt.Bcell_only.tsne.FeaturePlot.PBMC.png", width = 1200, height = 800)
 FeaturePlot(Apt, features = c("IL7R", "CCR7",
                               "S100A4",
                               "CD14", "LYZ",
							   "MS4A1","CD79A", "CD79B",
							   "CD8A",
							   "FCGR3A", "MS4A7",
							   "NKG7","GNLY",
							   "S100A9", 
							   "FCER1A", "CST3", 
							   "PPBP",
							   "ZNF90","MCL1","TYMS" , "HIST1H4C", "TUBA1B", "TUBB", "PCLAF"),reduction = "tsne", min.cutoff = "q9", ncol=6)
 dev.off()
 
 
 
 
  
 library(harmony, lib.loc = "/rsrch3/scratch/lym_myl_rsch/qcai1/R/library/4.0.0")
 #Apt <- subset(Apt,       subset = seurat_clusters %in% c(0:8,10,11,13,14,16) )
 #Apt <- NormalizeData(object = Apt, verbose = FALSE)
 #Apt <- ScaleData(Apt_Bcells_harmony,verbose = FALSE)
 Apt <- RunHarmony(Apt, "sample")
 Apt <- RunUMAP(Apt, reduction = "harmony", dims = 1:30)
 #Apt_Bcells_harmony <- RunTSNE( Apt_Bcells_harmony, reduction = "harmony", dims = 1:30)
 Apt <- FindNeighbors( Apt, reduction = "harmony", dims = 1:30)
 Apt <- FindClusters( Apt, resolution = 0.5)
saveRDS (Apt, file= "./Data/Object.PatientA&PBMC.Bcell_only.harmony.RDS")


 png("ZZZ.1020_Apt.Bcell_only.harmony.umap.png", width = 600, height = 600)
 DimPlot(Apt, reduction = "umap", label = TRUE)
 dev.off()

 png("ZZZ.1020_Apt.Bcell_only.harmony.umap.Finalname.png", width = 900, height = 600)
 DimPlot(Apt, reduction = "umap", label = TRUE, group.by= "Finalname")
 dev.off()

 png("ZZZ.1020.Apt.Bcell_only.harmony.umap.splitbyFinalname.groupbyPt.png", width = 1200, height = 800)
 DimPlot(Apt, reduction = "umap", label = TRUE, split.by= "Finalname",group.by= "species", ncol=8)
 dev.off()
 
 png("ZZZ.1020.Apt.Bcell_only.harmony.umap.splitbyFinalname.png", width = 1200, height = 800)
 DimPlot(Apt, reduction = "umap", label = TRUE, split.by= "Finalname", ncol=8)
 dev.off()
 
  
 png("ZZZ.1020.Apt.Bcell_only.harmony.umap.byPt.png", width = 600, height = 600)
 DimPlot(Apt, reduction = "umap", label = TRUE, group.by= "Pt")
 dev.off()
 
 
 png("ZZZ.1020_Apt.Bcell_only.harmony.umap.FeaturePlot.PBMC.png", width = 1200, height = 800)
 FeaturePlot(Apt, features = c("IL7R", "CCR7",
                               "S100A4",
                               "CD14", "LYZ",
							   "MS4A1","CD79A", "CD79B",
							   "CD8A",
							   "FCGR3A", "MS4A7",
							   "NKG7","GNLY",
							   "S100A9", 
							   "FCER1A", "CST3", 
							   "PPBP",
							   "ZNF90","MCL1","TYMS" , "HIST1H4C", "TUBA1B", "TUBB", "PCLAF"),reduction = "umap", min.cutoff = "q9", ncol=6)
 dev.off()
 







