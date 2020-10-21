
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




 #Apt <- readRDS("./Data/Object.PatientA&PBMC.RDS")
 #png("ZZZ.1020.Apt.tsne.byspecies.png", width = 600, height = 600)
 #DimPlot(Apt, reduction = "tsne", label = TRUE, group.by= "species")
 #dev.off()


 Objectall <- readRDS("./Data/Object.merged_updated.rds")
 Apt <- subset(Objectall, subset = Pt %in% c("A","P"))
 Apt <- subset(Apt, subset = species %in% "human")
 Apt <- NormalizeData(object = Apt, normalization.method = "LogNormalize", scale.factor = 10000)  #, verbose = FALSE)
 Apt <- FindVariableFeatures(object = Apt, selection.method = "vst", nfeatures = 2000)
 Apt <- ScaleData(Apt,verbose = FALSE)
 Apt <- RunPCA(Apt, features = VariableFeatures(object = Apt))
 Apt <- FindNeighbors(Apt, dims = 1:10)
 Apt <- FindClusters(Apt, resolution = 0.5)
 Apt <- RunUMAP(Apt, dims = 1:30)
 Apt <- RunTSNE(Apt, dims = 1:30)
saveRDS (Apt, file= "./Data/Object.PatientA&PBMC.RDS")

 png("ZZZ.1020_Apt.tsne.png", width = 600, height = 600)
 DimPlot(Apt, reduction = "tsne", label = TRUE)
 dev.off()

 png("ZZZ.1020_Apt.tsne.Finalname.png", width = 900, height = 600)
 DimPlot(Apt, reduction = "tsne", label = TRUE, group.by= "Finalname")
 dev.off()
 
 png("ZZZ.1020.Apt.tsne.byPt.png", width = 600, height = 600)
 DimPlot(Apt, reduction = "tsne", label = TRUE, group.by= "Pt")
 dev.off()
 
  png("ZZZ.1020.Apt.tsne.byspecies.png", width = 600, height = 600)
 DimPlot(Apt, reduction = "tsne", label = TRUE, group.by= "species")
 dev.off()


 png("ZZZ.1020.Apt.tsne.splitbyFinalname.groupbyPt.png", width = 1200, height = 800)
 DimPlot(Apt, reduction = "tsne", label = TRUE, split.by= "Finalname",group.by= "species", ncol=8)
 dev.off()
 
 png("ZZZ.1020.Apt.tsne.splitbyFinalname.png", width = 1200, height = 800)
 DimPlot(Apt, reduction = "tsne", label = TRUE, split.by= "Finalname", ncol=8)
 dev.off()
 
 
 png("ZZZ.1020_Apt.tsne.FeaturePlot.PBMC.png", width = 1200, height = 800)
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
 
 
 
 
 
 
 
 
 #Apt <- readRDS("./Data/Object.PatientA&PBMC.harmony.RDS")
 #png("ZZZ.1020_Apt.harmony.umap.species.png", width = 600, height = 600)
 #DimPlot(Apt, reduction = "umap", label = TRUE, group.by= "species")
 #dev.off()
 
 
 
 
library(harmony, lib.loc = "/rsrch3/scratch/lym_myl_rsch/qcai1/R/library/4.0.0")
 Apt<- RunHarmony( Apt, "sample")
 Apt<- ScaleData( Apt,verbose = FALSE)
 Apt<- RunUMAP( Apt, reduction = "harmony", dims = 1:30)
 Apt<- RunTSNE( Apt, reduction = "harmony", dims = 1:30)
 Apt<- FindNeighbors( Apt, reduction = "harmony", dims = 1:30)
 Apt<- FindClusters( Apt, resolution = 0.5)
saveRDS (Apt, file= "./Data/Object.PatientA&PBMC.harmony.RDS")

 png("ZZZ.1020_Apt.harmony.umap.png", width = 600, height = 600)
 DimPlot(Apt, reduction = "umap", label = TRUE)
 dev.off()
 
 png("ZZZ.1020.Apt.harmony.umap.byPt.png", width = 600, height = 600)
 DimPlot(Apt, reduction = "umap", label = TRUE, group.by= "Pt")
 dev.off()
 
 png("ZZZ.1020_Apt.harmony.umap.species.png", width = 600, height = 600)
 DimPlot(Apt, reduction = "umap", label = TRUE, group.by= "species")
 dev.off()
 
 
 png("ZZZ.1020_Apt.harmony.umap.png", width = 600, height = 600)
 DimPlot(Apt, reduction = "umap", label = TRUE)
 dev.off()

 png("ZZZ.1020.Apt.harmony.umap.splitbyFinalname.groupbyPt.png", width = 1200, height = 800)
 DimPlot(Apt, reduction = "umap", label = TRUE, split.by= "Finalname",group.by= "species", ncol=8)
 dev.off()
 
 png("ZZZ.1020.Apt.harmony.umap.splitbyFinalname.png", width = 1200, height = 800)
 DimPlot(Apt, reduction = "umap", label = TRUE, split.by= "Finalname", ncol=8)
 dev.off()
  
 png("ZZZ.1020_Apt.harmony.umap.FeaturePlot.PBMC.png", width = 1200, height = 800)
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
 


 png("ZZZ.1020.harmony_Apt.umap.png", width = 600, height = 600)
 DimPlot( Apt, reduction = "umap", label = TRUE)
 dev.off()
  
 png("ZZZ.1020.harmony_Apt.species.umap.png", width = 1200, height = 1200)
 DimPlot( Apt, reduction = "umap", label = TRUE, group.by= "species")
 dev.off()
 
 
 png("ZZZ.1020.harmony_Apt.umap.pt.png", width = 2400, height = 1200)
 DimPlot( Apt, reduction = "umap", label = TRUE, group.by= "Pt")
 dev.off()
 
 
 png("ZZZ.1020.harmony_Apt.umap.Finalname.png", width = 2400, height = 1200)
 DimPlot( Apt, reduction = "umap", label = TRUE, group.by= "Finalname")
 dev.off()
 
 
 png("ZZZ.1020.harmony_Apt.umap.splitbyFinalname.png", width = 2400, height = 2400)
 DimPlot( Apt, reduction = "umap", label = TRUE, split.by= "Finalname",group.by= "species", ncol=4)
 dev.off()
 
 
 
 png("ZZZ.1020.harmony_Apt.FeaturePlot.png", width = 2400, height = 2400)
 FeaturePlot( Apt, features = c("CD79A", "CD79B","NKG7", "GNLY","S100A9", "HBB", "ZNF90","MCL1",
                                    "HIST1H4C", "TUBA1B", "TUBB", "PCLAF","TYMS"),reduction = "umap", min.cutoff = "q9")
 dev.off()			
 
 
 
