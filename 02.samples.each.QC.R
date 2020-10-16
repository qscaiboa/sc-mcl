
 
 bsub -Is -q interactive -W 6:00 -M 200 -R rusage[mem=200] -n 1 /bin/bash
R

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
dataDir  =  "/rsrch3/scratch/lym_myl_rsch/qcai1/SC145/Data.RDS/"
figureDIR =  "/rsrch3/scratch/lym_myl_rsch/qcai1/SC145/Figures/"
workDIR = "/rsrch3/scratch/lym_myl_rsch/qcai1/SC145/10XMH/" 

samplePDX <- read_xlsx("/rsrch3/home/lym_myl_rsch/MCL_Lab/more_data_from_linghua/from_Shaojun/Vivian/scPDX/sampleInfo/PDX_scRNA_WES_info.xlsx", sheet = 1)
unique(as.character(samplePDX$scSeqName))
data.frame(samplePDX)

allSampleName = list.files(dataDir)
#allSampleName <- allSampleName[!allSampleName %in% "scripts"]
allSampleName  <- gsub(".RDS","",allSampleName)

initialSampleName = "S10"
allSampleName = allSampleName[!allSampleName %in% initialSampleName]


i <- initialSampleName
Object1 <- readRDS(paste0("./Data.RDS/", i , ".RDS"))
Object1@meta.data$Pt <-    samplePDX[samplePDX$scSeqName %in% i,]$Pt   ######################################
Object1@meta.data$sample <-    samplePDX[samplePDX$scSeqName %in% i,]$scSeqName
Object1@meta.data$Newname <-    samplePDX[samplePDX$scSeqName %in% i,]$Newname
Object1@meta.data$HumanMouse <-    samplePDX[samplePDX$scSeqName %in% i,]$HumanMouse
Object1@meta.data$PDXmodel <-    samplePDX[samplePDX$scSeqName %in% i,]$PDXmodel
Object1 <- NormalizeData(Object1)
Object1 <- FindVariableFeatures(Object1, selection.method = "vst", nfeatures= 2000)
Object1 <- ScaleData(Object1)

Object1 <- subset(Object1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)


Object.combined = c()
for ( i in allSampleName ) {
    print(i)
    j = i
    Object32 <- readRDS(paste0("./Data.RDS/", i , ".RDS"))
    
Object32@meta.data$Pt <-    samplePDX[samplePDX$scSeqName %in% i,]$Pt   ######################################
Object32@meta.data$sample <-    samplePDX[samplePDX$scSeqName %in% i,]$scSeqName
Object32@meta.data$Newname <-    samplePDX[samplePDX$scSeqName %in% i,]$Newname 
Object32@meta.data$HumanMouse <-    samplePDX[samplePDX$scSeqName %in% i,]$HumanMouse
Object32@meta.data$PDXmodel <-    samplePDX[samplePDX$scSeqName %in% i,]$PDXmodel
Object32[["percent.mt"]] <- PercentageFeatureSet(Object32, pattern = "^MT-")

Object32 <- NormalizeData(Object32)
Object32 <- FindVariableFeatures(Object32, selection.method = "vst", nfeatures= 2000)
Object32 <- ScaleData(Object32)
Object32 <- subset(Object32, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
Object32    
   Object.combined = c(Object.combined, Object32)
  }
  
  
  
  
  

 Objectall <- merge(Object1, y =   Object.combined, add.cell.ids =c(initialSampleName, allSampleName)  , project = "SC145")
 Objectall <- NormalizeData(object = Objectall, verbose = FALSE)
 Objectall <- FindVariableFeatures(object = Objectall, selection.method = "vst", nfeatures = 2000)
 Objectall <- RunPCA(Objectall, features = VariableFeatures(object = Objectall))
 saveRDS(Objectall, "./Data/Object1.CQS.rds")

library(harmony)
Objectall <- RunHarmony(Objectall, "sample")
Objectall <- ScaleData(Objectall,verbose = FALSE)
Objectall <- RunUMAP(Objectall, reduction = "harmony", dims = 1:30)
Objectall <- FindNeighbors(Objectall, reduction = "harmony", dims = 1:30)
Objectall <- FindClusters(Objectall, resolution = 0.5)
saveRDS(Objectall, "/data/exx/SC/SC145/Data/Hamed.CQS.rds")
 



for (samplename in allSampleName) {
Object2 <- readRDS(paste0("./Data.RDS/", samplename, ".RDS"))
 #print(table( Object2@meta.data$species))
a <- data.frame(t(as.matrix(table(Object2@meta.data$species))) )

if(!"mouse" %in% colnames(a))
 {
   a$mouse <- 0
 }

 print(samplename)
 print(a)
 
 

Object2[["percent.mt"]]<- PercentageFeatureSet(Object2, pattern = "MT-")
png(paste0(figureDIR,samplename,"_QC_before.png"), width = 920, height = 320)
print(VlnPlot(Object2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
dev.off()



Object2 <- subset(Object2, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 15)
png(paste0(figureDIR,samplename,"_QC_filtered.png"), width = 920, height = 320)
print(VlnPlot(Object2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
dev.off()


Object2 <- NormalizeData(Object2, normalization.method = "LogNormalize", scale.factor = 10000)
Object2 <- NormalizeData(Object2)
Object2 <- FindVariableFeatures(Object2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Object2)
Object2 <- ScaleData(Object2, features = all.genes)
Object2 <- RunPCA(Object2, features = VariableFeatures(object = Object2))


Object2 <- FindNeighbors(Object2, dims = 1:10)
Object2 <- FindClusters(Object2, resolution = 0.5)

Object2 <- RunUMAP(Object2, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
png(paste0(figureDIR,samplename,"_umap.png"), width = 480, height = 480)
print(DimPlot(Object2, reduction = "umap",label=T))
dev.off()
png(paste0(figureDIR,samplename,"_umap_s.png"), width = 480, height = 480)
print(DimPlot(Object2, reduction = "umap", label=T,group.by= "species"))
dev.off()


Object2 <- RunTSNE(Object2, dims = 1:10)


png(paste0(figureDIR,samplename,"_tsne.png"), width = 480, height = 480)
print(DimPlot(Object2, reduction = "tsne",label=T))
dev.off()
png(paste0(figureDIR,samplename,"_tsne_s.png"), width = 480, height = 480)
print(DimPlot(Object2, reduction = "tsne", label=T,group.by= "species"))
dev.off()

png(paste0(figureDIR,samplename,"_umap_CD79A.png"), width = 480, height = 480)
print(FeaturePlot(Object2, features = c("CD79A")))
dev.off()

png(paste0(figureDIR,samplename,"_tsne_CD79A.png"), width = 480, height = 480)
print(FeaturePlot(Object2, reduction = "tsne",features = c("CD79A")))
dev.off()

}

