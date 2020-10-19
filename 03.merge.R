#bsub -Is -q interactive -W 6:00 -M 200 -R rusage[mem=200] -n 1 /bin/bash
#R

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
samplePDX <- data.frame(samplePDX)

allSampleName = list.files(dataDir)
#allSampleName <- allSampleName[!allSampleName %in% "scripts"]
allSampleName  <- gsub(".RDS","",allSampleName)

initialSampleName = "PBMC-1"
allSampleName = allSampleName[!allSampleName %in% initialSampleName]


i <- initialSampleName
Object1 <- readRDS(paste0("./Data.RDS/", i , ".RDS"))

genelist <- rownames(Object1$RNA@counts)
length(genelist)
for ( i in allSampleName) {
    print(i)
    Object32 <- readRDS(paste0("./Data.RDS/", i , ".RDS"))
	genelist <- c(genelist, rownames(Object32$RNA@counts))
	#print(length(genelist))
	genelist <- unique(genelist)
	print(length(genelist))
	}


 all_genelist <- genelist
 length(all_genelist)
 
###################################################
 
 
allSampleName = list.files(dataDir)
#allSampleName <- allSampleName[!allSampleName %in% "scripts"]
allSampleName  <- gsub(".RDS","",allSampleName)

initialSampleName = "PBMC-1"
allSampleName = allSampleName[!allSampleName %in% initialSampleName]


i <- initialSampleName
Object1 <- readRDS(paste0("./Data.RDS/", i , ".RDS"))

 genelist_data <- rownames(Object1$RNA@counts)
 cell_list <- colnames(Object1$RNA@counts)
 cell_length <- ncol(Object1$RNA@counts)
 MM <- data.frame(Object1$RNA@counts)
 MM$AAA <- 100
 genelist_add <- all_genelist[!all_genelist  %in% genelist_data]
 m <- data.frame(matrix(0, nrow = length(genelist_add), ncol = (cell_length+1)))
 rownames(m) <-  genelist_add
 colnames(m) <-  c(cell_list, "AAA")
 
m$AAA <- 100 #set high number for CreateSeuratObject, will remove this later


  
colnames(m)  <- gsub(".1","",colnames(m))
colnames(m)  <- gsub("-1","",colnames(m))   
colnames(MM)  <- gsub(".1","",colnames(MM))
colnames(MM)  <- gsub("-1","",colnames(MM)) 

m2 <- rbind(MM, m)

m2 <- CreateSeuratObject(counts = m2, project = "temp")
m2 <- subset(x = m2, subset = colnames %in% "AAA")
m2 <- m2[,1:length(cell_list)]
m2@meta.data <- Object1@meta.data
Object1 <- m2

Object1@meta.data$Pt <-    samplePDX[samplePDX$scSeqName %in% i,]$Pt   ######################################
Object1@meta.data$sample <-    samplePDX[samplePDX$scSeqName %in% i,]$scSeqName
Object1@meta.data$Finalname <-    samplePDX[samplePDX$scSeqName %in% i,]$Finalname
Object1@meta.data$HumanMouse <-    samplePDX[samplePDX$scSeqName %in% i,]$HumanMouse
Object1@meta.data$PDXmodel <-    samplePDX[samplePDX$scSeqName %in% i,]$PDXmodel
Object1 <- NormalizeData(Object1)
Object1 <- FindVariableFeatures(Object1, selection.method = "vst", nfeatures= 2000)
Object1 <- ScaleData(Object1)
Object1[["percent.mt"]] <- PercentageFeatureSet(Object1, pattern = "MT-")
Object1 <- subset(Object1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
print(OBject1)
head(Object1@meta.data)
	rm(m)
	rm(MM)
	rm(m2)
###############################################

Object.combined = c()
for ( i in allSampleName ) {
    print(i)
    j = i

 Object32 <- readRDS(paste0("./Data.RDS/", i , ".RDS"))
 genelist_data <- rownames(Object32$RNA@counts)
 cell_list <- colnames(Object32$RNA@counts)
 cell_length <- ncol(Object32$RNA@counts)
 MM <- data.frame(Object32$RNA@counts)
 MM$AAA <- 100
 genelist_add <- all_genelist[!all_genelist  %in% genelist_data]
 m <- data.frame(matrix(0, nrow = length(genelist_add), ncol = (cell_length+1)))
 rownames(m) <-  genelist_add
 colnames(m) <-  c(cell_list, "AAA")

m$AAA <- 100 #set high number for CreateSeuratObject, will remove this later

colnames(m)  <- gsub(".1","",colnames(m))
colnames(m)  <- gsub("-1","",colnames(m))   
colnames(MM)  <- gsub(".1","",colnames(MM))
colnames(MM)  <- gsub("-1","",colnames(MM)) 

m2 <- rbind(MM, m)

m2 <- CreateSeuratObject(counts = m2, project = "temp")
 print("T1")	
m2 <- m2[,1:length(cell_list)]
m2@meta.data <- Object32@meta.data
Object32 <- m2
 print("T2")	
    
Object32@meta.data$Pt <-    samplePDX[samplePDX$scSeqName %in% i,]$Pt   ######################################
Object32@meta.data$sample <-    samplePDX[samplePDX$scSeqName %in% i,]$scSeqName
Object32@meta.data$Finalname <-    samplePDX[samplePDX$scSeqName %in% i,]$Finalname 
Object32@meta.data$HumanMouse <-    samplePDX[samplePDX$scSeqName %in% i,]$HumanMouse
Object32@meta.data$PDXmodel <-    samplePDX[samplePDX$scSeqName %in% i,]$PDXmodel
Object32[["percent.mt"]] <- PercentageFeatureSet(Object32, pattern = "^MT-")

Object32 <- NormalizeData(Object32)
Object32 <- FindVariableFeatures(Object32, selection.method = "vst", nfeatures= 2000)
Object32 <- ScaleData(Object32)
Object32 <- subset(Object32, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
print(Object32)
head(Object32@meta.data)   
   Object.combined = c(Object.combined, Object32)
	rm(m)
	rm(MM)
	rm(m2)
  }
  
  
  
  

 Objectall <- merge(Object1, y =   Object.combined, add.cell.ids =c(initialSampleName, allSampleName)  , project = "SC145")
 Objectall <- NormalizeData(object = Objectall, verbose = FALSE)
 Objectall <- FindVariableFeatures(object = Objectall, selection.method = "vst", nfeatures = 2000)
 Objectall <- ScaleData(Objectall,verbose = FALSE)
 Objectall <- RunPCA(Objectall, features = VariableFeatures(object = Objectall))
 Objectall <- RunUMAP(Objectall, dims = 1:30)
 Objectall <- RunTSNE(Objectall, dims = 1:30)
  saveRDS(Objectall, "./Data/Object.merged.rds")
  
 png("ZZZ.01.umap.png", width = 1200, height = 1200)
 DimPlot(Objectall, reduction = "umap", label = TRUE, group.by= "species")
 dev.off()
 
 
 png("ZZZ.01.umap.orig.ident.png", width = 2400, height = 1200)
 DimPlot(Objectall, reduction = "umap", label = TRUE, group.by= "orig.ident")
 dev.off()
 
 
 
 png("ZZZ.01.umap.pt.png", width = 2400, height = 1200)
 DimPlot(Objectall, reduction = "umap", label = TRUE, group.by= "Pt")
 dev.off()
 
 png("ZZZ.01.FeaturePlot.png", width = 2400, height = 2400)
 FeaturePlot(Objectall, features = c("CD79A", "CD79B","NKG7", "GNLY","S100A9", "HBB", "ZNF90","MCL1",
                                    "HIST1H4C", "TUBA1B", "TUBB", "PCLAF","TYMS"),reduction = "umap", min.cutoff = "q9")
dev.off()
 
 
 orig.ident
 
 saveRDS(Objectall, "./Data/Object.merged.rds")
