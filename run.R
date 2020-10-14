

library(Seurat)
library(ggplot2)
library(readr)
library(readxl)
library(Matrix)
library(tidyr)
library(dplyr)
library(readxl)
library(cowplot)
library(BUSpaRse)
library(tidyverse)
library(DropletUtils)
source("readPDX2.R")

dataDir  =  "/scratch/03988/qscai/hgmm/Data/"
figureDIR =  "/scratch/03988/qscai/hgmm/Figures/"
workDIR = "/scratch/03988/qscai/hgmm/sc3MAP/10Xhm/" 
workDIRh = "/scratch/03988/qscai/hgmm/sc3MAP/10Xhuman/"
workDIRm = "/scratch/03988/qscai/hgmm/sc3MAP/10Xmouse/"


allSampleName = list.files(workDIR)
allSampleName

for (i in allSampleName) {
print(i)
readPDX(samplename=i, hmDIR=workDIR , humanDIR=workDIRh, mouseDIR=workDIRm,  figureDIR= figureDIR, dataDir= dataDir)
list.files("./Data")
}



readPDX_human <- function(samplename, humanDIR, dataDir) {
  fileDIRh = paste(workDIRh, samplename, "/outs/filtered_feature_bc_matrix", sep="")
  Data_tmp_h<-Read10X(data.dir = fileDIRh )
  rownames(Data_tmp_h) <- toupper(rownames(Data_tmp_h))
  ObjectH<- CreateSeuratObject(counts = Data_tmp_h, project = samplename, min.cells = 3, min.features=200)
  ObjectH@meta.data$cellid <- rownames(ObjectH@meta.data)
  ObjectH@meta.data$species <- "human" 
  Object2 <- ObjectH

  saveRDS(Object2, file = paste0(dataDir,samplename,".RDS"))
  print(Object2)
  }




source("readPDX_human.R")
readPDX_human(samplename="PBMC-1", humanDIR=workDIRh,   dataDir= dataDir)
readPDX_human(samplename="PBMC-2", humanDIR=workDIRh,   dataDir= dataDir)
