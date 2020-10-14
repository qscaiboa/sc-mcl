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
source("readPDX.R")
dataDir  =  "/home/qcai1/Downloads/scCancer-master/Data/"
figureDIR =  "/home/qcai1/Downloads/scCancer-master/Figures/"
workDIR = "/home/qcai1/Downloads/scCancer-master/sc3MAP/10Xhm/" 
workDIRh = "/home/qcai1/Downloads/scCancer-master/sc3MAP/10Xhuman/"
workDIRm = "/home/qcai1/Downloads/scCancer-master/sc3MAP/10Xmouse/"
allSampleName = list.files(workDIR)

initialSampleName = "347795-Tumor"

readPDX(samplename=initialSampleName, hmDIR=workDIR , humanDIR=workDIRh, mouseDIR=workDIRm,  figureDIR= figureDIR, dataDir= dataDir)

