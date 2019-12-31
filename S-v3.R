#conda activate sc-tutorial-seurat-2
#R

library(Seurat)
library(cowplot)
library(ggplot2)
library(cowplot)

V0.data <- Read10X(data.dir = "/home/qcai1/Downloads/sc/mcl-sc/VJ_SC38/1135518-V0_analysis/outs/filtered_gene_bc_matrices/hg19/")
V0 <- CreateSeuratObject(counts = V0.data, project = "V0", min.cells = 3, min.features = 200)
V0
V0[["percent.mt"]] <- PercentageFeatureSet(V0, pattern = "^MT-")
VlnPlot(V0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
V0 <- subset(V0, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(V0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
V0 <- NormalizeData(V0)
V0 <- FindVariableFeatures(V0, selection.method = "vst", nfeatures= 2000)
V0 <- ScaleData(V0)
V0@meta.data$pt <- "V"
V0@meta.data$response <- "R"
V0@meta.data$tissue <- "BM"
V0@meta.data$sample <- "V0"
V0

V1.data <- Read10X(data.dir = "/home/qcai1/Downloads/sc/mcl-sc/DP_10XGenSC20/V1_aggr/outs/filtered_feature_bc_matrix")
V1 <- CreateSeuratObject(counts = V1.data, project = "V1", min.cells = 3, min.features = 200)
V1
V1[["percent.mt"]] <- PercentageFeatureSet(V1, pattern = "^MT-")
VlnPlot(V1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
V1 <- subset(V1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(V1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
V1 <- NormalizeData(V1)
V1 <- FindVariableFeatures(V1, selection.method = "vst", nfeatures= 2000)
V1 <- ScaleData(V1)
V1@meta.data$pt <- "V"
V1@meta.data$response <- "R"
V1@meta.data$tissue <- "BM"
V1@meta.data$sample <- "V1"
V1




V2.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/DP_10XGenSC20/V2-1138372_analysis/outs/filtered_gene_bc_matrices/hg19")
V2 <- CreateSeuratObject(counts = V2.data, project = "V2", min.cells = 3, min.features = 200)
V2
V2[["percent.mt"]] <- PercentageFeatureSet(V2, pattern = "^MT-")
VlnPlot(V2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
V2 <- subset(V2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(V2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
V2 <- NormalizeData(V2)
V2 <- FindVariableFeatures(V2, selection.method = "vst", nfeatures= 2000)
V2 <- ScaleData(V2)
V2@meta.data$pt <- "V"
V2@meta.data$response <- "R"
V2@meta.data$tissue <- "AP"
V2@meta.data$sample <- "V2"
V2



V3.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/DP_10XGenSC20/V3-1146484_analysis/outs/filtered_gene_bc_matrices/hg19/")
V3 <- CreateSeuratObject(counts = V3.data, project = "V3", min.cells = 3, min.features = 200)
V3
V3[["percent.mt"]] <- PercentageFeatureSet(V3, pattern = "^MT-")
VlnPlot(V3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
V3 <- subset(V3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(V3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
V3 <- NormalizeData(V3)
V3 <- FindVariableFeatures(V3, selection.method = "vst", nfeatures= 2000)
V3 <- ScaleData(V3)
V3@meta.data$pt <- "V"
V3@meta.data$response <- "R"
V3@meta.data$tissue <- "PB"
V3@meta.data$sample <- "V3"
V3



V4.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/DP_10XGenSC20/V4_aggr/outs/filtered_feature_bc_matrix")
V4 <- CreateSeuratObject(counts = V4.data, project = "V4", min.cells = 3, min.features = 200)
V4
V4[["percent.mt"]] <- PercentageFeatureSet(V4, pattern = "^MT-")
VlnPlot(V4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
V4 <- subset(V4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(V4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
V4 <- NormalizeData(V4)
V4 <- FindVariableFeatures(V4, selection.method = "vst", nfeatures= 2000)
V4 <- ScaleData(V4)
V4@meta.data$pt <- "V"
V4@meta.data$response <- "R"
V4@meta.data$tissue <- "PB"
V4@meta.data$sample <- "V4"
V4


V5.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ_SC38/1157356-V5_analysis/outs/filtered_gene_bc_matrices/hg19/")
V5 <- CreateSeuratObject(counts = V5.data, project = "V5", min.cells = 3, min.features = 200)
V5
V5[["percent.mt"]] <- PercentageFeatureSet(V5, pattern = "^MT-")
VlnPlot(V5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
V5 <- subset(V5, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(V5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
V5 <- NormalizeData(V5)
V5 <- FindVariableFeatures(V5, selection.method = "vst", nfeatures= 2000)
V5 <- ScaleData(V5)
V5@meta.data$pt <- "V"
V5@meta.data$response <- "R"
V5@meta.data$tissue <- "PB"
V5@meta.data$sample <- "V5"
V5


V6.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ_SC38/1165357-V6_analysis/outs/filtered_gene_bc_matrices/hg19/")
V6 <- CreateSeuratObject(counts = V6.data, project = "V6", min.cells = 3, min.features = 200)
V6
V6[["percent.mt"]] <- PercentageFeatureSet(V6, pattern = "^MT-")
VlnPlot(V6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
V6 <- subset(V6, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(V6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
V6 <- NormalizeData(V6)
V6 <- FindVariableFeatures(V6, selection.method = "vst", nfeatures= 2000)
V6 <- ScaleData(V6)
V6@meta.data$pt <- "V"
V6@meta.data$response <- "R"
V6@meta.data$tissue <- "PB"
V6@meta.data$sample <- "V6"
V6





C1.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ-RNA_10XGenSC68/1135981-C1_analysis/outs/filtered_gene_bc_matrices/hg19/")
C1 <- CreateSeuratObject(counts = C1.data, project = "C1", min.cells = 3, min.features = 200)
C1
C1[["percent.mt"]] <- PercentageFeatureSet(C1, pattern = "^MT-")
VlnPlot(C1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C1 <- subset(C1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(C1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C1 <- NormalizeData(C1)
C1 <- FindVariableFeatures(C1, selection.method = "vst", nfeatures= 2000)
C1 <- ScaleData(C1)
C1@meta.data$pt <- "C"
C1@meta.data$response <- "R"
C1@meta.data$tissue <- "BM"
C1@meta.data$sample <- "C1"
C1


C2.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ-RNA_10XGenSC68/1138530-C2_analysis/outs/filtered_gene_bc_matrices/hg19/")
C2 <- CreateSeuratObject(counts = C2.data, project = "C2", min.cells = 3, min.features = 200)
C2
C2[["percent.mt"]] <- PercentageFeatureSet(C2, pattern = "^MT-")
VlnPlot(C2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C2 <- subset(C2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(C2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C2 <- NormalizeData(C2)
C2 <- FindVariableFeatures(C2, selection.method = "vst", nfeatures= 2000)
C2 <- ScaleData(C2)
C2@meta.data$pt <- "C"
C2@meta.data$response <- "R"
C2@meta.data$tissue <- "AP"
C2@meta.data$sample <- "C2"
C2


C4.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ-RNA_10XGenSC68/1145995-C4_analysis/outs/filtered_gene_bc_matrices/hg19/")
C4 <- CreateSeuratObject(counts = C4.data, project = "C4", min.cells = 3, min.features = 200)
C4
C4[["percent.mt"]] <- PercentageFeatureSet(C4, pattern = "^MT-")
VlnPlot(C4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C4 <- subset(C4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(C4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C4 <- NormalizeData(C4)
C4 <- FindVariableFeatures(C4, selection.method = "vst", nfeatures= 2000)
C4 <- ScaleData(C4)
C4@meta.data$pt <- "C"
C4@meta.data$response <- "R"
C4@meta.data$tissue <- "PB"
C4@meta.data$sample <- "C4"
C4

C6.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ-RNA_10XGenSC68/1157929-C6_analysis/outs/filtered_gene_bc_matrices/hg19/")
C6 <- CreateSeuratObject(counts = C6.data, project = "C6", min.cells = 3, min.features = 200)
C6
C6[["percent.mt"]] <- PercentageFeatureSet(C6, pattern = "^MT-")
VlnPlot(C6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C6 <- subset(C6, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(C6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C6 <- NormalizeData(C6)
C6 <- FindVariableFeatures(C6, selection.method = "vst", nfeatures= 2000)
C6 <- ScaleData(C6)
C6@meta.data$pt <- "C"
C6@meta.data$response <- "R"
C6@meta.data$tissue <- "PB"
C6@meta.data$sample <- "C6"
C6


D1.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ-RNA2_10XGenSC68/1114173-D1_analysis/outs/filtered_gene_bc_matrices/hg19/")
D1 <- CreateSeuratObject(counts = D1.data, project = "D1", min.cells = 3, min.features = 200)
D1
D1[["percent.mt"]] <- PercentageFeatureSet(D1, pattern = "^MT-")
VlnPlot(D1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
D1 <- subset(D1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(D1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
D1 <- NormalizeData(D1)
D1 <- FindVariableFeatures(D1, selection.method = "vst", nfeatures= 2000)
D1 <- ScaleData(D1)
D1@meta.data$pt <- "D"
D1@meta.data$response <- "R"
D1@meta.data$tissue <- "PB"
D1@meta.data$sample <- "D1"
D1

D2.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ-RNA2_10XGenSC68/1115686-D2_analysis/outs/filtered_gene_bc_matrices/hg19/")
D2 <- CreateSeuratObject(counts = D2.data, project = "D2", min.cells = 3, min.features = 200)
D2
D2[["percent.mt"]] <- PercentageFeatureSet(D2, pattern = "^MT-")
VlnPlot(D2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
D2 <- subset(D2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(D2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
D2 <- NormalizeData(D2)
D2 <- FindVariableFeatures(D2, selection.method = "vst", nfeatures= 2000)
D2 <- ScaleData(D2)
D2@meta.data$pt <- "D"
D2@meta.data$response <- "R"
D2@meta.data$tissue <- "BM"
D2@meta.data$sample <- "D2"
D2


D4.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ-RNA2_10XGenSC68/1163458-D4_analysis/outs/filtered_gene_bc_matrices/hg19/")
D4 <- CreateSeuratObject(counts = D4.data, project = "D4", min.cells = 3, min.features = 200)
D4
D4[["percent.mt"]] <- PercentageFeatureSet(D4, pattern = "^MT-")
VlnPlot(D4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
D4 <- subset(D4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(D4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
D4 <- NormalizeData(D4)
D4 <- FindVariableFeatures(D4, selection.method = "vst", nfeatures= 2000)
D4 <- ScaleData(D4)
D4@meta.data$pt <- "D"
D4@meta.data$response <- "R"
D4@meta.data$tissue <- "BM"
D4@meta.data$sample <- "D4"
D4

D5.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ-RNA2_10XGenSC68/1163461-D5_analysis/outs/filtered_gene_bc_matrices/hg19/")
D5 <- CreateSeuratObject(counts = D5.data, project = "D5", min.cells = 3, min.features = 200)
D5
D5[["percent.mt"]] <- PercentageFeatureSet(D5, pattern = "^MT-")
VlnPlot(D5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
D5 <- subset(D5, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(D5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
D5 <- NormalizeData(D5)
D5 <- FindVariableFeatures(D5, selection.method = "vst", nfeatures= 2000)
D5 <- ScaleData(D5)
D5@meta.data$pt <- "D"
D5@meta.data$response <- "R"
D5@meta.data$tissue <- "PB"
D5@meta.data$sample <- "D5"
D5

B0.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ_SC38/964301-B0_analysis/outs/filtered_gene_bc_matrices/hg19/")
B0 <- CreateSeuratObject(counts = B0.data, project = "B0", min.cells = 3, min.features = 200)
B0
B0[["percent.mt"]] <- PercentageFeatureSet(B0, pattern = "^MT-")
VlnPlot(B0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
B0 <- subset(B0, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(B0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
B0 <- NormalizeData(B0)
B0 <- FindVariableFeatures(B0, selection.method = "vst", nfeatures= 2000)
B0 <- ScaleData(B0)
B0@meta.data$pt <- "B"
B0@meta.data$response <- "NR"
B0@meta.data$tissue <- "BM"
B0@meta.data$sample <- "B0"
B0

B1.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ_SC38/964303-B1_analysis/outs/filtered_gene_bc_matrices/hg19/")
B1 <- CreateSeuratObject(counts = B1.data, project = "B1", min.cells = 3, min.features = 200)
B1
B1[["percent.mt"]] <- PercentageFeatureSet(B1, pattern = "^MT-")
VlnPlot(B1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
B1 <- subset(B1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(B1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
B1 <- NormalizeData(B1)
B1 <- FindVariableFeatures(B1, selection.method = "vst", nfeatures= 2000)
B1 <- ScaleData(B1)
B1@meta.data$pt <- "B"
B1@meta.data$response <- "NR"
B1@meta.data$tissue <- "PB"
B1@meta.data$sample <- "B1"
B1

B4.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ_SC38/1106365-B4_analysis/outs/filtered_gene_bc_matrices/hg19/")
B4 <- CreateSeuratObject(counts = B4.data, project = "B4", min.cells = 3, min.features = 200)
B4
B4[["percent.mt"]] <- PercentageFeatureSet(B4, pattern = "^MT-")
VlnPlot(B4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
B4 <- subset(B4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(B4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
B4 <- NormalizeData(B4)
B4 <- FindVariableFeatures(B4, selection.method = "vst", nfeatures= 2000)
B4 <- ScaleData(B4)
B4@meta.data$pt <- "B"
B4@meta.data$response <- "NR"
B4@meta.data$tissue <- "AP"
B4@meta.data$sample <- "B4"
B4


E1.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ-RNA3_10XGenSC68/1008320-E1_analysis/outs/filtered_gene_bc_matrices/hg19/")
E1 <- CreateSeuratObject(counts = E1.data, project = "E1", min.cells = 3, min.features = 200)
E1
E1[["percent.mt"]] <- PercentageFeatureSet(E1, pattern = "^MT-")
VlnPlot(E1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E1 <- subset(E1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(E1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E1 <- NormalizeData(E1)
E1 <- FindVariableFeatures(E1, selection.method = "vst", nfeatures= 2000)
E1 <- ScaleData(E1)
E1@meta.data$pt <- "E"
E1@meta.data$response <- "NR"
E1@meta.data$tissue <- "PB"
E1@meta.data$sample <- "E1"
E1

E2.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ-RNA3_10XGenSC68/1024089-E2_analysis/outs/filtered_gene_bc_matrices/hg19/")
E2 <- CreateSeuratObject(counts = E2.data, project = "E2", min.cells = 3, min.features = 200)
E2
E2[["percent.mt"]] <- PercentageFeatureSet(E2, pattern = "^MT-")
VlnPlot(E2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E2 <- subset(E2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(E2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E2 <- NormalizeData(E2)
E2 <- FindVariableFeatures(E2, selection.method = "vst", nfeatures= 2000)
E2 <- ScaleData(E2)
E2@meta.data$pt <- "E"
E2@meta.data$response <- "NR"
E2@meta.data$tissue <- "PB"
E2@meta.data$sample <- "E2"
E2

E3.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ-RNA3_10XGenSC68/1034264-E3_analysis/outs/filtered_gene_bc_matrices/hg19/")
E3 <- CreateSeuratObject(counts = E3.data, project = "E3", min.cells = 3, min.features = 200)
E3
E3[["percent.mt"]] <- PercentageFeatureSet(E3, pattern = "^MT-")
VlnPlot(E3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E3 <- subset(E3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(E3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E3 <- NormalizeData(E3)
E3 <- FindVariableFeatures(E3, selection.method = "vst", nfeatures= 2000)
E3 <- ScaleData(E3)
E3@meta.data$pt <- "E"
E3@meta.data$response <- "NR"
E3@meta.data$tissue <- "PB"
E3@meta.data$sample <- "E3"
E3

N1.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ-RNA3_10XGenSC68/PBMC-1_analysis/outs/filtered_gene_bc_matrices/hg19/")
N1 <- CreateSeuratObject(counts = N1.data, project = "N1", min.cells = 3, min.features = 200)
N1
N1[["percent.mt"]] <- PercentageFeatureSet(N1, pattern = "^MT-")
VlnPlot(N1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
N1 <- subset(N1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(N1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
N1 <- NormalizeData(N1)
N1 <- FindVariableFeatures(N1, selection.method = "vst", nfeatures= 2000)
N1 <- ScaleData(N1)
N1@meta.data$pt <- "N"
N1@meta.data$response <- "NA"
N1@meta.data$tissue <- "PB"
N1@meta.data$sample <- "N1"
N1



N2.data <- Read10X(data.dir =  "/home/qcai1/Downloads/sc/mcl-sc/VJ-RNA3_10XGenSC68/PBMC-2_analysis/outs/filtered_gene_bc_matrices/hg19/")
N2 <- CreateSeuratObject(counts = N2.data, project = "N2", min.cells = 3, min.features = 200)
N2
N2[["percent.mt"]] <- PercentageFeatureSet(N2, pattern = "^MT-")
VlnPlot(N2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
N2 <- subset(N2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#VlnPlot(N2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
N2 <- NormalizeData(N2)
N2 <- FindVariableFeatures(N2, selection.method = "vst", nfeatures= 2000)
N2 <- ScaleData(N2)
N2@meta.data$pt <- "N"
N2@meta.data$response <- "NA"
N2@meta.data$tissue <- "PB"
N2@meta.data$sample <- "N2"
N2

pbmc.big <- merge(V0, y = c(V1, V2, V3, V4, V5, V6, N1, N2, C1, C2, C4, D1, D2, D4, D5, B0, B1, B4, E1, E2, E3), 
	add.cell.ids = c("V0", "V1", "V2", "V3", "V4", "V5", "V6", "N1", "N2", "C1", "C2", "C4", "D1", "D2", "D4", "D5", "B0", "B1", "B4", "E1", "E2", "E3"), project = "MCL-sc")

pancreas.list <- SplitObject(pbmc.big, split.by = "sample")
saveRDS(pancreas.list, file = "0.mcl.list.rds")



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
library(Seurat)
library(cowplot)
library(ggplot2)
pancreas.list  <- readRDS("0.mcl.list.rds")

for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
    pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
}

pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, dims = 1:30)
rm(pancreas.list)
saveRDS(pancreas.anchors, file= "1.mcl.integrated.rds")



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
library(Seurat)
library(cowplot)
library(ggplot2)
pancreas.anchors  <- readRDS("1.mcl.integrated.rds")


pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(pancreas.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)

saveRDS(pancreas.integrated, file= "2.mcl.integrated.rds")


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
library(Seurat)
library(cowplot)
library(ggplot2)
pancreas.integrated <- readRDS("2.mcl.integrated.rds")
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
pancreas.integrated <- FindNeighbors(pancreas.integrated, reduction = "pca", dims = 1:20)
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.5)
saveRDS(pancreas.integrated, file= "3.mcl.integrated.rds")

p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "pt")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "response")
p3 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "sample")
p4 <- DimPlot(pancreas.integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2, p3, p4)




pancreas.integrated <- RunTSNE(pancreas.integrated, reduction = "pca", dims = 1:30)

p5 <- DimPlot(pancreas.integrated, reduction = "tsne", group.by = "pt")
p6 <- DimPlot(pancreas.integrated, reduction = "tsne", group.by = "response")
p7 <- DimPlot(pancreas.integrated, reduction = "tsne", group.by = "sample")
p8 <- DimPlot(pancreas.integrated, reduction = "tsne", label = TRUE)
plot_grid(p5, p6, p7, p8)


saveRDS(pancreas.integrated, file= "3.mcl.integrated.rds")




library(Seurat)
library(cowplot)
library(ggplot2)
pancreas.integrated <- readRDS("3.mcl.integrated.rds")
# Visualization
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "pt")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "response")
p3 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "sample")
p4 <- DimPlot(pancreas.integrated, reduction = "umap", label = TRUE)

p5 <- DimPlot(pancreas.integrated, reduction = "tsne", group.by = "pt")
p6 <- DimPlot(pancreas.integrated, reduction = "tsne", group.by = "response")
p7 <- DimPlot(pancreas.integrated, reduction = "tsne", group.by = "sample")
p8 <- DimPlot(pancreas.integrated, reduction = "tsne", label = TRUE)


plot_grid(p1, p2, p3, p4)
plot_grid(p5, p6, p7, p8)


#################
slist <- SplitObject(pancreas.integrated, split.by = "sample")
sl <-  c("V0", "V1", "V2", "V3", "V4", "V5", "V6", "N1", "N2", "C1", "C2", "C4", "D1", "D2", "D4", "D5", "B0", "B1", "B4", "E1", "E2", "E3")
m <- c("$`0`", "$`1`", "$`2`", "$`3`", "$`4`", "$`5`", "$`6`", "$`7`", "$`8`", "$`9`", "$`10`")
for  ( i in 1:22){ 
	print(sl[i])
	for ( j in 1:11){ 
print(m[j], ncol(eval(parse(text =  paste("SplitObject(slist$", sl[i], ", split.by = 'ident')", m[j], sep = "") )))  /  ncol(eval(parse(text = paste("slist$", sl[i], sep ="") ))))
}
}
####################



DimPlot(pancreas.integrated, reduction = "tsne", split.by = "response", label = TRUE)
DimPlot(pancreas.integrated, reduction = "tsne", split.by = "sample", label = TRUE)
DimPlot(pancreas.integrated, reduction = "tsne", split.by = "pt", label = TRUE)



DimPlot(pancreas.integrated, reduction = "umap", split.by = "pt")
DimPlot(pancreas.integrated, label = TRUE)
DimPlot(pancreas.integrated, reduction = "umap", split.by = "pt", label = TRUE)
DimPlot(pancreas.integrated, reduction = "umap", split.by = "response", label = TRUE)


markers.0 <- FindConservedMarkers(pancreas.integrated, ident.1 = 0, grouping.var = "response", verbose = FALSE)
head(markers.0)

markers.1 <- FindConservedMarkers(pancreas.integrated, ident.1 = 1, grouping.var = "response", verbose = FALSE)
head(markers.1)
rownames(markers.1)[1:12]
"CD74" "HLA-DRB1" "HLA-DRA"  "HLA-DPA1" "HLA-DQA1" "HLA-DQB1" "HLA-DMA"  "CD79B"  "SPIB"     "HLA-DPB1" "HLA-DMB"  "BANK1"


FeaturePlot(pancreas.integrated, reduction = "tsne", features = rownames(markers.1)[1:12], min.cutoff = "q9")





pancreas.integrated <- RenameIdents(pancreas.integrated, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
    `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", `7` = "NK", `8` = "DC", `9` = "B Activated", 
    `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")





#####################

slist <- SplitObject(pancreas.integrated, split.by = "sample")

a <- c()

sl <- c('V0', 'V1')
eval(parse(text = paste("slist$", sl[1], sep ="") ))

ncol(eval(parse(text = paste("slist$", sl[1], sep ="") )))

SplitObject(  eval(parse(text = paste("slist$", sl[1], sep =""))), split.by = "ident")

m <- c("$`0`", "$`1`")


ncol(eval(parse(text =  paste("SplitObject(slist$V0, split.by = 'ident')", m[1], sep = "") )))  /  ncol(eval(parse(text = paste("slist$", sl[1], sep ="") )))
ncol(eval(parse(text =  paste("SplitObject(slist$", sl[1], ", split.by = 'ident')", m[1], sep = "") )))  /  ncol(eval(parse(text = paste("slist$", sl[1], sep ="") )))
###################


