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


Data_tmp<-Read10X(data.dir = "/home/qcai1/Downloads/scCancer-master/filtered_feature_bc_matrix_hm")
ObjectV1<- CreateSeuratObject(counts = Data_tmp, project = "hm", min.cells = 3, min.features=200)

res_mat <- as.matrix(ObjectV1@assays$RNA@data)

tot_counts <- Matrix::colSums(res_mat)
summary(tot_counts)

bc_rank <- barcodeRanks(res_mat)

qplot(bc_rank$total, bc_rank$rank, geom = "line") +
  geom_vline(xintercept = bc_rank@metadata$knee, color = "blue", linetype = 2) +
  geom_vline(xintercept = bc_rank@metadata$inflection, color = "green", linetype = 2) 


gene_species <- ifelse(str_detect(rownames(res_mat), "^mm10"), "mouse", "human")
mouse_inds <- gene_species == "mouse"
human_inds <- gene_species == "human"
# mark cells as mouse or human
cell_species <- tibble(n_mouse_umi = Matrix::colSums(res_mat[mouse_inds,]),
                       n_human_umi = Matrix::colSums(res_mat[human_inds,]),
                       tot_umi = Matrix::colSums(res_mat),
                       prop_mouse = n_mouse_umi / tot_umi,
                       prop_human = n_human_umi / tot_umi)



cell_species<- cell_species %>% 
  mutate(species = case_when(
    prop_mouse > 0.9 ~ "mouse",
    prop_human > 0.9 ~ "human",
    TRUE ~ "mixed"
  ))

table( cell_species$species)

cell_species$cellid <- colnames(res_mat) 
cell_species[, c("cellid","species")]

##########################
Data_tmp_h<-Read10X(data.dir = "/home/qcai1/Downloads/scCancer-master/filtered_feature_bc_matrix-h")
rownames(Data_tmp_h) <- toupper(rownames(Data_tmp_h))
ObjectH<- CreateSeuratObject(counts = Data_tmp_h, project = "huamn", min.cells = 3, min.features=200)

ObjectH@meta.data$cellid <- rownames(ObjectH@meta.data)
tmp.meta <- merge(ObjectH@meta.data, cell_species[, c("cellid","species")], by= "cellid")
rownames(tmp.meta)  <- tmp.meta$cellid
ObjectH@meta.data <- tmp.meta

ObjectH_h <- subset(ObjectH, subset = species %in% "human")
ObjectH_h[["percent.mt"]] <- PercentageFeatureSet(ObjectH_h, pattern = "MT-")
png("human_cell_to_human_genome.png", width = 960, height = 320)
VlnPlot(ObjectH_h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

png("mouse_cell_to_human_genome.png", width = 960, height = 320)
ObjectH_m <- subset(ObjectH, subset = species %in% "mouse")
ObjectH_m[["percent.mt"]] <- PercentageFeatureSet(ObjectH_h, pattern = "MT-")
VlnPlot(ObjectH_m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#####################################
Data_tmp_m<-Read10X(data.dir = "/home/qcai1/Downloads/scCancer-master/filtered_feature_bc_matrix-m")
rownames(Data_tmp_m) <- toupper(rownames(Data_tmp_m))
ObjectM<- CreateSeuratObject(counts = Data_tmp_m, project = "mouse", min.cells = 3, min.features=200)

ObjectM@meta.data$cellid <- rownames(ObjectM@meta.data)

tmp.meta <- merge(ObjectM@meta.data, cell_species[, c("cellid","species")], by= "cellid")
rownames(tmp.meta)  <- tmp.meta$cellid

ObjectM@meta.data <- tmp.meta

ObjectM_m <- subset(ObjectM, subset = species %in% "mouse")

ObjectM_m[["percent.mt"]]<- PercentageFeatureSet(ObjectM_m, pattern = "MT-")
png("mouse_cell_to_mouse_genome.png", width = 960, height = 320)
VlnPlot(ObjectM_m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


##################
Object1 <-merge(ObjectH_h, ObjectM_m)
Object2<- CreateSeuratObject(counts = as.matrix(Object1@assays$RNA@data), project = "sample_name", min.cells = 3, min.features=200)
Object2@meta.data$species <- Object1@meta.data$species 
Object2[["percent.mt"]]<- PercentageFeatureSet(Object2, pattern = "MT-")
VlnPlot(Object2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(Object2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Object2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


Object2 <- subset(Object2, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 15)
VlnPlot(Object2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Object2 <- NormalizeData(Object2, normalization.method = "LogNormalize", scale.factor = 10000)
Object2 <- NormalizeData(Object2)
Object2 <- FindVariableFeatures(Object2, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Object2), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Object2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(Object2)
Object2 <- ScaleData(Object2, features = all.genes)
Object2 <- RunPCA(Object2, features = VariableFeatures(object = Object2))
print(Object2[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Object2, dims = 1:2, reduction = "pca")

Object2 <- FindNeighbors(Object2, dims = 1:10)
Object2 <- FindClusters(Object2, resolution = 0.5)
Object2 <- RunUMAP(Object2, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Object2, reduction = "umap",label=T)
DimPlot(Object2, reduction = "umap", label=T,group.by= "species")

Object2 <- RunTSNE(Object2, dims = 1:10)
DimPlot(Object2, reduction = "tsne",label=T)
DimPlot(Object2, reduction = "tsne", label=T,group.by= "species")

FeaturePlot(Object2, features = c("CD79A"))
FeaturePlot(Object2, reduction = "tsne",features = c("CD79A"))


#########################################################################










plot1 <- FeatureScatter(Object1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Object1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


Object1 <- subset(Object1, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 15)
VlnPlot(Object1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Object1 <- NormalizeData(Object1, normalization.method = "LogNormalize", scale.factor = 10000)
Object1 <- NormalizeData(Object1)
Object1 <- FindVariableFeatures(Object1, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Object1), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Object1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(Object1)
Object1 <- ScaleData(Object1, features = all.genes)
Object1 <- RunPCA(Object1, features = VariableFeatures(object = Object1))
print(Object1[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Object1, dims = 1:2, reduction = "pca")

Object1 <- FindNeighbors(Object1, dims = 1:10)
Object1 <- FindClusters(Object1, resolution = 0.5)
Object1 <- RunUMAP(Object1, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Object1, reduction = "umap",label=T)
DimPlot(Object1, reduction = "umap", label=T,group.by= "species")

Object1 <- RunTSNE(Object1, dims = 1:10)
DimPlot(Object1, reduction = "tsne",label=T)
DimPlot(Object1, reduction = "tsne", label=T,group.by= "species")

FeaturePlot(Object1, features = c("CD79A"))
FeaturePlot(Object1, reduction = "tsne",features = c("CD79A"))


