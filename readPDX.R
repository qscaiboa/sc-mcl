readPDX <- function(samplename, hmDIR, humanDIR, mouseDIR, figureDIR, dataDir) {

  fileDIR_hm = paste(hmDIR, samplename, "/outs/filtered_feature_bc_matrix", sep="")  #S10
  Data_tmp<-Read10X(data.dir = fileDIR_hm)
  ObjectV1 <- CreateSeuratObject(counts = Data_tmp, project = paste0(samplename,"_hm"), min.cells = 3, min.features=200)
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
print(samplename)
print(table( cell_species$species))
  cell_species$cellid <- colnames(res_mat)
  human_cellid <- cell_species[cell_species$species %in% "human",]$cellid
  mouse_cellid <- cell_species[cell_species$species %in% "mouse",]$cellid

print(nrow(data.frame(table( cell_species$species))))


  ##########################
  fileDIRh = paste(workDIRh, samplename, "/outs/filtered_feature_bc_matrix", sep="")
  Data_tmp_h<-Read10X(data.dir = fileDIRh )
  rownames(Data_tmp_h) <- toupper(rownames(Data_tmp_h))
  ObjectH<- CreateSeuratObject(counts = Data_tmp_h, project = samplename, min.cells = 3, min.features=200)

  ObjectH@meta.data$cellid <- rownames(ObjectH@meta.data)


  ObjectH_h <- subset(ObjectH, cellid %in% human_cellid)
  ObjectH_h@meta.data$species <- "human"
  ObjectH_h[["percent.mt"]] <- PercentageFeatureSet(ObjectH_h, pattern = "MT-")
  png(paste0(figureDIR,samplename,"_human_cell_to_human_genome.png"), width = 920, height = 320)
  VlnPlot(ObjectH_h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  dev.off()

  #
  # ObjectH_m <- subset(ObjectH, cellid %in% mouse_cellid)
  # ObjectH_m
  # ObjectH_m[["percent.mt"]] <- PercentageFeatureSet(ObjectH_h, pattern = "MT-")
  # png(paste0(figureDIR,samplename,"_mouse_cell_to_human_genome.png"), width = 960, height = 320)
  # VlnPlot(ObjectH_m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # dev.off()
  #######################################


a <- data.frame(t(as.matrix(table( cell_species$species))) )
if(!"mixed" %in% colnames(a))
 {
   a$mixed <- 0
 }

if(!"mouse" %in% colnames(a))
 {
   a$mouse <- c("mouse",0)
 }


print(a)
if  (a$mouse> 1) {

  fileDIRm = paste(workDIRm, samplename, "/outs/filtered_feature_bc_matrix", sep="")
  Data_tmp_m<-Read10X(data.dir = fileDIRh )
  rownames(Data_tmp_m) <- toupper(rownames(Data_tmp_m))
  ObjectM<- CreateSeuratObject(counts = Data_tmp_m, project = samplename, min.cells = 3, min.features=200)

  ObjectM@meta.data$cellid <- rownames(ObjectM@meta.data)


  ObjectM_m <- subset(ObjectH, cellid %in% mouse_cellid)
  ObjectM_m@meta.data$species <- "mouse"
  ObjectM_m[["percent.mt"]] <- PercentageFeatureSet(ObjectM_m, pattern = "MT-")
  png(paste0(figureDIR,samplename,"mouse_cell_to_mouse_genome.png"), width = 920, height = 320)
    VlnPlot(ObjectM_m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  dev.off()


  # ObjectM_h <- subset(ObjectH, cellid %in% human_cellid)
  # ObjectM_h
  # ObjectM_h[["percent.mt"]] <- PercentageFeatureSet(ObjectM_h, pattern = "MT-")
  # png("huamn_cell_to_mouse_genome.png", width = 960, height = 320)
  # VlnPlot(ObjectM_h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # dev.off()
  #########################

  Object1 <-merge(ObjectH_h, ObjectM_m)
  Object2<- CreateSeuratObject(counts = as.matrix(Object1@assays$RNA@data), project = samplename, min.cells = 3, min.features=200)
  Object2
  mdata <- rbind(ObjectH_h@meta.data, ObjectM_m@meta.data)
  dim(mdata)

  Object2@meta.data <- mdata
}

else {
Object2 <- ObjectH_h

}


print(head(Object2@meta.data))

  Object2[["percent.mt"]]<- PercentageFeatureSet(Object2, pattern = "MT-")
  png(paste0(figureDIR,samplename,"_QC.png"), width = 320, height = 320)
  VlnPlot(Object2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  dev.off()



  saveRDS(Object2, file = paste0(dataDir,samplename,".RDS"))
  print(Object2)
  }
