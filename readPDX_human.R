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


