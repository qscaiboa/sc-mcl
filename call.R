
library(BUSpaRse)
library(tidyverse)
library(DropletUtils)
library(Matrix)
library(Seurat)
source("barcodeRanks2.R")

res_mat1 <- read_count_output("hgmm_genecount",name = "genes", tcc = FALSE)
res_mat2 <- read_count_output("genecount",name = "genes", tcc = FALSE)

dim(res_mat1)
dim(res_mat2)
res_mat <- cbind(res_mat1, res_mat2)
dim(res_mat)

tot_counts <- Matrix::colSums(res_mat)
summary(tot_counts)

bc_rank <- barcodeRanks2(res_mat)
bc_rank
bc_rank@metadata





# Filter the matrix
res_mat <- res_mat[, tot_counts > bc_rank@metadata$inflection]
dim(res_mat)
 


gene_species <- ifelse(str_detect(rownames(res_mat), "^ENSMUSG"), "mouse", "human")
mouse_inds <- gene_species == "mouse"
human_inds <- gene_species == "human"
# mark cells as mouse or human
cell_species <- tibble(n_mouse_umi = Matrix::colSums(res_mat[mouse_inds,]),
                       n_human_umi = Matrix::colSums(res_mat[human_inds,]),
                       tot_umi = Matrix::colSums(res_mat),
                       prop_mouse = n_mouse_umi / tot_umi,
                       prop_human = n_human_umi / tot_umi)
		   
   


cell_species<-  cell_species %>% mutate(species = case_when(
    prop_mouse > 0.9 ~ "mouse",
    prop_human > 0.9 ~ "human",
    TRUE ~ "mixed"
  ))

table(cell_species$ground_truth, cell_species$species)


cell_species<- cell_species %>% 
  mutate(ground_truth = case_when(
    str_detect(colnames(res_mat), "^AAAAAAAA") ~ "mouse",
    str_detect(colnames(res_mat), "^TTTTTTTT") ~ "human",
    TRUE ~ "unknown"
  )) %>%
  mutate(species = case_when(
    prop_mouse > 0.9 ~ "mouse",
    prop_human > 0.9 ~ "human",
    TRUE ~ "mixed"
  ))

table(cell_species$ground_truth, cell_species$species)

cell_species[cell_species$ground_truth = "unknown",]


