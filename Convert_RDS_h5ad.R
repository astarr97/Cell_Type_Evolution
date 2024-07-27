library(Seurat)
library(SeuratObject)
library(SeuratDisk)

### USED TO CONVERT RDS TO H5AD FOR USE IN SCANPY ####

s <- readRDS("Human_M1_10xV3_Matrix.RDS")
s <- CreateSeuratObject(s)
SaveH5Seurat(s, "Human_M1_10xV3_Matrix.h5seurat")
Convert("Human_M1_10xV3_Matrix.h5seurat", dest="Human_M1_10xV3_Matrix.h5ad")

s <- readRDS("Marmoset_M1_10xV3_Matrix.RDS")
s <- CreateSeuratObject(s)
SaveH5Seurat(s, "Marmoset_M1_10xV3_Matrix.h5seurat")
Convert("Marmoset_M1_10xV3_Matrix.h5seurat", dest="Marmoset_M1_10xV3_Matrix.h5ad")

s <- readRDS("Mouse_M1_10xV3_Matrix.RDS")
s <- CreateSeuratObject(s)
SaveH5Seurat(s, "Mouse_M1_10xV3_Matrix.h5seurat")
Convert("Mouse_M1_10xV3_Matrix.h5seurat", dest="Mouse_M1_10xV3_Matrix.h5ad")

s <- readRDS("MTG_Data/rhesus_SCT_UMI_expression_matrix.RDS")
#s <- CreateSeuratObject(s)
SaveH5Seurat(s, "rhesus_SCT_UMI_expression_matrix.h5seurat")
Convert("rhesus_SCT_UMI_expression_matrix.h5seurat", dest="rhesus_SCT_UMI_expression_matrix.h5ad")
s <- 0
s <- readRDS("MTG_Data/Rhesus_Master_metadata_for_plots_and_sharing_12_16_21.RDS")
write.table(s, "Rhesus_Master_metadata.txt", sep = "\t", row.names=FALSE, quote = FALSE)

s <- readRDS("MTG_Data/chimp_SCT_UMI_expression_matrix.RDS")
#s <- CreateSeuratObject(s)
SaveH5Seurat(s, "chimp_SCT_UMI_expression_matrix.h5seurat")
Convert("chimp_SCT_UMI_expression_matrix.h5seurat", dest="chimp_SCT_UMI_expression_matrix.h5ad")
s <- 0
s <- readRDS("MTG_Data/Chimp_Master_metadata_for_plots_and_sharing_12_16_21.RDS")
write.table(s, "Chimp_Master_metadata.txt", sep = "\t", row.names=FALSE, quote = FALSE)

s <- readRDS("MTG_Data/gorilla_SCT_UMI_expression_matrix.RDS")
#s <- CreateSeuratObject(s)
SaveH5Seurat(s, "gorilla_SCT_UMI_expression_matrix.h5seurat")
Convert("gorilla_SCT_UMI_expression_matrix.h5seurat", dest="gorilla_SCT_UMI_expression_matrix.h5ad")
s <- 0
s <- readRDS("MTG_Data/Gorilla_Master_metadata_for_plots_and_sharing_12_16_21.RDS")
write.table(s, "Gorilla_Master_metadata.txt", sep = "\t", row.names=FALSE, quote = FALSE)

s <- readRDS("MTG_Data/human_SCT_UMI_expression_matrix.RDS")
#s <- CreateSeuratObject(s)
SaveH5Seurat(s, "human_SCT_UMI_expression_matrix.h5seurat")
Convert("human_SCT_UMI_expression_matrix.h5seurat", dest="human_SCT_UMI_expression_matrix.h5ad")
s <- 0
s <- readRDS("MTG_Data/Human_Master_metadata_for_plots_and_sharing_12_16_21.RDS")
write.table(s, "Human_Master_metadata.txt", sep = "\t", row.names=FALSE, quote = FALSE)

s <- readRDS("MTG_Data/marmoset_SCT_UMI_expression_matrix.RDS")
#s <- CreateSeuratObject(s)
SaveH5Seurat(s, "marmoset_SCT_UMI_expression_matrix.h5seurat")
Convert("marmoset_SCT_UMI_expression_matrix.h5seurat", dest="marmoset_SCT_UMI_expression_matrix.h5ad")
s <- 0
s <- readRDS("MTG_Data/Marmoset_Master_metadata_for_plots_and_sharing_12_16_21.RDS")
write.table(s, "Marmoset_Master_metadata.txt", sep = "\t", row.names=FALSE, quote = FALSE)
