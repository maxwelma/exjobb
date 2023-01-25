
#Load data from files
source("functions.R")
library(Seurat)

#speed up sct transform with glmgampoi
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("glmGamPoi")

#Direcotry to matrixm feature and barcode files
matrix_dir ="c:\\Users\\matil\\OneDrive\\Dokument\\Exjobb\\" #path to sample folders
db<-F #sample without doublets present, doublets = False

#Load data from matrix to seurat object
alldata <- LoadAllData(matrix_dir,c("Sample_1","Sample_2","Sample_3","Sample_4"),db)

#Filter data to remove bad cells
filterdata <- FilterData(alldata)

#normalize and produce anchorset for integration
anchors <- ProduceAnchorsSCT(filterdata)

#integrate data using anchorset
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
#cluster data and produce umap vizualisation (due to the stochatic nature of umap and FindClusters the vizualizations differ between computers and versions of R)
combined.sct <- ClusterDataSCT(combined.sct,0.9,0.8) #save 90% of variance, set cluster resolution to 0.8

#uncomment to save file for future use or load file
#saveRDS(combined.sct, file = "combined_data_sct.RDS") 
combined.sct <- readRDS("combined_data_sct.RDS")

