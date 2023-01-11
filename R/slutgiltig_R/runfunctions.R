
#Load data from files
source("functions.R")
library(Seurat)

#speed up sct transform with glmgampoi
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("glmGamPoi")

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

#uncomment to save file for future use
#saveRDS(combined.sct, file = "combined_data_sct.RDS") 


#Extras 
#calculate counts of events per sample, for comparison
counts.by.ident <- data.frame(matrix(nrow=0,ncol=4))
colnames(counts.by.ident) <- c("Sample_1","Sample_2","Sample_3","Sample_4")
combined.sct.1 <- subset(combined.sct,orig.ident=="Sample_1")
combined.sct.2 <- subset(combined.sct,orig.ident=="Sample_2")
combined.sct.3 <- subset(combined.sct,orig.ident=="Sample_3")
combined.sct.4 <- subset(combined.sct,orig.ident=="Sample_4")
total_counts<-lapply(c(combined.sct.1,combined.sct.2,combined.sct.3,combined.sct.4), cellsSample) 

#Annotate cells and add to metadata
combined.sct <- annotateClusters(combined.sct)
DimPlot(combined.sct,group.by = "celltype")

#plot stacked bars of fractions of cell types, post annotation
png(file2,width=800, height=2000)
print(plotStackedBars(combined.sct))
dev.off()

