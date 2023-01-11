source("functions.R")
library(Seurat)
library(ggplot2)
require(gridExtra)
library(patchwork)
#samp <- Seurat::ReadMtx(mtx="c:\\Users\\matil\\OneDrive\\Dokument\\Exjobb\\Sample_1\\no_doublets\\matrix_no_doublets.mtx",cells="c:\\Users\\matil\\OneDrive\\Dokument\\Exjobb\\Sample_1\\no_doublets\\barcodes_no_doublets.tsv",features="c:\\Users\\matil\\OneDrive\\Dokument\\Exjobb\\Sample_1\\no_doublets\\features.tsv")


matrix_dir ="c:\\Users\\matil\\OneDrive\\Dokument\\Exjobb\\"
db<-F
alldata <- LoadAllData(matrix_dir,c("Sample_1","Sample_2","Sample_3","Sample_4"),db) #sample without doublets present
alldata_with_dbs <- LoadAllData(matrix_dir,c("Sample_1","Sample_2","Sample_3","Sample_4"),T) #sample with doublets present

filterdata <- FilterData(alldata)

#alldata <- PercentageFeatureSet(alldata, "^MT", col.name = "percent_mito") #final values fall between 0 - 100

'%nin%' <- Negate('%in%')
dbs <- WhichCells(alldata_with_dbs)[WhichCells(alldata_with_dbs)%nin%WhichCells(alldata)]
dbs <- subset(alldata_with_dbs, cells = dbs)

#violin plot of all 3
feats <- c("nFeature_RNA", "nCount_RNA","percent_mito","percent_hb")

VlnPlot(alldata_with_dbs, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 4)+
  NoLegend()+plot_annotation("Before filtering and doublet removal")& 
  theme(plot.title = element_text(hjust = 0.5))

VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 4)+
  NoLegend()+plot_annotation("After doublet removal")& 
  theme(plot.title = element_text(hjust = 0.5))

VlnPlot(filterdata, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 4)+
  NoLegend()+ plot_annotation("After filtering")& 
  theme(plot.title = element_text(hjust = 0.5))



#violin plot of mito
feats <- c("percent_mito")
VlnPlot(alldata,   group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 1) + NoLegend() #all mitochondrial counts fall below ~1%, no need for filtering?
#Histograms of mitochondiral percent
hi
#Scatterplot of count to feature relationship
plot1 <- FeatureScatter(alldata_with_dbs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(filterdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

grid.arrange(plot1, plot2,plot3, ncol=3)

plot + geom_hline(yintercept=300, linetype="dashed", color = "red", size=0.5) #use lower threshold to not loose to much biological signas
st(alldata@meta.data[["percent_mito"]],breaks=100)
abline(v=0.025, col="blue")

#violin plot of hb
feats <- c("percent_hb", "percent_mito")
VlnPlot(alldata,   group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 1) + NoLegend() #all mitochondrial counts fall below ~1%, no need for filtering?
#Histograms of hb percent
hist(alldata@meta.data[["percent_hb"]],breaks=100)
abline(v=0.025, col="blue")

#Histograms of RNA counts
par(mfrow=c(1, 3))
hist(alldata_with_dbs@meta.data[["nCount_RNA"]],breaks = 100,main="Counts, data with doublets")
hist(alldata@meta.data[["nCount_RNA"]],breaks = 100,main="Counts, data after doublet removal")
hist(filterdata@meta.data[["nCount_RNA"]],breaks = 100,main="Counts, data after filtering")


hgA <- hist(alldata_with_dbs@meta.data[["nCount_RNA"]],breaks = 200, plot = FALSE) # Save first histogram data
hgB <- hist(alldata@meta.data[["nCount_RNA"]],breaks = 200, plot = FALSE) # Save 2nd histogram data
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
plot(hgB, col = c1) # Plot 1st histogram using a transparent color
plot(hgB, col = c2, add = TRUE) # Add 2nd histogram using different color


abline(v=1100, col="blue")
#Counts up to 4000
hist(alldata@meta.data[["nCount_RNA"]][alldata@meta.data[["nCount_RNA"]] < 4000],breaks = 100,col=c1) #small peak between 500-100, set threshold at 1100 ?
abline(v=1100, col="blue")

#Histograms of feature count
#Histograms of RNA counts
par(mfrow=c(1, 3))
hist(alldata_with_dbs@meta.data[["nFeature_RNA"]],breaks = 80,main="Feature count, data with doublets")
hist(alldata@meta.data[["nFeature_RNA"]],breaks = 80,main="Feature count, data after doublet removal")
hist(filterdata@meta.data[["nFeature_RNA"]],breaks = 80,main="Feature count, data after filtering")



hist(alldata@meta.data[["nFeature_RNA"]],breaks = 100,col=c1) #filtered at 300
abline(v=650, col="blue")
hist(alldata@meta.data[["nFeature_RNA"]][alldata@meta.data[["nFeature_RNA"]] < 1000],breaks = 100) 


