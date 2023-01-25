#appendix data


## compare results of clustering algorithms
combined.sct.l1 <- FindClusters(combined.sct,verbose = TRUE, algorithm=1)
plot1 <- DimPlot(combined.sct.l1, reduction = "umap_30", label=T)+ NoLegend()+ labs(title="Louvain algorithm (Default)")


combined.sct.l2 <- FindClusters(combined.sct,verbose = TRUE,  algorithm=2)
plot2 <- DimPlot(combined.sct.l2, reduction = "umap_30", label=T)+ NoLegend()+ labs(title="Louvain algorithm with multilevel refinement")


combined.sct.slm <- FindClusters(combined.sct,verbose = TRUE,  algorithm=3)
plot3 <- DimPlot(combined.sct.slm, reduction = "umap_30", label=T)+ NoLegend()+ labs(title="Smart local moving algorithm")


grid.arrange(plot1,plot2,plot3, ncol=3,nrow=1)

## Compare umaps
######## testing min.dist
#plot0 <- DimPlot(combined.sct, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()+ labs(title="Default settings")

combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:27, min.dists = 0.001, reduction.name = "umap_0.001")
plot01<- DimPlot(combined.sct, reduction = "umap_0.001", label = TRUE, repel = TRUE) + NoLegend()+ labs(title="Min dists = 0.001")

combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:27, min.dists = 0.15, reduction.name = "umap_0.15")
plot1<- DimPlot(combined.sct, reduction = "umap_0.15", label = TRUE, repel = TRUE) + NoLegend()+ labs(title="Min dists = 0.15")

combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:27 , reduction.name = "umap_0.3")
plot02<- DimPlot(combined.sct, reduction = "umap_0.3", label = TRUE, repel = TRUE) + NoLegend()+ labs(title="Min dists = 0.3 (default)")

#combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:27 , reduction.name = "umap_30_1", umap.method= 'umap-learn', metric = "correlation")

combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:27,min.dists = 0.45 , reduction.name = "umap_0.45")
plot2 <- DimPlot(combined.sct, reduction = "umap_0.45", label = TRUE, repel = TRUE) + NoLegend()+ labs(title="Neighbors=0.45")

combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:27,min.dists = 0.5, reduction.name = "umap_0.5")
plot3 <- DimPlot(combined.sct, reduction = "umap_0.5", label = TRUE, repel = TRUE) + NoLegend()+ labs(title="Neighbors=0.5")

require(gridExtra)

grid.arrange(plot01,plot1, plot02,plot2,plot3, ncol=3,nrow=2)


######## testing n.neighbours
#plot0 <- DimPlot(combined.sct, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()+ labs(title="Default settings")

combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:27, n.neighbors = 5, reduction.name = "umap_5")
plot01<- DimPlot(combined.sct, reduction = "umap_5", label = TRUE, repel = TRUE) + NoLegend()+ labs(title="Neighbors=5")

combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:27, n.neighbors = 25, reduction.name = "umap_25")
plot1<- DimPlot(combined.sct, reduction = "umap_25", label = TRUE, repel = TRUE) + NoLegend()+ labs(title="Neighbors=25")

combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:27 , reduction.name = "umap")
plot02<- DimPlot(combined.sct, reduction = "umap_30", label = TRUE, repel = TRUE) + NoLegend()+ labs(title="Neighbors=30 (default)")

#combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:27 , reduction.name = "umap_30_1", umap.method= 'umap-learn', metric = "correlation")

combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:27, n.neighbors = 50 , reduction.name = "umap_50")
plot2 <- DimPlot(combined.sct, reduction = "umap_50", label = TRUE, repel = TRUE) + NoLegend()+ labs(title="Neighbors=50")

combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:27, n.neighbors = 75, reduction.name = "umap_75")
plot3 <- DimPlot(combined.sct, reduction = "umap_75", label = TRUE, repel = TRUE) + NoLegend()+ labs(title="Neighbors=75")

combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:27, n.neighbors = 100, reduction.name = "umap_100")
plot4 <- DimPlot(combined.sct, reduction = "umap_100", label = TRUE, repel = TRUE) + NoLegend()+ labs(title="Neighbors=100")

require(gridExtra)

grid.arrange(plot01,plot1, plot02,plot2,plot3, plot4, ncol=3,nrow=2)

grid.arrange(plot0,plot02, ncol=2)

#fractions gdT
combined.sct.1.t1 <- subset(subset(combined.sct,orig.ident=="Sample_1"),idents = c(0,1,23,25))
combined.sct.2.t1 <- subset(subset(combined.sct,orig.ident=="Sample_2"),idents = c(0,1,23,25))
combined.sct.3.t1 <- subset(subset(combined.sct,orig.ident=="Sample_3"),idents = c(0,1,23,25))
combined.sct.4.t1 <- subset(subset(combined.sct,orig.ident=="Sample_4"),idents = c(0,1,23,25))

#fractions abT
combined.sct.1.t2 <- subset(subset(combined.sct,orig.ident=="Sample_1"),idents = c(3,8,11,12,13,14,15,16,24,30))
combined.sct.2.t2 <- subset(subset(combined.sct,orig.ident=="Sample_2"),idents = c(3,8,11,12,13,14,15,16,24,30))
combined.sct.3.t2 <- subset(subset(combined.sct,orig.ident=="Sample_3"),idents = c(3,8,11,12,13,14,15,16,24,30))
combined.sct.4.t2 <- subset(subset(combined.sct,orig.ident=="Sample_4"),idents = c(3,8,11,12,13,14,15,16,24,30))

#Fractions
combined.sct.mono <- ReclusterSCT(combined.sct,sub.idents =c(5,6))
combined.sct.1.het <- subset(subset(combined.sct.mono,orig.ident=="Sample_1"),idents = 2)
combined.sct.2.het <- subset(subset(combined.sct.mono,orig.ident=="Sample_2"),idents = 2)
combined.sct.3.het <- subset(subset(combined.sct.mono,orig.ident=="Sample_3"),idents = 2)
combined.sct.4.het <- subset(subset(combined.sct.mono,orig.ident=="Sample_4"),idents = 2)

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


#Extra stuff
clust.17 <- ReclusterSCT(combined.sct,sub.idents = 17)




