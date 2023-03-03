#create figures
library(ggplot2)
source("markers.R")

#stor figur
colours <- c(rep("#de8a00ff",3),rep("#00bf8aff",6),rep("#00b3f0ff",1),rep("#ff63b0ff",2),rep("#c77affff",1))

stora.fig <- c("LOC396098","SOX5","PAX5","CD3E","CD3D","TARP","CD4","CD8A","CD8BP","MMR1L4","ITGA2B","ITGB3","HBBA")

DotPlot(combined.sct,features = rev(stora.fig), cluster.idents = T)+coord_flip()+
  theme(axis.text.y = element_text(colour=rev(colours)))

#shah markers
#B-cells, DC, T
colours <- c(rep("#de8a00ff",9),rep("#00bf8aff",8),rep("#c77affff",9))

DotPlot(combined.sct,features = rev(shah_markers), cluster.idents = T)+coord_flip()+
  theme(axis.text.y = element_text(colour=rev(colours)))

#proliferating
combined.sct.prol <- subset(combined.sct,idents = 22)
DimPlot(combined.sct.prol,label=T, label.box=T, repel=T)+NoLegend()+ggtitle("Putative T-cells")+
  theme(plot.title = element_text(hjust = 0.5))

combined.sct.prol.re <- ReclusterSCT(combined.sct, 0.9, 22,0.8)
DimPlot(combined.sct.prol.re,label=T, label.box=T, repel=T)+NoLegend()+ggtitle("Putative proliferating cells")+
  theme(plot.title = element_text(hjust = 0.5))

#T-cell figure
colours <- c(rep("black",9),rep("pink",3),rep("green",8),rep("purple",8),rep("blue",11),rep("black",7))

dat <- data.frame(name=t.cell.figure.221216,value=colours)
t.clusts <- c(0,	1,	23,	25,	29,	4,	7,	3,	8,	11,	12,	13,	14,	15,	16,	24,	30)
combined.sct.t <- subset(combined.sct,idents = t.clusts)
levels(combined.sct.t) <- t.clusts

DimPlot(combined.sct.t,label=T, label.box=T, repel=T)+NoLegend()+ggtitle("Putative T-cells")+
  theme(plot.title = element_text(hjust = 0.5))
DotPlot(combined.sct.t,features = rev(t.cell.figure.221216))+coord_flip()+
  theme(axis.text.y = element_text(colour=rev(colours)))

t.cell.4 <- ReclusterSCT(combined.sct, 0.9, c(4),0.3)


DimPlot(t.cell.4,label=T, label.box=T, repel=T)+NoLegend()+ggtitle("Putative cytolytic cells, Re-clustered")+
  theme(plot.title = element_text(hjust = 0.5))
DotPlot(t.cell.4,features = rev(t.cell.figure.221216))+coord_flip()+
  theme(axis.text.y = element_text(colour=rev(colours)))


#Monocyte figure
colours <- c(rep("black",5),rep("yellow2",6),rep("green",6),rep("black",16),rep("purple",3))

dat <- data.frame(name=mono.figure.221230,value=colours)
mono.clusts <- c(5,6)
combined.sct.mono <- subset(combined.sct,idents = mono.clusts)
levels(combined.sct.mono) <- mono.clusts

DimPlot(combined.sct.mono,label=T, label.box=T, repel=T)+NoLegend()+ggtitle("Putative Monocytes")+
  theme(plot.title = element_text(hjust = 0.5))

DotPlot(combined.sct,features = rev(mono.figure.221230),idents=c(5,6),scale=F)+coord_flip()+
  theme(axis.text.y = element_text(colour=rev(colours)))

#reclustered monocytes
colours <- c(rep("black",5),rep("orange",6),rep("green",6),rep("black",16),rep("purple",3),rep("#00bf8aff",8),rep("#00bf8bff",6))

combined.sct.mono.re <- ReclusterSCT(combined.sct, 0.9, c(5,6),0.8)

DimPlot(combined.sct.mono.re,label=T, label.box=T, repel=T)+NoLegend()+ggtitle("Putative Monocytes, Reclustered")+
  theme(plot.title = element_text(hjust = 0.5))

DotPlot(combined.sct.mono.re,features = rev(mono.figure.221230))+coord_flip()+
  theme(axis.text.y = element_text(colour=rev(colours)))

#B cell figure
colours <- c(rep("black",10),rep("orange",5),rep("green",2),rep("purple",3))
colours <- c(rep("black",1),rep("blue",7),rep("green",4),rep("black",5),rep("red",11),rep("black",3),rep("purple",4))


b.clusts <- c(2,10,17,18,19,20,21,26)
combined.sct.b <- subset(combined.sct,idents = b.clusts)
levels(combined.sct.b) <- b.clusts

bcellfeats <- c("EBF1", "PAX5", "IGLL1", "CXCR4","TNFSF13B","HMGB1", "LAMP3", "CD79B","JCHAIN", "PRDM1", "TNFRSF13C", "LOC396098", "CD74", "BLB1", "SOX5")

DimPlot(combined.sct.b,label=T, label.box=T, repel=T)+NoLegend()+ggtitle("Putative B-cells")+
  theme(plot.title = element_text(hjust = 0.5))

DotPlot(combined.sct.b,features = rev( b.cell.figure.230112))+coord_flip()+
  theme(axis.text.y = element_text(colour=rev(colours)))

#cluster 22
DotPlot(combined.sct,features = rev(prol),scale=F)+coord_flip()+
  theme(axis.text.y = element_text(colour=rev(colours)))

