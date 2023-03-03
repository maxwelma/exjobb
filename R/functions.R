# Data filtering of the doublet free data from python
#based on tutorial: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Calculate_QC
#Seurat is required

#load data from samples directory
#Give argument for whether or not doublets are included in the data or removed
#True = doublets are present
#False = Doublets are removed 
#If no argument is given doublets are assumed to be present (default)
LoadData <- function(matrix_dir,projectn,doublets)
{
  if(missing(doublets)) {
    doublets <- TRUE
  }
  #load data
  #one sample
  if(doublets==FALSE){
    barcode.path <- paste0(matrix_dir, "barcodes_no_doublets.tsv") #remove no_doublets in final version
    features.path <- paste0(matrix_dir, "features.tsv")
    matrix.path <- paste0(matrix_dir, "matrix_no_doublets.mtx")
  }
  if(doublets==TRUE){
    barcode.path <- paste0(matrix_dir, "barcodes.tsv") #remove no_doublets in final version
    features.path <- paste0(matrix_dir, "features.tsv")
    matrix.path <- paste0(matrix_dir, "matrix.mtx")
  }
  
  samp <- Seurat::ReadMtx(mtx=matrix.path,cells=barcode.path,features=features.path)
  sdata <- CreateSeuratObject(samp, project=projectn)
  #Calculate matadata
  #automatically includes
  #nFeature_RNA - The number of genes detected in each cell. Low feature count can indicate dead or dying cells
  #nCount_RNA - The total number of molecules detected within a cell. High RNA count can indicate doublet
  #percentage_mito - The percentage of reads per cell that correspont to mitochondira. High percentage can idicate dead/dying cells
  #percentage of mitochondiral, ribosomal, hemoglobin genes, and thrombocyte genes
  sdata <- PercentageFeatureSet(sdata, "^MT", col.name = "percent_mito")
  sdata <- PercentageFeatureSet(sdata, "^RP[SL]", col.name = "percent_ribo")
  sdata <- PercentageFeatureSet(sdata, "^HB[^(P)]", col.name = "percent_hb")
  sdata <- PercentageFeatureSet(sdata, "PECAM1|PF4", col.name = "percent_plat")
  
  return(sdata) 
}


#thresholds based on quality_control.R
#Filter out cells that falls below certain thresholds
FilterData <- function(alldata, cellcount, uniquefeaturecount, percentmito, percenthb)
{
  if(missing(cellcount)) {
    cellcount <- 3
  }
  if(missing(uniquefeaturecount)){
    uniquefeaturecount <- 300
  }
  if(missing(percentmito)){
    percentmito <- 20
  }
  if(missing(percenthb)){
    percenthb <- 5
  }
  
  #Filter based on quality metrics
  # filter genes that are expressed  in less that 3 cells
  # filter cells that have unique feature counts less than 300 
  selected_c <- WhichCells(alldata, expression = nFeature_RNA > uniquefeaturecount)
  selected_f <- rownames(alldata)[Matrix::rowSums(alldata) > cellcount]
  
  data.filt <- subset(alldata, features = selected_f, cells = selected_c)
  #dim(data.filt)
  
  # filter cells that have >5% red blood cell reads and < 20% mitochondrial genes
  selected_mito <- WhichCells(data.filt, expression = percent_mito < percentmito)
  selected_hb <- WhichCells(data.filt, expression = percent_hb < percenthb)
  
  # filter cells that have 
  
  #subset the object to only keep those cells
  data.filt <- subset(data.filt, cells = selected_mito)
  data.filt <- subset(data.filt, cells = selected_hb)
  
  #Inspect filtered data
  #dim(data.filt)
  #table(data.filt$orig.ident)
  return(data.filt)
}


#Load all datasets into one object
LoadAllData <- function(matrix.dir,dirs.list,doublets){
  samples.list <- list()
  if(missing(doublets)) {
    print("no doublet argument given, assuming doublets present")
    doublets <- TRUE
    for (i in 1:length(dirs.list)) {
      samples.list[[i]] <- LoadData(paste0(paste0(matrix_dir,dirs.list[[i]]),"\\"),projectn=dirs.list[[i]],doublets) #remove no doublets in final v
    }
  }
  if(doublets==FALSE){
    for (i in 1:length(dirs.list)) {
      samples.list[[i]] <- LoadData(paste0(paste0(matrix_dir,dirs.list[[i]]),"\\no_doublets\\"),projectn=dirs.list[[i]],doublets) #remove no doublets in final v
    }
  }else{
    for (i in 1:length(dirs.list)) {
      samples.list[[i]] <- LoadData(paste0(paste0(matrix_dir,dirs.list[[i]]),"\\"),projectn=dirs.list[[i]],doublets) #remove no doublets in final v
    }
  }
  #merge data into one set
  alldata <- merge( samples.list[[1]],  samples.list[-c(1)], add.cell.ids = dirs.list)
  return(alldata)
}


#Normalize filtered data using sct and prep data for integration by producing integration anchors
ProduceAnchorsSCT<- function(object){
  samples.list <- SplitObject(object=object, split.by = "orig.ident")
  
  samples.list <-lapply(X = samples.list, FUN = SCTransform, return.only.var.genes=F, method = "glmGamPoi")
  
  features <- SelectIntegrationFeatures(object.list = samples.list, nfeatures = 3000)
  samples.list <- PrepSCTIntegration(object.list = samples.list, anchor.features = features)
  samples.anchors <- FindIntegrationAnchors(object.list = samples.list, normalization.method = "SCT",
                                            anchor.features = features)
  
  saveRDS(samples.anchors, file = "anchors_sct.RDS") 
  return(samples.anchors)
}


#cluster data using pca and print umap
ClusterDataSCT <- function(object,variance,resolution_clusts){
  
  if(missing(resolution_clusts)) {
    resolution_clusts <- 0.8
  } 
  if(missing(variance)) {
    variance <- 0.9
  } 
  
  if ("integrated"%in%Assays(object)){
    DefaultAssay(object) <- "integrated"
  }else{
    assay <- DefaultAssay(object)
    print(paste0("This object does not include assay: integrated. Continuing with assay: ",assay))
  }
  
  object <- RunPCA(object, verbose = TRUE)
  
  PoV <- object@reductions[["pca"]]@stdev^2/sum(object@reductions[["pca"]]@stdev^2)
  cumulative.PoV <- cumsum(PoV)
  dims.var <- min(which(cumulative.PoV > variance))
  print(paste("PCA dimensions: ", dims.var))
  
  object <- FindNeighbors(object, dims = 1:dims.var, verbose = FALSE)
  object <- FindClusters(object, verbose = TRUE, resolution = resolution_clusts)
  
  object <- RunUMAP(object, reduction="pca", dims = 1:dims.var, verbose = FALSE)
  
  if ("SCT"%in%Assays(object)){
    DefaultAssay(object) <- "SCT"
  }else{
    assay <- DefaultAssay(object)
    print(paste0("This object does not include assay: SCT. Continuing with assay: ",assay))
  }
  print(DimPlot(object, reduction = "umap", label = TRUE, repel = TRUE)+ggplot2::ggtitle(label = "UMAP on PCA")+ NoLegend())
  return(object)
}

#cluster data subset based on cluster identities
ReclusterSCT <- function(object, variance, sub.idents,resolution_clusts){
  object.subset<- subset(object, idents=sub.idents)
  object.subset <- ClusterDataSCT(object.subset,variance,resolution_clusts)
  return(object.subset)
}

###Functions related to producing excel outputs

#Return ncbi descriptions for gene list, requires internet
getDescritions <- function(genes,l){
  search=rep(NA, l)
  j=0
  for (k in genes){
    j=j+1
    print(k)
    #    if (grepl( "LOC", k, fixed = TRUE) || grepl( "J6367", k, fixed = TRUE)){
    #      search[j]=k
    #    }else{
    sterm = paste0(paste0( "Homo sapiens[ORGN] OR Homo sapiens[ALL] AND ", k), " [Gene]")
    sid=entrez_search(db="gene", term=sterm,retmax=1)[["ids"]]
    
    if (length(sid) == 0){
      sterm = paste0(paste0("Gallus gallus[ORGN] OR Gallus gallus[ALL] AND ", k), " [Gene]")
      sid=entrez_search(db="gene", term=sterm,retmax=1)[["ids"]]
    }
    
    if (length(sid) == 0){
      search[j]=k
    }else{
      search[j]=paste0(paste0(k," - "),paste(entrez_summary(db="gene",id=sid)[c(2,3,17)],collapse = " - "))
      #print(search[j])
    }
    #    }
  }
  
  #print(search)
  return(search)
}



#print the resulting markers to 
PrintMarkers <- function(orig.object,object,markers, filename,featuresdotplot){
  if(file.exists(filename)==T){
    ## delete the file if already existing
    unlink(filename)
    print("File is being overwriten")
  }
  
  
  wb <- createWorkbook()
  #strating sheet for file, includes dotplots and violinplots of relevant genes (featuresdotplot), dimplot of relevant clusters
  sheet  <- createSheet(wb, sheetName=as.character("Clustered cells"))
  
  file <- paste0(filename,"_whole_dotplot.png")
  print(file)
  png(file,width=1000, height=1000)
  print(DotPlot(object,features =  featuresdotplot,scale=F)+coord_flip())
  dev.off()
  addPicture(file, sheet,startColumn = 1)
  
  file <- paste0(filename,"_whole_dotplot_scaled.png")
  print(file)
  png(file,width=1000, height=1000)
  print(DotPlot(object,features =  featuresdotplot,scale=T)+coord_flip())
  dev.off()
  addPicture(file, sheet,startRow = 55)
  
  file <- paste0(filename,"_whole_vln.png")
  png(file,width=1000, height=1000)
  print(VlnPlot(object,features =  featuresdotplot,assay = ))
  dev.off()
  addPicture(file, sheet,startColumn = 18, startRow = 55)
  
  file <- paste0(filename,"_whole_dimplot.png")
  png(file,width=1000, height=1000)
  print(DimPlot(orig.object,cells.highlight = Cells(object),label=T))
  dev.off()
  addPicture(file, sheet,startColumn = 18)
  
  
  file <- paste0(filename,"_subset_dimplot.png")
  png(file,width=1000, height=1000)
  print(DimPlot(object,label=T))
  dev.off()
  addPicture(file, sheet,startColumn = 36)
  
  #insert clustree picture if it exists in clustree/subs
  dir <-  paste("clustree/subs/",strsplit(filename,"/")[[1]][2],"/",strsplit(filename,"/")[[1]][2],"_clustree_stability.png",sep="")
  if (exists(dir)){
    addPicture(dir, sheet, startColumn=54)
  }
  
  #check if the markers df has a "cluster" column
  if(length(markers$cluster)>0){
    #loop over all clusters in subset
    for (i in unique(markers$cluster)){
      #subset markers from cluster
      markers_sub <- subset(markers,cluster==i)
      markers_sub<- markers_sub[order(-markers_sub$avg_log2FC),]
      
      #check if markers df has a "gene" column, if not use rownames 
      #(if markers has several clusters it will name the ronames gene.1, gene.2, etc. Then the gene column is needed)
      
      if ("gene"%in%names(markers_sub)){
        genes <- head(subset(markers_sub,avg_log2FC >0.25)$gene,n=100)
      }else{
        genes <- head(rownames(subset(markers_sub,avg_log2FC >0.25)),n=100)
      }
      
      #Uncomment below to add gene descriptions from ncbi
      #descs <- getDescritions(genes, length(rownames(markers_sub)))
      #markers_sub["description"] <- descs
      
      #creare
      sheet  <- createSheet(wb, sheetName=as.character(paste0("Cluster",i)))
      cs <- CellStyle(wb, alignment = Alignment(wrapText = TRUE))
      dfColIndex <- rep(list(cs), ncol(markers_sub)) 
      names(dfColIndex) <- ncol(markers_sub)
      
      setColumnWidth(sheet, colIndex=ncol(markers_sub)+1, colWidth=100)
      addDataFrame(markers_sub,sheet=sheet,colStyle= dfColIndex)
      
      #add go terms 
      sheet_go  <- createSheet(wb, sheetName=as.character(paste0("GO_cluster",i)))
      gostres <- gost(query =genes, organism = "ggallus") 
      if (length(gostres)==0){
        res <- data.frame(0) 
      }else{
        res <- gostres$result[order(gostres$result$p_value),]
      }
      addDataFrame(res,sheet=sheet_go)
      
      file = paste0(paste0(paste0(filename,"_dimplot_subcluster_"),i),".png")
      png(file,width=800, height=800)
      print(DimPlot(object,cells.highlight = WhichCells(object,idents = i)))
      dev.off()
      addPicture(file, sheet,startRow=2, startColumn=45)
      
      file2 = paste0(paste0(paste0(filename,"_dotplot_subcluster_"),i),".png")
      png(file2,width=800, height=2000)
      print(DotPlot(object,features =  rev(genes),scale=T)+coord_flip())
      dev.off()
      addPicture(file2, sheet,startRow=2, startColumn=30)
      
      file2 = paste0(paste0(paste0(filename,"_dotplot_subcluster_"),i),".png")
      png(file2,width=800, height=2000)
      print(DotPlot(object,features =  rev(genes),scale=F)+coord_flip())
      dev.off()
      addPicture(file2, sheet,startRow=2, startColumn=15)
      
      saveWorkbook(wb, paste0(filename,".xlsx"))
    }
  }else{
    markers_sub<-markers
    markers_sub<- markers_sub[order(-markers_sub$avg_log2FC),]
    
    if ("gene"%in%names(markers_sub)){
      genes <- head(subset(markers_sub,avg_log2FC >0.25)$gene,n=100)
    }else{
      genes <- head(rownames(subset(markers_sub,avg_log2FC >0.5)),n=100)
    }
    
    descs <- getDescritions(genes, length(rownames(markers_sub)))
    
    markers_sub["description"] <- descs
    
    sheet  <- createSheet(wb, sheetName=as.character("All"))
    
    cs <- CellStyle(wb, alignment = Alignment(wrapText = TRUE))
    dfColIndex <- rep(list(cs), ncol(markers_sub)) 
    names(dfColIndex) <- ncol(markers_sub)
    
    setColumnWidth(sheet, colIndex=ncol(markers_sub)+1, colWidth=100)
    addDataFrame(markers_sub,sheet=sheet,colStyle= dfColIndex)
    
    #add kegg go terms
    sheet_go  <- createSheet(wb, sheetName=as.character("GO_All"))
    gostres <- gost(query =genes, organism = "ggallus") 
    if (length(gostres)==0){
      res <- data.frame(0) 
    }else{
      res <- gostres$result[order(gostres$result$p_value),]
    }
    addDataFrame(res,sheet=sheet_go)
    
    file = paste0(paste0(filename,"_dimplot_subcluster_"),".png")
    png(file,width=800, height=800)
    print(DimPlot(object,cells.highlight = Cells(object)))
    dev.off()
    addPicture(file, sheet,startRow=2, startColumn=45)
    
    
    file2 = paste0(paste0(filename,"_dotplot_subcluster_"),".png")
    png(file2,width=800, height=2000)
    print(DotPlot(object,features =  rev(genes),scale=T)+coord_flip())
    dev.off()
    addPicture(file2, sheet,startRow=2, startColumn=30)
    
    
    file2 = paste0(paste0(filename,"_dotplot_subcluster_"),".png")
    png(file2,width=800, height=2000)
    print(DotPlot(object,features =  rev(genes),scale=F)+coord_flip())
    dev.off()
    addPicture(file2, sheet,startRow=2, startColumn=15)
    
    saveWorkbook(wb, paste0(filename,".xlsx"))
  }
}



#print the resulting markers to 
PrintMarkersRaw <- function(orig.object,object,markers, filename){
  if(file.exists(filename)==T){
    ## delete the file if already existing
    unlink(filename)
    print("File is being overwriten")
  }
  
  
  wb <- createWorkbook()

  #check if the markers df has a "cluster" column
  if(length(markers$cluster)>0){
    #loop over all clusters in subset
    for (i in unique(markers$cluster)){
      #subset markers from cluster
      markers_sub <- subset(markers,cluster==i)
      markers_sub<- markers_sub[order(-markers_sub$avg_log2FC),]
      
      #check if markers df has a "gene" column, if not use rownames 
      #(if markers has several clusters it will name the ronames gene.1, gene.2, etc. Then the gene column is needed)
      
      if ("gene"%in%names(markers_sub)){
        genes <- head(subset(markers_sub,avg_log2FC >0.25)$gene,n=100)
      }else{
        genes <- head(rownames(subset(markers_sub,avg_log2FC >0.25)),n=100)
      }
      
      #Uncomment below to add gene descriptions from ncbi
      #descs <- getDescritions(genes, length(rownames(markers_sub)))
      #markers_sub["description"] <- descs
      
      #creare
      sheet  <- createSheet(wb, sheetName=as.character(paste0("Cluster",i)))
      addDataFrame(markers_sub,sheet=sheet,colStyle= dfColIndex)
      
      #add go terms 
      sheet_go  <- createSheet(wb, sheetName=as.character(paste0("GO_cluster",i)))
      gostres <- gost(query =genes, organism = "ggallus") 
      if (length(gostres)==0){
        res <- data.frame(0) 
      }else{
        res <- gostres$result[order(gostres$result$p_value),]
      }
      addDataFrame(res,sheet=sheet_go)
      
      saveWorkbook(wb, paste0(filename,".xlsx"))
    }
  }else{
    markers_sub<-markers
    markers_sub<- markers_sub[order(-markers_sub$avg_log2FC),]
    
    if ("gene"%in%names(markers_sub)){
      genes <- head(subset(markers_sub,avg_log2FC >0.25)$gene,n=100)
    }else{
      genes <- head(rownames(subset(markers_sub,avg_log2FC >0.5)),n=100)
    }
    
    sheet  <- createSheet(wb, sheetName=as.character("All"))
    addDataFrame(markers_sub,sheet=sheet,colStyle= dfColIndex)
    
    #add kegg go terms
    sheet_go  <- createSheet(wb, sheetName=as.character("GO_All"))
    gostres <- gost(query =genes, organism = "ggallus") 
    if (length(gostres)==0){
      res <- data.frame(0) 
    }else{
      res <- gostres$result[order(gostres$result$p_value),]
    }
    addDataFrame(res,sheet=sheet_go)
    
    saveWorkbook(wb, paste0(filename,".xlsx"))
  }
}

#produce and save markers for object subset
Markers2xlsx <- function(object, clusters, directory, filename, resolution_clusts,featuresdotplot,print_plot){
  if (missing(print_plot)){
    print_plot = T
  }
  #markers between clusters within subset 
  clust <- ReclusterSCT(object = object , 0.9, clusters,resolution_clusts)
  clust<- PrepSCTFindMarkers(clust)
  cluster.markers.in <- FindAllMarkers(clust,assay = "SCT",only.pos = T)
  
  if(dir.exists(directory)==F){
    dir.create(directory, recursive = TRUE) 
  }
  
  file <- paste0(directory,paste0(filename,"_1_in"))
  if (print_plot == T){
    PrintMarkers(orig.object=object,object=clust,markers = cluster.markers.in, filename=file,featuresdotplot=featuresdotplot)
  }else{
    PrintMarkersRaw(orig.object=object,object=clust,markers = cluster.markers.in, filename=file)
  }
  #markers between clusters in subset and all other cells
  #initiate empty df
  cluster.markers.out<-setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj"))
  
  #loop over subclusters and find out markers
  for (i in levels(clust)){
    marks<- FindMarkers(object,ident.1 = WhichCells(clust,idents = i),assay = "SCT",only.pos = T)
    marks["cluster"] <- rep(i,nrow(marks))
    marks["gene"] <- rownames(marks) 
    cluster.markers.out <-rbind(cluster.markers.out,marks)
  }
  
  file <- paste0(directory,paste0(filename,"_2_out"))
  if (print_plot == T){
    PrintMarkers(orig.object=object, object=clust,markers = cluster.markers.out, filename=file,featuresdotplot=featuresdotplot)
  }else{
    PrintMarkersRaw(orig.object=object, object=clust,markers = cluster.markers.out, filename=file)
  }
  return(list(cluster.markers.in,cluster.markers.out))
}


#annotate clusters manually
annotateClusters <- function(object,levelmap){
  #give a map of what clusters belong to what cell type, if not given use the map below
  if(missing(levelmap)){
    #make dataframe that is the length of clusters in object
    #set values in col 2 to manually annotated cell type
    levelmap <- data.frame(matrix(ncol = 2,nrow = length(levels(object))))
    levelmap[,1] <- levels(object)
    levelmap[c(0,1,23,25)+1,2] <- "g/dT-cells"
    levelmap[c(3,8,11,12,13,14,15,16,24,30)+1,2] <- "CD4T-cells"
    levelmap[c(4)+1,2] <- "Cytolytic-cells"
    levelmap[c(7,29)+1,2] <- "CD8T-cells"
    levelmap[c(27)+1,2] <- "RBC"
    levelmap[c(28)+1,2] <- "Basophils"
    levelmap[c(9)+1,2] <- "Thrombocytes"
    levelmap[c(5,6)+1,2] <- "Monocytes"
    levelmap[c(2,10,18,20,21,26)+1,2] <- "B-cells"
    levelmap[c(19)+1,2] <- "SOX5B-cells"
    levelmap[c(22)+1,2] <- "Proliferating"
    levelmap[c(17)+1,2] <- "ProB-cells"
    rownames(levelmap)<-levelmap$X1
  }
  #get cluster assignemt of all cells
  int.clusters <-data.frame(object@meta.data$seurat_clusters)
  
  #return cellt type for a given cluster nr
  setCluster <- function(n){
    celltype <- as.character(levelmap[n,][2])
    return(celltype)
  }
  
  #make list of celltype per cell
  #Add celtype to metadata
  cell_types<- as.character(lapply(as.character(int.clusters[,1]), setCluster))
  #match name of cell to celltype
  names(cell_types) <- colnames(x = object)
  
  object<-AddMetaData(object, factor(cell_types), col.name = "celltype")
  return(object)
}

#Calculate counts of cells that express >0.2 of a gene
cellsSample  <- function(object){
  counts.by.ident<-c()
  counts.by.ident["TRBV6-5"]<- length(WhichCells(object = object, expression = `TRBV6-5` >0.2, slot = 'data'))
  counts.by.ident["TARP"]<- length(WhichCells(object = object, expression = TARP >0.2, slot = 'data'))
  counts.by.ident["LOC112530376"]<- length(WhichCells(object = object, expression = LOC121112191 >0.2, slot = 'data'))
  counts.by.ident["CD4"]<- length(WhichCells(object = object, expression = CD4 >0.2, slot = 'data'))
  counts.by.ident["CD8A"]<- length(WhichCells(object = object, expression = CD8A >0.2, slot = 'data'))
  counts.by.ident["CD8BP"]<- length(WhichCells(object = object, expression = CD8BP>0.2, slot = 'data'))
  counts.by.ident["MMR1L4"]<- length(WhichCells(object = object, expression = MMR1L4 >0.2, slot = 'data'))
  counts.by.ident["LOC396098"]<- length(WhichCells(object = object, expression = LOC396098 >0.2, slot = 'data'))
  counts.by.ident["ITGA2B"]<- length(WhichCells(object = object, expression = ITGA2B >0.2, slot = 'data'))
  counts.by.ident["ITGB3"]<- length(WhichCells(object = object, expression = ITGB3 >0.2, slot = 'data'))
  return(counts.by.ident)
}

#Plot fractions of celltypes, requires manual annotation of cells
plotStackedBars <- function(object){
  cluster_freq.table <- table(object$orig.ident, object$celltype) %>% melt()
  colnames(cluster_freq.table) <- c("Sample", "Celltype", "Fraction_of_all_cells")
  
  cluster_freq.table$celltype <- as.factor(cluster_freq.table$Celltype)
  levs <- levels(cluster_freq.table$celltype )
  cluster_freq.table$celltype <- droplevels(cluster_freq.table$Celltype)
  
  ggplot(cluster_freq.table, aes(
    x = Sample,
    y = Fraction_of_all_cells,
    fill = factor(Celltype, levels = levs%>%rev())
  )) +
    geom_bar(stat = "identity", position = "fill", width = 0.9666) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip() +
    ylab(label = "Fraction of all cells") +
    scale_x_discrete(limits = rev(levels(cluster_freq.table$Subject))) +
    labs(fill = "celltype")+
    ggtitle("Fractions of cells per sample")+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text = element_text(size = 12))
}


