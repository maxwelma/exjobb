#Recluster and produce excel files for different subsets and resolutions
source("functions.R")
source("markers.R")
#all
all.markers.both <- Markers2xlsx(combined.sct, c(0:30), "reclustering_GO_new_res/whole/08/", "whole_go_08", 0.8) #- samma
all.markers.both <- Markers2xlsx(combined.sct, c(0:30), "reclustering_GO_new_res/whole/02/", "whole_go_02", 0.2) #- samma

#cluster 28
clust28.markers.both <- Markers2xlsx(combined.sct, c(28), "reclustering_GO_new_res/clust28/08/", "clust22_GO_0.8", 0.8, featuresdotplot =all_subset_marks ) #- samma

#active cells
clust22.markers.both <- Markers2xlsx(combined.sct, c(22), "reclustering_GO_new_res/clust22/08/", "clust22_GO_0.8", 0.8,featuresdotplot=unique(c(tcellfeats,bcellfeats,monocytefeats,all_subset_marks))) #- samma
clust22.markers.both <- Markers2xlsx(combined.sct, c(22), "reclustering_GO_new_res/clust22/03/", "clust22_GO_0.3", 0.3, featuresdotplot=c(tcellfeats,bcellfeats))
clust17.markers.both <- Markers2xlsx(combined.sct, c(17), "reclustering_GO_new_res/clust17/06/descs/", "clust17_GO_06",0.6,featuresdotplot =all_subset_marks)
clust28.markers.both <- Markers2xlsx(combined.sct, c(28), "reclustering_GO_new_res/clust28/03/", "clust28_GO_0.3", 0.3, featuresdotplot=bcellfeats)

#likely T-cells
tcells1.markers.both <- Markers2xlsx(combined.sct, c(0,1,23,25), "reclustering_GO_new_res/tcells1/0.3/", "tcells1_GO_0.3",0.3,featuresdotplot=tcellfeats)
tcells2.markers.both <- Markers2xlsx(combined.sct, c(3,8,11,12,13,14,15,16,24,30), "reclustering_GO_new_res/tcells2/03/", "tcells2_GO_03", 0.3,featuresdotplot=tcellfeats)#- samma
tcells2.markers.both <- Markers2xlsx(combined.sct, c(3,8,11,12,13,14,15,16,24,30), "reclustering_GO_new_res/tcells2/01/", "tcells2_GO_01", 0.1,featuresdotplot=tcellfeats)#- samma
tcells2.markers.both <- Markers2xlsx(combined.sct, c(3,8,11,12,13,14,15,16,24,30), "reclustering_GO_new_res/tcells2/02/", "tcells2_GO_02", 0.2,featuresdotplot=tcellfeats)#- samma


tcells3.markers.both <- Markers2xlsx(combined.sct, c(7,29), "reclustering_GO_new_res/tcells3/", "tcells3_GO_0.2", 0.2)# - samma

tcells_clust4.markers.both <- Markers2xlsx(combined.sct, c(4), "reclustering_GO_new_res/tcells_clust4/", "tcells_clust4_GO", 0.3)

tcells_all.markers.both <- Markers2xlsx(combined.sct, c(0,1,23,25,3,8,11,12,13,14,15,16,24,30,7,29,4), "reclustering_GO_new_res/tcells2_all/03/", "tcells_all_GO_03", 0.3,featuresdotplot=tcellfeats)#- samma


#likely b cells
bcells1.markers.both <- Markers2xlsx(combined.sct, c(2,10,17,18,19,20,21,26), "reclustering_GO_new_res/bcells1/03_all_clusters/", "bcells1_GO_0.3", 0.3,featuresdotplot=bcellfeats)

clust19.markers.both <- Markers2xlsx(combined.sct, c(19), "reclustering_GO_new_res/clust19/03/descs/", "clust19_GO_03", 0.3,featuresdotplot=all_subset_marks) #- samma

#likely monocytes
Mono.markers.both <- Markers2xlsx(combined.sct, c(5,6), "reclustering_GO_new_res/monocytes/0.4/", "monocytes_GO_0.4",0.4,featuresdotplot=monocytefeats)
Mono.markers.both <- Markers2xlsx(combined.sct, c(5,6), "reclustering_GO_new_res/monocytes/0.6/", "monocytes_GO_0.6",0.6,featuresdotplot=monocytefeats)
Mono.markers.both <- Markers2xlsx(combined.sct, c(5,6), "reclustering_GO_new_res/monocytes/0.7/", "monocytes_GO_0.7",0.7,featuresdotplot=monocytefeats)

#likely thrombocytes
clust9.markers.both <- Markers2xlsx(combined.sct, c(9), "reclustering_GO_new_res/clust9/", "clust9_GO_0.3",0.3)
