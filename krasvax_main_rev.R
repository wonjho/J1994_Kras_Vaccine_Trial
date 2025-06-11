#Won Jin Ho
#Last updated 3/14/2024
#R version 4.0.2 (2020-06-22)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 18.04.5 LTS

##Environment Setup

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
workd<-getwd()
configd<-paste0(workd,'/Config')


library(readxl)
library(reshape2)
library(Hmisc)
library(stringr)
library(matrixStats)
library(dplyr)
library(scales)

library(flowCore)
library(FlowSOM)
library(umap)

library(pals)
library(RColorBrewer)
library(ggplot2)
library(ggridges)
library(ggrepel)

library(ComplexHeatmap)
library(pheatmap)
library(circlize)
library(corrplot)

library(edgeR)
library(ggpubr)
library(tidyr)


####READ and CLUSTER FUNCTIONS####
returnfcs <- function(FDR_cutoff=.05,
                      metaDataFile='~/',
                      panelDataFile='~/',
                      dataDirectory='~/',
                      shape_timepoint=NULL,
                      color_timepoint=NULL){
  #This function generates an fcs file, subtype_markers, colors and shapes for clustering 
  require(scales);
  require(readxl);
  require(dplyr);
  require(flowCore)
  ##directory and metadatafile checking
  if(!dir.exists(dataDirectory)) {stop('ERR: cannot find data directory')}
  if(!file.exists(metaDataFile)) {stop('ERR: cannot find metadata.xlsx or .csv file')}
  ##readin metadata and clean
  ifelse(grepl(metaDataFile,pattern='.xls'),
         md <- read_excel(metaDataFile),
         md <- read.csv(metaDataFile,header = TRUE))#must be in xls format or csv
  md$condition <- factor(md$condition)
  md$stim <- factor(md$stim)
  md$patient_id <- factor(md$patient_id)
  md$batch <- factor(md$batch)
  md$timepoint <- factor(md$timepoint)
  md$response <- factor(md$response)
  md$OS <- factor(md$OS)
  rownames(md) = md$sample_id;md$sample_id <- md$sample_id
  #Make sure all files in metadata present in datadirectory
  if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])){
    print(paste('ERR: not all filenames in metadata present in data folder - missing',md$file_name[!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])],'Subsetting...'))
    md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])),]
  }
  ##Define shapes for timepoint
  if(is.null(shape_timepoint)){shape_timepoint <- c(0:25)[1:length(levels(md$timepoint))]}#can specify as long as number is same
  if(length(shape_timepoint)!=length(levels(md$timepoint))){stop(paste0('ERR no. shapes specified is less than no. of timepoint (',length(levels(md$timepoint)),')'))}
  names(shape_timepoint) <- levels(md$timepoint)
  ## Define colors for the timepoint
  if(is.null(color_timepoint)){color_timepoint <- hue_pal()(length(levels(md$timepoint)))}#can specify as long as number is same
  if(length(color_timepoint)!=length(levels(md$timepoint))){stop(paste0('ERR no. shapes specified is less than no. of timepoint (',length(levels(md$timepoint)),')'))}
  ## read fcs
  fcs_raw <- read.flowSet(md$file_name, path = dataDirectory, transformation = FALSE, truncate_max_range = FALSE)
  sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
  panel <- read_excel(panelDataFile)
  head(data.frame(panel))
  ## Replace problematic characters
  panel$Metal <- gsub('-', '_', panel$Metal)
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  panel_fcs$desc <- gsub('-', '_', panel_fcs$desc)
  panel_fcs$desc[is.na(panel_fcs$desc)] <- paste0('NA_',which(is.na(panel_fcs$desc))) #was labelled 'todo'(mapping based on isotope for now just getting rid of NA keeping rownum) 
  # use panel$Antigen to fix description in panel_fcs
  # use metal+isotope as mapping between panel from xlsx and panel from the fcs files
  rownames(panel_fcs) = panel_fcs$name
  panel_fcs[paste0(panel$Metal,panel$Isotope,'Di'),2] <- panel$Antigen
  ## Replace paramater data in flowSet
  pData(parameters(fcs_raw[[1]])) <- panel_fcs
  ## Define variables indicating marker types
  subtype_markers <- panel$Antigen[panel$Subtype == 1]
  functional_markers <- panel$Antigen[panel$Functional == 1]
  if(!all(subtype_markers %in% panel_fcs$desc)){stop('ERR: Not all subtype_markers in panel_fcs$desc (isotopes)')}
  if(!all(functional_markers %in% panel_fcs$desc)){stop('ERR: Not all functional_markers in panel_fcs$desc (isotopes)')}
  ## arcsinh transformation and column subsetting
  fcs <- fsApply(fcs_raw, function(x, cofactor = 5){
    colnames(x) <- panel_fcs$desc
    expr <- exprs(x)
    expr <- asinh(expr[, union(subtype_markers,functional_markers)] / cofactor)
    exprs(x) <- expr
    x
  })
  return(list('fcs'=fcs,
              'subtype_markers'=subtype_markers,
              'functional_markers'=functional_markers,
              'shape_timepoint'=shape_timepoint,
              'color_timepoint'=color_timepoint,
              'sample_ids'=sample_ids,
              'meta_data'=md))
}

clusterfcs <- function(fcs=output$fcs,
                       subtype_markers = output$subtype_markers,
                       seed=1234,plottitle='consensus_plots',
                       numclusters=40){
  ## Cell population identification with FlowSOM and ConsensusClusterPlus
  require(dplyr);require(FlowSOM);require(ConsensusClusterPlus)
  set.seed(seed)
  som <- ReadInput(fcs, transform = FALSE, scale = FALSE) %>% BuildSOM(colsToUse = subtype_markers)
  ## Get the cell clustering into 100 SOM codes
  cell_clustering_som <- som$map$mapping[,1]
  ## Metaclustering into numclusters with ConsensusClusterPlus
  codes <- som$map$codes
  mc <- ConsensusClusterPlus(t(codes), maxK = numclusters, reps = 100,
                             pItem = 0.9, pFeature = 1, title = plottitle, 
                             plot = "png", clusterAlg = "hc", 
                             innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = 1234)
  
  ## Get cluster ids for each cell
  code_clustering <- mc[[numclusters]]$consensusClass#metaclusters consensus
  cell_clustering <- code_clustering[cell_clustering_som]#cell clustering from som
  return(list('code_clustering'=code_clustering,'cell_clustering'=cell_clustering,'metaclusters'=mc))
}


####CLUSTER HEATMAP FUNCTIONS ####
plot_clustering_heatmap_wrapper <- function(fcs, cell_clustering, nclusters=40,
                                            color_clusters=clustercolors2, cluster_merging = NULL, 
                                            subtype_markers,
                                            clusterMergeFile=NULL){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.05, 0.95))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(list(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(list(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  color_heat <- colorRampPalette(brewer.pal(n = 9, name = "YlOrBr"))(100)
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_mean$cell_clustering, " (", clustering_prop ,
                       "%)")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  ## Annotation for the merged clusters
  
  if(!is.null(clusterMergeFile)){
    ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Merged <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Merged <- color_clusters2
  }
  
  p <- pheatmap(expr_heat, color = rev(colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(100)), 
                cluster_cols = T,
                cluster_rows = T, 
                labels_row = labels_row,
                #scale="column",
                display_numbers = FALSE, number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8
  )
  print('Colors:')
  print(color_clusters)
  return(p)
}

plot_clustering_heatmap_wrapper2 <- function(fcs, cell_clustering, nclusters=40,
                                             color_clusters=clustercolors,
                                             subtype_markers){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.05, 0.95))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  
  ## Calculate the mean expression##################################################
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(list(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(list(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  #labels_row <- expr01_mean$cell_clustering
  
  labels_row <- paste0(expr01_mean$cell_clustering, " ")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  p <- pheatmap(expr_heat, 
                color = rev(colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(100)), 
                cluster_cols = T,
                cluster_rows = F, 
                labels_row = labels_row,
                #scale="column",
                display_numbers = F, 
                number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8,
                border_color = NA,
                annotation_legend = F
  )
  print('Colors:')
  print(color_clusters)
  return(p)
}





#======================
#     RUNNING DATA
#======================

####DATA LOADING AND CLUSTERING####

#from a previous save

output<-readRDS("backup_output.rds")
output2<-readRDS("backup_output2r.rds")     
                     
#read and cluster -- skip if loading prior save

output <- returnfcs(metaDataFile = paste0(configd,"/metadata.xlsx"), #output includes positive (dynabead) controls
                    panelDataFile = paste0(configd,"/panel.xlsx"),
                    dataDirectory = paste0(workd,'/Data'))
output[8:10] <- clusterfcs(fcs = output$fcs,
                           subtype_markers = output$subtype_markers, 
                           numclusters = 50)
names(output)[8:10] <- c('code_clustering','cell_clustering','metaclusters')


output2 <- returnfcs(metaDataFile = paste0(configd,"/metadata3.xlsx"), #output2 excludes positive (dynabead) controls
                    panelDataFile = paste0(configd,"/panel.xlsx"),
                    dataDirectory = paste0(workd,'/Data'))
output2[8:10] <- clusterfcs(fcs = output2$fcs,
                            subtype_markers = output2$subtype_markers, 
                            numclusters = 50)
names(output2)[8:10] <- c('code_clustering','cell_clustering','metaclusters')


#load cell cluster annotations

clusterMergeFile = paste0(configd,'/merged.xlsx')
cluster_merging <- read_xlsx(clusterMergeFile)

clusterMergeFile2 = paste0(configd,'/merged3r_final.xlsx')
cluster_merging2 <- read_xlsx(clusterMergeFile2)

#set up factor levels and colors

clusterlevels = c(1:50)
samplevels<-output$meta_data$sample_id
stimlevels<-c("UN","FLC","V","A","R","C","D","D13","WT","POS")
timelevels=c("C1","C3","C8","C16","C18")
ptlevels=levels(factor(output$meta_data$patient_id))
clustercolors <- as.character(c(cols25(n=25),alphabet(n=19),alphabet2(6)))
stimcolors<-cols25(length(stimlevels))
timecolors<-brewer.set2(length(timelevels))
ptcolors<-as.vector(polychrome(length(ptlevels)))

clusterlevels2 = c(1:50)
samplevels2<-output2$meta_data$sample_id
stimlevels2<-c("UN","FLC","V","A","R","C","D","D13","WT")
timelevels2=c("C1","C3","C8","C16","C18")
ptlevels2=levels(factor(output2$meta_data$patient_id))
clustercolors2 <- as.character(c(cols25(n=25),alphabet(n=19),alphabet2(6)))
stimcolors2<-c("#CCCCCC","#666666","#FB8808","#D84420","#690002","#7F32BD","#5770FF","#87BEFF","black")
timecolors2<-brewer.set2(length(timelevels2))
ptcolors2<-as.vector(polychrome(length(ptlevels2)))

responselevels <- c("NA","short","long")

#match up annotations

mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)
cell_clustering1m <- cluster_merging$new_cluster[mm1]
output$cell_clustering1m <- cell_clustering1m

mm2 <- match(output2$cell_clustering, cluster_merging2$original_cluster)
cell_clustering1m <- cluster_merging2$new_cluster[mm2]
cell_clustering1md <- cluster_merging2$detailed_cluster[mm2]
output2$cell_clustering1m <- cell_clustering1m #broad
output2$cell_clustering1md <- cell_clustering1md #detailed

mergedlevels2 <- c("CD4","CD8","DNT","NKT","NK","B","Mono")
mergedlevels3d <- c("CD4 N","CD4 Tcm_I","CD4 Tcm_II","CD4 Tcm_Acv","CD4 Tem_I","CD4 Tem_II","CD4 Tem_Acv",
                    "CD8 N","CD8 Tem","CD8 Tem_Acv","CD8 Teff_I","CD8 Teff_II","CD8 Teff_Acv",
                    "DNT","Mono","Mono_CD16","NK","NKT","B","B mem")

cd8clusters <- c("CD8 N","CD8 Tem","CD8 Tem_Acv","CD8 Teff_I","CD8 Teff_II","CD8 Teff_Acv")
cd4clusters <- c("CD4 N","CD4 Tcm_I","CD4 Tcm_II","CD4 Tcm_Acv","CD4 Tem_I","CD4 Tem_II","CD4 Tem_Acv")

#save clustered output

saveRDS(output, file = "backup_output.rds")
saveRDS(output2, file = "backup_output2r.rds")

####DIAGNOSTIC PLOTS####

#spot check - number of cells per sample

cell_table <- table(output$sample_ids)
color_batch <-c("red","blue","purple","cyan","yellow","orange")
ggdf <- data.frame(sample_id = names(cell_table), cell_counts = as.numeric(cell_table))
mm <- match(ggdf$sample_id, output$meta_data$sample_id)
ggdf$patient_id <- output$meta_data$patient_id[mm]
ggdf$timepoint <- output$meta_data$timepoint[mm]
ggdf$batch <- output$meta_data$batch[mm]

ggdf <- ggdf[ggdf$patient_id %nin% c("P003","P004"),] #get rid of CRC for this run

ggp<- ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = batch)) + 
  geom_bar(stat = 'identity') + 
  theme_bw() +
  theme(axis.text.x = element_blank()) +  
  scale_fill_manual(values = color_batch, drop = FALSE) + 
  scale_x_discrete(drop = FALSE)
pdf("plot_cellcounts.pdf",width=12,height=6);ggp;dev.off()

#extract expression data and normalize to 0 - 1 based on 1st and 99th percentile

expr <- fsApply(output$fcs, exprs)
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1])) #scaling 0-1
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

#generate per-sample level MDS plots

expr_mean_sample_tbl <- data.frame(sample_id = output$sample_ids, expr01) %>%
  group_by(sample_id) %>%  summarize_all(funs(mean))

expr_mean_sample <- t(expr_mean_sample_tbl[, -1])
colnames(expr_mean_sample) <- expr_mean_sample_tbl$sample_id

mds <- plotMDS(expr_mean_sample, plot = FALSE)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                   sample_id = colnames(expr_mean_sample))
mm <- match(ggdf$sample_id, output$meta_data$sample_id)
ggdf$timepoint <- output$meta_data$timepoint[mm]
ggdf$batch <- factor(output$meta_data$batch[mm])
ggdf$patient_id <- factor(output$meta_data$patient_id[mm])
ggdf$stim <- factor(output$meta_data$stim[mm], levels=stimlevels)

expr_mean_sample_nonPOS<-expr_mean_sample[,output$meta_data[output$meta_data$stim!="POS",]$sample_id]

mds2 <- plotMDS(expr_mean_sample_nonPOS, plot = F)
ggdf2 <- data.frame(MDS1 = mds2$x, MDS2 = mds2$y,
                    sample_id = colnames(expr_mean_sample_nonPOS))
mm2 <- match(ggdf2$sample_id, output$meta_data$sample_id)
ggdf2$timepoint <- output$meta_data$timepoint[mm2]
ggdf2$batch <- factor(output$meta_data$batch[mm2])
ggdf2$patient_id <- factor(output$meta_data$patient_id[mm2])
ggdf2$stim <- factor(output$meta_data$stim[mm2], levels=stimlevels)

ggdf2 <- ggdf2[ggdf2$patient_id %nin% c("P003","P004"),] #get rid of CRC for this run

mdsggp <- ggplot(ggdf, aes(x = MDS1, y = MDS2, color = stim, shape = timepoint)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_bw() +
  theme(axis.text = element_text(color="black")) +
  #ylim(-0.1, 0.08) +
  scale_shape_manual(values = c(2,0,1,3,4,5,6,7)) +
  coord_fixed()
mdsggp2 <- ggplot(ggdf2, aes(x = MDS1, y = MDS2, color = patient_id, shape = batch)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_bw() +
  theme(axis.text = element_text(color="black")) +
  #ylim(-0.1, 0.08) +
  scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8)) +
  coord_fixed()

pdf("plot_mds.pdf", width=4, height=5);mdsggp;mdsggp2;dev.off()

#for cluster-specific mean expression heatmaps using output2 detailed clusters

markerlist_plot <- c('CCR7', 'CD45RA', 'CD45RO', 'CD137', 'CD25', 'CD28', 'OX40', 'IFNg', 'IL2', 'TNFa', 'GZMB', 'KI67', 'PD1', 'CTLA4', 'LAG3', 'TIM3', 'CXCR3', 'CCR5')
markerlist_plot2 <- c('CD137', 'CD25', 'CD28', 'OX40', 'IFNg', 'IL2', 'TNFa', 'GZMB', 'KI67', 'PD1', 'CTLA4', 'LAG3', 'TIM3', 'CXCR3', 'CCR5','CD45RO','CD45RA')

expr <- fsApply(output2$fcs, exprs);expr <-expr[,markerlist_plot]
rng <- colQuantiles(expr, probs = c(0.05, 0.95))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,markerlist_plot]

expr_mean <- data.frame(expr, cell_clustering = output2$cell_clustering1md, check.names = FALSE) %>%
  group_by(cell_clustering) %>% summarize_all(list(mean))
expr01_mean <- data.frame(expr01, cell_clustering = output2$cell_clustering1md, check.names = FALSE) %>%
  group_by(cell_clustering) %>% summarize_all(list(mean))

expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
rownames(expr_heat) <- expr01_mean$cell_clustering

p1 <- pheatmap(expr_heat[cd4clusters,markerlist_plot2], 
               main="CD4",
               color = rev(colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(100)), 
               cluster_cols = T,
               cluster_rows = F, 
               #labels_row = labels_row,
               scale="column",
               display_numbers = F, 
               number_color = "black",
               fontsize = 9, fontsize_number = 6,  
               #legend_breaks = legend_breaks,
               #annotation_row = annotation_row, 
               #annotation_colors = annotation_colors,
               cellwidth = 8,
               cellheight = 8,
               border_color = NA,
               annotation_legend = F
)

p2 <- pheatmap(expr_heat[cd8clusters,markerlist_plot2], 
               main="CD8",
               color = rev(colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(100)), 
               cluster_cols = T,
               cluster_rows = F, 
               #labels_row = labels_row,
               scale="column",
               display_numbers = F, 
               number_color = "black",
               fontsize = 9, fontsize_number = 6,  
               #legend_breaks = legend_breaks,
               #annotation_row = annotation_row, 
               #annotation_colors = annotation_colors,
               cellwidth = 8,
               cellheight = 8,
               border_color = NA,
               annotation_legend = F
)

dev.off()
pdf('clusteringheatmap_cd4_cd8_columnscaled_sub_revisedJan2025.pdf')
print(p1)
plot.new()
print(p2)
dev.off()

#for sample-level mean expression heatmap

mm <- match(colnames(expr_mean_sample[,samplevels]), output$meta_data$sample_id)
annotation_col <- data.frame(timepoint = factor(output$meta_data$timepoint[mm]),
                             patient_id = factor(output$meta_data$patient_id[mm]),
                             stim = factor(output$meta_data$stim[mm], levels = stimlevels),
                             row.names = colnames(expr_mean_sample[,samplevels]))

annotation_colors <- list(condition = clustercolors[1:length(levels(annotation_col$timepoint))],
                          stim = c(UN=stimcolors[1],FLC=stimcolors[2],V=stimcolors[3],A=stimcolors[4],R=stimcolors[5],
                                   C=stimcolors[6],D=stimcolors[7],D13=stimcolors[8],POS=stimcolors[9],WT=stimcolors[10]))

color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)

meanht<- pheatmap(expr_mean_sample[,samplevels], color = color, display_numbers = F,
                  #number_color = "black", fontsize_number = 5, 
                  cluster_cols = T,
                  show_colnames = F,
                  border_color = NA,
                  annotation_col = annotation_col, main = '',
                  annotation_colors = annotation_colors, clustering_method = "average")
dev.off()
pdf("plot_meanheatmap.pdf", width=12, height=7);meanht;dev.off()

#metacluster/annotated cluster-level heatmaps

p<-plot_clustering_heatmap_wrapper(fcs=output$fcs,
                                color_clusters = c(stepped2(20),stepped3(20),rev(cubicl(20))),
                                cell_clustering = output$cell_clustering, 
                                subtype_markers = output$subtype_markers,
                                clusterMergeFile = clusterMergeFile)
dev.off()
pdf("clusteringheatmap.pdf",width=8, height=11);p;dev.off()


p<-plot_clustering_heatmap_wrapper(fcs=output2$fcs,
                                   color_clusters = c(stepped2(20),stepped3(20),rev(cubicl(20))),
                                   cell_clustering = output2$cell_clustering, 
                                   subtype_markers = output2$subtype_markers,
                                   clusterMergeFile = clusterMergeFile2)
dev.off()
pdf("clusteringheatmap2r.pdf",width=8, height=11);p;dev.off()

p<- plot_clustering_heatmap_wrapper2(fcs=output2$fcs,
                                     color_clusters = clustercolors2,
                                     cell_clustering = factor(output2$cell_clustering1m, levels=mergedlevels2), 
                                     subtype_markers=output2$subtype_markers)
dev.off()
pdf("clusteringheatmap2r.pdf",width=8, height=11);p;dev.off()

p_d<- plot_clustering_heatmap_wrapper2(fcs=output2$fcs,
                                 color_clusters = clustercolors2,
                                 cell_clustering = factor(output2$cell_clustering1md, levels=mergedlevels3d), 
                                 subtype_markers=output2$subtype_markers)
dev.off()
pdf("clusteringheatmap2r_detailed.pdf",width=8, height=11);p_d;dev.off()

####DIFFERENTIAL PLOTS####

#set up count and prop matrices
counts_table <- table(output$cell_clustering1m, output$sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

write.csv(counts, file='table_counts.csv')
write.csv(props, file='table_props.csv')

counts_table2 <- table(output2$cell_clustering1m, output2$sample_ids)
props_table2 <- t(t(counts_table2) / colSums(counts_table2)) * 100
counts2 <- as.data.frame.matrix(counts_table2)
props2 <- as.data.frame.matrix(props_table2)

counts_table2d <- table(output2$cell_clustering1md, output2$sample_ids)
props_table2d <- t(t(counts_table2d) / colSums(counts_table2d)) * 100
counts2d <- as.data.frame.matrix(counts_table2d)
props2d <- as.data.frame.matrix(props_table2d)

write.csv(counts2, file='table_counts2.csv')
write.csv(props2, file='table_props2.csv')

clustermatching <- data.frame(unique(cluster_merging2[,c("new_cluster","detailed_cluster")]))
clustermatching <- with(clustermatching, clustermatching[order(detailed_cluster),])

CD4anno <- clustermatching[clustermatching$new_cluster=="CD4",]$detailed_cluster
counts_table2d_CD4 <- table(output2$cell_clustering1md[output2$cell_clustering1md %in% CD4anno], output2$sample_ids[output2$cell_clustering1md %in% CD4anno])
props_table2d_CD4 <- t(t(counts_table2d_CD4) / colSums(counts_table2d_CD4)) * 100
counts2d_CD4 <- as.data.frame.matrix(counts_table2d_CD4)
props2d_CD4 <- as.data.frame.matrix(props_table2d_CD4)

CD8anno <- clustermatching[clustermatching$new_cluster=="CD8",]$detailed_cluster
counts_table2d_CD8 <- table(output2$cell_clustering1md[output2$cell_clustering1md %in% CD8anno], output2$sample_ids[output2$cell_clustering1md %in% CD8anno])
props_table2d_CD8 <- t(t(counts_table2d_CD8) / colSums(counts_table2d_CD8)) * 100
counts2d_CD8 <- as.data.frame.matrix(counts_table2d_CD8)
props2d_CD8 <- as.data.frame.matrix(props_table2d_CD8)

Tanno <- clustermatching[clustermatching$new_cluster %in% c("CD8","CD4"),]$detailed_cluster
counts_table2d_T <- table(output2$cell_clustering1md[output2$cell_clustering1md %in% Tanno], output2$sample_ids[output2$cell_clustering1md %in% Tanno])
props_table2d_T <- t(t(counts_table2d_T) / colSums(counts_table2d_T)) * 100
counts2d_T <- as.data.frame.matrix(counts_table2d_T)
props2d_T <- as.data.frame.matrix(props_table2d_T)

#heatmaps of proportional abundance

  #to exclude CRC patients
P003levels<-!str_detect(samplevels2, "P003")
P004levels<-!str_detect(samplevels2, "P004")
samplevels2_pdac <- samplevels2[P003levels & P004levels]

  #for the run without dynabead stim and CRC patients:

annotation_colors2=list(stim = c(UN=stimcolors2[1],FLC=stimcolors2[2],V=stimcolors2[3],A=stimcolors2[4],R=stimcolors2[5],
                                 C=stimcolors2[6],D=stimcolors2[7],D13=stimcolors2[8],WT=stimcolors2[9]),
                        timepoint = c(C1=timecolors[1],C3=timecolors[2],C8=timecolors[3],C16=timecolors[4],C18=timecolors[5]),
                        patient_id = c(P001=ptcolors[1],P002=ptcolors[2],#P003=ptcolors[3],P004=ptcolors[4],
                                       P005=ptcolors[5],
                                       P006=ptcolors[6],P007=ptcolors[7],P008=ptcolors[8],P010=ptcolors[9],P012=ptcolors[10],P013=ptcolors[11],
                                       P009=ptcolors[12],P014=ptcolors[13],P016=ptcolors[14]))

  #all patients:
column_ha <- HeatmapAnnotation(patient_id=factor(output2$meta_data$patient_id,levels=ptlevels2),
                               stim=factor(output2$meta_data$stim,levels=stimlevels2),
                               timepoint=output2$meta_data$timepoint,
                               col=annotation_colors2)

  #for without CRC patients:
column_ha <- HeatmapAnnotation(patient_id=factor(ggdf2$patient_id,levels=ptlevels2),
                               stim=factor(ggdf2$stim,levels=stimlevels2),
                               timepoint=ggdf2$timepoint,
                               col=annotation_colors2)

pdf('plot_heatmap_props.pdf',width=10, height=8)
Heatmap(name="global", as.data.frame.matrix(props2[,samplevels2_pdac]),
        top_annotation = column_ha)
Heatmap(name="global_detailed", as.data.frame.matrix(props2d[,samplevels2_pdac]),
        top_annotation = column_ha)
Heatmap(name="CD4", as.data.frame.matrix(props2d_CD4[,samplevels2_pdac]),
        top_annotation = column_ha)
Heatmap(name="CD8", as.data.frame.matrix(props2d_CD8[,samplevels2_pdac]),
        top_annotation = column_ha)
dev.off()

#correlation among abundances (no CRC)
corr_propCD4<-cor(t(props2d_CD4[,samplevels2_pdac]))
corr_propCD8<-cor(t(props2d_CD8[,samplevels2_pdac]))
corr_propT<-cor(t(props2d_T[,samplevels2_pdac]))

pdf("plot_corr_prop.pdf", width=4, height=4)
par(oma=c(0,0,1,0), xpd=NA)
corrplot(corr_propT, col = rev(brewer.rdbu(100)), tl.col = "black", title = "All T",
         method = 'circle', order = 'hclust', hclust.method = 'ward.D', addrect = 4)
corrplot(corr_propCD4, col = rev(brewer.rdbu(100)), tl.col = "black", title = "CD4",
         method = 'circle', order = 'hclust', hclust.method = 'ward.D', addrect = 2)
corrplot(corr_propCD8, col = rev(brewer.rdbu(100)), tl.col = "black", title = "CD8",
         method = 'circle', order = 'hclust', hclust.method = 'ward.D', addrect = 2)
dev.off()

#set up the data frame for proportional plotting of all major T cell types only at the unstim condition

ggdf <- melt(data.frame(cluster = rownames(props2d),props2d, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels2)
ggdf$cluster <- factor(ggdf$cluster, levels=mergedlevels3d)
ggdf$timepoint <- factor(output2$meta_data$timepoint[match(ggdf$sample_id,output2$meta_data$sample_id)], levels=timelevels2)
ggdf$patient_id <- factor(output2$meta_data$patient_id[match(ggdf$sample_id,output2$meta_data$sample_id)])
ggdf$stim <- factor(output2$meta_data$stim[match(ggdf$sample_id,output2$meta_data$sample_id)], levels=stimlevels2)
ggdf$batch <- output2$meta_data$batch[match(ggdf$sample_id,output2$meta_data$sample_id)]

#remove CRC patients and keep only UN

ggdf<-ggdf[ggdf$patient_id %nin% c("P003","P004"),]
ggdf<-ggdf[ggdf$stim == "UN",]
ggdf$patient_id <- factor(ggdf$patient_id)
ptlevels_ggdf <- levels(ggdf$patient_id)

pdf('plot_proportion_line_unstim.pdf', width=8, height=8)
lp <- ggplot(ggdf, aes(x = timepoint, y = proportion, group = patient_id, color=patient_id)) +
  geom_line()+
  geom_point()+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(linewidth=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  xlab("Timepoint") +
  ylab("% of CD45 Cells")+
  facet_wrap(~cluster, ncol=4, scales="free")+
  #facet_grid(rows=vars(patient_id), cols= vars(cluster), scales="free")+
  scale_color_manual(values = ptcolors2) +
  #scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8,9,10)) +
  #scale_y_continuous(expand = c(0,.2))+
  guides(fill=guide_legend(ncol=1))
print(lp)
dev.off()

pdf('plot_proportion_line_unstim_bypt.pdf', width=8, height=8)
for(i in 1:length(ptlevels_ggdf)){
  ggdfsubset<-ggdf[ggdf$patient_id==ptlevels_ggdf[i],]
  lp <- ggplot(ggdfsubset, aes(x = timepoint, y = proportion, group = stim, shape=patient_id)) +
    ggtitle(ptlevels_ggdf[i])+
    geom_line()+
    geom_point()+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
          axis.text.y = element_text(color="black"),
          axis.ticks = element_line(linewidth=0.25),
          strip.text=element_text(size=8),
          strip.background = element_rect(fill=NA, color=NA),
          legend.title = element_blank(),
          legend.key.size = unit(.75,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white"),
    ) +
    xlab("Timepoint") +
    ylab("% of CD4 Cells")+
    facet_wrap(~cluster, ncol=4, scales="free")+
    #facet_grid(rows=vars(patient_id), cols= vars(cluster), scales="free")+
    scale_color_manual(values = stimcolors2[c(2:8)]) +
    scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8,9,10)) +
    #scale_y_continuous(expand = c(0,.2))+
    guides(fill=guide_legend(ncol=1))
  print(lp)
}
dev.off()

#set up the data frame for proportional plotting of CD4 T cells

ggdf <- melt(data.frame(cluster = rownames(props2d_CD4),props2d_CD4, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels2)
ggdf$cluster <- factor(ggdf$cluster, levels=mergedlevels3d)
ggdf$timepoint <- factor(output2$meta_data$timepoint[match(ggdf$sample_id,output2$meta_data$sample_id)], levels=timelevels2)
ggdf$patient_id <- factor(output2$meta_data$patient_id[match(ggdf$sample_id,output2$meta_data$sample_id)])
ggdf$stim <- factor(output2$meta_data$stim[match(ggdf$sample_id,output2$meta_data$sample_id)], levels=stimlevels2)
ggdf$batch <- output2$meta_data$batch[match(ggdf$sample_id,output2$meta_data$sample_id)]

#remove CRC patients
ggdf<-ggdf[ggdf$patient_id %nin% c("P003","P004"),]
ggdf<-ggdf[ggdf$stim %nin% c("UN","WT"),]
ggdf$patient_id <- factor(ggdf$patient_id)
ptlevels_ggdf <- levels(ggdf$patient_id)

#normalized ggdf
ggdf_n <- ggdf
stimn<-length(unique(ggdf_n$stim))

ggdf_n$comb <- paste(ggdf_n$patient_id, ggdf_n$timepoint, ggdf_n$cluster, "s", ggdf_n$stim, sep="_")
ggdf_n <-arrange(ggdf_n, comb)
ggdf_n_FLC <- ggdf_n[ggdf_n$stim=="FLC",]
ggdf_n_FLC <- ggdf_n_FLC %>% slice(rep(1:n(), each=stimn))
ggdf_n$comb2 <- ggdf_n_FLC$comb
ggdf_n$norm <- ggdf_n_FLC$proportion
ggdf_n$propn <- ggdf_n$proportion/ggdf_n$norm

write.csv(ggdf_n, 'table_CD4prop_norm.csv')

#also plot baseline differences of proportions normalized to FLC
ggdf_n$response <- factor(output2$meta_data$response[match(ggdf_n$sample_id,output2$meta_data$sample_id)], levels=responselevels)
ggdf_n_sub <- ggdf_n[ggdf_n$stim!="FLC",]
ggdf_n_sub <- ggdf_n_sub[is.finite(ggdf_n_sub$propn),]

ggdf_n_baseline <- ggdf_n_sub[ggdf_n_sub$timepoint=="C1",]
ggp_res_c1 <- ggplot(ggdf_n_baseline, aes(x=response, y=propn))+
  facet_wrap(~cluster, ncol=5, scales="free")+
  geom_boxplot(outlier.color=NA)+
  geom_jitter(aes(color=patient_id, shape=stim), width=0.2)+
  scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8,9,10))+
  stat_compare_means(label = "p.format",
                     label.x = 1.4,
                     label.y.npc = "top",
                     hide.ns = F,
                     size=1.5)
ggdf_n_c3 <- ggdf_n_sub[ggdf_n_sub$timepoint=="C3",]
ggp_res_c3 <- ggplot(ggdf_n_c3, aes(x=response, y=propn))+
  facet_wrap(~cluster, ncol=5, scales="free")+
  geom_boxplot(outlier.color=NA)+
  geom_jitter(aes(color=patient_id, shape=stim), width=0.2)+
  scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8,9,10))+
  stat_compare_means(label = "p.format",
                     label.x = 1.4,
                     label.y.npc = "top",
                     hide.ns = F,
                     size=1.5)
pdf('plot_proportion_box_res_CD4.pdf', width=10, height=10);ggp_res_c1;ggp_res_c3;dev.off()

pdf('plot_proportion_line_CD4.pdf', width=6, height=6)
for(i in 1:length(ptlevels_ggdf)){
  ggdfsubset<-ggdf[ggdf$patient_id==ptlevels_ggdf[i],]
  lp <- ggplot(ggdfsubset, aes(x = timepoint, y = proportion, color = stim, group = stim, shape=patient_id)) +
    ggtitle(ptlevels_ggdf[i])+
    geom_line()+
    geom_point()+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
          axis.text.y = element_text(color="black"),
          axis.ticks = element_line(linewidth=0.25),
          strip.text=element_text(size=8),
          strip.background = element_rect(fill=NA, color=NA),
          legend.title = element_blank(),
          legend.key.size = unit(.75,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white"),
    ) +
    xlab("Timepoint") +
    ylab("% of CD4 Cells")+
    facet_wrap(~cluster, ncol=3, scales="free")+
    #facet_grid(rows=vars(patient_id), cols= vars(cluster), scales="free")+
    scale_color_manual(values = stimcolors2[c(2:8)]) +
    scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8,9,10)) +
    #scale_y_continuous(expand = c(0,.2))+
    guides(fill=guide_legend(ncol=1))
  print(lp)
}
dev.off()

pdf('plot_proportion_line_CD4_norm.pdf', width=6, height=6)
for(i in 1:length(ptlevels_ggdf)){
  ggdfsubset<-ggdf_n[ggdf_n$patient_id==ptlevels_ggdf[i],]
  lp <- ggplot(ggdfsubset, aes(x = timepoint, y = propn, color = stim, group = stim, shape=patient_id)) +
    ggtitle(ptlevels_ggdf[i])+
    geom_line()+
    geom_point()+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
          axis.text.y = element_text(color="black"),
          axis.ticks = element_line(linewidth=0.25),
          strip.text=element_text(size=8),
          strip.background = element_rect(fill=NA, color=NA),
          legend.title = element_blank(),
          legend.key.size = unit(.75,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white"),
    ) +
    xlab("Timepoint") +
    ylab("% of CD4 Cells (normalized by FLC)")+
    facet_wrap(~cluster, ncol=3, scales="free")+
    #facet_grid(rows=vars(patient_id), cols= vars(cluster), scales="free")+
    scale_color_manual(values = stimcolors2[c(2:8)]) +
    scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8,9,10)) +
    #scale_y_continuous(expand = c(0,.2))+
    guides(fill=guide_legend(ncol=1))
  print(lp)
}
dev.off()

#set up the data frame for proportional plotting of CD8 T cells
ggdf <- melt(data.frame(cluster = rownames(props2d_CD8),props2d_CD8, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels2)
ggdf$cluster <- factor(ggdf$cluster, levels=mergedlevels3d)
ggdf$timepoint <- factor(output2$meta_data$timepoint[match(ggdf$sample_id,output2$meta_data$sample_id)], levels=timelevels2)
ggdf$patient_id <- factor(output2$meta_data$patient_id[match(ggdf$sample_id,output2$meta_data$sample_id)])
ggdf$stim <- factor(output2$meta_data$stim[match(ggdf$sample_id,output2$meta_data$sample_id)], levels=stimlevels2)
ggdf$batch <- output2$meta_data$batch[match(ggdf$sample_id,output2$meta_data$sample_id)]

#remove CRC patients
ggdf<-ggdf[ggdf$patient_id %nin% c("P003","P004"),]
ggdf<-ggdf[ggdf$stim %nin% c("UN","WT"),]
ggdf$patient_id <- factor(ggdf$patient_id)
ptlevels_ggdf <- levels(ggdf$patient_id)

#normalized ggdf
ggdf_n <- ggdf
stimn<-length(unique(ggdf_n$stim))

ggdf_n$comb <- paste(ggdf_n$patient_id, ggdf_n$timepoint, ggdf_n$cluster, "s", ggdf_n$stim, sep="_")
ggdf_n <-arrange(ggdf_n, comb)
ggdf_n_FLC <- ggdf_n[ggdf_n$stim=="FLC",]
ggdf_n_FLC <- ggdf_n_FLC %>% slice(rep(1:n(), each=stimn))
ggdf_n$comb2 <- ggdf_n_FLC$comb
ggdf_n$norm <- ggdf_n_FLC$proportion
ggdf_n$propn <- ggdf_n$proportion/ggdf_n$norm

write.csv(ggdf_n, 'table_CD8prop_norm.csv')

ggdf_n$response <- factor(output2$meta_data$response[match(ggdf_n$sample_id,output2$meta_data$sample_id)], levels=responselevels)
ggdf_n_sub <- ggdf_n[ggdf_n$stim!="FLC",]
ggdf_n_sub <- ggdf_n_sub[is.finite(ggdf_n_sub$propn),]

ggdf_n_baseline <- ggdf_n_sub[ggdf_n_sub$timepoint=="C1",]
ggp_res_c1 <- ggplot(ggdf_n_baseline, aes(x=response, y=propn))+
  facet_wrap(~cluster, ncol=5, scales="free")+
  geom_boxplot(outlier.color=NA)+
  geom_jitter(aes(color=patient_id, shape=stim), width=0.2)+
  scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8,9,10))+
  stat_compare_means(label = "p.format",
                     label.x = 1.4,
                     label.y.npc = "top",
                     hide.ns = F,
                     size=1.5)  
ggdf_n_c3 <- ggdf_n_sub[ggdf_n_sub$timepoint=="C3",]
ggp_res_c3 <- ggplot(ggdf_n_c3, aes(x=response, y=propn))+
  facet_wrap(~cluster, ncol=5, scales="free")+
  geom_boxplot(outlier.color=NA)+
  geom_jitter(aes(color=patient_id, shape=stim), width=0.2)+
  scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8,9,10))+
  stat_compare_means(label = "p.format",
                     label.x = 1.4,
                     label.y.npc = "top",
                     hide.ns = F,
                     size=1.5)  
pdf('plot_proportion_box_res_CD8.pdf', width=10, height=10);ggp_res_c1;ggp_res_c3;dev.off()

pdf('plot_proportion_line_CD8.pdf', width=6, height=6)
for(i in 1:length(ptlevels_ggdf)){
  ggdfsubset<-ggdf[ggdf$patient_id==ptlevels_ggdf[i],]
  lp <- ggplot(ggdfsubset, aes(x = timepoint, y = proportion, color = stim, group = stim, shape=patient_id)) +
    ggtitle(ptlevels_ggdf[i])+
    geom_line()+
    geom_point()+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
          axis.text.y = element_text(color="black"),
          axis.ticks = element_line(linewidth=0.25),
          strip.text=element_text(size=8),
          strip.background = element_rect(fill=NA, color=NA),
          legend.title = element_blank(),
          legend.key.size = unit(.75,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white"),
    ) +
    xlab("Timepoint") +
    ylab("% of CD8 Cells")+
    facet_wrap(~cluster, ncol=3, scales="free")+
    #facet_grid(rows=vars(patient_id), cols= vars(cluster), scales="free")+
    scale_color_manual(values = stimcolors2[c(2:8)]) +
    scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8,9,10)) +
    #scale_y_continuous(expand = c(0,.2))+
    guides(fill=guide_legend(ncol=1))
  print(lp)
}
dev.off()

pdf('plot_proportion_line_CD8_norm.pdf', width=6, height=6)
for(i in 1:length(ptlevels_ggdf)){
  ggdfsubset<-ggdf_n[ggdf_n$patient_id==ptlevels_ggdf[i],]
  lp <- ggplot(ggdfsubset, aes(x = timepoint, y = propn, color = stim, group = stim, shape=patient_id)) +
    ggtitle(ptlevels_ggdf[i])+
    geom_line()+
    geom_point()+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
          axis.text.y = element_text(color="black"),
          axis.ticks = element_line(linewidth=0.25),
          strip.text=element_text(size=8),
          strip.background = element_rect(fill=NA, color=NA),
          legend.title = element_blank(),
          legend.key.size = unit(.75,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white"),
    ) +
    xlab("Timepoint") +
    ylab("% of CD8 Cells (normalized by FLC)")+
    facet_wrap(~cluster, ncol=3, scales="free")+
    #facet_grid(rows=vars(patient_id), cols= vars(cluster), scales="free")+
    scale_color_manual(values = stimcolors2[c(2:8)]) +
    scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8,9,10)) +
    #scale_y_continuous(expand = c(0,.2))+
    guides(fill=guide_legend(ncol=1))
  print(lp)
}
dev.off()

####FUNCTIONAL MARKER####

markerlist = output2$functional_markers

exprtbl <- data.frame(fsApply(output2$fcs,exprs)[, markerlist], 
                      sample_id = output2$sample_ids, 
                      cluster = output2$cell_clustering1m)

#corrplots

exprtbl_CD4 <- exprtbl[exprtbl$cluster=="CD4",]
corr_CD4 <- cor(exprtbl_CD4[,markerlist])

exprtbl_CD8 <- exprtbl[exprtbl$cluster=="CD8",]
corr_CD8 <- cor(exprtbl_CD8[,markerlist])

exprtbl_T <- exprtbl[exprtbl$cluster %in% c("CD8","CD4"),]
corr_T <- cor(exprtbl_T[,markerlist])

pdf("plot_corr_expr.pdf", width=4, height=4)
par(oma=c(0,0,1,0), xpd=NA)
corrplot(corr_T, col = rev(brewer.rdbu(100)), tl.col = "black", title = "All T",
         method = 'circle', order = 'hclust', hclust.method = 'ward.D2', addrect = 4)
corrplot(corr_CD4, col = rev(brewer.rdbu(100)), tl.col = "black", title = "CD4",
         method = 'circle', order = 'hclust', hclust.method = 'ward.D2', addrect = 4)
corrplot(corr_CD8, col = rev(brewer.rdbu(100)), tl.col = "black", title = "CD8",
         method = 'circle', order = 'hclust', hclust.method = 'ward.D2', addrect = 4)
dev.off()

#boxplots at a per patient basis- based on detailed clusters

markerlist_plot = c("TNFa","IL2","IFNg","KI67","CD137","GZMB","CD25","PD1", "CTLA4", "LAG3", "TIM3", "TIGIT")

exprtbl_mn <- data.frame(fsApply(output2$fcs,exprs)[, markerlist_plot], 
                         sample_id = output2$sample_ids, 
                         cluster = output2$cell_clustering1md) %>%
  group_by(sample_id, cluster) %>%
  summarize_all(funs(median))

ggdf2_mn<-melt(exprtbl_mn, id.var=c("cluster","sample_id"))
ggdf2_mn$cluster <- factor(ggdf2_mn$cluster, levels=mergedlevels3d)
ggdf2_mn$sample_id <- factor(ggdf2_mn$sample_id, levels=unique(ggdf2_mn$sample_id))
ggdf2_mn$condition <- factor(output2$meta_data$condition[match(ggdf2_mn$sample_id,output2$meta_data$sample_id)])
ggdf2_mn$timepoint <- factor(output2$meta_data$timepoint[match(ggdf2_mn$sample_id,output2$meta_data$sample_id)],levels = timelevels2)
ggdf2_mn$patient_id <- factor(output2$meta_data$patient_id[match(ggdf2_mn$sample_id,output2$meta_data$sample_id)])
ggdf2_mn$stim <- factor(output2$meta_data$stim[match(ggdf2_mn$sample_id,output2$meta_data$sample_id)],levels = stimlevels2)
ggdf2_mn$response <- factor(output2$meta_data$response[match(ggdf2_mn$sample_id,output2$meta_data$sample_id)],levels = responselevels)

ggdf2_mn<-ggdf2_mn[ggdf2_mn$patient_id %nin% c("P003","P004"),]
ggdf2_mn<-ggdf2_mn[ggdf2_mn$stim %nin% c("UN","WT"),]
ggdf2_mn$patient_id <- factor(ggdf2_mn$patient_id)
ggdf2_mn_sub <- ggdf2_mn[ggdf2_mn$stim!="FLC",]

#cd8 cells

ggdf2_mn_sub_cd8 <- ggdf2_mn_sub[ggdf2_mn_sub$cluster %in% cd8clusters,]

#choose timepoint here

#ggdf2_mn_sub_cd8 <- ggdf2_mn_sub_cd8[ggdf2_mn_sub_cd8$timepoint == "C1",] 
ggdf2_mn_sub_cd8 <- ggdf2_mn_sub_cd8[ggdf2_mn_sub_cd8$timepoint == "C3",]

pdf("plot_functional_box_CD8_C3_med.pdf",width=6,height=6)
for (i in 1:length(markerlist_plot)){
  ggdfsubset<-ggdf2_mn_sub_cd8[ggdf2_mn_sub_cd8$variable==markerlist_plot[i],]
  ggp <- ggplot(ggdfsubset, aes(x=response, y=value))+
    ggtitle(markerlist_plot[i])+
    facet_wrap(~cluster, scales='free')+
    geom_boxplot(outlier.color = NA)+
    geom_jitter(aes(color=patient_id, shape=stim), width=0.2)+
    scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8,9,10))+
    #scale_fill_manual(values = stimcolors2[c(2:8)]) +
    theme(axis.text.x = element_text(angle = 90, vjust=.5, hjust =1, color="black", size=8),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(linewidth=0.25),
          axis.line = element_line(linewidth=0.25),
          axis.text = element_text(color="black"),
          axis.title.y = element_text(size=10, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=10),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(1,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white"))+
    stat_compare_means(label = "p.format",
                       label.x = 1.4,
                       label.y.npc = "top",
                       hide.ns = F,
                       size=1.5)  
  print(ggp)
}
dev.off()

#cd4 cells

ggdf2_mn_sub_cd4 <- ggdf2_mn_sub[ggdf2_mn_sub$cluster %in% cd4clusters,]

#choose timepoint here

ggdf2_mn_sub_cd4 <- ggdf2_mn_sub_cd4[ggdf2_mn_sub_cd4$timepoint == "C1",]
ggdf2_mn_sub_cd4 <- ggdf2_mn_sub_cd4[ggdf2_mn_sub_cd4$timepoint == "C3",]

pdf("plot_functional_box_CD4_C3_med.pdf",width=6,height=6)
for (i in 1:length(markerlist_plot)){
  ggdfsubset<-ggdf2_mn_sub_cd4[ggdf2_mn_sub_cd4$variable==markerlist_plot[i],]
  ggp <- ggplot(ggdfsubset, aes(x=response, y=value))+
    ggtitle(markerlist_plot[i])+
    facet_wrap(~cluster, scales='free')+
    geom_boxplot(outlier.color = NA)+
    geom_jitter(aes(color=patient_id, shape=stim), width=0.2)+
    scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8,9,10))+
    #scale_fill_manual(values = stimcolors2[c(2:8)]) +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(linewidth=0.25),
          axis.line = element_line(linewidth=0.25),
          axis.text = element_text(color="black"),
          axis.title.y = element_text(size=10, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=10),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(1,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white"))+
    stat_compare_means(label = "p.format",
                       label.x = 1.4,
                       label.y.npc = "top",
                       hide.ns = F,
                       size=1.5)  
  print(ggp)
}
dev.off()

#boxplots at a per patient basis and also averaged across stims- based on detailed clusters

exprtbl_mn <- data.frame(fsApply(output2$fcs,exprs)[, markerlist_plot], 
                         sample_id = output2$sample_ids, 
                         cluster = output2$cell_clustering1md)
exprtbl_mn$cluster <- factor(exprtbl_mn$cluster, levels=mergedlevels3d)
exprtbl_mn$sample_id <- factor(exprtbl_mn$sample_id, levels=unique(exprtbl_mn$sample_id))
exprtbl_mn$condition <- factor(output2$meta_data$condition[match(exprtbl_mn$sample_id,output2$meta_data$sample_id)])
exprtbl_mn$timepoint <- factor(output2$meta_data$timepoint[match(exprtbl_mn$sample_id,output2$meta_data$sample_id)],levels = timelevels2)
exprtbl_mn$patient_id <- factor(output2$meta_data$patient_id[match(exprtbl_mn$sample_id,output2$meta_data$sample_id)])
exprtbl_mn$stim <- factor(output2$meta_data$stim[match(exprtbl_mn$sample_id,output2$meta_data$sample_id)],levels = stimlevels2)
exprtbl_mn$response <- factor(output2$meta_data$response[match(exprtbl_mn$sample_id,output2$meta_data$sample_id)],levels = responselevels)

exprtbl_mn<-exprtbl_mn[exprtbl_mn$patient_id %nin% c("P003","P004"),]
exprtbl_mn<-exprtbl_mn[exprtbl_mn$stim %nin% c("UN","WT"),]
exprtbl_mn<-exprtbl_mn[exprtbl_mn$stim!="FLC",]
exprtbl_mn$patient_id <- factor(exprtbl_mn$patient_id)
exprtbl_mn<-exprtbl_mn[,c(markerlist_plot, "cluster","timepoint","patient_id")] %>%
  group_by(patient_id, timepoint, cluster) %>%
  summarize_all(funs(median))

ggdf2_mn<-melt(exprtbl_mn, id.var=c("cluster","timepoint","patient_id"))
ggdf2_mn$cluster <- factor(ggdf2_mn$cluster, levels=mergedlevels3d)
ggdf2_mn$patient_id <- factor(ggdf2_mn$patient_id, levels=unique(ggdf2_mn$patient_id))
ggdf2_mn$response <- factor(output2$meta_data$response[match(ggdf2_mn$patient_id,output2$meta_data$patient_id)],levels = responselevels)

ggdf2_mn_sub <- ggdf2_mn

#cd8 cells

ggdf2_mn_sub_cd8 <- ggdf2_mn_sub[ggdf2_mn_sub$cluster %in% cd8clusters,]

#choose timepoint here

ggdf2_mn_sub_cd8 <- ggdf2_mn_sub_cd8[ggdf2_mn_sub_cd8$timepoint == "C1",]
ggdf2_mn_sub_cd8 <- ggdf2_mn_sub_cd8[ggdf2_mn_sub_cd8$timepoint == "C3",]

pdf("plot_functional_box_CD8_C3_med_simplified.pdf",width=4.5,height=6)
for (i in 1:length(markerlist_plot)){
  ggdfsubset<-ggdf2_mn_sub_cd8[ggdf2_mn_sub_cd8$variable==markerlist_plot[i],]
  ggp <- ggplot(ggdfsubset, aes(x=response, y=value))+
    ggtitle(markerlist_plot[i])+
    facet_wrap(~cluster, scales='free')+
    geom_boxplot(outlier.color = NA, aes(fill=response), alpha=0.75)+
    geom_jitter(width=0, size=1)+
    #scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8,9,10))+
    scale_fill_manual(values = c(muted('blue'),muted('red'))) +
    theme(axis.text.x = element_text(angle = 45, vjust=.5, hjust=.5, color="black", size=8),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(linewidth=0.25),
          axis.line = element_line(linewidth=0.25),
          axis.text = element_text(color="black"),
          axis.title.y = element_text(size=8, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=8),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(1,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white"))
  print(ggp)
}
dev.off()

#cd4 cells

ggdf2_mn_sub_cd4 <- ggdf2_mn_sub[ggdf2_mn_sub$cluster %in% cd4clusters,]

#choose timepoint here

ggdf2_mn_sub_cd4 <- ggdf2_mn_sub_cd4[ggdf2_mn_sub_cd4$timepoint == "C1",]
ggdf2_mn_sub_cd4 <- ggdf2_mn_sub_cd4[ggdf2_mn_sub_cd4$timepoint == "C3",]

pdf("plot_functional_box_CD4_C3_med_simplified.pdf",width=4.5,height=6)
for (i in 1:length(markerlist_plot)){
  ggdfsubset<-ggdf2_mn_sub_cd4[ggdf2_mn_sub_cd4$variable==markerlist_plot[i],]
  ggp <- ggplot(ggdfsubset, aes(x=response, y=value))+
    ggtitle(markerlist_plot[i])+
    facet_wrap(~cluster, scales='free')+
    geom_boxplot(outlier.color = NA, aes(fill=response), alpha=0.75)+
    geom_jitter(width=0, size=1)+
    #scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8,9,10))+
    scale_fill_manual(values = c(muted('blue'),muted('red'))) +
    theme(axis.text.x = element_text(angle = 45, vjust=.5, hjust=.5, color="black", size=8),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(linewidth=0.25),
          axis.line = element_line(linewidth=0.25),
          axis.text = element_text(color="black"),
          axis.title.y = element_text(size=8, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=8),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(1,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white"))
  print(ggp)
}
dev.off()

#export the p values for functional comparisons


ggdf2_mn_sub_cd8 <- ggdf2_mn_sub[ggdf2_mn_sub$cluster %in% cd8clusters,]
ggdf2_mn_sub_cd8 <- ggdf2_mn_sub_cd8[ggdf2_mn_sub_cd8$timepoint == "C1",]
ggdf2_mn_sub_cd8_c1_p <- ggdf2_mn_sub_cd8 %>% 
  dplyr::select(cluster, patient_id, variable, response, value) %>% 
  group_by(response, cluster, variable) %>% 
  summarise(value=list(value), .groups = 'drop') %>%
  spread(response, value) %>%
  group_by(cluster, variable) %>%
  mutate(p_value = wilcox.test(unlist(short),unlist(long))$p.value,
         diffSL = mean(unlist(short))-mean(unlist(long))) #short - long
ggdf2_mn_sub_cd8_c1_p$padj <- p.adjust(ggdf2_mn_sub_cd8_c1_p$p_value, method="BH")

ggdf2_mn_sub_cd8 <- ggdf2_mn_sub[ggdf2_mn_sub$cluster %in% cd8clusters,]
ggdf2_mn_sub_cd8 <- ggdf2_mn_sub_cd8[ggdf2_mn_sub_cd8$timepoint == "C3",]
ggdf2_mn_sub_cd8_c3_p <- ggdf2_mn_sub_cd8 %>% 
  dplyr::select(cluster, patient_id, variable, response, value) %>% 
  group_by(response, cluster, variable) %>% 
  summarise(value=list(value), .groups = 'drop') %>%
  spread(response, value) %>%
  group_by(cluster, variable) %>%
  mutate(p_value = wilcox.test(unlist(short),unlist(long))$p.value,
         diffSL = mean(unlist(short))-mean(unlist(long))) #short - long
ggdf2_mn_sub_cd8_c3_p$padj <- p.adjust(ggdf2_mn_sub_cd8_c3_p$p_value, method="BH")

write.csv(ggdf2_mn_sub_cd8_c1_p[,c('cluster','variable','p_value','padj','diffSL')],'CD8_C1_functional_SvL.csv')
write.csv(ggdf2_mn_sub_cd8_c3_p[,c('cluster','variable','p_value','padj','diffSL')],'CD8_C3_functional_SvL.csv')


ggdf2_mn_sub_cd4 <- ggdf2_mn_sub[ggdf2_mn_sub$cluster %in% cd4clusters,]
ggdf2_mn_sub_cd4 <- ggdf2_mn_sub_cd4[ggdf2_mn_sub_cd4$timepoint == "C1",]
ggdf2_mn_sub_cd4_c1_p <- ggdf2_mn_sub_cd4 %>% 
  dplyr::select(cluster, patient_id, variable, response, value) %>% 
  group_by(response, cluster, variable) %>% 
  summarise(value=list(value), .groups = 'drop') %>%
  spread(response, value) %>%
  group_by(cluster, variable) %>%
  mutate(p_value = wilcox.test(unlist(short),unlist(long))$p.value,
         diffSL = mean(unlist(short))-mean(unlist(long))) #short - long
ggdf2_mn_sub_cd4_c1_p$padj <- p.adjust(ggdf2_mn_sub_cd4_c1_p$p_value, method="BH")

ggdf2_mn_sub_cd4 <- ggdf2_mn_sub[ggdf2_mn_sub$cluster %in% cd4clusters,]
ggdf2_mn_sub_cd4 <- ggdf2_mn_sub_cd4[ggdf2_mn_sub_cd4$timepoint == "C3",]
ggdf2_mn_sub_cd4_c3_p <- ggdf2_mn_sub_cd4 %>% 
  dplyr::select(cluster, patient_id, variable, response, value) %>% 
  group_by(response, cluster, variable) %>% 
  summarise(value=list(value), .groups = 'drop') %>%
  spread(response, value) %>%
  group_by(cluster, variable) %>%
  mutate(p_value = wilcox.test(unlist(short),unlist(long))$p.value,
         diffSL = mean(unlist(short))-mean(unlist(long))) #short - long
ggdf2_mn_sub_cd4_c3_p$padj <- p.adjust(ggdf2_mn_sub_cd4_c3_p$p_value, method="BH")

write.csv(ggdf2_mn_sub_cd4_c1_p[,c('cluster','variable','p_value','padj','diffSL')],'CD4_C1_functional_SvL.csv')
write.csv(ggdf2_mn_sub_cd4_c3_p[,c('cluster','variable','p_value','padj','diffSL')],'CD4_C3_functional_SvL.csv')


