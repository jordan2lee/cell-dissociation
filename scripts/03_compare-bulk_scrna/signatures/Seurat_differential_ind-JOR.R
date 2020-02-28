
args = commandArgs(trailingOnly=TRUE)
#subject_id <- as.character(args[1])
subject_id<- 'TEST'

library(ggplot2)
library(dplyr)
library(Seurat) # install.packages('devtools') --> devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(openxlsx) # install.packages("openxlsx", dependencies=TRUE)
library(reticulate)
library(cowplot) # install.packages("cowplot")
library(RColorBrewer)

################ Manual
indir ="/Users/leejor/Ellrott_Lab/04_CSBC/cell-dissociation/raw-data/TESTDATA--filtered_gene_bc_matrices/hg19"
outdir='/Users/leejor/Ellrott_Lab/04_CSBC/cell-dissociation/data/03_compare-bulk_scrna/signatures'
#########################
print(R.version.string)

# Declaring path where raw data is and creating new directory as output
dir_input_sc <- indir
dir.create(file.path(paste(outdir,"MechEnz/scRNAseq", sep='/'),subject_id), recursive=TRUE)
dir_output <- file.path(paste(outdir,"MechEnz/scRNAseq", sep='/'),subject_id)
setwd(dir_output)



############################################# SKIP FOR NOW UNTIL HAVE SEVERAL SAMPLE FILES
##info <- read.xlsx(file.path(dir_input,"Data.xlsx"), sheet = 1)
##sample_dirs <- list.files(dir_input_sc,pattern = paste("cellranger_"))
##sample_ids <- sample_dirs[grep(subject_id,sample_dirs)]
############################################# SKIP FOR NOW UNTIL HAVE SEVERAL SAMPLE FILES


#######################
# Create Seurat obj
#######################
# Loading 10X Genomics data
ctrl.data <- Read10X(data.dir = indir)
#ctrl.data <- Read10X(data.dir = file.path(dir_input_sc,sample_ids[2],"outs/filtered_feature_bc_matrix"))

# Loading gene expression counts into Seurat object for both mechanical and enzymatic (separately)
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "DIGESTION_MEC", min.cells=3, min.features=200)
ctrl$stim <- "MECHANICAL"


#######################
# QC and select cells
#######################
# Calculating % mitochondrial genes in each cell ("percent.mt" ~ indicator of dead or broken cells) for both mechanical and enzymatic
ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^MT-") # QC metrics are stored in ctrl@meta.data. 

# Showing metrics for the first 5 cells for both mechanical and enzymatic
head(ctrl@meta.data,5)

# Checking the number of cells for both mechanical and enzymatic
length(Idents(ctrl))

# Visualizing QC metrics using violin plots for both mechanical and enzymatic
VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(dir_output,"violin_mech.png"),width = 12, height =10)

# Visualizing some feature-feature relationships for both mechanical and enzymatic
plot1 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
ggsave(file.path(dir_output,"FeatureScatter_mech.png"),width = 12, height =10)

###################### Manual
ctrl_min_ft = 200
ctrl_max_ft = 2500
ctrl_max_mt = 5
######################
# Removing cells based on violin and feature plots
ctrl <- subset(ctrl, subset= nFeature_RNA > ctrl_min_ft & nFeature_RNA < ctrl_max_ft & percent.mt < ctrl_max_mt)
length(Idents(ctrl))


#######################
# Normalize
#######################
# Normalizing the feature expression measurements for each cell by the total expression multiplied by a scale factor (10,000 by default), and log-transform the result. 
ctrl <- NormalizeData(ctrl, normalization.method="LogNormalize", scale.factor=10000) # Normalized values are stored in pbmc[["RNA"]]@data.


#######################
# Optional: ft selection
#######################
# Calculating a subset of features that exhibit high cell-to-cell variation in the dataset 
ctrl <- FindVariableFeatures(ctrl, selection.method="vst", nfeatures=2000)

# Showing top 10 variable genes and plotting variable features with and without labels
top10 <- head(VariableFeatures(ctrl), 10)
top10
plot1 <- VariableFeaturePlot(ctrl)
plot1
LabelPoints(plot=plot1, points=top10, repel=TRUE)
ggsave(file.path(dir_output,"VariableFeaturePlot_10_mech.png"),width = 12, height =10)


#######################
# Scale - impacts heatmap interpretability
#######################
# Not used downstream: Linear transformation (scaling): mean and variace expression across cells is 0 and 1 respectively
ctrl_scaled <- ScaleData(ctrl, features=VariableFeatures(ctrl))
write.table(ctrl_scaled@assays$RNA@scale.data,"scaled.data_mech.txt",sep="\t",row.names=TRUE,col.names=NA)

#During scaling, you can additionally regress out any effects arising from technical variation. eg. nUMI for variation in sequencing depth. Any variable added in the metadata field can be used for this purpose
#ctrl <- ScaleData(object = ctrl, vars.to.regress = c("nUMI""))

# ############################################
# # Integrating both mechanical and enzymatic
# ############################################
# immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
# immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
# DefaultAssay(immune.combined) <- "integrated"
# 
# # Running the standard workflow for visualization and clustering
# immune.combined <- ScaleData(immune.combined, verbose = FALSE)
# immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE) #Linear dimension reduction
# 
# # Performing t-SNE and Clustering
# immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
# immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
# immune.combined <- FindClusters(immune.combined, resolution = 0.5)
# 
# # Exporting the proportion of cells for each cluster
# prop <- prop.table(table(Idents(immune.combined)))
# write.table(prop,file.path(dir_output,"prop.table.txt"),sep="\t",row.names=FALSE)
# prop_mech <- table(immune.combined@meta.data$seurat_clusters[which(immune.combined@meta.data$stim == "MECHANICAL")])
# write.table(prop_mech,file.path(dir_output,"mech.table.txt"),sep="\t",row.names=FALSE)
# prop_enz <- table(immune.combined@meta.data$seurat_clusters[which(immune.combined@meta.data$stim == "ENZYMATIC")])
# write.table(prop_enz,file.path(dir_output,"enz.table.txt"),sep="\t",row.names=FALSE)
# 
# # Visualizing both conditions
# p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
# p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
# plot_grid(p1, p2)
# ggsave(file.path(dir_output,"plot_grid.png"),width = 12, height =10)
# 
# DimPlot(immune.combined, reduction = "umap", split.by = "stim")
# ggsave(file.path(dir_output,"umap_split.png"),width = 12, height =10)
# 
# 
# # Identifying markers for every cluster compared to all remaining cells, report only the positive ones
# immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# immune.combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% print(n = Inf)
# 
# Markers_clusters <- data.frame(immune.combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% print(n = Inf))
# write.table(Markers_clusters,file.path(dir_output,"Markers_clusters.txt"),sep="\t",row.names=FALSE,col.names=TRUE)
# Markers_clusters_top2 <- data.frame(immune.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) %>% print(n = Inf))
# write.table(Markers_clusters_top2,file.path(dir_output,"Markers_clusters_top2.txt"),sep="\t",row.names=FALSE,col.names=TRUE)
# write.table(immune.combined.markers,file.path(dir_output,"Markers_clusters_all.txt"),sep="\t",row.names=FALSE,col.names=TRUE)
# 
# # Visualizing top 10 markers
# top10 <- immune.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# DoHeatmap(immune.combined, features = top10$gene, size=2.5, hjust = 0.5, angle = 90) + NoLegend() +
# 	theme(axis.text.y = element_text(size = 7))
# ggsave(file.path(dir_output,"DoHeatmap_topmarkers.png"),width = 12, height =12)
# 
# # Saving workspace
# saveRDS(immune.combined, file = file.path(dir_output,paste0(subject_id,".rds")))
