
args = commandArgs(trailingOnly=TRUE)
subject_id <- as.character(args[1])


library(ggplot2)
library(dplyr)
library(Seurat)
library(openxlsx)
library(reticulate)
library(cowplot)
library(RColorBrewer)
#use_python("/share/software/user/open/python/3.6.1/bin/python3")
#source("/oak/stanford/groups/andrewg/users/aespin/Project1/Scripts/Run_CIBERSORT_dir.R")
#source("/oak/stanford/groups/andrewg/users/aespin/Project1/Scripts/Run_CIBERSORT_signature2_dir.R")
#source("/oak/stanford/groups/andrewg/users/aespin/Scripts/CSx/CIBERSORTxFractions.R")

################ Manual
#indir ="/oak/stanford/groups/andrewg/users/zinaida/u54_data/cellranger_aligned"
#outdir='/oak/stanford/groups/andrewg/users/aespin'

indir ="../../raw-data/TESTDATA--filtered_gene_bc_matrices/hg19"
outdir='../../data/02_process-scrna'
#########################

# Declaring path where raw data is and creating new directory as output
dir_input_sc <- indir
dir.create(file.path(paste(outdir,"MechEnz/scRNAseq", sep='/'),subject_id))
dir_output <- file.path("/oak/stanford/groups/andrewg/users/aespin/MechEnz/scRNAseq",subject_id)
setwd(dir_output)
dir_input_sc <- "../../raw-data/TESTDATA--filtered_gene_bc_matrices/hg19"


#info <- read.xlsx(file.path(dir_input,"Data.xlsx"), sheet = 1)
sample_dirs <- list.files(dir_input_sc,pattern = paste("cellranger_"))
sample_ids <- sample_dirs[grep(subject_id,sample_dirs)]


# Loading 10X Genomics data
ctrl.data <- Read10X(data.dir = file.path(dir_input_sc,sample_ids[2],"outs/filtered_feature_bc_matrix"))
stim.data <- Read10X(data.dir = file.path(dir_input_sc,sample_ids[1],"outs/filtered_feature_bc_matrix"))

# Loading gene expression counts into Seurat object for both mechanical and enzymatic (separately)
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "DIGESTION_MEC", min.cells=3, min.features=200)
ctrl$stim <- "MECHANICAL"

stim <- CreateSeuratObject(counts = stim.data, project = "DIGESTION_ENZ", min.cells=3, min.features=200)
stim$stim <- "ENZYMATIC"

# Calculating % mitochondrial genes in each cell ("percent.mt" ~ indicator of dead or broken cells) for both mechanical and enzymatic
ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^MT-") # QC metrics are stored in ctrl@meta.data. 
stim[["percent.mt"]] <- PercentageFeatureSet(stim, pattern = "^MT-")

# Showing metrics for the first 5 cells for both mechanical and enzymatic
head(ctrl@meta.data,5)
head(stim@meta.data,5)

# Checking the number of cells for both mechanical and enzymatic
length(Idents(ctrl))
length(Idents(stim))

# Visualizing QC metrics using violin plots for both mechanical and enzymatic
VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(dir_output,"violin_mech.png"),width = 12, height =10)
VlnPlot(stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(dir_output,"violin_enz.png"),width = 12, height =10)

# Visualizing some feature-feature relationships for both mechanical and enzymatic
plot1 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
ggsave(file.path(dir_output,"FeatureScatter_mech.png"),width = 12, height =10)
plot1 <- FeatureScatter(stim, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(stim, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
ggsave(file.path(dir_output,"FeatureScatter_enz.png"),width = 12, height =10)

# Removing cells based on violin and feature plots
ctrl <- subset(ctrl, subset= nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
length(Idents(ctrl))
stim <- subset(stim, subset= nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
length(Idents(stim))

# Normalizing the feature expression measurements for each cell by the total expression multiplied by a scale factor (10,000 by default), and log-transform the result. 
ctrl <- NormalizeData(ctrl, normalization.method="LogNormalize", scale.factor=10000) # Normalized values are stored in pbmc[["RNA"]]@data.
stim <- NormalizeData(stim, normalization.method="LogNormalize", scale.factor=10000) # Normalized values are stored in pbmc[["RNA"]]@data.

# Calculating a subset of features that exhibit high cell-to-cell variation in the dataset 
ctrl <- FindVariableFeatures(ctrl, selection.method="vst", nfeatures=2000)
stim <- FindVariableFeatures(stim, selection.method="vst", nfeatures=2000)

# Showing top 10 variable genes and plotting variable features with and without labels
top10 <- head(VariableFeatures(ctrl), 10)
top10
plot1 <- VariableFeaturePlot(ctrl)
plot1
LabelPoints(plot=plot1, points=top10, repel=TRUE)
ggsave(file.path(dir_output,"VariableFeaturePlot_10_mech.png"),width = 12, height =10)
top10 <- head(VariableFeatures(stim), 10)
top10
plot1 <- VariableFeaturePlot(stim)
plot1
LabelPoints(plot=plot1, points=top10, repel=TRUE)
ggsave(file.path(dir_output,"VariableFeaturePlot_10_enz.png"),width = 12, height =10)

# Not used downstream: Linear transformation (scaling): mean and variace expression across cells is 0 and 1 respectively
ctrl_scaled <- ScaleData(ctrl, features=VariableFeatures(ctrl))
stim_scaled <- ScaleData(stim, features=VariableFeatures(ctrl))
write.table(ctrl_scaled@assays$RNA@scale.data,"scaled.data_mech.txt",sep="\t",row.names=TRUE,col.names=NA)
write.table(stim_scaled@assays$RNA@scale.data,"scaled.data_enz.txt",sep="\t",row.names=TRUE,col.names=NA)

#During scaling, you can additionally regress out any effects arising from technical variation. eg. nUMI for variation in sequencing depth. Any variable added in the metadata field can be used for this purpose
#ctrl <- ScaleData(object = ctrl, vars.to.regress = c("nUMI""))

# Integrating both mechanical and enzymatic
immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
DefaultAssay(immune.combined) <- "integrated"

# Running the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE) #Linear dimension reduction

# Performing t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Exporting the proportion of cells for each cluster
prop <- prop.table(table(Idents(immune.combined)))
write.table(prop,file.path(dir_output,"prop.table.txt"),sep="\t",row.names=FALSE)
prop_mech <- table(immune.combined@meta.data$seurat_clusters[which(immune.combined@meta.data$stim == "MECHANICAL")])
write.table(prop_mech,file.path(dir_output,"mech.table.txt"),sep="\t",row.names=FALSE)
prop_enz <- table(immune.combined@meta.data$seurat_clusters[which(immune.combined@meta.data$stim == "ENZYMATIC")])
write.table(prop_enz,file.path(dir_output,"enz.table.txt"),sep="\t",row.names=FALSE)

# Visualizing both conditions
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
ggsave(file.path(dir_output,"plot_grid.png"),width = 12, height =10)

DimPlot(immune.combined, reduction = "umap", split.by = "stim")
ggsave(file.path(dir_output,"umap_split.png"),width = 12, height =10)


# Identifying markers for every cluster compared to all remaining cells, report only the positive ones
immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
immune.combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% print(n = Inf)

Markers_clusters <- data.frame(immune.combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% print(n = Inf))
write.table(Markers_clusters,file.path(dir_output,"Markers_clusters.txt"),sep="\t",row.names=FALSE,col.names=TRUE)
Markers_clusters_top2 <- data.frame(immune.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) %>% print(n = Inf))
write.table(Markers_clusters_top2,file.path(dir_output,"Markers_clusters_top2.txt"),sep="\t",row.names=FALSE,col.names=TRUE)
write.table(immune.combined.markers,file.path(dir_output,"Markers_clusters_all.txt"),sep="\t",row.names=FALSE,col.names=TRUE)

# Visualizing top 10 markers
top10 <- immune.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(immune.combined, features = top10$gene, size=2.5, hjust = 0.5, angle = 90) + NoLegend() +
	theme(axis.text.y = element_text(size = 7))
ggsave(file.path(dir_output,"DoHeatmap_topmarkers.png"),width = 12, height =12)

# Saving workspace
saveRDS(immune.combined, file = file.path(dir_output,paste0(subject_id,".rds")))

# Running CIBERSORT
mean_raw.data_all <- c()
for (i in 0:(length(names(table(Idents(immune.combined))))-1)){
raw.data <- as.matrix(GetAssayData(immune.combined)[, WhichCells(immune.combined, ident = i)])
dir.create(file.path(dir_output,i))
#try(Run_CIBERSORT(raw.data,file.path(dir_output,i)))
#try(Run_CIBERSORT_signature2(raw.data,file.path(dir_output,i)))
mean_raw.data <- rowMeans(raw.data,na.rm = TRUE)
mean_raw.data_all <- cbind(mean_raw.data_all,mean_raw.data)
}
colnames(mean_raw.data_all) <- names(table(Idents(immune.combined)))
dir.create(file.path(dir_output,"Cib_mean"))
try(Run_CIBERSORT(mean_raw.data_all,file.path(dir_output,"Cib_mean")))
try(Run_CIBERSORT_signature2(mean_raw.data_all,file.path(dir_output,"Cib_mean")))


# Plotting conserved cell type markers
DefaultAssay(immune.combined) <- "RNA"

#markers <- FindConservedMarkers(immune.combined, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
#FeaturePlot(immune.combined, features = c(), min.cutoff = "q9")
#ggsave(file.path(dir_output,"markers.png"),width = 12, height =10)
#immune.combined <- RenameIdents(immune.combined, `0` = "Monocytes",...)
#DotPlot(immune.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + RotatedAxis()
#ggsave(file.path(dir_output,"Markers_conditions.png"),width = 12, height =10)


# Identifying differential expressed genes across conditions
for (i in 0:(length(names(table(Idents(immune.combined))))-1)){
cluster0.cells <- subset(immune.combined, idents = i)
Idents(cluster0.cells) <- "stim"
avg.cluster0.cells <- log1p(AverageExpression(cluster0.cells, verbose = FALSE)$RNA)
avg.cluster0.cells$gene <- rownames(avg.cluster0.cells)
	if(length(levels(Idents(cluster0.cells))) > 1){
	avg.cluster0.cells$Diff <- avg.cluster0.cells[,1] - avg.cluster0.cells[,2]
	avg.cluster0.cells <- avg.cluster0.cells[order(-abs(avg.cluster0.cells$Diff)),]
	write.table(avg.cluster0.cells,file.path(dir_output,paste0("avg.cluster",i,".cells.txt")),sep="\t",row.names=FALSE,col.names=TRUE)
	}else{
	write.table(avg.cluster0.cells,file.path(dir_output,paste0("avg.cluster",i,".cells.txt")),sep="\t",row.names=FALSE,col.names=TRUE)
	}
}

#genes.to.label = rownames(avg.cluster0.cells)[1:10]
#p1 <- ggplot(avg.cluster0.cells, aes(MECHANICAL, ENZYMATIC)) + geom_point() + ggtitle("Cluster 0")
#p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
#p2 <- ggplot(avg.cluster9.cells, aes(MECHANICAL, ENZYMATIC)) + geom_point() + ggtitle("Cluster 9")
#p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
#plot_grid(p1, p2)
#ggsave(file.path(dir_output,"Scatter_top10_cluster0.png"),width = 12, height =10)


#FeaturePlot(immune.combined, features = rownames(avg.cluster0.cells_2)[1:3], split.by = "stim", max.cutoff = 3, 
#    cols = c("grey", "red"))
#ggsave(file.path(dir_output,"FeaturePlot_top3_cluster0_up.png"),width = 12, height =10)

#plots <- VlnPlot(immune.combined, features = rownames(avg.cluster0.cells_2)[1:3], split.by = "stim", group.by = "celltype", 
#    pt.size = 0, combine = FALSE)
#CombinePlots(plots = plots, ncol = 1)
#ggsave(file.path(dir_output,"VlnPlot_top3_cluster0_up.png"),width = 12, height =10)




