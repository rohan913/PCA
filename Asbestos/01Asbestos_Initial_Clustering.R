#Top ----
#Analysis of the single cell RNA-seq data from asbestos-induced pulmonary fibrosis in mice.
#SC33_multirun - Naive mouse, SC14 - TiO2, SC15 - Asbestos
#Authors: Nikita Joshi, Rohan Verma, Scott Budinger, Alexander Misharin

#load libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
#SC33 control, SC14 Titanium Dioxide, SC15 Asbestos
#Initial data processing, PCA and clustering ----
#Plan: Create Merged Seurat Object without CCA to Identify Clusters 
#Load data and create aggregated Seurat object::
SC33.data <- Read10X(data.dir = "01_data/SC33_multirun/filtered_gene_bc_matrices/mm10/")
SC33 <- CreateSeuratObject(raw.data = SC33.data, min.cells = 3, min.genes = 200, project = "SC33")
SC33@cell.names <- paste("SC33", SC33@cell.names, sep = "_")

#Modify cell names so they are all unique between samples 
colnames(x = SC33@raw.data) <- paste("SC33", colnames(x = SC33@raw.data), sep = "_")
rownames(x = SC33@meta.data) <- paste("SC33", rownames(x = SC33@meta.data), sep = "_")
head(SC33@cell.names)
tail(SC33@cell.names)
SC14.data <- Read10X(data.dir = "01_data/SC14/filtered_gene_bc_matrices/mm10/")
lung <- AddSamples(object = SC33, new.data = SC14.data, add.cell.id = "SC14")
SC15.data <- Read10X(data.dir = "01_data/SC15/filtered_gene_bc_matrices/mm10/")
lung <- AddSamples(object = lung, new.data = SC15.data, add.cell.id = "SC15")
head(lung@cell.names)
tail(lung@cell.names)

lung #27998 genes across 24618 samples.

#Normalize and find variable genes
#Confirm that there are 958 variable genes
lung <- NormalizeData(object = lung, normalization.method = "LogNormalize", scale.factor = 10000) #Log normalize
lung <- FindVariableGenes(object = lung, mean.function = ExpMean, dispersion.function = LogVMR) #Find variable genes
length(x = lung@var.genes) #958
hv.genes <- head(rownames(lung@hvg.info), 958)

#Data inspection and basic filtering
mito.genes <- grep(pattern = "^mt-", x = rownames(x = lung@data), value = TRUE)
percent.mito <- Matrix::colSums(lung@raw.data[mito.genes, ])/Matrix::colSums(lung@raw.data)
lung <- AddMetaData(object = lung, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = lung, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = lung, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = lung, gene1 = "nUMI", gene2 = "nGene")

#In the beginning, skip filtering and keep all cells, irrespecive of %mito and number of genes
#lung <- FilterCells(object = lung, subset.names = c("nGene", "percent.mito"), low.thresholds = c(300, -Inf), high.thresholds = c(4000, 0.1))

#Scale and regress
lung <- ScaleData(object = lung, display.progress = T,vars.to.regress = c("nUMI"))

#Perform PCA
lung <- RunPCA(object = lung, pc.genes = lung@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 100)
PCAPlot(object = lung, dim.1 = 1, dim.2 = 2)
PCAPlot(object = lung, dim.1 = 1, dim.2 = 3)
#lung <- ProjectPCA(object = lung, do.print = FALSE)
PCHeatmap(object = lung, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = lung, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = lung, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = lung, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = lung, pc.use = 91:100, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = lung, num.pc = 100)
#lung <- FindClusters(object = lung, reduction.type = "pca", dims.use = 1:44, resolution = 1, print.output = F, save.SNN = TRUE, force.recalc = T)
lung <- FindClusters(object = lung, reduction.type = "pca", dims.use = 1:75, 
                     resolution = 3, print.output = F, save.SNN = TRUE, force.recalc = T, n.start = 100)
PrintFindClustersParams(object = lung)
lung <- RunTSNE(object = lung, dims.use = 1:75, do.fast = TRUE, check_duplicates = FALSE, perplexity = 30)
TSNEPlot(object = lung, do.label = T)
TSNEPlot(object = lung, do.label = T, group.by = "orig.ident")
save(lung, file = "lung.Robj")
load("lung.Robj")

TSNEPlot(object = lung, do.label = T)
lung <-ValidateSpecificClusters(lung, cluster1 = 0, cluster2 = 2, top.genes = 50) #B cells
TSNEPlot(object = lung, do.label = T)
lung <-ValidateSpecificClusters(lung, cluster1 = 1, cluster2 = 7, top.genes = 50) #AT2
TSNEPlot(object = lung, do.label = T)
lung <-ValidateSpecificClusters(lung, cluster1 = 19, cluster2 = 7, top.genes = 50) #AT2
TSNEPlot(object = lung, do.label = T)
lung <-ValidateSpecificClusters(lung, cluster1 = 14, cluster2 = 7, top.genes = 50) #AT2
TSNEPlot(object = lung, do.label = T)
lung <-ValidateSpecificClusters(lung, cluster1 = 12, cluster2 = 4, top.genes = 50) #AM
TSNEPlot(object = lung, do.label = T) 

lung.markers <- FindAllMarkers(object = lung, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.5, max.cells.per.ident = 200)
write.csv(lung.markers, "lung.markers_PC1_75_res3.csv")

FeaturePlot(lung, features.plot = c("Sftpc", "Ager","Foxj1", "Scgb3a1"))

#Initial list of clusters based on examination of markers----
#C2 - B cells
#C3 - Classcial monocytes
#C4 - Alveolar macrophages
#C5 - Classcial monocytes
#C6 - T cells
#C7 - AT2 cells
#C8 - NK cells
#C9 - Endothelial cells (Cldn5, Kdr, Kitl, Igfbp7)
#C10 - Non-classical monocytes
#C11 - Neutrophils
#C13 - Endothelial cells (Cd93, Ptprb, Aqp1)
#C15 - T cells
#C16 - Peribronchial interstitial macrophages
#C17 - DC2
#C18 - AT1 cells
#C20 - DC1
#C21 - T cells
#C22 - T cells
#C23 - Ciliated cells
#C24 - Fibroblasts
#C25 - Neutrophils
#C26 - Club cells
#C27 - Ciliated cells
#C28 - pDC
#C29 - Doublets: Neutrophils and B cells
#C30 - Endothelial cells (C13)
#C31 - Endothelial cells (C13)
#C32 - Cycling AMs
#C33 - Mesothelial cells
#C34 - AT2 cells
#C35 - Smooth muscle cells
#C36 - AT2 cells
#C37 - Ccr7+ DC
#C38 - Doublets: AT2 and Neutrophils
#C39 - Doublets: B cells and monocytes
#C40 - Cell cycle/Damaged DC1
#C41 - Cell cycle/Damaged T cells
#C42 - Mast cells/basophils
#C43 - Lymphatic progenitors
#C44 - Doublets: AM/AT2
#C45 - Megakaryocytes
#C46 - Doublets: AM and Neutrophils
#C47 - Doublets: AM and NK cells

#cleanup/remove doublets 
lung <- SubsetData(lung, ident.remove = c(29, 38, 39, 40, 41, 44, 46, 47))
TSNEPlot(object = lung, do.label = T) 
#Validating Clusters
lung <-ValidateSpecificClusters(lung, cluster1 = 30, cluster2 = 13, top.genes = 50) #Endothelial
TSNEPlot(object = lung, do.label = T) 
lung <-ValidateSpecificClusters(lung, cluster1 = 31, cluster2 = 13, top.genes = 50) #Endothelial
TSNEPlot(object = lung, do.label = T) 
lung <-ValidateSpecificClusters(lung, cluster1 = 27, cluster2 = 23, top.genes = 20) #Ciliated
TSNEPlot(object = lung, do.label = T) 
lung <-ValidateSpecificClusters(lung, cluster1 = 25, cluster2 = 11, top.genes = 50) #Neutrophils
TSNEPlot(object = lung, do.label = T) 
lung <-ValidateSpecificClusters(lung, cluster1 = 32, cluster2 = 4, top.genes = 50) #AM
TSNEPlot(object = lung, do.label = T) 
lung <-ValidateSpecificClusters(lung, cluster1 = 36, cluster2 = 7, top.genes = 50) #AT2
TSNEPlot(object = lung, do.label = T) 
lung <-ValidateSpecificClusters(lung, cluster1 = 34, cluster2 = 7, top.genes = 50) #AT2
TSNEPlot(object = lung, do.label = T) 

save(lung, file = "lung.clean.Robj")

#To vizualize top X genes
topgenes <- lung.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)
cluster.averages <- AverageExpression(object = lung, return.seurat = TRUE, show.progress = T)
DoHeatmap(object = cluster.averages, genes.use = topgenes$gene, group.label.rot = TRUE, group.cex = 0)

#Barplot to show that all clusters are composed of cells coming from all 3 libraries

load('~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/lung.clean.Robj')
lung@metadata.cluster <- lung@ident
ggplot(lung@metadata,aes(x=cluster,fill=orig.ident))+
  geom_bar(position='fill')+
  ggtitle('Initial Cluster Proportions') +
  labs(x='cluster',y='proportion',fill='library.id')+
  scale_fill_manual(values=c("#00BA38","#F8766D","#619CFF"))