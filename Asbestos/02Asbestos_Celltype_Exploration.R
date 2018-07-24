##### Explore specific cell types #####

#Explore AT2 cells ----
AT2 <- SubsetData(lung, ident.use = c(7, 34, 36))
TSNEPlot(AT2, do.label = T)
TSNEPlot(AT2, do.label = T, group.by = "orig.ident")
AT2 <- FindVariableGenes(object = AT2, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = AT2@var.genes)
AT2 <- ScaleData(object = AT2, vars.to.regress = c("nUMI", "percent.mito"))
AT2 <- RunPCA(object = AT2, pc.genes = AT2@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 50)
PCAPlot(object = AT2, dim.1 = 1, dim.2 = 2)
PCAPlot(object = AT2, dim.1 = 1, dim.2 = 3)
PCAPlot(object = AT2, dim.1 = 1, dim.2 = 3, group.by = "orig.ident")
AT2 <- ProjectPCA(object = AT2, do.print = FALSE)
PCHeatmap(object = AT2, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT2, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT2, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT2, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = AT2, num.pc = 50)
AT2 <- FindClusters(object = AT2, reduction.type = "pca", dims.use = 1:15, 
                    resolution = 0.3, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = AT2)
AT2 <- RunTSNE(object = AT2, dims.use = 1:12, do.fast = TRUE, check_duplicates = FALSE)
TSNEPlot(object = AT2, do.label = T)
TSNEPlot(object = AT2, do.label = T, group.by = "orig.ident")


FeaturePlot(AT2, features.plot = c("Ccr2", "Cd3e","Nkg7", "Cd79a"))
AT2.markers <- FindAllMarkers(object = AT2, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.5, max.cells.per.ident = 200)
#Clusters 3, 4 and 5 contain B, monocytes and T/NK cell doublets, let's remove them. 
AT2 <- SubsetData(AT2, ident.use = c(0,1,2))
TSNEPlot(AT2, do.label = T)
VlnPlot(AT2, features.plot = "nGene")

AT2 <- FindVariableGenes(object = AT2, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = AT2@var.genes)
AT2 <- ScaleData(object = AT2, vars.to.regress = c("nUMI"))
AT2 <- RunPCA(object = AT2, pc.genes = AT2@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 20)
PCAPlot(object = AT2, dim.1 = 1, dim.2 = 2)
PCAPlot(object = AT2, dim.1 = 1, dim.2 = 3)
AT2 <- ProjectPCA(object = AT2, do.print = FALSE)
PCHeatmap(object = AT2, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT2, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT2, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT2, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = AT2, num.pc = 20)
AT2 <- FindClusters(object = AT2, reduction.type = "pca", dims.use = 1:8, 
                    resolution = 0.4, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = AT2)
AT2 <- RunTSNE(object = AT2, dims.use = 1:8, do.fast = TRUE, check_duplicates = FALSE, perplexity = 20)
TSNEPlot(object = AT2, do.label = T)
TSNEPlot(object = AT2, do.label = F, group.by = "orig.ident")

save(AT2, file = "AT2.Robj")

AT2.markers <- FindAllMarkers(object = AT2, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.25, test.use = "wilcox")
write.csv(AT2.markers, "AT2.markers_PC1_8_res0.4wilcox.csv")
top10 <- AT2.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = AT2, genes.use = top10$gene, slim.col.label = TRUE, remove.key = F)

VlnPlot(AT2, features.plot = c("Retnla", "Chil3"), group.by = "orig.ident")
#We might do better if we simply run DEGs between the samples. 

#Vizualize origin of cells across the clusters
gSC33=grep('SC33',AT2@meta.data$orig.ident)
gSC14=grep('SC14',AT2@meta.data$orig.ident)
gSC15=grep('SC15',AT2@meta.data$orig.ident)
a=TSNEPlot(AT2,cells.use = rownames((AT2@meta.data)[gSC33,]), do.label=T,no.legend=T,do.return=T)
b=TSNEPlot(AT2,cells.use = rownames((AT2@meta.data)[gSC14,]), do.label=T,no.legend=T,do.return=T)
c=TSNEPlot(AT2,cells.use = rownames((AT2@meta.data)[gSC15,]), do.label=T,no.legend=T,do.return=T)
plot_grid(a,b,c, ncol = 1)

barplotdata = AT2@meta.data 
barplotdata$ident = as.character(AT2@ident)
toplotdata <- barplotdata %>% group_by(ident,orig.ident) %>% count()
ggplot(data=toplotdata,mapping = aes(x=ident,y=n,fill=orig.ident)) + geom_col(inherit.aes = T) 
sums <- toplotdata %>% group_by(ident) %>% summarise(Sum=sum(n))
sums <- rep(sums$Sum,each=3)
toplotdata$sums <- sums
toplotdata$proportion <- (toplotdata$n/toplotdata$sums)
ggplot(data=toplotdata,mapping = aes(x=ident,y=proportion,fill=orig.ident)) + geom_col(inherit.aes = T) +
  theme(axis.text = element_text(angle=90,hjust = 1))


#Explore AT1 cells ----
AT1 <- SubsetData(lung, ident.use = c(18))
TSNEPlot(AT1, do.label = T)
TSNEPlot(AT1, do.label = T, group.by = "orig.ident")
AT1 <- FindVariableGenes(object = AT1, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = AT1@var.genes)
AT1 <- ScaleData(object = AT1, vars.to.regress = c("nUMI", "percent.mito"))
AT1 <- RunPCA(object = AT1, pc.genes = AT1@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 50)
PCAPlot(object = AT1, dim.1 = 1, dim.2 = 2)
PCAPlot(object = AT1, dim.1 = 1, dim.2 = 3)
PCAPlot(object = AT1, dim.1 = 1, dim.2 = 3, group.by = "orig.ident")
AT1 <- ProjectPCA(object = AT1, do.print = FALSE)
PCHeatmap(object = AT1, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT1, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT1, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT1, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = AT1, num.pc = 50)
AT1 <- FindClusters(object = AT1, reduction.type = "pca", dims.use = 1:9, 
                    resolution = 0.3, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = AT1)
AT1 <- RunTSNE(object = AT1, dims.use = 1:9, do.fast = TRUE, check_duplicates = FALSE)
TSNEPlot(object = AT1, do.label = T)
TSNEPlot(object = AT1, do.label = T, group.by = "orig.ident")


FeaturePlot(AT1, features.plot = c("Ccr2", "Cd3e","Nkg7", "Cd79a"))
AT1.markers <- FindAllMarkers(object = AT1, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.5, max.cells.per.ident = 200)
#Clusters 1 and 2 contain B, monocytes and T/NK cell doublets, let's remove them. 
AT1 <- SubsetData(AT1, ident.use = c(0))
TSNEPlot(AT1, do.label = T)

AT1 # 27998 genes across 397 samples.

AT1 <- FindVariableGenes(object = AT1, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = AT1@var.genes)
AT1 <- ScaleData(object = AT1, vars.to.regress = c("nUMI"))
AT1 <- RunPCA(object = AT1, pc.genes = AT1@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 20)
PCAPlot(object = AT1, dim.1 = 1, dim.2 = 2)
PCAPlot(object = AT1, dim.1 = 1, dim.2 = 3)
AT1 <- ProjectPCA(object = AT1, do.print = FALSE)
PCHeatmap(object = AT1, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT1, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT1, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT1, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = AT1, num.pc = 20)
AT1 <- FindClusters(object = AT1, reduction.type = "pca", dims.use = 1:8, 
                    resolution = 0.3, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = AT1)
AT1 <- RunTSNE(object = AT1, dims.use = 1:8, do.fast = TRUE, check_duplicates = FALSE, perplexity = 20)
TSNEPlot(object = AT1, do.label = T)
TSNEPlot(object = AT1, do.label = F, group.by = "orig.ident")
#No clear clustering by the group/exposure.

save(AT1, file = "AT1.Robj")

#Explore Fibroblasts ----
Fibroblasts <- SubsetData(lung, ident.use = c(24))
TSNEPlot(Fibroblasts, do.label = T)
TSNEPlot(Fibroblasts, do.label = T, group.by = "orig.ident")
Fibroblasts <- FindVariableGenes(object = Fibroblasts, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = Fibroblasts@var.genes)
Fibroblasts <- ScaleData(object = Fibroblasts, vars.to.regress = c("nUMI", "percent.mito"))
Fibroblasts <- RunPCA(object = Fibroblasts, pc.genes = Fibroblasts@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 50)
PCAPlot(object = Fibroblasts, dim.1 = 1, dim.2 = 2)
PCAPlot(object = Fibroblasts, dim.1 = 1, dim.2 = 3)
PCAPlot(object = Fibroblasts, dim.1 = 1, dim.2 = 3, group.by = "orig.ident")
Fibroblasts <- ProjectPCA(object = Fibroblasts, do.print = FALSE)
PCHeatmap(object = Fibroblasts, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Fibroblasts, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Fibroblasts, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Fibroblasts, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = Fibroblasts, num.pc = 50)
Fibroblasts <- FindClusters(object = Fibroblasts, reduction.type = "pca", dims.use = 1:6, 
                    resolution = 0.3, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = Fibroblasts)
Fibroblasts <- RunTSNE(object = Fibroblasts, dims.use = 1:6, do.fast = TRUE, check_duplicates = FALSE)
TSNEPlot(object = Fibroblasts, do.label = T)
TSNEPlot(object = Fibroblasts, do.label = T, group.by = "orig.ident")


FeaturePlot(Fibroblasts, features.plot = c("Ccr2", "Cx3cr1","Nkg7", "Cd79a"))
Fibroblasts.markers <- FindAllMarkers(object = Fibroblasts, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.5, max.cells.per.ident = 200)
#Cluster 2 contains B cells and monocytes, let's remove them. 
Fibroblasts <- SubsetData(Fibroblasts, ident.use = c(0,1))
TSNEPlot(Fibroblasts, do.label = T)
VlnPlot(Fibroblasts, features.plot = "nGene")

Fibroblasts <- FindVariableGenes(object = Fibroblasts, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = Fibroblasts@var.genes)
Fibroblasts <- ScaleData(object = Fibroblasts, vars.to.regress = c("nUMI"))
Fibroblasts <- RunPCA(object = Fibroblasts, pc.genes = Fibroblasts@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 20)
PCAPlot(object = Fibroblasts, dim.1 = 1, dim.2 = 2)
PCAPlot(object = Fibroblasts, dim.1 = 1, dim.2 = 3)
PCAPlot(object = Fibroblasts, dim.1 = 1, dim.2 = 3, group.by = "orig.ident")
Fibroblasts <- ProjectPCA(object = Fibroblasts, do.print = FALSE)
PCHeatmap(object = Fibroblasts, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Fibroblasts, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Fibroblasts, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Fibroblasts, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = Fibroblasts, num.pc = 20)
Fibroblasts <- FindClusters(object = Fibroblasts, reduction.type = "pca", dims.use = 1:6, 
                    resolution = 0.3, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = Fibroblasts)
Fibroblasts <- RunTSNE(object = Fibroblasts, dims.use = 1:6, do.fast = TRUE, check_duplicates = FALSE, perplexity = 30)
TSNEPlot(object = Fibroblasts, do.label = T)
TSNEPlot(object = Fibroblasts, do.label = F, group.by = "orig.ident")

save(Fibroblasts, file = "Fibroblasts.Robj")

Fibroblasts.markers <- FindAllMarkers(object = Fibroblasts, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.25, test.use = "wilcox")
write.csv(Fibroblasts.markers, "Fibroblasts.markers_PC1_6_res0.3wilcox.csv")
top10 <- Fibroblasts.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = Fibroblasts, genes.use = top10$gene, slim.col.label = TRUE, remove.key = F)

VlnPlot(Fibroblasts, features.plot = c("Igfbp4", "Col3a1"), group.by = "orig.ident")
VlnPlot(Fibroblasts, features.plot = c("Igfbp4", "Col3a1"))
#We might do better if we simply run DEGs between the samples. 

#Vizualize origin of cells across the clusters
gSC33=grep('SC33',Fibroblasts@meta.data$orig.ident)
gSC14=grep('SC14',Fibroblasts@meta.data$orig.ident)
gSC15=grep('SC15',Fibroblasts@meta.data$orig.ident)
a=TSNEPlot(Fibroblasts,cells.use = rownames((Fibroblasts@meta.data)[gSC33,]), do.label=T,no.legend=T,do.return=T)
b=TSNEPlot(Fibroblasts,cells.use = rownames((Fibroblasts@meta.data)[gSC14,]), do.label=T,no.legend=T,do.return=T)
c=TSNEPlot(Fibroblasts,cells.use = rownames((Fibroblasts@meta.data)[gSC15,]), do.label=T,no.legend=T,do.return=T)
plot_grid(a,b,c, ncol = 1)

barplotdata = Fibroblasts@meta.data 
barplotdata$ident = as.character(Fibroblasts@ident)
toplotdata <- barplotdata %>% group_by(ident,orig.ident) %>% count()
ggplot(data=toplotdata,mapping = aes(x=ident,y=n,fill=orig.ident)) + geom_col(inherit.aes = T) 
sums <- toplotdata %>% group_by(ident) %>% summarise(Sum=sum(n))
sums <- rep(sums$Sum,each=3)
toplotdata$sums <- sums
toplotdata$proportion <- (toplotdata$n/toplotdata$sums)
ggplot(data=toplotdata,mapping = aes(x=ident,y=proportion,fill=orig.ident)) + geom_col(inherit.aes = T) +
  theme(axis.text = element_text(angle=90,hjust = 1))




#Explore peribronchial IM ----
IM <- SubsetData(lung, ident.use = c(16))
TSNEPlot(IM, do.label = T)
TSNEPlot(IM, do.label = T, group.by = "orig.ident")
IM <- FindVariableGenes(object = IM, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = IM@var.genes)
IM <- ScaleData(object = IM, vars.to.regress = c("nUMI", "percent.mito"))
IM <- RunPCA(object = IM, pc.genes = IM@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 50)
PCAPlot(object = IM, dim.1 = 1, dim.2 = 2)
PCAPlot(object = IM, dim.1 = 1, dim.2 = 3)
PCAPlot(object = IM, dim.1 = 1, dim.2 = 3, group.by = "orig.ident")
IM <- ProjectPCA(object = IM, do.print = FALSE)
PCHeatmap(object = IM, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = IM, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = IM, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = IM, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = IM, num.pc = 50)
IM <- FindClusters(object = IM, reduction.type = "pca", dims.use = 1:9, 
                    resolution = 0.7, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = IM)
IM <- RunTSNE(object = IM, dims.use = 1:9, do.fast = TRUE, check_duplicates = FALSE)
TSNEPlot(object = IM, do.label = T)
TSNEPlot(object = IM, do.label = T, group.by = "orig.ident")


FeaturePlot(IM, features.plot = c("Pdgfb", "Cd3e","Nkg7", "Cd79a"))
IM.markers <- FindAllMarkers(object = IM, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.5, max.cells.per.ident = 200)
#Clusters 3, 4 and 5 contain B, T/NK cell doublets, let's remove them. 
IM <- SubsetData(IM, ident.use = c(0,1,2))
TSNEPlot(IM, do.label = T)
VlnPlot(IM, features.plot = "nGene")

IM <- FindVariableGenes(object = IM, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = IM@var.genes)
IM <- ScaleData(object = IM, vars.to.regress = c("nUMI"))
IM <- RunPCA(object = IM, pc.genes = IM@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 20)
PCAPlot(object = IM, dim.1 = 1, dim.2 = 2)
PCAPlot(object = IM, dim.1 = 1, dim.2 = 3)
IM <- ProjectPCA(object = IM, do.print = FALSE)
PCHeatmap(object = IM, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = IM, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = IM, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = IM, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = IM, num.pc = 20)
IM <- FindClusters(object = IM, reduction.type = "pca", dims.use = 1:4, 
                    resolution = 0.4, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = IM)
IM <- RunTSNE(object = IM, dims.use = 1:4, do.fast = TRUE, check_duplicates = FALSE, perplexity = 50)
TSNEPlot(object = IM, do.label = T)
TSNEPlot(object = IM, do.label = F, group.by = "orig.ident")
FeaturePlot(IM, features.plot = "Chil3")

VlnPlot(IM, features.plot = c("Fn1", "Chil3", "Retnla"))
VlnPlot(IM, features.plot = c("Retnla", "Fn1"),  group.by = "orig.ident")

save(IM, file = "IM.Robj")

IM.markers <- FindAllMarkers(object = IM, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.25, test.use = "wilcox")
write.csv(IM.markers, "IM.markers_PC1_4_res0.4wilcox.csv")
top10 <- IM.markers %>% group_by(cluster) %>% top_n(15, avg_logFC)
DoHeatmap(object = IM, genes.use = top10$gene, slim.col.label = TRUE, remove.key = F)

VlnPlot(IM, features.plot = c("Retnla", "Chil3"), group.by = "orig.ident")
#We might do better if we simply run DEGs between the samples. 

#Vizualize origin of cells across the clusters
gSC33=grep('SC33',IM@meta.data$orig.ident)
gSC14=grep('SC14',IM@meta.data$orig.ident)
gSC15=grep('SC15',IM@meta.data$orig.ident)
a=TSNEPlot(IM,cells.use = rownames((IM@meta.data)[gSC33,]), do.label=T,no.legend=T,do.return=T)
b=TSNEPlot(IM,cells.use = rownames((IM@meta.data)[gSC14,]), do.label=T,no.legend=T,do.return=T)
c=TSNEPlot(IM,cells.use = rownames((IM@meta.data)[gSC15,]), do.label=T,no.legend=T,do.return=T)
plot_grid(a,b,c, ncol = 1)

barplotdata = IM@meta.data 
barplotdata$ident = as.character(IM@ident)
toplotdata <- barplotdata %>% group_by(ident,orig.ident) %>% count()
ggplot(data=toplotdata,mapping = aes(x=ident,y=n,fill=orig.ident)) + geom_col(inherit.aes = T) 
sums <- toplotdata %>% group_by(ident) %>% summarise(Sum=sum(n))
sums <- rep(sums$Sum,each=3)
toplotdata$sums <- sums
toplotdata$proportion <- (toplotdata$n/toplotdata$sums)
ggplot(data=toplotdata,mapping = aes(x=ident,y=proportion,fill=orig.ident)) + geom_col(inherit.aes = T) +
  theme(axis.text = element_text(angle=90,hjust = 1))


#Explore Monocytes ----
Monocytes <- SubsetData(lung, ident.use = c(3,5,10))
TSNEPlot(Monocytes, do.label = T)
TSNEPlot(Monocytes, do.label = T, group.by = "orig.ident")
Monocytes <- FindVariableGenes(object = Monocytes, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = Monocytes@var.genes)
Monocytes <- ScaleData(object = Monocytes, vars.to.regress = c("nUMI", "percent.mito"))
Monocytes <- RunPCA(object = Monocytes, pc.genes = Monocytes@var.genes, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5, pcs.compute = 50)
PCAPlot(object = Monocytes, dim.1 = 1, dim.2 = 2)
PCAPlot(object = Monocytes, dim.1 = 1, dim.2 = 3)
PCAPlot(object = Monocytes, dim.1 = 1, dim.2 = 3, group.by = "orig.ident")
Monocytes <- ProjectPCA(object = Monocytes, do.print = FALSE)
PCHeatmap(object = Monocytes, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Monocytes, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Monocytes, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Monocytes, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = Monocytes, num.pc = 50)
Monocytes <- FindClusters(object = Monocytes, reduction.type = "pca", dims.use = 1:14, 
                   resolution = 0.7, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = Monocytes)
Monocytes <- RunTSNE(object = Monocytes, dims.use = 1:14, do.fast = TRUE, check_duplicates = FALSE)
TSNEPlot(object = Monocytes, do.label = T)
TSNEPlot(object = Monocytes, do.label = T, group.by = "orig.ident")


FeaturePlot(Monocytes, features.plot = c("Pdgfb", "Cd3e","Nkg7", "Cd79a"))
Monocytes.markers <- FindAllMarkers(object = Monocytes, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.5, max.cells.per.ident = 200)
#Clusters 3, 4 and 5 contain B, T/NK cell doublets, let's remove them. 
Monocytes <- SubsetData(Monocytes, ident.remove = c(7))
TSNEPlot(Monocytes, do.label = T)
VlnPlot(Monocytes, features.plot = "nGene")

Monocytes <- FindVariableGenes(object = Monocytes, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = Monocytes@var.genes)
Monocytes <- ScaleData(object = Monocytes, vars.to.regress = c("nUMI", "orig.ident"))
Monocytes <- RunPCA(object = Monocytes, pc.genes = Monocytes@var.genes, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5, pcs.compute = 20)
PCAPlot(object = Monocytes, dim.1 = 1, dim.2 = 2)
PCAPlot(object = Monocytes, dim.1 = 1, dim.2 = 3)
Monocytes <- ProjectPCA(object = Monocytes, do.print = FALSE)
PCHeatmap(object = Monocytes, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Monocytes, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Monocytes, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = Monocytes, num.pc = 20)
Monocytes <- FindClusters(object = Monocytes, reduction.type = "pca", dims.use = 1:4, 
                   resolution = 0.4, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = Monocytes)
Monocytes <- RunTSNE(object = Monocytes, dims.use = 1:4, do.fast = TRUE, check_duplicates = FALSE, perplexity = 50)
TSNEPlot(object = Monocytes, do.label = T)
TSNEPlot(object = Monocytes, do.label = F, group.by = "orig.ident")
FeaturePlot(Monocytes, features.plot = "Ccr2")

VlnPlot(Monocytes, features.plot = c("Fn1", "Chil3", "Retnla"))
VlnPlot(Monocytes, features.plot = c("Retnla", "Fn1"),  group.by = "orig.ident")

save(Monocytes, file = "Monocytes.Robj")

Monocytes.markers <- FindAllMarkers(object = Monocytes, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.25, test.use = "wilcox")
write.csv(Monocytes.markers, "Monocytes.markers_PC1_4_res0.4wilcox.csv")
top10 <- Monocytes.markers %>% group_by(cluster) %>% top_n(15, avg_logFC)
DoHeatmap(object = Monocytes, genes.use = top10$gene, slMonocytes.col.label = TRUE, remove.key = F)

VlnPlot(Monocytes, features.plot = c("Ccr2"), group.by = "orig.ident")
#We might do better if we sMonocytesply run DEGs between the samples. 

#Vizualize origin of cells across the clusters
gSC33=grep('SC33',Monocytes@meta.data$orig.ident)
gSC14=grep('SC14',Monocytes@meta.data$orig.ident)
gSC15=grep('SC15',Monocytes@meta.data$orig.ident)
a=TSNEPlot(Monocytes,cells.use = rownames((Monocytes@meta.data)[gSC33,]), do.label=T,no.legend=T,do.return=T)
b=TSNEPlot(Monocytes,cells.use = rownames((Monocytes@meta.data)[gSC14,]), do.label=T,no.legend=T,do.return=T)
c=TSNEPlot(Monocytes,cells.use = rownames((Monocytes@meta.data)[gSC15,]), do.label=T,no.legend=T,do.return=T)
plot_grid(a,b,c, ncol = 1)

barplotdata = Monocytes@meta.data 
barplotdata$ident = as.character(Monocytes@ident)
toplotdata <- barplotdata %>% group_by(ident,orig.ident) %>% count()
ggplot(data=toplotdata,mapping = aes(x=ident,y=n,fill=orig.ident)) + geom_col(inherit.aes = T) 
sums <- toplotdata %>% group_by(ident) %>% summarise(Sum=sum(n))
sums <- rep(sums$Sum,each=3)
toplotdata$sums <- sums
toplotdata$proportion <- (toplotdata$n/toplotdata$sums)
ggplot(data=toplotdata,mapping = aes(x=ident,y=proportion,fill=orig.ident)) + geom_col(inherit.aes = T) +
  theme(axis.text = element_text(angle=90,hjust = 1))




#Explore Endothelial1 cells ----
Endothelial1 <- SubsetData(lung, ident.use = c(9))
TSNEPlot(Endothelial1, do.label = T)
TSNEPlot(Endothelial1, do.label = T, group.by = "orig.ident")
Endothelial1 <- FindVariableGenes(object = Endothelial1, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = Endothelial1@var.genes)
Endothelial1 <- ScaleData(object = Endothelial1, vars.to.regress = c("nUMI", "percent.mito"))
Endothelial1 <- RunPCA(object = Endothelial1, pc.genes = Endothelial1@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 50)
PCAPlot(object = Endothelial1, dim.1 = 1, dim.2 = 2)
PCAPlot(object = Endothelial1, dim.1 = 1, dim.2 = 3)
PCAPlot(object = Endothelial1, dim.1 = 1, dim.2 = 3, group.by = "orig.ident")
Endothelial1 <- ProjectPCA(object = Endothelial1, do.print = FALSE)
PCHeatmap(object = Endothelial1, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Endothelial1, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Endothelial1, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Endothelial1, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = Endothelial1, num.pc = 50)
Endothelial1 <- FindClusters(object = Endothelial1, reduction.type = "pca", dims.use = 1:8, 
                    resolution = 0.3, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = Endothelial1)
Endothelial1 <- RunTSNE(object = Endothelial1, dims.use = 1:8, do.fast = TRUE, check_duplicates = FALSE)
TSNEPlot(object = Endothelial1, do.label = T)
TSNEPlot(object = Endothelial1, do.label = T, group.by = "orig.ident")


FeaturePlot(Endothelial1, features.plot = c("Ccr2", "Cd3e","Nkg7", "Cd79a"))
Endothelial1.markers <- FindAllMarkers(object = Endothelial1, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.5, max.cells.per.ident = 200)
#Clusters 3, 4 and 5 contain B, monocytes and T/NK cell doublets, let's remove them. 
Endothelial1 <- SubsetData(Endothelial1, ident.use = c(0,1))
TSNEPlot(Endothelial1, do.label = T)
VlnPlot(Endothelial1, features.plot = "nGene")

Endothelial1 <- FindVariableGenes(object = Endothelial1, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = Endothelial1@var.genes)
Endothelial1 <- ScaleData(object = Endothelial1, vars.to.regress = c("nUMI"))
Endothelial1 <- RunPCA(object = Endothelial1, pc.genes = Endothelial1@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 20)
PCAPlot(object = Endothelial1, dim.1 = 1, dim.2 = 2)
PCAPlot(object = Endothelial1, dim.1 = 1, dim.2 = 3)
Endothelial1 <- ProjectPCA(object = Endothelial1, do.print = FALSE)
PCHeatmap(object = Endothelial1, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Endothelial1, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Endothelial1, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Endothelial1, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = Endothelial1, num.pc = 20)
Endothelial1 <- FindClusters(object = Endothelial1, reduction.type = "pca", dims.use = 1:5, 
                    resolution = 0.4, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = Endothelial1)
Endothelial1 <- RunTSNE(object = Endothelial1, dims.use = 1:5, do.fast = TRUE, check_duplicates = FALSE, perplexity = 20)
TSNEPlot(object = Endothelial1, do.label = T)
TSNEPlot(object = Endothelial1, do.label = F, group.by = "orig.ident")

save(Endothelial1, file = "Endothelial1.Robj")

Endothelial1.markers <- FindAllMarkers(object = Endothelial1, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.25, test.use = "wilcox")
write.csv(Endothelial1.markers, "Endothelial1.markers_PC1_5_res0.4wilcox.csv")
top10 <- Endothelial1.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = Endothelial1, genes.use = top10$gene, slim.col.label = TRUE, remove.key = F)

VlnPlot(Endothelial1, features.plot = c("Retnla", "Chil3"), group.by = "orig.ident")
#We might do better if we simply run DEGs between the samples. 

#Vizualize origin of cells across the clusters
gSC33=grep('SC33',Endothelial1@meta.data$orig.ident)
gSC14=grep('SC14',Endothelial1@meta.data$orig.ident)
gSC15=grep('SC15',Endothelial1@meta.data$orig.ident)
a=TSNEPlot(Endothelial1,cells.use = rownames((Endothelial1@meta.data)[gSC33,]), do.label=T,no.legend=T,do.return=T)
b=TSNEPlot(Endothelial1,cells.use = rownames((Endothelial1@meta.data)[gSC14,]), do.label=T,no.legend=T,do.return=T)
c=TSNEPlot(Endothelial1,cells.use = rownames((Endothelial1@meta.data)[gSC15,]), do.label=T,no.legend=T,do.return=T)
plot_grid(a,b,c, ncol = 1)

barplotdata = Endothelial1@meta.data 
barplotdata$ident = as.character(Endothelial1@ident)
toplotdata <- barplotdata %>% group_by(ident,orig.ident) %>% count()
ggplot(data=toplotdata,mapping = aes(x=ident,y=n,fill=orig.ident)) + geom_col(inherit.aes = T) 
sums <- toplotdata %>% group_by(ident) %>% summarise(Sum=sum(n))
sums <- rep(sums$Sum,each=3)
toplotdata$sums <- sums
toplotdata$proportion <- (toplotdata$n/toplotdata$sums)
ggplot(data=toplotdata,mapping = aes(x=ident,y=proportion,fill=orig.ident)) + geom_col(inherit.aes = T) +
  theme(axis.text = element_text(angle=90,hjust = 1))


#Explore AMs
AM <- SubsetData(lung, ident.use = c(4))
TSNEPlot(AM, do.label = T)
TSNEPlot(AM, do.label = T, group.by = "orig.ident")
AM <- FindVariableGenes(object = AM, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = AM@var.genes)
AM <- ScaleData(object = AM, vars.to.regress = c("nUMI", "percent.mito"))
AM <- RunPCA(object = AM, pc.genes = AM@var.genes, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5, pcs.compute = 50)
PCAPlot(object = AM, dim.1 = 1, dim.2 = 2)
PCAPlot(object = AM, dim.1 = 1, dim.2 = 3)
PCAPlot(object = AM, dim.1 = 1, dim.2 = 3, group.by = "orig.ident")
AM <- ProjectPCA(object = AM, do.print = FALSE)
# PCHeatmap(object = AM, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
#           label.columns = FALSE, use.full = FALSE)
# PCHeatmap(object = AM, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
#           label.columns = FALSE, use.full = FALSE)
# PCHeatmap(object = AM, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
#           label.columns = FALSE, use.full = FALSE)
# PCHeatmap(object = AM, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
#           label.columns = FALSE, use.full = FALSE)
#PCElbowPlot(object = AM, num.pc = 50)
AM <- FindClusters(object = AM, reduction.type = "pca", dims.use = 1:7, 
                   resolution = 0.7, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = AM)
AM <- RunTSNE(object = AM, dims.use = 1:7, do.fast = TRUE, check_duplicates = FALSE)
TSNEPlot(object = AM, do.label = T)
TSNEPlot(object = AM, do.label = T, group.by = "orig.ident")


FeaturePlot(AM, features.plot = c("Pdgfb", "Cd3e","Nkg7", "Cd79a"))
AM.markers <- FindAllMarkers(object = AM, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.5, max.cells.per.ident = 200)
#Cluster 5 contain B cell doublets
AM <- SubsetData(AM, ident.use = c(0:4,6))
TSNEPlot(AM, do.label = T)
VlnPlot(AM, features.plot = "nGene")

AM <- FindVariableGenes(object = AM, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = AM@var.genes)
AM <- ScaleData(object = AM, vars.to.regress = c("nUMI"))
AM <- RunPCA(object = AM, pc.genes = AM@var.genes, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5, pcs.compute = 20)
PCAPlot(object = AM, dim.1 = 1, dim.2 = 2)
PCAPlot(object = AM, dim.1 = 1, dim.2 = 3)
AM <- ProjectPCA(object = AM, do.print = FALSE)
PCHeatmap(object = AM, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
           label.columns = FALSE, use.full = FALSE)
# PCHeatmap(object = AM, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
#           label.columns = FALSE, use.full = FALSE)
# PCHeatmap(object = AM, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
#           label.columns = FALSE, use.full = FALSE)
# PCHeatmap(object = AM, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
#           label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = AM, num.pc = 20)
AM <- FindClusters(object = AM, reduction.type = "pca", dims.use = 1:3, 
                   resolution = 0.15, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = AM)
AM <- RunTSNE(object = AM, dims.use = 1:3, do.fast = TRUE, check_duplicates = FALSE, perplexity = 50)
TSNEPlot(object = AM, do.label = T)
TSNEPlot(object = AM, do.label = F, group.by = "orig.ident")
FeaturePlot(AM, features.plot = "Chil3")

VlnPlot(AM, features.plot = c("Fn1", "Chil3", "Retnla"))
VlnPlot(AM, features.plot = c("Retnla", "Fn1"),  group.by = "orig.ident")

save(AM, file = "~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/AM.Robj")

AM.markers <- FindAllMarkers(object = AM, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.25, test.use = "wilcox")
write.csv(AM.markers, "~/Box/Asbestos_SC/03_Analysis_V2/02_Markers/AM.markers_PC1_4_res0.4wilcox.csv")
top10 <- AM.markers %>% group_by(cluster) %>% top_n(15, avg_logFC)
DoHeatmap(object = AM, genes.use = top10$gene, slim.col.label = TRUE, remove.key = F)

VlnPlot(AM, features.plot = c("Retnla", "Chil3"), group.by = "orig.ident")

#Classical Monocytes
CM <- SubsetData(lung,ident.use = c(3,5,10))
TSNEPlot(CM, do.label = T)
TSNEPlot(CM, do.label = T, group.by = "orig.ident")
CM <- FindVariableGenes(object = CM, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = CM@var.genes)
CM <- ScaleData(object = CM, vars.to.regress = c("orig.ident"))
CM <- RunPCA(object = CM, pc.genes = CM@var.genes, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5, pcs.compute = 50)
PCAPlot(object = CM, dim.1 = 1, dim.2 = 2)
PCAPlot(object = CM, dim.1 = 1, dim.2 = 3)
PCAPlot(object = CM, dim.1 = 1, dim.2 = 3, group.by = "orig.ident")
CM <- ProjectPCA(object = CM, do.print = FALSE)
PCHeatmap(object = CM, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
  label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = CM, num.pc = 50)
CM <- FindClusters(object = CM, reduction.type = "pca", dims.use = 1:3, 
                   resolution = 0.1, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParCMs(object = CM)
CM <- RunTSNE(object = CM, dims.use = 1:3, do.fast = TRUE, check_duplicates = FALSE)
TSNEPlot(object = CM, do.label = T)
TSNEPlot(object = CM, do.label = T, group.by = "orig.ident")
FeaturePlot(CM,c("Ccr2","Sell"))
CM <- SubsetData(CM,ident.use = 0)
CM <- ScaleData(object = CM, vars.to.regress = c("nUMI","percent.mito"))
CM <- RunPCA(object = CM, pc.genes = CM@var.genes, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5, pcs.compute = 50)

CM <- ProjectPCA(object = CM, do.print = FALSE)
PCHeatmap(object = CM, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = CM, num.pc = 50)

CM <- FindClusters(object = CM, reduction.type = "pca", dims.use = 1:3, 
                   resolution = 0.2, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = CM)
CM <- RunTSNE(object = CM, dims.use = 1:3, do.fast = TRUE, check_duplicates = FALSE, perplexity = 50)

TSNEPlot(object = CM, do.label = T)
TSNEPlot(object = CM, do.label = F, group.by = "orig.ident")
FeaturePlot(CM, features.plot = "Chil3")

VlnPlot(CM, features.plot = c("Fn1", "Chil3", "Retnla"))
VlnPlot(CM, features.plot = c("Retnla", "Fn1"),  group.by = "orig.ident")

save(CM, file = "~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/CM.Robj")

CM.markers <- FindAllMarkers(object = CM, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.25, test.use = "wilcox")
write.csv(CM.markers, "~/Box/Asbestos_SC/03_Analysis_V2/02_Markers/CM.markers_PC1_4_res0.4wilcox.csv")