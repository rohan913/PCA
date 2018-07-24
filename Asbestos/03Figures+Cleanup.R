#Overall#####
load('~/Box/Asbestos_SC/03_Analysis_V2/00_Asbestos/lung.clean.Robj')
lung=SubsetData(lung,ident.remove = 44)#rm  AM/AT2 doublets
current.cluster.ids <- c(2:11,13,15,16:18,20:24,26,28,32:37,42,43,45)
new.cluster.ids <- c('Bcell','C Mon','AM','C Mon','Tcell','AT2','NKcell','Endo(Cldn5, Kdr, Kitl, Igfbp7)',
                     'NC Mon','Neu','Endo(Cd93, Ptprb, Aqp1)','Tcell','IM','DC2','AT1','DC1','Tcell',
                     'Tcell','Ciliated cells','Fibroblasts','Club','pDC','Cycling AMs','Meso','AT2'
                      ,'Smooth Musc','AT2 cells','Ccr7+ DC','Mast cells/basophils','Lprogen','Megakar')
lung@ident <- plyr::mapvalues(x = lung@ident, from = current.cluster.ids, to = new.cluster.ids)
# CM=SubsetData(lung,ident.use = 'C Mon')
# NCM=SubsetData(lung,ident.use = 'NC Mon')
# save(CM,file='~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/CM.Robj')
# save(NCM,file='~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/NCM.Robj')

#recreate overall tsne for supplement 1
SC33=SubsetData(lung,cells.use = rownames(lung@meta.data[grep('SC33',lung@meta.data$orig.ident),]))
SC14=SubsetData(lung,cells.use = rownames(lung@meta.data[grep('SC14',lung@meta.data$orig.ident),]))
SC15=SubsetData(lung,cells.use = rownames(lung@meta.data[grep('SC15',lung@meta.data$orig.ident),]))
pdf('~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Figure 1/SupplementalFig1_individualTSNEs_Vln.pdf')
TSNEPlot(SC33,do.label = T,no.legend=T,do.return=T) + ggtitle(paste0('Naive (',nrow(SC33@meta.data),') cells'))
VlnPlot(object = SC33, features.plot = c("nGene"),x.lab.rot=T)
VlnPlot(object = SC33, features.plot = c("nUMI"),x.lab.rot=T)
VlnPlot(object = SC33, features.plot = c("percent.mito"),x.lab.rot=T)
TSNEPlot(SC14,do.label = T,no.legend=T,do.return=T) + ggtitle(paste0('TiO2 (',nrow(SC14@meta.data),') cells'))
VlnPlot(object = SC14, features.plot = c("nGene"),x.lab.rot=T)
VlnPlot(object = SC14, features.plot = c("nUMI"),x.lab.rot=T)
VlnPlot(object = SC14, features.plot = c("percent.mito"),x.lab.rot=T)
TSNEPlot(SC15,do.label = T,no.legend=T,do.return=T) + ggtitle(paste0('Asbestos (',nrow(SC15@meta.data),') cells'))
VlnPlot(object = SC15, features.plot = c("nGene"),x.lab.rot=T)
VlnPlot(object = SC15, features.plot = c("nUMI"),x.lab.rot=T)
VlnPlot(object = SC15, features.plot = c("percent.mito"),x.lab.rot=T)
dev.off()
colors=c("#619CFF","#00BA38","#F8766D")

pdf('~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Figure 1/Fig1b_mergedTSNE.pdf')
TSNEPlot(lung,do.label = T,no.legend=T,do.return=T) + ggtitle(paste0('Merged (',nrow(lung@meta.data),') cells'))
TSNEPlot(lung,do.label = F,no.legend=F,do.return=T,group.by='orig.ident') +scale_color_manual(values = colors) + ggtitle(paste0('Merged (',nrow(lung@meta.data),') cells'))
dev.off()

#regenerate fig 1c barplot overall
barplotdata = lung@meta.data 
barplotdata$ident = as.character(lung@ident)
toplotdata <- barplotdata %>% group_by(ident,orig.ident) %>% count()
sums <- toplotdata %>% group_by(ident) %>% summarise(Sum=sum(n))
sums <- rep(sums$Sum,each=3)
#megakar has only SC33andSC14
sums=sums[-c(55)]
toplotdata$sums <- sums
toplotdata$proportion <- (toplotdata$n/toplotdata$sums)
colors=c("#619CFF","#00BA38","#F8766D")
pdf('~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Figure 1/Fig1c_barplotforcombined.pdf')
ggplot(data=toplotdata,mapping = aes(x=ident,y=proportion,fill=orig.ident)) + geom_col(inherit.aes = T) +
  theme(axis.text = element_text(angle=90,hjust = 1))+scale_fill_manual(values=colors)
dev.off()

pdf("~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/SupplementalFig2_MergedVln.pdf")
VlnPlot(object = lung, features.plot = c("nGene"),x.lab.rot=T)
VlnPlot(object = lung, features.plot = c("nUMI"),x.lab.rot=T)
VlnPlot(object = lung, features.plot = c("percent.mito"),x.lab.rot=T)
dev.off()

#regenfig 1d heatmap overall
All.markers2 <- FindAllMarkers(lung)
write.csv(All.markers2,file="~/Box/Asbestos_SC/03_Analysis_V2/02_Markers/All.markers.2_cellnames.csv")
All.markers2 <- read.csv("~/Box/Asbestos_SC/03_Analysis_V2/02_Markers/All.markers.2_cellnames.csv",row.names = 1)
top5 <- All.markers2 %>% group_by(cluster) %>% top_n(3, avg_logFC)
smallALl <- SubsetData(lung, max.cells.per.ident=100)
pdf("~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Fig1d_heatmapoftop3genes.pdf")
DoHeatmap(
  object = lung,
  genes.use = top5$gene, group.cex = 8,
  slim.col.label = TRUE,cells.use = rownames(smallALl@meta.data),
  remove.key = TRUE,group.label.rot = T, cex.row = 3
)
dev.off()


#Figure 2######
#Barplots
Barplot <- function(x){
  load(paste0('~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/',x,'.Robj'))
  barplotdata=get(x)
  barplotdata=barplotdata@meta.data
  barplotdata$ident = as.character(get(x)@ident)
  toplotdata <- barplotdata %>% group_by(ident,orig.ident) %>% count()
  sums <- toplotdata %>% group_by(ident) %>% summarise(Sum=sum(n))
  sums <- rep(sums$Sum,each=3)
  toplotdata$sums <- sums
  toplotdata$proportion <- (toplotdata$n/toplotdata$sums)
  colors=c("#619CFF","#00BA38","#F8766D")
  ggplot(data=toplotdata,mapping = aes(x=ident,y=proportion,fill=orig.ident)) + geom_col(inherit.aes = T) +
    theme(axis.text = element_text(angle=90,hjust = 1))+scale_fill_manual(values=colors)
}
pdf('~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Figure 2/Fig2g_AT2barplot.pdf')
Barplot('AT2')
dev.off()

pdf('~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Figure 2/Fig2h_IMbarplot.pdf')
Barplot('IM')
dev.off()

pdf('~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Figure 2/Fig2i_CMbarplot.pdf')
Barplot('CM')
dev.off()

pdf('~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Figure 2/Fig2k_Fibbarplot.pdf')
Barplot('Fibroblasts')
dev.off()

#manually edit for CM
barplotdata = CM@meta.data 
barplotdata$ident = as.character(CM@ident)
toplotdata <- barplotdata %>% group_by(ident,orig.ident) %>% count()
sums <- toplotdata %>% group_by(ident) %>% summarise(Sum=sum(n))
sums <- rep(sums$Sum,each=3)
#last group has only SC33andSC14
sums=sums[-c(9)]
toplotdata$sums <- sums
toplotdata$proportion <- (toplotdata$n/toplotdata$sums)
colors=c("#619CFF","#00BA38","#F8766D")

pdf('~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Figure 2/Fig2i1barplot.pdf')
ggplot(data=toplotdata,mapping = aes(x=ident,y=proportion,fill=orig.ident)) + geom_col(inherit.aes = T) +
  theme(axis.text = element_text(angle=90,hjust = 1))+scale_fill_manual(values=colors)
dev.off()

#endo1 is Endo(Cldn5, Kdr, Kitl, Igfbp7)
barplotdata = Endothelial1@meta.data 
barplotdata$ident = as.character(Endothelial1@ident)
toplotdata <- barplotdata %>% group_by(ident,orig.ident) %>% count()
sums <- toplotdata %>% group_by(ident) %>% summarise(Sum=sum(n))
sums <- rep(sums$Sum,each=3)
#megakar has only SC33andSC14
sums=sums[-c(1)]
toplotdata$sums <- sums
toplotdata$proportion <- (toplotdata$n/toplotdata$sums)
colors=c("#619CFF","#00BA38","#F8766D")

pdf('~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Figure 2/Fig2k_Endo1barplot.pdf')
ggplot(data=toplotdata,mapping = aes(x=ident,y=proportion,fill=orig.ident)) + geom_col(inherit.aes = T) +
  theme(axis.text = element_text(angle=90,hjust = 1))+scale_fill_manual(values=colors)
dev.off()


#figure 2 TSNEs/Heatmap
Fig2TSNE <- function(x){
  y <- get(x)
  TSNEPlot(y,do.label = T, no.legend=T,do.return=T)+ggtitle(paste0(x,'(',nrow(y@meta.data),') cells'))
}

cells=c('AT2', 'IM', 'CM', 'Endothelial1', 'Fibroblasts')
pdf('~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Figure 2/Fig2_abcef_TSNEs(AT2.IM.CM.Endo.Fib).pdf')
for (i in cells){print(Fig2TSNE(i))}
dev.off()

Fig2Heatmap <- function(x){
  require(dplyr)
  y<-get(x);a<-FindAllMarkers(y)
  top10 <- a %>% group_by(cluster) %>% top_n(10, avg_logFC)
  DoHeatmap(
    object = y,
    genes.use = top10$gene, 
    slim.col.label = TRUE, 
    remove.key = TRUE,group.label.rot = T,title = paste0(x,'Heatmaptop10')
  )
}
pdf('~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Figure 2/Fig2_mnoqr_Heatmaps(AT2.IM.CM.Endo.Fib).pdf')
for (i in cells){print(Fig2Heatmap(i))}
dev.off()


#Figure 3####
pdf("~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Figure 3/Fig3a.pdf")
TSNEPlot(AM,do.label = T, no.legend=T,do.return=T)+ggtitle(paste0('AM (',nrow(AM@meta.data),') cells'))
TSNEPlot(AM,group.by='orig.ident',do.label=F,colors.use = c("#619CFF","#00BA38","#F8766D"))
dev.off()


barplotdata = AM@meta.data 
barplotdata$ident = as.character(AM@ident)
toplotdata <- barplotdata %>% group_by(ident,orig.ident) %>% count()


barplot=ggplot(data=toplotdata,mapping = aes(x=ident,y=n,fill=orig.ident)) + geom_col(inherit.aes = T) 
colors=c("#619CFF","#00BA38","#F8766D")
sums <- toplotdata %>% group_by(ident) %>% summarise(Sum=sum(n))
sums <- rep(sums$Sum,each=3)
toplotdata$sums <- sums
toplotdata$proportion <- (toplotdata$n/toplotdata$sums)

pdf("~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Fig3b_AMbarplot.pdf")
ggplot(data=toplotdata,mapping = aes(x=ident,y=proportion,fill=orig.ident)) + geom_col(inherit.aes = T) +
  theme(axis.text = element_text(angle=90,hjust = 1))+scale_fill_manual(values=colors)
dev.off()
top10 <- AM.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
pdf("~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Fig3c_AMHeatmap.pdf")
DoHeatmap(
  object = AM,
  genes.use = top10$gene, 
  slim.col.label = TRUE, 
  remove.key = TRUE,group.label.rot = T
)
dev.off()

pdf("~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Figure 3/Fig3e.pdf")
VlnPlot(AM, c("Mertk","Mrc1",'Fcgr1','Adgre1','Axl'))
DoHeatmap(AM,genes.use = c("Mertk","Mrc1",'Fcgr1','Adgre1','Axl'),use.scaled = T,slim.col.label = T,remove.key = T)
DotPlot(AM,c("Mertk","Mrc1",'Fcgr1','Adgre1','Axl'),plot.legend = T)
#FeatureHeatmap(AM2,features.plot = c("Mertk","Mrc1",'Fcgr1','Adgre1','Axl'))
dev.off()

pdf("~/Box/Asbestos_SC/04_Manuscript/01_Figures_Tables/Figure 3/Fig3f.pdf")
VlnPlot(AM, c("Siglecf","Car4",'Chil3','Marco','Lpl','Il18','Csf1r','Csf1','Csf2','Itgam','Itgax','Cx3cr1'))
DoHeatmap(AM,genes.use = c("Siglecf","Car4",'Chil3','Marco','Lpl','Il18','Csf1r','Csf1','Csf2','Itgam','Itgax','Cx3cr1'),use.scaled = T,slim.col.label = T,remove.key = T)
DotPlot(AM,c("Siglecf","Car4",'Chil3','Marco','Lpl','Il18','Csf1r','Csf1','Csf2','Itgam','Itgax','Cx3cr1'))
FeatureHeatmap(AM, c('Csf2'))
dev.off()


#reshape/cleanup AM TSNE
#a <- c('Mmp12', 'Retnla','Chil3','Pdgfa','Mertk','Siglecf')
a <- c('Mertk','Siglecf','Mmp12','Retnla','Chil3','Pdgfa','Pdgfb')
pdf('~/Desktop/AM_Vln_Featureplot.pdf')
VlnPlot(AM,a)
# for (i in a){
#   print(VlnPlot(AM, features.plot = c(i)))
# }
# for (i in a){
#   FeaturePlot(AM, features.plot = c(i))
# }
FeaturePlot(AM,a)
dev.off()

#testing
load('~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/AM.Robj')
# AM <- RunTSNE(object = AM, dims.use = 1:3, do.fast = TRUE, check_duplicates = FALSE, perplexity = 90,seed.use = 400)
# TSNEPlot(object = AM, do.label = T)
# TSNEPlot(object = AM, do.label = F, group.by = "orig.ident")

load('~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/AT1.Robj')
AT1@meta.data$cluster <- AT1@ident
pdf('~/Desktop/AT1_TSNE_Barplot.pdf')
TSNEPlot(object = AT1, do.label = T,no.legend=T)
ggplot(AT1@meta.data,aes(x=cluster,fill=orig.ident))+geom_bar(position='fill')+
  ggtitle(paste('AT1','Cluster Proportions'))+labs(x='cluster',y='proportion',fill='library.id')+
  scale_fill_manual(values=c("#00BA38","#F8766D","#619CFF"))
dev.off()