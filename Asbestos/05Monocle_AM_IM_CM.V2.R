#Analysis in Monocle 


#Load in objects and narrow down to populations of interest 
load("~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/IM.Robj")
load("~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/CM.Robj")
load("~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/AM.Robj")

#AM <- SubsetData(AM, cells.use = rownames(AM@meta.data)[which(AM@meta.data$orig.ident=='SC15')])
IM <- SubsetData(IM, ident.use=2)

#add initial ID to metadata
AM@meta.data$initial.id=paste0(AM@meta.data$orig.ident,"AM")#~2 clusters
IM@meta.data$initial.id=paste0(IM@meta.data$orig.ident,"IM")#1 cluster
CM@meta.data$initial.id=paste0(CM@meta.data$orig.ident,"CM")#3clusters
#merge objects normalize and scale before using as cell dataset in monocle
AM_IM_CM <- MergeSeurat(AM,IM,do.normalize = F)
AM_IM_CM <- MergeSeurat(AM_IM_CM,CM,do.normalize = F)
AM_IM_CM <- NormalizeData(AM_IM_CM)
AM_IM_CM <- ScaleData(AM_IM_CM)
rm(AM);rm(IM);rm(CM)
# monocle
library(monocle)
AM_IM_CM <- importCDS(AM_IM_CM,import_all = F)

##estimate size factors
AM_IM_CM <- estimateSizeFactors(AM_IM_CM)
AM_IM_CM <- estimateDispersions(AM_IM_CM)


#######Quality control (http://pklab.med.harvard.edu/scw2014/monocle_tutorial_oedsymbol.html)##
#plotting density of expressed genes in each library
library(reshape2)
L <- log10(exprs(AM_IM_CM)+1)
L[L==0] <- NA
melted.dens.df <- melt(t(scale(t(L))))
qplot(value, geom = 'density', data = melted.dens.df) + stat_function(fun = dnorm, size = 0.5, color = 'red') + xlab('Standardized log(Expression)') + ylab('Density')
AM_IM_CM <- detectGenes(AM_IM_CM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(AM_IM_CM),
                                    num_cells_expressed >= 10))
L <- log(exprs(AM_IM_CM[expressed_genes,]))

melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")

pData(AM_IM_CM)$Total_mRNAs <- Matrix::colSums(exprs(AM_IM_CM))

AM_IM_CM <- AM_IM_CM[,pData(AM_IM_CM)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(AM_IM_CM)$Total_mRNAs)) +
                     2*sd(log10(pData(AM_IM_CM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(AM_IM_CM)$Total_mRNAs)) -
                     2*sd(log10(pData(AM_IM_CM)$Total_mRNAs)))

qplot(Total_mRNAs, data = pData(AM_IM_CM),color=orig.ident,  geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

qplot(Total_mRNAs, data = pData(AM_IM_CM),color=initial.id,  geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)



# monocle differential expression test
DE_genes <- differentialGeneTest(AM_IM_CM,fullModelFormulaStr = '~Cluster',
                                 cores = 6)# not sure how many cores but i should try 1 to see
dim(DE_genes)

###
library(dplyr)
#select genes for ordering trajectory
#my_ordering_genes <- row.names (subset(DE_genes, qval < 0.01))
my_ordering_genes <- DE_genes %>% top_n(-200,qval) %>% select(gene_short_name)#DE_genes use pval cutoff or something around here
AM_IM_CM <- setOrderingFilter(AM_IM_CM, ordering_genes = my_ordering_genes)
AM_IM_CM<- reduceDimension(AM_IM_CM, method = 'DDRTree')
plot_ordering_genes(AM_IM_CM)
AM_IM_CM <- orderCells(AM_IM_CM, reverse = F) #can specify root state

save(AM_IM_CM,file="~/Box/Asbestos_SC/03_Analysis_V2/04_Monocle/AM_IM_CM.V2.Rdata")
pdf("~/Box/Asbestos_SC/03_Analysis_V2/04_Monocle/figsV2.pdf")
#plot1
plot_cell_trajectory(AM_IM_CM, color_by = "initial.id",show_branch_points = T)+facet_wrap(~initial.id)
plot_cell_trajectory(AM_IM_CM, color_by = "State")
plot_cell_trajectory(AM_IM_CM, color_by = "res.0.15")
#AM_IM_CM <- orderCells(AM_IM_CM, root_state = 7)
dev.off()
#plot2
plot_cell_trajectory(AM_IM_CM, color_by = "State")