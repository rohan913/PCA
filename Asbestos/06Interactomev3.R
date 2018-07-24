library(readxl)
library(tidyverse)
library(assertthat)
library(gtools)
###load R objects to be used####
load('~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/AM.Robj')
load('~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/AT2.Robj')
load('~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/IM.Robj')
load('~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/Fibroblasts.Robj')
load('~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/CM.Robj')
AM_0 <- SubsetData(AM, ident.use = 0);AM_1 <- SubsetData(AM, ident.use = 1);AM_2 <- SubsetData(AM, ident.use = 2)
IM_0 <- SubsetData(IM, ident.use = 0);IM_1 <- SubsetData(IM, ident.use = 1);IM_2 <- SubsetData(IM, ident.use = 2)
rm(AM);rm(IM)
cell_pops <- ls()
#extract raw data from Robj cells of interest and assign to global env.
for (i in cell_pops){assign(i,get(i)@raw.data[,colnames(get(i)@raw.data) %in% rownames(get(i)@meta.data)])};rm(i)

###read in list of ligand receptor interactions from cellreports paper####
setwd('~/Box/Asbestos_SC/03_Analysis_V2/05_Interactome')
LR_pairs <- read_xlsx('mmc4.xlsx',skip = 4,col_names = T)
ligand_genes <- LR_pairs$ligand_symbol
receptor_genes <- LR_pairs$receptor_symbol

#create list of dataframes of ligand and receptor gene expression for each cell type
#do this by subsetting each celltype obj in global env. to just the rows containing ligand or receptor genes
#that are already known
raw_ligand <- list()
raw_receptor <- list()
for (i in cell_pops){
  raw_ligand[[i]] <- as.data.frame(as.matrix(get(i)))[ligand_genes,]
  raw_receptor[[i]] <- as.data.frame(as.matrix(get(i)))[receptor_genes,]
}
#combine dfs from raw ligand and name each cell
raw_ligand1 <- do.call(cbind,raw_ligand)
raw_receptor1 <-do.call(cbind,raw_receptor)

### function to test different resolutions (pct cell expressing) for interactions detected####
#default 20%
whole=function(jj=.2){xx=(jj)
#function to check if expressed (count >0) in jj percent of cells
expressed <- function(x, cell_pops) {
  assert_that(length(x) == length(cell_pops))
  tapply(x, cell_pops, function(y) mean(y > 0) > xx)
}

#set population names for each cell to be used to count up which cells have expressed ligand
pops <- character()
for(i in cell_pops){pops <- c(pops,rep(i,ncol(get(i))))}
#check if ligand or receptor expressed - ligands and receptors duplicated for each unique L-R combination
ligand_expressed <- apply(raw_ligand1, 1, expressed, cell_pops=pops)
receptor_expressed <- apply(raw_receptor1, 1, expressed, cell_pops=pops)

#NAs mean our dataset did not have those interactions = FALSE (just at cutoff or not present)
ligand_expressed[is.na(ligand_expressed)] <- FALSE
receptor_expressed[is.na(receptor_expressed)] <- FALSE
#logical matrix of interactions of all possible ligand and receptor cells
#below take cross product of the two logical matricies to find whether each L-R pair is expressed
interactions <- tcrossprod(ligand_expressed, receptor_expressed)

#gather each ligand receptor value
idat <- rownames_to_column(as.data.frame(interactions), "ligand") %>%
  gather("receptor", 'value', -ligand)  
#reformat dataframe sorting by rows with least interactions to most
#also refactor and plot
lev <- names(sort(rowMeans(interactions)))
idat$ligand <- factor(idat$ligand, levels=lev)
idat$receptor <- factor(idat$receptor, levels=lev)
#plot just for quick visualization 
gi <- ggplot(idat, aes(x=ligand, y=receptor)) +
  geom_tile(aes(fill=value)) +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_distiller(palette=("Greens"),direction = 1) + theme_bw(base_size=16) +
  theme(axis.text.x=element_text(angle=90))
print(gi+ ggtitle(paste('res',xx)))

#get all permutations of celltypes interacting as ligand and receptor
celltypes <- rownames(ligand_expressed)
perms <- permutations(n=length(celltypes),r=2,v=celltypes,repeats.allowed = T)
perms <- as.data.frame(perms,stringsAsFactors=F)
colnames(perms) <- c('ligand','receptor')
#columns 3 to end of current df become all possible LR molecule pairs 
for (i in 1:nrow(LR_pairs)){
  name <- paste0(LR_pairs$ligand_symbol[i],'_',LR_pairs$receptor_symbol[i])
  perms[,name] <- logical(nrow(perms))
}
#check if each combination of molecules is detected for each combination of cells
#for this we go into ligand expressed and rec expressed logical matrices and check for our celltype of interest
for (i in 3:ncol(perms)){
  genes <- strsplit(colnames(perms)[i],'_')
  liggene <- genes[[1]][1]; recgene <- genes[[1]][2]
  for (j in 1:nrow(perms)){
    liglog <- ligand_expressed[perms$ligand[j],grep(liggene,colnames(ligand_expressed))]#first grep will work
    reclog <- receptor_expressed[perms$receptor[j],grep(recgene,colnames(receptor_expressed))]#bc
    if(length(liglog)==0){liglog=FALSE};if(length(reclog)==0){reclog=FALSE}
    if (liglog[1] & reclog[1]){perms[j,i] <- TRUE}#duplicated genes have the same value TT or FF so only need one
  }
}
return(perms)
}

#resolutions we wanted to test
resns <- c(.01,.05,.1,.15,.2,.25,.3,.35,.4)

for (i in resns){
  permtbl <- whole(jj=i)
  write.csv(permtbl,file = paste0('~/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Tbl.res.',i,'.csv'))
}

#venn diagrams for AMs
Venn_lists <- function (i){
  permtable=read.csv(file = paste0('~/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Tbl.res.',i,'.csv'),row.names = 1,stringsAsFactors = F)
  Ams <- c('AM_0','AM_1','AM_2')
  Fib <- 'Fibroblasts'
  both <- c(Ams,Fib)
  #subset permutation table created earlier to cell pops of interest
  rows=intersect(which(permtable$ligand %in% both),which(permtable$receptor %in% both))
  rows=unique(rows)
  permtable <- permtable[rows,]
  comps <- list()#list of all interactions for each with fibs and seeing which overlap
  for (k in (Ams)){#first line uncommented = AM as ligand Fib as receptor -- selecting these rows only
    toFib <- unique(c(intersect(grep(k,permtable$ligand),grep('Fibroblast',permtable$receptor)),intersect(grep('Fibroblast',permtable$ligand),grep(k,permtable$receptor))))
    #toFib <-c(intersect(grep(k,permtable$ligand),grep('Fibroblast',permtable$receptor)))
    #toFib <- c(intersect(grep('Fibroblast',permtable$ligand),grep(k,permtable$receptor)))
    #total <- sum(rowSums(permtable[toFib,3:ncol(permtable)]))
    names=c()
    for (j in toFib){#for each row selected above get colnames of each TRUE value (this gives us the LR pairs that are expressed for each AM subset to Fib)
      names <- c(names,which(permtable[j,]==TRUE))
    }
    names=colnames(permtable)[unique(names)]
    comps[[k]]=(names)#comps[[k]] becomes comps[['AM_0']], _1, _2... ie which LR interaction is detected for AM_0/1/2 to Fib -> list of 3 vectors 
  }
  require(VennDiagram)
  require(gridExtra)
  #create Venndiagram
  Venn <- venn.diagram(x=comps,filename = NULL,output=T,height = 480,width = 480,
                       compression = 'lzw',resolution = 300,lwd=2,lty='blank',imagetype = 'png',
                       fill=c('blue','green','red'),cex=1,fontface='bold', fontfamily='sans',
                       cat.cex=.6,cat.fontface='bold',
                       cat.fontfamily='sans',
                       rotation=1)
  grid.draw(Venn)#save from here or list file name in above
  #using list of three chr vectors created earlier get lists of overlaps 0/1,0/2,1/2,0/1/2,0,1,2
  ov012 <- intersect(intersect(comps[["AM_0"]],comps[["AM_1"]]),comps[["AM_2"]])
  ov01 <- setdiff(intersect(comps[["AM_0"]],comps[["AM_1"]]),ov012)
  ov02 <- setdiff(intersect(comps[["AM_0"]],comps[["AM_2"]]),ov012)
  ov12 <- setdiff(intersect(comps[["AM_1"]],comps[["AM_2"]]),ov012)
  just0 <- setdiff(comps[["AM_0"]],c(ov01,ov02,ov012))
  just1 <- setdiff(comps[["AM_1"]],c(ov01,ov12,ov012))
  just2 <- setdiff(comps[["AM_2"]],c(ov02,ov12,ov012))
  #create overlap list and pad to same length for converting to dataframe and saving as csv file
  maxi <- max(length(comps[[1]]),length(comps[[2]]),length(comps[[3]]))
  x.pad <- function(x){
    pad <- maxi-length(x)
    return(x <- c(x,rep("X",(pad))))
  }
  topad <- list(ov012,ov01,ov02,ov12,just0,just1,just2)
  for (q in 1:length(topad)){topad[[q]]=x.pad(topad[[q]])}
  overlaps <- do.call(cbind,topad)
  colnames(overlaps) <- c('ov012','ov01','ov02','ov12','just0','just1','just2')
  overlaps <- as.data.frame(overlaps)
  return(list('Venn'=Venn,'overlaps'=overlaps))
}


for (i in resns){
  resulting <- Venn_lists(i)
  pdf(paste0('/Users/rvp277/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Venn AM <> Fib/Resn.',i,'venn.pdf'))
  grid.draw(resulting[["Venn"]])
  dev.off()
  write.csv(resulting[["overlaps"]],file=paste0('/Users/rvp277/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Venn AM <> Fib/Resn.',i,'overlaps.csv'))
}


#### Fib Autocrine ####
FibAuto <- list()
#for each resolution:
#get the table of LR pairs at resolution requested
#select Fib-Fib row (51) and 
#pick which LR pairs are detected (==TRUE)
for (i in resns){
  conn <- paste0('/Users/rvp277/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Tbl.res.',i,'.csv')
  tbl <- read.csv(conn, row.names = 1,stringsAsFactors = F)
  tbl <- tbl[51,c(3:ncol(tbl))]
  tbl <- colnames(tbl)[which(tbl==T)]
  name <- paste0('res.',i)
  FibAuto[[name]] <- tbl
}

#save as one dataframe
FibAuto <- do.call(rbind2,FibAuto)
for (i in 1:length(FibAuto)){
  write.csv(FibAuto[[i]],file=paste0('/Users/rvp277/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Fib Auto/FibAuto',names(FibAuto)[i],'.csv'))
}

### AM Auto ####
#same process for AMs but now it is AM0_0 0_1 0_2 0_1_2 
#similar method to the AM-Ligand -> Fib-Receptor Venn Diagram and overlap table creation
Venn_lists_AM <- function (i){
  permtable=read.csv(file = paste0('~/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Tbl.res.',i,'.csv'),row.names = 1)
  Ams <- c('AM_0','AM_1','AM_2')
  both <- c(Ams)
  rows=c(which(permtable$ligand %in% both),which(permtable$receptor %in% both))
  rows=unique(rows)
  permtable <- permtable[rows,]
  
  comps <- list()#all interactions for each with self and seeing which overlap
  for (k in (Ams)){
    toself <- unique(c(intersect(grep(k,permtable$ligand),grep(k,permtable$receptor))))
    names=c()
    for (j in toself){
      names <- c(names,which(permtable[j,]==TRUE))
    }
    names=colnames(permtable)[unique(names)]
    comps[[k]]=(names)
  }
  require(VennDiagram)
  require(gridExtra)
  Venn <- venn.diagram(x=comps,filename = NULL,output=T,height = 480,width = 480,
                       compression = 'lzw',resolution = 300,lwd=2,lty='blank',imagetype = 'png',
                       fill=c('blue','green','red'),cex=1,fontface='bold', fontfamily='sans',
                       cat.cex=.6,cat.fontface='bold',
                       cat.fontfamily='sans',
                       rotation=1)
  grid.draw(Venn)#save from here or list file name in above
  ov012 <- intersect(intersect(comps[["AM_0"]],comps[["AM_1"]]),comps[["AM_2"]])
  ov01 <- setdiff(intersect(comps[["AM_0"]],comps[["AM_1"]]),ov012)
  ov02 <- setdiff(intersect(comps[["AM_0"]],comps[["AM_2"]]),ov012)
  ov12 <- setdiff(intersect(comps[["AM_1"]],comps[["AM_2"]]),ov012)
  just0 <- setdiff(comps[["AM_0"]],c(ov01,ov02,ov012))
  just1 <- setdiff(comps[["AM_1"]],c(ov01,ov12,ov012))
  just2 <- setdiff(comps[["AM_2"]],c(ov02,ov12,ov012))
  maxi <- max(length(comps[[1]]),length(comps[[2]]),length(comps[[3]]))
  x.pad <- function(x){
    pad <- maxi-length(x)
    return(x <- c(x,rep("X",(pad))))
  }
  topad <- list(ov012,ov01,ov02,ov12,just0,just1,just2)
  for (q in 1:length(topad)){topad[[q]]=x.pad(topad[[q]])}
  overlaps <- do.call(cbind,topad)
  colnames(overlaps) <- c('ov012','ov01','ov02','ov12','just0','just1','just2')
  overlaps <- as.data.frame(overlaps)
  return(list('Venn'=Venn,'overlaps'=overlaps))
}

#testing at all resolutions of interest
for (i in resns){
  resulting <- Venn_lists_AM(i)
  pdf(paste0('/Users/rvp277/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Venn AM Auto/Resn.',i,'venn.pdf'))
  grid.draw(resulting[["Venn"]])
  dev.off()
  write.csv(resulting[["overlaps"]],file=paste0('/Users/rvp277/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Venn AM Auto/Resn.',i,'overlaps.csv'))
}


###resolution images####
pdf('~/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTest.pdf')
whole(jj=.01);whole(jj=.05);whole(jj=.1);whole(jj=.15);whole()
whole(jj=.25);whole(jj=.3);whole(jj=.35);whole(jj=.4)
dev.off()
###AT2 complete####
resns <- c(.01,.05,.1,.15,.2,.25,.3,.35,.4)

#### AT2 Autocrine 
#for each resolution:
#get the table of LR pairs at resolution requested
#select AT2-AT2 row (31) and 
#pick which LR pairs are detected (==TRUE)
AT2Auto <- list()
for (i in resns){
  conn <- paste0('/Users/rvp277/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Tbl.res.',i,'.csv')
  tbl <- read.csv(conn, row.names = 1,stringsAsFactors = F)
  tbl <- tbl[31,c(3:ncol(tbl))]#31 is AT2 <> AT2
  tbl <- colnames(tbl)[which(tbl==T)]
  name <- paste0('res.',i)
  AT2Auto[[name]] <- tbl
}

for (i in 1:length(AT2Auto)){
  write.csv(AT2Auto[[i]],file=paste0('/Users/rvp277/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/AT2 Autocrine/AT2Auto',names(AT2Auto)[i],'.csv'))
}


###AT2 AM instead of Fib

AMs <- c('AM_0','AM_1','AM_2')
AT2s <- c('AT2')

#venn diagrams for AT2s (copied function changed commented sections before running twice)
Venn_lists <- function (i){
  permtable=read.csv(file = paste0('~/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Tbl.res.',i,'.csv'),row.names = 1,stringsAsFactors = F)
  Ams <- c('AM_0','AM_1','AM_2')
  Fib <- 'AT2'
  both <- c(Ams,Fib)
  rows=intersect(which(permtable$ligand %in% both),which(permtable$receptor %in% both))
  rows=unique(rows)
  permtable <- permtable[rows,]
  
  comps <- list()#all interactions for each with fibs and seeing which overlap
  for (k in (Ams)){
    #toFib <-c(intersect(grep(k,permtable$ligand),grep('AT2',permtable$receptor)))#for first set AML AT2R
    toFib <- c(intersect(grep('AT2',permtable$ligand),grep(k,permtable$receptor)))#for second set AMR AT2L
    #total <- sum(rowSums(permtable[toFib,3:ncol(permtable)]))
    names=c()
    for (j in toFib){
      names <- c(names,which(permtable[j,]==TRUE))
    }
    names=colnames(permtable)[unique(names)]
    comps[[k]]=(names)
  }
  require(VennDiagram)
  require(gridExtra)
  Venn <- venn.diagram(x=comps,filename = NULL,output=T,height = 480,width = 480,
                       compression = 'lzw',resolution = 300,lwd=2,lty='blank',imagetype = 'png',
                       fill=c('blue','green','red'),cex=1,fontface='bold', fontfamily='sans',
                       cat.cex=.6,cat.fontface='bold',
                       cat.fontfamily='sans',
                       rotation=1)
  grid.draw(Venn)#save from here or list file name in above
  ov012 <- intersect(intersect(comps[["AM_0"]],comps[["AM_1"]]),comps[["AM_2"]])
  ov01 <- setdiff(intersect(comps[["AM_0"]],comps[["AM_1"]]),ov012)
  ov02 <- setdiff(intersect(comps[["AM_0"]],comps[["AM_2"]]),ov012)
  ov12 <- setdiff(intersect(comps[["AM_1"]],comps[["AM_2"]]),ov012)
  just0 <- setdiff(comps[["AM_0"]],c(ov01,ov02,ov012))
  just1 <- setdiff(comps[["AM_1"]],c(ov01,ov12,ov012))
  just2 <- setdiff(comps[["AM_2"]],c(ov02,ov12,ov012))
  maxi <- max(length(comps[[1]]),length(comps[[2]]),length(comps[[3]]))
  x.pad <- function(x){
    pad <- maxi-length(x)
    return(x <- c(x,rep("X",(pad))))
  }
  topad <- list(ov012,ov01,ov02,ov12,just0,just1,just2)
  for (q in 1:length(topad)){topad[[q]]=x.pad(topad[[q]])}
  overlaps <- do.call(cbind,topad)
  colnames(overlaps) <- c('ov012','ov01','ov02','ov12','just0','just1','just2')
  overlaps <- as.data.frame(overlaps)
  return(list('Venn'=Venn,'overlaps'=overlaps))
}


for (i in resns){
  resulting <- Venn_lists(i)
  pdf(paste0('/Users/rvp277/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Venn AM-L > AT2-R/Resn.',i,'venn.pdf'))
  grid.draw(resulting[["Venn"]])
  dev.off()
  write.csv(resulting[["overlaps"]],file=paste0('/Users/rvp277/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Venn AM-L > AT2-R/Resn.',i,'overlaps.csv'))
}

#second set 
for (i in resns){
  resulting <- Venn_lists(i)
  pdf(paste0('/Users/rvp277/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Venn AM-R < AT2-L/Resn.',i,'venn.pdf'))
  grid.draw(resulting[["Venn"]])
  dev.off()
  write.csv(resulting[["overlaps"]],file=paste0('/Users/rvp277/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Venn AM-R < AT2-L/Resn.',i,'overlaps.csv'))
}


#last part AT2 > Fib and Fib < AT2 at different resolutions
for (i in resns){
  permtable=read.csv(file = paste0('~/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Tbl.res.',i,'.csv'),row.names = 1,stringsAsFactors = F)
  rows=intersect(which(permtable$ligand %in% c('AT2')),which(permtable$receptor %in% c('Fibroblasts')))
  permtable <- permtable[rows,]
  result <- permtable[,3:ncol(permtable)]
  result <- names(result)[which(result==TRUE)]
  write.csv(result,file=paste0('/Users/rvp277/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/AT2-L > Fib-R/Resn.',i,'overlaps.csv'))
}


for (i in resns){
  permtable=read.csv(file = paste0('~/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/Tbl.res.',i,'.csv'),row.names = 1,stringsAsFactors = F)
  rows=intersect(which(permtable$ligand %in% c('Fibroblasts')),which(permtable$receptor %in% c('AT2')))
  permtable <- permtable[rows,]
  result <- permtable[,3:ncol(permtable)]
  result <- names(result)[which(result==TRUE)]
  write.csv(result,file=paste0('/Users/rvp277/Box/Asbestos_SC/03_Analysis_V2/05_Interactome/ResolutionTests/AT2-R < Fib-L/Resn.',i,'overlaps.csv'))
}