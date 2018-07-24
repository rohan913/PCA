GO <- function(x){
  require(Seurat)
  load(paste0('~/Box/Asbestos_SC/03_Analysis_V2/01_Robj/',x,'.Robj'))
  #get background genes
  background=FindAllMarkers(get(x), thresh.use = 0, min.pct = .1,min.diff.pct = -Inf)
  background=unique(background$gene)
  background=noquote(background)
  if (dir.exists(paste0('~/Box/Asbestos_SC/03_Analysis_V2/03_GO/',x))){setwd(paste0('~/Box/Asbestos_SC/03_Analysis_V2/03_GO/',x))
  }else{dir.create(paste0('~/Box/Asbestos_SC/03_Analysis_V2/03_GO/',x));setwd(paste0('~/Box/Asbestos_SC/03_Analysis_V2/03_GO/',x))}
  write.table(background,file=paste0(x,'_','GObackground.txt'),row.names = F,col.names=F,quote = F, sep="\n")
  #get cluster DE genes
  DE=FindAllMarkers(get(x),only.pos = T)
  DE$cluster <- as.character(DE$cluster)
  for(i in 0:max(unique(DE$cluster))){write.table(DE[DE$cluster==i,]['gene'],file=paste0('C',i,x,'_GOpositivetargets.txt'),row.names = F,col.names=FALSE,quote = F, sep="\n")}
}

GO('AM');GO('Fibroblasts');GO('AT2')