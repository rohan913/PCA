####Quick WGCNA intitial testing for Nasal Epithelium RNAseq Data ####

#following WGCNA tutorial/keeping varnames same: 

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#Read in the data
femData = read.csv("~/Downloads/ne2counts.csv",row.names = 1)
femData = femData[3:ncol(femData)]
colnames(femData) = sapply(colnames(femData), function(x) {paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x)))})
#outliers found (looked at %alignment QCs and hierarchical clustering) removed here
remove = c("T074","T203","T284")
femData = femData[-c(match(remove,colnames(femData)))]

# Take a quick look at what is in the data set:
dim(femData)

datExpr0 = as.data.frame(t(femData))
names(datExpr0) = rownames(femData)
rownames(datExpr0) = colnames(femData)

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(18,12)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
#dev.off()

#got here and saw T074/T203 as outlier and checked in PC3 which is only 4% variance 
###- could exlude if running agian but will not for now

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

traitData = read.csv("~/Downloads/NE_ClinicalTraits.csv");
#custom age groups (just for testing purposes)
ranges = c(20,30,40,50,60,70)
Intervals = findInterval(traitData$Age,ranges)
traitData$Age = Intervals
traitData$Age = paste0(ranges[traitData$Age],"s")

dim(traitData)
names(traitData)

# Form a data frame analogous to expression data that will hold the clinical traits.

femaleSamples = rownames(datExpr0);
traitRows = match(femaleSamples, traitData$ID);
datTraits = traitData[traitRows,-1]
rownames(datTraits) = traitData[traitRows,1]
#J0012 changed from J012 in clinical trial sheet
collectGarbage()


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
#all datTraits must be numeric
datTraits$Age = as.numeric(substr(datTraits$Age,1,nchar(datTraits$Age)-1))
for (i in 1:nrow(datTraits)){
  if (datTraits$Gender[i]=="M"){datTraits$Gender[i]=0}
  else(datTraits$Gender[i]=1)
  if (datTraits$Smoker[i]=="no"){datTraits$Smoker[i]=0}
  else(datTraits$Smoker[i]=1)
  if (datTraits$Group[i]=="Young Healthy"){datTraits$Group[i]=0}
  if (datTraits$Group[i]=="Old Healthy"){datTraits$Group[i]=1}
  if (datTraits$Group[i]=="ILA"){datTraits$Group[i]=2}
}

datTraits2=data.matrix(datTraits)
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits2,signed = F)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = colnames(datTraits2), 
                    main = "Sample dendrogram and trait heatmap")

save(datExpr0, datTraits2, file = "FemaleLiver-01-dataInput.RData")

lnames = load(file = "FemaleLiver-01-dataInput.RData")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
datExpr = datExpr0
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
datExpr = datExpr0
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

datExpr = sapply(datExpr0,as.double)
datExpr = data.matrix(datExpr)

net = blockwiseModules(datExpr, power = 7,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "FemaleLiver-02-networkConstruction-auto.RData")
#tested network construction using "FemaleLiver-02-networkConstruction-auto.RData"