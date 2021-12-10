library(WGCNA)
library(flashClust)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# load in vsd transformed counts, filtered to remove DEGs with metamorphosis
da<-read.table("vsd-metamorph.txt",header=T)
head(da)
colnames(da)= c ("26.1", "26.2", "26.3", "26.4","26.5","26.6",
	"26.1sym", "26.2sym", "26.3sym","26.4sym","26.5sym","26.6sym",
	"32.1", "32.2", "32.3", "32.4","32.5","32.6",
	"32.1sym", "32.2sym", "32.3sym","32.4sym","32.5sym","32.6sym",
	"72_26.1", "72_26.2", "72_26.3", "72_26.4","72_26.5","72_26.6",
	"72_26.1sym", "72_26.2sym","72_26.3sym","72_26.4sym","72_26.5sym","72_26.6sym", 
	"72_32.1", "72_32.2", "72_32.3", "72_32.4","72_32.5","72_32.6",
	"72_32.1sym", "72_32.2sym","72_32.3sym","72_32.4sym","72_32.5sym","72_32.6sym")

names(da)
datExp=as.data.frame(t(da[,-c(1)]))
names(datExp) = da$gene
rownames(datExp)=names(da)[-c(1)] 

# filter out samples/genes with too much missing data
gsg= goodSamplesGenes(datExp, verbose=3)
gsg$allOK

# hierchical clustering of expression data
sampleTree = hclust(dist(datExp), method = "average");
sizeGrWindow(12,9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree, main="Sample", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)

#remove outlier sample if necessary
abline(h=15, col="red")
clust=cutreeStatic(sampleTree, cutHeight=15, minSize=10)
table(clust)
keepSamples= (clust==1)
datExp= datExper0[keepSamples,]
nGenes=ncol(datExp)
nSamples= nrow(datExp)

# Choose a set of soft-thresholding powers
powers = c(c(6:10), seq(from = 12, to=24, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExp, powerVector = powers, verbose = 5, networkType="signed")
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
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

#Single Block analysis
#set power and construct matrix of adjacency
softPower = 16;
adj = adjacency(datExp, power = softPower, type="signed", corFnc="bicor");
# Turn adjacency into topological overlap0
dissTOM = TOMdist(adj)

# Call the hierarchical clustering function
geneTree=flashClust(as.dist(dissTOM), method="average")

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04);

# Set minimum gene set size
minModuleSize = 15;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExp, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-bicor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

# set threshold for dissimilarity
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExp, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
MET=orderMEs(MEs)
plotEigengeneNetworks(MET,"",marDendro=c(0,5,1,4.5),marHeatmap=c(2,5,1,2),cex.lab=0.8, xLabelsAngle=90)

#load in trait data
traitData=read.csv("sample_info.csv")
dim(traitData)
names(traitData)

Samples=rownames(datExp)
traitRows=match(Samples, traitData$Samples)
datTraits=traitData[traitRows,-1]
rownames(datTraits)=traitData[traitRows,1]

sampleTree2=flashClust(dist(datExp),method="average")
traitColors=numbers2colors(datTraits, signed=FALSE)
plotDendroAndColors(sampleTree2, traitColors, groupLabels=names(datTraits),main="Samples dendrogram and trait heatmap")

# Define numbers of genes and samples
nGenes = ncol(datExp);
nSamples = nrow(datExp);
# Recalculate MEs with color labels
moduleTraitCor = bicor(MEs, datTraits, use = "p");
write.csv(moduleTraitCor, file="moduleTraitCor.csv")

# Calculate Student asymptotic p-value for given correlations:
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#????
mtp<-data.frame(cbind("r" = moduleTraitCor, "p" = moduleTraitPvalue))
names(mtp)<-c("r","p")
mtp
rownames(mtp)<-modNames
mtp
mtps<-mtp[mtp$p<0.05,]
mtps

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 10, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))

# GS and MM calculations
# Define variable weight containing the weight column of datTrait
met_Y = as.data.frame(datTraits$Met_Y);
names(met_Y) = "met_Y"
GS.met_Y=as.numeric(bicor(datExp, datTraits$Met_Y, use="p"))
GS.met_YColor=numbers2colors(GS.met_Y, signed=T)

temp= as.data.frame(datTraits$Temperature);
names(temp) = "temp"
GS.temp=as.numeric(bicor(datExp, datTraits$Temperature, use="p"))
GS.tempColor=numbers2colors(GS.temp, signed=T)

time= as.data.frame(datTraits$Time);
names(time) = "time"
GS.time=as.numeric(bicor(datExp, datTraits$Time, use="p"))
GS.timeColor=numbers2colors(GS.time, signed=T)

inf= as.data.frame(datTraits$Infection);
names(inf) = "infection"
GS.inf=as.numeric(bicor(datExp, datTraits$Infection, use="p"))
GS.infColor=numbers2colors(GS.inf, signed=T)

plotDendroAndColors(geneTree, cbind(moduleColors,GS.tempColor,GS.infColor,GS.timeColor,GS.met_YColor),
c("Merged modules", "GS.temperature", "GS.infection", "GS.time", "GS.metamorphisis"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(bicor(datExp, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(bicor(datExp, temp, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(temp), sep="");
names(GSPvalue) = paste("p.GS.", names(met_Y), sep="");

module=c("purple")
column=match(module,modNames)
moduleGenes=moduleColors==module

sizeGrWindow(7,7)
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
	abs(geneTraitSignificance[moduleGenes,1]),
	xlab=paste("Module Membership in", module, "module"),
	ylab="Gene Significance for Colonization",
	main=paste("Module membership vs. gene significance\n"),
	cex.main=1.2, cex.lab=1.2, cex.axis=1.2, col=module)
	
cmd1=cmdscale(as.dist(dissTOM),2)
plot(cmd1,col=moduleColors, main="MDS plot", xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")
	
#write to file with annotation	
annot=read.csv(file="all.annotation.ADIG.csv")	
dim(annot)
names(annot)
probes=names(datExp)
probes2annot=match(probes,annot$Gene_Name)
sum(is.na(probes2annot))

geneInfo0= data.frame(Gene_Name = probes, geneSymbol= annot$Gene.Symbol[probes2annot],
	KEGG= annot$KEGG[probes2annot],
	GOterms=annot$GO_terms[probes2annot],
	moduleColors=moduleColors,
	geneTraitSignificance,
	GSPvalue)
names(geneInfo0)
modOrder=order(-abs(cor(MEs, met_Y, use="p")))

for (mod in 1:ncol(geneModuleMembership))
{
	oldNames=names(geneInfo0)
	geneInfo0=data.frame(geneInfo0,geneModuleMembership[,modOrder[mod]],
		MMPvalue[,modOrder[mod]]);
	names(geneInfo0)= c(oldNames, paste("MM", modNames[modOrder[mod]], sep= ""), paste("p.MM", modNames[modOrder[mod]], sep=""))
}

geneOrder=order(geneInfo0$moduleColor, -abs(geneInfo0$GS.met_Y))
geneInfo=geneInfo0[geneOrder,]
write.csv(geneInfo, file="geneInfo.csv")

# export for cytoscape
#select module
modules= c("blue")
inModule=is.finite(match(moduleColors, modules))
modProbes=probes[inModule]
modGenes=annot$Gene.Symbol[match(modProbes,annot$Gene_Name)]
#TOM=TOMsimilarityFromExpr(datExp,power=9)
modTOM=TOM[inModule, inModule]
dimnames(modTOM)= list(modProbes, modProbes)
cyt=exportNetworkToCytoscape(modTOM[top,top],
	edgeFile=paste("CytoscapeInput-edges-",paste(modules, collapse="-"), "top.txt", sep=""),
	nodeFile=paste("CytoscapeInput-nodes-",paste(modules, collapse="-"), "top.txt", sep=""),
	weighted=TRUE,
	threshold=0.02,
	nodeNames=modProbes,
	altNodeNames=modGenes,
	nodeAttr=moduleColors[inModule])
	
#export to VisANT
nTop=30
IMConn=softConnectivity(datExp[,modProbes])
top=(rank(-IMConn)<= nTop)

vis=exportNetworkToVisANT(modTOM,
	file=paste("VisANTInput-", modules, "top30.txt", sep=""),
	weighted=TRUE,
	threshold=0.02,
	probeToGene=NULL)

#heatmap of genes	
plotTOM=dissTOM^7
diag(plotTOM)=NA
TOMplot(plotTOM, geneTree, moduleColors, main="Network heatmap, all genes")
	
#barplot of expression	
sizeGrWindow(8,10)
which.module="darkgrey"
ME2=MEs[,paste("ME",which.module,sep="")]

par(mar=c(5,4.2,2,0.7))
barplot(ME2, col=which.module, main="", cex.main=1.5, ylab="eigengene expression")

gs1=as.numeric(cor(datTraits$Me_num,datExp,use="p"))
geneSig=abs(gs1)
ModuleSig=tapply(geneSig, moduleColors,mean, na.rm=T)
ModuleSig
sizeGrWindow(8,7)
par(mfrow=c(1,1), mar=c(3,3,3,3))
plotModuleSignificance(geneSig, moduleColors, cex.main=1)

plotMEpairs(MEs, y=datTraits$Time)

# save data to file
save(MEs,moduleLabels,moduleColors,geneTree, TOM, dissTOM, cmd1, file="stepByStep-networkConstruction.RData")