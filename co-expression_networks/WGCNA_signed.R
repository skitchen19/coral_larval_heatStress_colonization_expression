library(WGCNA)
library(flashClust)
library(dplyr)
library(stringr)
options(stringsAsFactors = FALSE)

da<-read.csv("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/2022_STAR_FC_countTable_analysis/vst_star_noMetamorph.csv",header=T)
head(da)
colnames(da)= c ("gene","26.1", "26.2", "26.3", "26.4","26.5","26.6","26.1sym", "26.2sym", "26.3sym","26.4sym","26.5sym","26.6sym",
                 "32.1", "32.2", "32.3", "32.4","32.5","32.6","32.1sym", "32.2sym", "32.3sym","32.4sym","32.5sym","32.6sym",
                 "72_26.1", "72_26.2", "72_26.3", "72_26.4","72_26.5","72_26.6",
                 "72_26.1sym", "72_26.2sym","72_26.3sym","72_26.4sym","72_26.5sym","72_26.6sym",
                 "72_32.1", "72_32.2", "72_32.3", "72_32.4","72_32.5","72_32.6",
                 "72_32.1sym", "72_32.2sym","72_32.3sym","72_32.4sym","72_32.5sym","72_32.6sym")

names(da)
datExp=as.data.frame(t(da[,-c(1)]))
names(datExp) = da$gene
rownames(datExp)=names(da)[-c(1)] 

gsg= goodSamplesGenes(datExp, verbose=3)
gsg$allOK
gsg$goodSamples

sampleTree = hclust(dist(datExp), method = "average");
sizeGrWindow(12,9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree, main="Sample", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)
nGenes=ncol(datExp)
nSamples= nrow(datExp)

# Choose a set of soft-thresholding powers
powers = c(c(6:10), seq(from = 12, to=24, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExp, powerVector = powers, verbose = 5, networkType="signed")

# Plot the results:
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

net <- blockwiseModules(datExp, power = 16,maxBlockSize=14000,
                        corType = "bicor", # use robust correlation
                        networkType = "signed", minModuleSize = 15,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = F, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM",
                        verbose = 3)
table(net$colors)

# Convert labels to colors for plotting
moduleColors = net$colors
geneTree = net$dendrograms[[1]]
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(geneTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = net$MEs
MET=orderMEs(MEs)
plotEigengeneNetworks(MET,"",marDendro=c(0,5,1,4.5),marHeatmap=c(2,5,1,2),cex.lab=0.8, xLabelsAngle=90)

#load in trait data
traitData=read.csv("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/2022_STAR_FC_countTable_analysis/WGCNA_repeat/sample_info.csv")
dim(traitData)
names(traitData)

Samples=rownames(datExp)
traitRows=match(Samples, traitData$Samples)
datTraits=traitData[traitRows,-1]
rownames(datTraits)=traitData[traitRows,1]

sampleTree2=flashClust(dist(datExp),method="average")
traitColors=numbers2colors(datTraits, signed=TRUE)
plotDendroAndColors(sampleTree2, traitColors, groupLabels=names(datTraits),main="Samples dendrogram and trait heatmap")

# Define numbers of genes and samples
nGenes = ncol(datExp)
nSamples = nrow(datExp)

moduleTraitCor = bicor(MEs, datTraits, use = "p");
#write.csv(moduleTraitCor, file="moduleTraitCor.csv")

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#write.csv(moduleTraitPvalue, file="moduleTraitPvalue.csv")

# change module names to numbers
modNames = paste0("M", seq(from = 1, to=length(names(net$MEs)), by=1))

modName_color<-as.data.frame(unlist(cbind(modNames, "colName"=colnames(net$MEs))))
modName_color$colName2<-gsub('ME', '', modName_color$colName)
num_col<-as.data.frame(table(net$colors))
modName_color$num<-num_col$Freq[match(modName_color$colName2,as.character(num_col$Var1))]

par(mfrow=c(1,2))
# barplot of size of modules
#library("RColorBrewer")
library("viridis")

n_colors<- length(modNames)
palette <- rev(turbo(n = n_colors))
#palette<-sample(c("#130702", "#5C240A", "#82320d", "#f25d18", "#f8a077", "#0B5383", "#48acf0", "#b4ddf9", "#daeefc", "#ffffff",
#            "#8AA1E5","#7EC5F5","#128FE2","#FCD0BB","#F57F48","#F5E663","#F4A23E","#A08236","#365DD3","#F9DB8F",
#            "#FCEDC7","#8D402A","#1E3888", "#886C5F", "grey30","grey70"),replace =F)

modName_color$palette<-rev(palette)

par(mar = c(5, 2, 3, 2))
par(xaxs = "i")
barplot(-as.numeric(rev(modName_color$num)),horiz=TRUE,axes=FALSE, xlim=c(-3200,0), col=palette, 
        main = paste("Module Size"),space = c(0.1,0.1,0.1,0.1))
axis(1, at= c(0, -500, -1000, -1500, -2000, -2500, -3000, -3500), labels=seq(from = 0, to=3500, by=500))

# Will display correlations and their p-values
textMatrix = ifelse(moduleTraitPvalue< 0.05, paste0(signif(moduleTraitCor, 2), " (p=",
                                                   signif(moduleTraitPvalue, 1), ")"), "")
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(5, 2, 3, 2))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor[,1:3],
xLabels = names(datTraits[1:3]),
yLabels = modNames,
ySymbols = modNames,
colorLabels = FALSE,
colors = colorRampPalette(c("#2B4C8C", "white","#B85B14"), bias=1.05)(100),
textMatrix = textMatrix[,1:3],
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
xpd=T,
main = paste("Module Trait Correlation"))


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

# convert color to numbers
geneColor<-as.data.frame(moduleColors)
geneColor$moduleID<-modName_color[match(geneColor$moduleColors,modName_color$colName2),1]

geneSig<-cbind(geneColor,GS.inf,GS.temp,GS.time,GS.met_Y)
write.csv(geneSig, "gene_correlation_traits.csv")

# module membership
geneModuleMembership = as.data.frame(bicor(datExp, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste(modNames)
names(MMPvalue) = paste("p.", modNames, sep="")

write.csv(geneModuleMembership, "geneModuleMembership.csv")
write.csv(MMPvalue, "geneModulePvalue.csv")

# signed kME values (same as above values)
skme<-signedKME(datExp, MEs,outputColumnName="KME", corFnc="bicor")
write.csv(skme, "signed_kme.csv")

# correlated genes with traits
geneTraitSignificance = as.data.frame(bicor(datExp, temp, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(temp), sep="");
names(GSPvalue) = paste("p.GS.", names(met_Y), sep="");

# plot module membership
module=c("M13")
column=match(module,modNames)
moduleGenes=geneColor$moduleID==module
moduleCol= modName_color[match(module,modName_color$modNames),3]

par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
	abs(geneTraitSignificance[moduleGenes,1]),
	xlab=paste("Module Membership in", module, "module"),
	ylab="Gene Significance for Colonization",
	main=paste("Module membership vs. gene significance\n"),
	cex.main=1.2, cex.lab=1.2, cex.axis=1.2, col=moduleCol)

	
#cmd1=cmdscale(as.dist(dissTOM),2)
#plot(cmd1,col=moduleColors, main="MDS plot", xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")
	
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

write.csv(MEs,"eigengeneValues_Modules.csv")

#barplot of expression
par(mfrow=c(4,1), mar=c(0.3,5.5,3,2))

which.module=1
which.color=modName_color$palette[which.module]
ME2=MEs[,which.module]
#par(mar=c(5,4.2,2,0.7))
barplot(ME2, col=which.color, main="", cex.main=1.5, ylab="eigengene expression", ylim(-0.3,0.3))
#text(x = 52, y = 0.2, label = paste("M",which.module), cex=2)
mtext(paste0("M",which.module), side = 2, adj = 1,line = 5, cex=1.5, las=1)
# if reloading
bar_dat<-read.csv("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/2022_STAR_FC_countTable_analysis/WGCNA_repeat/eigengeneValues_Modules.csv")

palette <- c("#FCEDC7","#82320d","#F4A23E", "#FCD0BB","#0B5383","#365DD3","#130702", 
             "#A08236","#8D402A","#b4ddf9","#128FE2","grey70","#f8a077","#F9DB8F",
             "#8AA1E5","#7EC5F5","#F57F48", "#5C240A","grey30","#F5E663","#48acf0",
             "black","#886C5F","#daeefc", "#f25d18")

par(mfrow=c(8,1), mar=c(2,10,2,1))
which.module=
which.color=palette[which.module]
bar_dat2=bar_dat[,which.module+1]
barplot(bar_dat2, col=which.color, main="", cex.main=1.5, ylab="", ylim=c(-0.3,0.3))
mtext(paste0("M",which.module), side = 2, adj = 0.3,line = 5, cex=1.5, las=1)

# module significance bar plots
gs1=as.numeric(cor(datTraits$Met_Y,datExp,use="p"))
geneSig=abs(gs1)
ModuleSig=tapply(geneSig, moduleColors,mean, na.rm=T)
ModuleSig

par(mfrow=c(1,1), mar=c(3,3,3,3))
plotModuleSignificance(geneSig, moduleColors, cex.main=1)

#t-test for modules

colnames(MEs)<-modNames

test_ID<- numeric(0)
pvalues.d1_26<- numeric(0)
pvalues.d1_32<- numeric(0)
pvalues.d3_26<- numeric(0)
pvalues.d3_32<- numeric(0)
result <- list()

for (i in 1:26) {
  test_id<-print(paste0("M",i))
  pvalues.d1_26<-t.test(MEs[1:6,i],MEs[7:12,i])$p.value
  pvalues.d1_32<-t.test(MEs[13:18,i],MEs[19:24,i])$p.value
  pvalues.d3_26<-t.test(MEs[25:30,i],MEs[31:36,i])$p.value
  pvalues.d3_32<-t.test(MEs[37:42,i],MEs[43:48,i])$p.value
  result[[i]] <- cbind(test_id,pvalues.d1_26,pvalues.d1_32,pvalues.d3_26,pvalues.d3_32)
}

res<-do.call("rbind", lapply(result, data.frame))
res
