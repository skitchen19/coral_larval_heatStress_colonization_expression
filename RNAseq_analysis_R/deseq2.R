#read in file, header false
cd= read.delim ("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/2014/Adig_2014/all_data/ALL_gene_counts.tab",
                head=T, row.names=1)

# assign names to columns
colnames(cd)= c ("27.1", "27.2", "27.3", "27.4","27.5","27.6","27.1sym", "27.2sym", 
"27.3sym","27.4sym","27.5sym","26.7sym", "32.1", "32.2", "32.3", "32.4","32.5","32.6","32.1sym", 
"32.2sym", "32.3sym","32.4sym","32.5sym","32.6sym", "72_27.1", "72_27.2", "72_27.3", "72_27.4",
"72_27.5","72_27.6","72_27.1sym", "72_27.2sym","72_27.3sym","72_27.4sym","72_27.5sym","72_27.6sym", 
"72_32.1", "72_32.2", "72_32.3", "72_32.4","72_32.5","72_32.6","72_32.1sym", "72_32.2sym","72_32.3sym",
"72_32.4sym","72_32.5sym","72_32.6sym")

# number of samples
icolno<-48

# coverage filter
cd<-cd[rowSums(cd)>=icolno,] 

# create data frame (table) of the sample information
colInfo <- data.frame( 
 row.names= colnames(cd),
 temp= c( "27", "27", "27", "27", "27", "27", "27", "27", "27", "27", "27", "27", "32", "32", "32", "32", "32", "32","32", "32", "32", "32", "32", "32","27", "27", "27", "27", "27", "27", "27", "27", "27", "27", "27", "27", "32", "32", "32", "32", "32", "32","32", "32", "32", "32", "32", "32"),
 inf= c("a", "a", "a","a", "a", "a", "s", "s", "s", "s", "s", "s","a", "a", "a", "a", "a", "a","s", "s", "s","s", "s", "s", "a", "a", "a","a", "a", "a", "s", "s", "s", "s", "s", "s","a", "a", "a", "a", "a", "a","s", "s", "s","s", "s", "s"),
 time= c("24","24","24","24","24","24","24","24","24","24","24","24","24","24","24","24","24","24","24","24","24","24","24","24", "72", "72", "72", "72", "72", "72", "72", "72", "72", "72", "72", "72", "72", "72", "72", "72", "72", "72", "72", "72", "72", "72", "72", "72"),
metamorph=c("N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","Y", "N","N","Y", "N","Y","Y", "Y","N","N","Y", "N", "Y", "N","N","N","Y", "Y", "Y","N","Y", "Y")
#metamorph=c("0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","5", "0","0","10", "0","10","40", "85","0","0","50", "0", "50", "0","0","0","80", "60", "30","0","5", "85")
 )

# look at the generated table
colInfo
#changes metamorph from factor to numeric
#colInfo$metamorph<-as.numeric(as.character(colInfo$metamorph))

#open DESeq2
library(DESeq2)

#construct a DESeq data set, last factor will be the one compared unless using a contrast argument, set to blind dispersion estimation (not using design formula)
ddsMF<-DESeqDataSetFromMatrix(countData=cd,
	colData= colInfo,
	design = ~ temp + inf + time + metamorph + temp:inf + temp:time + inf:time + temp:inf:time)

#dataframe information 
colData(ddsMF)

#set aposymbiotic as the control condition
ddsMF$inf<-relevel(ddsMF$inf, "a")
#set 26 degrees to the control temperature
ddsMF$temp<-relevel(ddsMF$temp, "27")
#set 24 hours to the control time
ddsMF$time<-relevel(ddsMF$time, "24")

#run DESeq as Wald's test (combines estimateSizeFactors, estimateDispersion, nbinomWaldTest)
ddsMF<-DESeq(ddsMF, betaPrior=FALSE)

# size factor added
colData(ddsMF)

#information on model tested by DESeq, should be standard
attr(ddsMF, "modelMatrixType") 
attr(ddsMF, "modelMatrix")

#create result table of the 3-way interaction
res<-results(ddsMF)

#move small adjusted pvalue to top of result table
resOrdered<-res[order(res$padj),]
head(resOrdered,15)

#summarize table of 3-way interaction
summary(resOrdered)

#subset the significant value to less than 0.01
resSig<-subset(resOrdered, padj <0.05)

#write to .csv file
#write.csv(as.data.frame(resSig), file= "3-way_interaction.csv")
#write.csv(as.data.frame(res), file= "3-way_interaction_ALL.csv")

ncount<-log2(counts(ddsMF, normalized=T)+1)
#write.table(ncount,"adig_norm_counts.txt", sep="\t")

#give names of result columns for each variable, what is available for building contrast
resultsNames(ddsMF)

#temp*infection interaction, 24h
resMFtempXinf<-results(ddsMF, name="temp32.infs")
summary(resMFtempXinf)
resSigtempXinf<-subset(resMFtempXinf, padj <0.05)
#write.csv(as.data.frame(resSigtempXinf), file= "tempXinf_interaction.csv")
#write.csv(as.data.frame(resMFtempXinf), file= "tempXinf_interaction_ALL.csv")

#time*infection interaction
resMFtimeXinf<-results(ddsMF, name="infs.time72")
summary(resMFtimeXinf)
resSigtimeXinf<-subset(resMFtimeXinf, padj <0.05)
#write.csv(as.data.frame(resSigtimeXinf), file= "timeXinf_interaction.csv")
#write.csv(as.data.frame(resMFtimeXinf), file= "timeXinf_interaction_ALL.csv")

#temp*time interaction
resMFtempXtime<-results(ddsMF, name= "temp32.time72")
summary(resMFtempXtime)
resSigtempXtime<-subset(resMFtempXtime, padj <0.05)
#write.csv(as.data.frame(resSigtempXtime), file= "tempXtime_interaction.csv")
#write.csv(as.data.frame(resMFtempXtime), file= "tempXtime_interaction_ALL.csv")

#time effect at 26deg, aposymbiotic
resMFtime<-results(ddsMF, contrast=c("time","72","24"))
summary(resMFtime)
resMFtimeOrder<-resMFtime[order(resMFtime$padj),]
head(resMFtimeOrder)
resSigtime<-subset(resMFtime, padj <0.05)
#write.csv(as.data.frame(resSigtime), file= "time_mainEffect.csv")
#write.csv(as.data.frame(resMFtime), file= "time_mainEffect_ALL.csv")

#temperature effect in aposymbiotic larvae, time 24
resMFtemp<-results(ddsMF, contrast=c("temp", "32", "27"))
summary(resMFtemp)
resMFtempOrder<-resMFtemp[order(resMFtemp$padj),]
head(resMFtempOrder)
resSigtemp<-subset(resMFtemp, padj <0.05)
#write.csv(as.data.frame(resSigtemp), file= "temp_mainEffect.csv")
#write.csv(as.data.frame(resMFtemp), file= "temp_mainEffect_ALL.csv")

#infection effect at 26deg, time 24
resMFinf<-results(ddsMF, contrast=c("inf", "s", "a"))
summary(resMFinf)
resMFinfOrder<-resMFinf[order(resMFinf$padj),]
head(resMFinfOrder)
resSiginf<-subset(resMFinf, padj <0.05)
#write.csv(as.data.frame(resSiginf), file= "inf_mainEffect.csv")
#write.csv(as.data.frame(resMFinf), file= "inf_mainEffect_ALL.csv")

#metamorph effect at 26deg, aposymbiotic
resMFmet<-results(ddsMF, name="metamorph_Y_vs_N") 
summary(resMFmet)
resMFmetOrder<-resMFmet[order(resMFmet$padj),]
head(resMFmetOrder)
resSigmet<-subset(resMFmet, padj <0.05)
#write.csv(as.data.frame(resSigmet), file= "metamorph_mainEffect.csv")
#write.csv(as.data.frame(resMFmet), file= "metamorph_mainEffect_ALL.csv")

##################
# relevel to 72 hrs for temp x colonization interaction

#set 72 hours to the control time
ddsMF$time<-relevel(ddsMF$time, "72")

#run DESeq as Wald's test (combines estimateSizeFactors, estimateDispersion, nbinomWaldTest)
ddsMF<-DESeq(ddsMF, betaPrior=FALSE)

#temp*infection interaction, 24h
resMFtempXinf72<-results(ddsMF, name="temp32.infs")
summary(resMFtempXinf72)
resSigtempXinf72<-subset(resMFtempXinf72, padj <0.05)
#write.csv(as.data.frame(resSigtempXinf), file= "tempXinf72hr_interaction.csv")
#write.csv(as.data.frame(resMFtempXinf), file= "tempXinf72hr_interaction_ALL.csv")


###################################################
#extract transformed values, rlog=regularized log, vsd=variance stablizing dispersion (Better when libraries are comparable size)
rld<-rlog(ddsMF, blind=FALSE) ##better! blind used when many differences are expected
rlogMat<-assay(rld)
#write.csv(as.data.frame(rlogMat), file= "rlog.csv")

vsd<-varianceStabilizingTransformation(ddsMF, blind=FALSE)
vstMat<-assay(vsd)
#write.csv(as.data.frame(vstMat), file= "vsd.csv")

#sample distances
sampleDists<-dist(t(assay(rld)))

library("RColorBrewer")
sampleDistMatrix<-as.matrix(sampleDists)
rownames(sampleDistMatrix)<-paste(rld$temp, rld$metamorph, rld$inf, rld$time, sep="-")
hc<-hclust(sampleDists)
heatmap.2(sampleDistMatrix,Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col=redblue(48-1), margins=c(2,10), labCol=FALSE)
#dev.copy(png,"sampleDis_rldF.png")
#dev.off()

library("PoiClaClu") #calculates Poisson Distance
poisd<-PoissonDistance(t(counts(ddsMF)))
samplePoisMatrix<-as.matrix(poisd$dd)
rownames(samplePoisMatrix)<-paste(rld$temp, rld$metamorph, rld$inf, rld$time, sep="-")
hc<-hclust(poisd$dd)
heatmap.2(samplePoisMatrix,Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col=redblue(48-1), margins=c(2,10), labCol=FALSE)

########
#setting a grouping variable- compare specific groups as combination of variables
#dds$group<-factor(paste0(dds$temp, dds$inf, dds$time, dds$metamorph))
ddsMF$group<-factor(paste0(ddsMF$temp, ddsMF$inf, ddsMF$time))
design(ddsMF) <-~ group
ddsMF<- DESeq(ddsMF)
resultsNames(ddsMF)

#subset 27apo to 32apo, 24hr
temper<-results(ddsMF, contrast=c("group", "32a24", "27a24"))
summary(temper)
tempOrder<-temper[order(temper$padj),]
head(tempOrder)
tempSig<-subset(temper, padj <0.05)
tmpr<-as.data.frame(temper)
tmpr$id<-row.names(tmpr)
#write.csv(as.data.frame(tempSig), file= "tempONLY24_mainEffect.csv")
#write.csv(as.data.frame(temper), file= "tempONLY24_mainEffect_ALL.csv")

#subset 27sym to 32sym, 24hr
temperS<-results(ddsMF, contrast=c("group", "32s24", "27s24"))
summary(temperS)
tempSOrder<-temperS[order(temperS$padj),]
head(tempSOrder)
tempSigS<-subset(temperS, padj <0.05)
tmprS<-as.data.frame(temperS)
tmprS$id<-row.names(tmprS)

#subset 27apo to 32apo, 72hr
temp72<-results(ddsMF, contrast=c("group", "32a72", "27a72"))
summary(temp72)
temp72Order<-temp72[order(temp72$padj),]
head(temp72Order)
temp72Sig<-subset(temp72, padj <0.05)
#write.csv(as.data.frame(temp72Sig), file= "tempONLY72_mainEffect.csv")
#write.csv(as.data.frame(temp72), file= "tempONLY72_mainEffect_ALL.csv")
tmpr72<-as.data.frame(temp72)
tmpr72$id<-row.names(temp72)

#subset 27sym to 32sym, 72hr
temp72S<-results(ddsMF, contrast=c("group", "32s72", "27s72"))
summary(temp72S)
temp72SOrder<-temp72[order(temp72$padj),]
head(temp72SOrder)
temp72SigS<-subset(temp72S, padj <0.05)
tmpr72S<-as.data.frame(temp72S)
tmpr72S$id<-row.names(temp72S)

#subset 27apo to 27sym,24hr
apo<-results(ddsMF, contrast=c("group", "27s24", "27a24"))
summary(apo)
apoOrder<-apo[order(apo$padj),]
head(apoOrder)
#apoSig<-subset(apo, padj <0.05)
#write.csv(as.data.frame(apoSig), file= "infONLY24_mainEffect.csv")
#write.csv(as.data.frame(apo), file= "infONLY24_mainEffect_ALL.csv")
L24<-as.data.frame(apo)
L24$id<-row.names(apo)

#subset 27apo to 27sym,72 hr
apo72<-results(ddsMF, contrast=c("group", "27s72", "27a72"))
summary(apo72)
apo72Order<-apo72[order(apo72$padj),]
head(apo72Order)
apo72Sig<-subset(apo72, padj <0.05)
#write.csv(as.data.frame(apo72Sig), file= "infONLY72_mainEffect.csv")
#write.csv(as.data.frame(apo72), file= "infONLY72_mainEffect_ALL.csv")
L72<-as.data.frame(apo72)
L72$id<-row.names(apo72)

#subset 32apo to 32sym,24 hr
inter24<-results(ddsMF, contrast=c("group", "32s24", "32a24"))
summary(inter24)
inter24Order<-inter24[order(inter24$padj),]
head(inter24Order)
inter24Sig<-subset(inter24, padj <0.05)
#write.csv(as.data.frame(inter24Sig), file= "tempXinf24_mainEffect.csv")
#write.csv(as.data.frame(inter24), file= "tempXinf24_mainEffect_ALL.csv")
H24<-as.data.frame(inter24)
H24$id<-row.names(inter24)

#subset 32apo to 32sym,72 hr
inter72<-results(ddsMF, contrast=c("group", "32s72", "32a72"))
summary(inter72)
inter72Order<-inter72[order(inter72$padj),]
head(inter72Order)
inter72Sig<-subset(inter72, padj <0.05)
H72_sig<-as.data.frame(inter72Sig)
H72_sig$id<-rownames(H72_sig)
#write.csv(as.data.frame(inter72Sig), file= "tempXinf72_mainEffect.csv")
#write.csv(as.data.frame(inter72), file= "tempXinf72_mainEffect_ALL.csv")
H72<-as.data.frame(inter72)
H72$id<-row.names(inter72)

# delta apo to sym
delta_tmp24<-cbind(Name=row.names(tmpr),"Day 1"=(abs(tmprS$log2FoldChange)-abs(tmpr$log2FoldChange)))
head(delta_tmp24)

delta_tmp72<-cbind(Name=row.names(temp72),"Day 3"=(abs(temp72S$log2FoldChange)-abs(temp72$log2FoldChange)))
head(delta_tmp72)

# load in WGCNA genes
mod<-read.csv("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/wgcna_gene_list.csv",
              check.names=FALSE, header=T, 
              na.strings=c("", "NA"), stringsAsFactors=FALSE, sep=",",quote="")
library(dplyr)
library(reshape2)
library(ggplot2)

tab<-mod %>%
  left_join(as.data.frame(delta_tmp24), by=c("Name")) %>%
  left_join(as.data.frame(delta_tmp72), by="Name")

tab2<-melt(tab, id=c("Name","Module_Color","Module_Number","significant"))
head(tab2)
tab2$Module_Color<-as.factor(tab2$Module_Color)
tab2$Module_Number<-factor(tab2$Module_Number, levels=c("25","24","21","15","14","13","12"))
tab2$significant<-as.factor(tab2$significant)
tab2$value<-as.numeric(paste(tab2$value))

tab3<-subset(tab2,significant == "yes")
res <- dcast(tab3, variable ~ Module_Color, mean)
res

ggplot(tab3,aes(x=Module_Number,y=value,color=Module_Number))+ geom_violin(trim=FALSE)+
  geom_point(size=1.5,alpha=0.2)+ theme_classic()+theme(legend.position="none")+
  stat_summary(fun.y=mean, geom="point", shape=18, size=4)+
  coord_flip()+facet_grid(.~variable)+ylab("Symbiotic log2 FC (32C/27C) - Aposymbiotic log2 FC (32C/27C)")+
  xlab("")+scale_color_manual(values=rev(c("orange","lightcoral","maroon","orangered3","darkgrey","darkmagenta","plum2"))) +
  scale_fill_manual(values=rev(c("orange","lightcoral","maroon","orangered3","darkgrey","darkmagenta","plum2"))) +
  geom_hline(yintercept = 0, linetype="dashed")

# delta temp by symbiotic state
delta_sym_d1<-cbind(Name=row.names(apo),"Day 1"=(abs(H24$log2FoldChange)-abs(L24$log2FoldChange)))
head(delta_sym_d1)

delta_sym_d3<-cbind(Name=row.names(apo72),"Day 3"=(abs(H72$log2FoldChange)-abs(L72$log2FoldChange)))
head(delta_sym_d3)

tab4<-mod %>%
  left_join(as.data.frame(delta_sym_d1), by=c("Name")) %>%
  left_join(as.data.frame(delta_sym_d3), by="Name")

head(tab4)

str(tab4)

tab5<-melt(tab4, id=c("Name","Module_Color","Module_Number","significant"))
head(tab5)
tab5$Module_Color<-as.factor(tab5$Module_Color)
tab5$Module_Number<-factor(tab5$Module_Number, levels=c("25","24","21","15","14","13","12"))
tab5$significant<-as.factor(tab5$significant)
tab5$value<-as.numeric(paste(tab5$value))

tab6<-subset(tab5,significant == "yes")
res2 <- dcast(tab6, variable ~ Module_Color, mean)
res2

ggplot(tab6,aes(x=Module_Number,y=value,color=Module_Number))+ geom_violin(trim=FALSE)+
  geom_point(size=1.5,alpha=0.2)+ theme_classic()+theme(legend.position="none")+
  stat_summary(fun.y=mean, geom="point", shape=18, size=4)+
  coord_flip()+facet_grid(.~variable)+ylab("32 C log2 FC (C/A) - 27 C log2 FC (C/A)")+
  xlab("")+scale_color_manual(values=rev(c("orange","lightcoral","maroon","orangered3","darkgrey","darkmagenta","plum2"))) +
  scale_fill_manual(values=rev(c("orange","lightcoral","maroon","orangered3","darkgrey","darkmagenta","plum2"))) +
  geom_hline(yintercept = 0, linetype="dashed")

###################################################
library("gplots")

#single gene analysis
plotCounts(ddsMF, "aug_v2a.15940", pch=c(21,21,21),intgroup=c("time","temp","inf"), bg=c("orange", "black","red"))

# plot HUB genes for the WGCNA networks
library("ggplot2")
par(mfrow=c(1,3))

#Orange (M12)
data<- plotCounts(ddsMF, "aug_v2a.08636", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p1<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.08636: Nephrocystin-3")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "orange"))

data<- plotCounts(ddsMF, "aug_v2a.08637", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p2<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.08637")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "orange"))

data<- plotCounts(ddsMF, "aug_v2a.01392", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p3<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.01392: Nephrocystin-3")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "orange"))

data<- plotCounts(ddsMF, "aug_v2a.23317", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p4<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.23317: Nephrocystin-3")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "orange"))

data<- plotCounts(ddsMF, "aug_v2a.03979", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p5<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.03979: transcription factor MafK-like")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "orange"))

data<- plotCounts(ddsMF, "aug_v2a.07724", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p6<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.07724: Tenascin-X ")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "orange"))

data<- plotCounts(ddsMF, "aug_v2a.08638", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=FALSE) #requires ggpolot2
p7<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.08638: uncharacterized protein")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "orange"))

data<- plotCounts(ddsMF, "aug_v2a.16588", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p8<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.16588: uncharacterized protein")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "orange"))


library(gridExtra)
grid.arrange(p1, p2, p3,p4,
                  p5,p6,p7,p8,ncol=2)

## Maroon (M14)
data<- plotCounts(ddsMF, "aug_v2a.20600", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p1<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.20600: WD repeat-containing protein 81")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "maroon"))

data<- plotCounts(ddsMF, "aug_v2a.08698", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p2<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.08698: WD repeat-containing protein 81")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "maroon"))

data<- plotCounts(ddsMF, "aug_v2a.01190", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p3<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.01190: transcription elongation regulator 1")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "maroon"))

data<- plotCounts(ddsMF, "aug_v2a.13193", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p4<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.13193: Zinc finger B-box domain-containing protein 1")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "maroon"))

data<- plotCounts(ddsMF, "aug_v2a.10807", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p5<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.10807: ubiquitin carboxyl-terminal hydrolase 7")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "maroon"))

grid.arrange(p1, p2, p3,p4,p5,ncol=2)

## Orangered3 (M15)
data<- plotCounts(ddsMF, "aug_v2a.23405", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p1<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.23405: uncharacterized protein")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "orangered3"))

data<- plotCounts(ddsMF, "aug_v2a.00902", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p2<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.00902: protein argonaute-2")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "orangered3"))

data<- plotCounts(ddsMF, "aug_v2a.04020", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p3<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.04020: Tubulin-specific chaperone E")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "orangered3"))

data<- plotCounts(ddsMF, "aug_v2a.18493", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p4<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.18493: E3 ubiquitin-protein ligase MIB2")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "orangered3"))

data<- plotCounts(ddsMF, "aug_v2a.13046", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p5<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.13046: Ras-related C3 botulinum toxin substrate 1 ")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "orangered3"))

grid.arrange(p1, p2, p3,p4,p5,ncol=2)


## DarkGrey (M21)
data<- plotCounts(ddsMF, "aug_v2a.20705", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p1<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.20705: cyclic AMP-responsive element-binding protein 3")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "darkgrey"))

data<- plotCounts(ddsMF, "aug_v2a.19366", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p2<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.19366: dnaJ homolog subfamily C ")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "darkgrey"))

data<- plotCounts(ddsMF, "aug_v2a.13869", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p3<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.13869: TNF receptor-associated factor 3")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "darkgrey"))

data<- plotCounts(ddsMF, "aug_v2a.14164", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p4<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.14164: regulator of nonsense transcripts 3B")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "darkgrey"))

data<- plotCounts(ddsMF, "aug_v2a.16692", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p5<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.16692: JmjC domain-containing protein 7 ")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "darkgrey"))

data<- plotCounts(ddsMF, "aug_v2a.19197", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p6<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.19197")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "darkgrey"))

data<- plotCounts(ddsMF, "aug_v2a.03511", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p7<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.03511: probable protein disulfide-isomerase A6")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "darkgrey"))

data<- plotCounts(ddsMF, "aug_v2a.07010", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p8<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.07010:tumor necrosis factor receptor superfamily member 16")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "darkgrey"))

grid.arrange(p1, p2, p3,p4,p5,p6,p7,p8,ncol=2)

## lightcoral (M13)
data<- plotCounts(ddsMF, "aug_v2a.01195", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p1<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.01195: adhesion G protein-coupled receptor A3")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "lightcoral"))

data<- plotCounts(ddsMF, "aug_v2a.16170", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p2<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.16170: E3 ubiquitin-protein ligase UBR5 ")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "lightcoral"))

data<- plotCounts(ddsMF, "aug_v2a.11031", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p3<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.11031: protein virilizer")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "lightcoral"))

data<- plotCounts(ddsMF, "aug_v2a.08445", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p4<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.08445: Putative glutamate synthase")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "lightcoral"))

data<- plotCounts(ddsMF, "aug_v2a.10011", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p5<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.10011: ATP-binding cassette sub-family A member 2")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "lightcoral"))

data<- plotCounts(ddsMF, "aug_v2a.16292", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p6<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.16292:vacuolar protein sorting-associated protein 13B")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "lightcoral"))

data<- plotCounts(ddsMF, "aug_v2a.00957", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p7<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.00957: regulator of G-protein signaling 12")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "lightcoral"))

data<- plotCounts(ddsMF, "aug_v2a.05906", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p8<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank()) + 
  ggtitle("aug_v2a.05906: piezo-type mechanosensitive ion channel component 1")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "lightcoral"))


grid.arrange(p1, p2, p3,p4,p5,p6,p7,p8,ncol=2)

## darkmagenta (M24)
data<- plotCounts(ddsMF, "aug_v2a.10329", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p1<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.10329: Sjoegren syndrome nuclear autoantigen 1")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "darkmagenta"))

data<- plotCounts(ddsMF, "aug_v2a.11120", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p2<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.11120: uncharacterized protein LOC107342077")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "darkmagenta"))

data<- plotCounts(ddsMF, "aug_v2a.14382", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p3<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.14382: endothelial differentiation-related factor 1")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "darkmagenta"))

data<- plotCounts(ddsMF, "aug_v2a.16098", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p4<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.16098: nuclear cap-binding protein subunit 3")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "darkmagenta"))

grid.arrange(p1, p2,p3,p4,ncol=2)

## plum2 (M25)
data<- plotCounts(ddsMF, "aug_v2a.01212", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p1<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.01212: uncharacterized protein LOC107352417")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "plum2"))

data<- plotCounts(ddsMF, "aug_v2a.05318", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p2<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.05318: transmembrane emp24 domain-containing protein 3")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "plum2"))

data<- plotCounts(ddsMF, "aug_v2a.04654", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p3<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.04654")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "plum2"))

data<- plotCounts(ddsMF, "aug_v2a.02108", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p4<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.02108: achaete-scute homolog 3")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "plum2"))

data<- plotCounts(ddsMF, "aug_v2a.13958", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p5<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.13958: Probable RNA-dependent RNA polymerase 1")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "plum2"))

data<- plotCounts(ddsMF, "aug_v2a.11316", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p6<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.11316: short coiled-coil protein B")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "plum2"))

data<- plotCounts(ddsMF, "aug_v2a.02038", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p7<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    legend.position = "none") + 
  ggtitle("aug_v2a.02038: superoxide dismutase [Cu-Zn]")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "plum2"))

data<- plotCounts(ddsMF, "aug_v2a.14754", intgroup=c("inf","time", "temp"), 
                  returnData=TRUE, transform=TRUE) #requires ggpolot2
p8<-ggplot(data, aes(x=time, y=log2(count), color=inf, fill=temp)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_dodge(0.7))+
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank()) + 
  ggtitle("aug_v2a.14754: Probable RNA-dependent RNA polymerase 1")+ 
  scale_color_manual(values=c("black", "orange3"))+
  scale_fill_manual(values=c("white", "plum2"))

grid.arrange(p1, p2, p3,p4,p5,p6,p7,p8,ncol=2)

###############################
#principle component analysis
plotPCA(rld, intgroup=c("inf", "temp"))

library("genefilter")
library(vegan)

# plotPCA uses top 500 genes by default
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(500, 
        length(rv)))]
pca <- prcomp(t(assay(rld)[select, ]))

# all genes
#pca <- prcomp(t(assay(rld)))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
 
intgroup.df <- as.data.frame(colData(ddsMF)[, c("temp","metamorph"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
        intgroup.df, names = colnames(ddsMF))

ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group",  shape="group")) + 
        geom_point(size = 4,alpha = 0.5) + scale_shape_manual(values=c(17,17,19,19)) + 
	  scale_color_manual(values=c("gray50","hotpink","black","deeppink3"))+ xlab(paste0("PC1: ", round(percentVar[1] * 
        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
        100), "% variance")) + stat_ellipse(type="t", linetype=2) + theme_bw()+xlim(-50,125)+ylim(-20,20)+
	 theme(text = element_text(size=20))+theme(legend.title=element_blank()) 

# dataframe for permANOVA
intgroup.df <- as.data.frame(colData(ddsMF)[, c("temp","metamorph", "inf","time"), drop = FALSE])
e <- data.frame(intgroup.df, names = colnames(ddsMF))

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(t(assay(rld)) ~ temp + inf + time + metamorph + temp:inf + temp:time+inf:time, data = e, permutations = 10000, method="eu")

#removed metamorph DEGs
rldS<-read.csv("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/2014/RNASeq/All_counts/rlog_subset.csv", head=T, row.names=1)
rv <- rowVars(as.matrix(rldS))
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))] 
pca <- prcomp(t((rldS)[select, ]))

#pca <- prcomp(t(as.matrix(rldS)))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(ddsMF)[, c("temp","metamorph"), drop = FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, intgroup.df, names = colnames(ddsMF))
ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group", shape="group")) + 
        geom_point(size = 4,alpha = 0.5) + scale_shape_manual(values=c(17,17,19,19,18,18,16,16)) +
        xlab(paste0("PC1: ", round(percentVar[1] *100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
        100), "% variance")) + stat_ellipse(type="t", linetype=2) + theme_bw()+
       scale_color_manual(values=c("gray50","hotpink","black","deeppink3"))+xlim(-50,125)+ ylim(-20,20)+
  theme(text = element_text(size=20))+theme(legend.title=element_blank()) 

# PerMANOVA - partitioning the euclidean distance matrix by species
adonis2(t((rldS)) ~ temp + inf + time + metamorph + temp:inf + temp:time + inf:time, data = e, permutations = 10000, method="eu")

########################################################
#correlation graphs of fold change

c<-read.csv("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Figure 3/Correlation_Figure/cor_input.csv", row.names=1, head=T)
head(c)

cor(c$CL.1d_log2,c$AH.1d_log2, use="pairwise.complete.obs", method="pearson")
cor.test(c$CL.1d_log2,c$AH.1d_log2,na.action=omit,use="pairwise.complete.obs", method="pearson")
cor.test(c$CL.3d_log2,c$AH.3d_log2,na.action=omit,use="pairwise.complete.obs", method="pearson")
cor.test(c$CL.1d_log2,c$AH.3d_log2,na.action=omit,use="pairwise.complete.obs", method="pearson")
cor.test(c$CL.3d_log2,c$AH.1d_log2,na.action=omit,use="pairwise.complete.obs", method="pearson")

library(car)
library(scales)
par(mar=c(4,6,4,2))
plot(c$CL.3d_log2,c$AH.3d_log2,type="n",
 xlab=expression("SC (log"[2]*" fold-change)"), ylab=expression("AH (log"[2]*" fold-change)"),
 xlim= c(-2,2), ylim=c(-2,2),font.axis=2,cex.lab=1.5, cex.axis=1.2, las=1, font=2)

#polygon(x = c(-2.2,-2.2,2.2,2.2), y=c(-0.6,0.6,0.6,-0.6) ,col=alpha("darkgrey", 0.5))
#polygon(x = c(-0.6,-0.6,0.6,0.6), y=c(-2.2,2.2,2.2,-2.2) ,col=alpha("darkgrey", 0.5))

points(c$CL.1d_log2,c$AH.1d_log2, col="black",bg=alpha("black",0.1), pch=21, cex=1.1)
points(c$CL.3d_log2,c$AH.3d_log2, col="darkgrey",bg=alpha("darkgrey",0.1), pch=24, cex=1)

l<-lm(c$AH.3d_log2~c$CL.3d_log2,na.action=na.omit)
summary(l)
l2<-lm(c$AH.1d_log2~c$CL.1d_log2,na.action=na.omit)
summary(l2)

abline(lm(c$AH.3d_log2~c$CL.3d_log2,na.action=na.omit),col="darkgrey", lty=1,lwd=2)
abline(lm(c$AH.1d_log2~c$CL.1d_log2, na.action=na.omit),col="black", lty=4,lwd=2)
abline(h=0.6,col="black", lty=1)
abline(h=-0.6,col="black", lty=1)
abline(v=-0.6,col="black", lty=1)
abline(v=0.6,col="black", lty=1)

polygon(x = c(-2.2,-2.2,2.2,2.2), y=c(-0.6,0.6,0.6,-0.6) ,col=alpha("grey", 0.5))
polygon(x = c(-0.6,-0.6,0.6,0.6), y=c(-2.2,2.2,2.2,-2.2) ,col=alpha("grey", 0.5))

legend(0.2,2.6, inset=.01,bty="n", cex=1.2,
  	c("Day 1","Day 3"), pch=c(21,24), col=c("black","darkgrey"),pt.bg=c(alpha("black",0.1),
      alpha("darkgrey",0.1)),pt.cex=c(1.6,1.4), text.font=2, horiz=T,xpd = TRUE)

###############################################
#heat maps of expression data for Figure 5
library("genefilter")
library("gplots")
library("RColorBrewer")
library(ComplexHeatmap)
library("tidyr")
library("reshape2")

fig5<-read.table("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/figure5_heatmap_genes2.txt", 
                 header = TRUE, stringsAsFactors=FALSE, sep="\t")

rld<-read.csv("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/2014/RNASeq/All_counts/rlog.csv", head=T, row.names=1)
head(rld)

#average of replicates
#rlogMat<-as.data.frame(assay(rld)) %>%
rlogMat<-as.data.frame(rld) %>%
  mutate(ID=rownames(.)) %>%
  rowwise() %>%
  mutate(avg_27=mean(c_across(c(1:6))), avg_27sym=mean(c_across(c(7:12))),
          avg_32=mean(c_across(c(13:18))), avg_32sym=mean(c_across(c(19:24))),
         avg_72_27=mean(c_across(c(25:30))), avg_72_27sym=mean(c_across(c(31:36))),
         avg_72_32=mean(c_across(c(37:42))), avg_72_32sym=mean(c_across(c(43:48)))) %>%
  select(ID,avg_27,avg_27sym,avg_32,avg_32sym,avg_72_27,avg_72_27sym,avg_72_32,avg_72_32sym)

head(rlogMat)

fig5_mat<-rlogMat[rlogMat$ID %in% fig5$gene_id, ]

fig5_mat2<-fig5_mat %>%
  left_join(fig5, by=c("ID"="gene_id"))

mat<-fig5_mat[,2:9]
mat<-as.data.frame(t(scale(t(mat))))
mat$pathway <-fig5_mat2$pathway2

fig5_mat3<-mat %>%
 separate(pathway, into = c("pathway", "direction"),sep = "_(?!.*_)") %>%
  dplyr::mutate_if(is.numeric, ~ case_when(direction=="up" & !pathway==c("cell_cycle","autophagy") ~ . * -1, TRUE ~ .)) %>%
  dplyr::mutate_if(is.numeric, ~ case_when(direction=="down" & pathway=="cell_cycle" ~ . * -1, TRUE ~ .)) %>%
  dplyr::mutate_if(is.numeric, ~ case_when(direction=="down" & pathway=="autophagy" ~ . * -1, TRUE ~ .)) %>%
  dplyr::select(-direction) %>%
  dplyr::group_by(pathway) %>%
  dplyr::summarize(across(everything(), median)) %>%
  dplyr::ungroup() 

fig5_mat4<-melt(fig5_mat3, id=c("pathway")) %>%
  mutate(day=ifelse(grepl("72",variable), "day3", "day1"))

fig5_mat4$pathway<-factor(fig5_mat4$pathway, levels=rev(c("ox_stress", "er_stress", "dna_damage","apoptosis-extrinsic", "apoptosis-intrinsic",
                                                         "apoptosis-both","autophagy","immune-recognition","immune-phagosome","inflammation_pathway",
                                                         "cell_cycle","cell_senescence")))

library(viridis)
library(cowplot)
library(ggplot2)
#ggplot(fig5_mat4, aes(x = variable, y=pathway, size = value, fill=value))+
#  geom_point(color="black",pch=21) + theme_bw()+ ylab("")+xlab("")+
#  scale_size_continuous(limits = c(-1.5,1.5), range = c(1,8), breaks = c(-1.5,0,1.5))+
#  scale_fill_gradientn(colors=rev(brewer.pal(n = 11, name = "RdYlBu")),breaks = c(-1,0,1),limits = c(-1.5,1.5)) +
#  #scale_fill_gradient2(low = "#3953A3", mid="white",high = "#EC2024",limits = c(-1.6,1.6))+
#  #scale_fill_viridis(option="magma")+
#  theme(panel.grid.major  = element_line(colour = "white", size = 0.2))+
#  theme(panel.grid.minor  = element_line(colour = "white", size = 0.5))+
#  facet_grid(~ day, scales="free") +theme(axis.text.x=element_text(angle=60,hjust=1))

ggplot(fig5_mat4, aes(x = variable, y=pathway, fill=value))+
  geom_tile(color="black",pch=21) + theme_bw()+ ylab("")+xlab("")+
  scale_fill_gradientn(colors=rev(brewer.pal(n = 11, name = "RdYlBu")),breaks = c(-1,0,1),limits = c(-1.5,1.5)) +
  #scale_fill_gradient2(low = "#3953A3", mid="white",high = "#EC2024",limits = c(-1.6,1.6))+
  #scale_fill_viridis(option="magma")+
  theme(panel.grid.major  = element_line(colour = "white", size = 0.2))+
  theme(panel.grid.minor  = element_line(colour = "white", size = 0.5))+
  facet_grid(~ day, scales="free") +theme(axis.text.x=element_text(angle=60,hjust=1))

ggplot(fig5_mat4, aes(x = variable, y=pathway, size = value, alpha=value))+
  scale_size_continuous(limits = c(-1.5,1.5), range = c(1,8), breaks = c(-1.5,0,1.5))+
  geom_point(color="black",fill="black",pch=21) + theme_bw()+ ylab("")+xlab("")+
  scale_fill_gradientn(colors=rev(brewer.pal(n = 11, name = "RdYlBu")),breaks = c(-1,0,1),limits = c(-1.5,1.5)) +
  theme(panel.grid.major  = element_line(colour = "white", size = 0.2))+
  theme(panel.grid.minor  = element_line(colour = "white", size = 0.5))+
  facet_grid(~ day, scales="free") +theme(axis.text.x=element_text(angle=60,hjust=1))

# plot phenotype data in same way as summarized expression
pheno<-read.table("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Figure 5/phenotype_data.txt", 
                  header = TRUE, stringsAsFactors=FALSE, sep="\t")

pheno$treatment<-factor(pheno$treatment, levels=c("AC", "SC", "AH","SH"))

par(mar=c(0,0,03,0))
a<-ggplot(pheno %>% filter(variable=="colonization"), aes(x = treatment, y=variable, fill=value))+
  #geom_point(color="black",pch=21) + theme_bw()+ ylab("")+xlab("")+
  geom_tile(color="black",pch=21) + theme_bw()+ ylab("")+xlab("")+
  #scale_size_continuous(limits = c(0,100),range=c(1.5,15),breaks=c(0,50,100))+
  #scale_fill_viridis(option="magma",limits = c(0,100))+
  scale_fill_gradientn(colors=rev(brewer.pal(n = 11, name = "RdYlBu")),limits = c(0,100)) +
  facet_grid(~ Day, scales="free") +theme(axis.text.x=element_blank(),plot.margin = unit(c(0,0.1,0,0.5), "cm"),
                                          axis.ticks = element_blank(),panel.border = element_blank())

b<-ggplot(pheno %>% filter(variable=="sym_density"), aes(x = treatment, y=variable, fill=value))+
  #geom_point(color="black",pch=21) + theme_bw()+ ylab("")+xlab("")+
  geom_tile(color="black",pch=21) + theme_bw()+ ylab("")+xlab("")+
  #scale_size_continuous(limits = c(0,110),range=c(1.5,15),breaks=c(0,50,100))+
  #scale_fill_viridis(option="magma",limits = c(0,110))+
  scale_fill_gradientn(colors=rev(brewer.pal(n = 11, name = "RdYlBu")),limits = c(0,110)) +
  facet_grid(~ Day, scales="free") +theme(axis.text.x=element_blank(),
                                          strip.text.x = element_text(size=0),strip.background = element_blank(),
                                          plot.margin = unit(c(-0.8,0.1,0.0,0.5), "cm"),axis.ticks = element_blank(),
                                          panel.border = element_blank())

c<-ggplot(pheno %>% filter(variable=="survival"), aes(x = treatment, y=variable, fill=value))+
  #geom_point(color="black",pch=21) + theme_bw()+ ylab("")+xlab("")+
  geom_tile(color="black",pch=21) + theme_bw()+ ylab("")+xlab("")+
  #scale_size_continuous(limits = c(80,100),range=c(12,15), breaks=c(80,100))+
  #scale_fill_viridis(option="magma",limits = c(80,100))+
  scale_fill_gradientn(colors=rev(brewer.pal(n = 11, name = "RdYlBu")),limits = c(0,100)) +
  facet_grid(~ Day, scales="free") +theme(axis.text.x=element_blank(),
                                          strip.text.x = element_text(size=0),strip.background = element_blank(),
                                          plot.margin = unit(c(-0.8,0.1,0.0,0.5), "cm"),
                                          axis.ticks = element_blank(),panel.border = element_blank())

d<-ggplot(pheno %>% filter(variable=="metamorph"), aes(x = treatment, y=variable, fill=value))+
  #geom_point(color="black",pch=21) + theme_bw()+ ylab("")+xlab("")+
  geom_tile(color="black",pch=21) + theme_bw()+ ylab("")+xlab("")+
  #scale_size_continuous(limits = c(0,60), range=c(1.5,9), breaks=c(0,30,60))+
  #scale_fill_viridis(option="magma",limits = c(0,100))+
  scale_fill_gradientn(colors=rev(brewer.pal(n = 11, name = "RdYlBu")),limits = c(0,100)) +
  facet_grid(~ Day, scales="free") +
  theme(axis.text.x=element_text(angle=60,hjust=1),strip.text.x = element_text(size=0),
        strip.background = element_blank(),plot.margin = unit(c(-0.8,0.1,0.0,0.5), "cm"),
        panel.border = element_blank())

legend <- get_legend(
  # create some space to the left of the legend
  a +  guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "right")
)

pdf("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Figure 5/bubble_plot_pheno.pdf",  width=8, height=4.42)
plot_grid(a,b,c,d,ncol=1, align="v", axis="b")
dev.off()


#cell cycle/DNA damage
fig5_mat_cellcycle<-fig5_mat2 %>% filter(grepl("cell cycle",pathway))

mat<-fig5_mat_cellcycle[,2:9]

ht_opt(heatmap_column_names_gp = gpar(fontface = "italic"), 
       heatmap_column_title_gp = gpar(fontsize = 6),
       heatmap_row_title_gp = gpar(fontsize = 6),
       legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE
)

scaled_mat = t(scale(t(mat)))
#labels<-fig5_mat_cellcycle$gene_symbol
labels<-fig5_mat_cellcycle$annotation
path<-fig5_mat_cellcycle$pathway
cn <- colnames(scaled_mat)
#colors <- colorRampPalette(c("blue","white", "red"), bias=1)(100)
#colors<-colorRampPalette(rocket(12), bias=1)(100)
hm<-Heatmap(scaled_mat, name = " ", col=rev(brewer.pal(n = 11, name = "RdYlBu")), border = TRUE,
            column_names_gp = gpar(fontsize = 8),cluster_columns = FALSE, 
            show_column_names =TRUE, row_names_gp = gpar(fontsize = 6),
            cluster_rows = TRUE, show_row_names = TRUE,
            row_title_rot = 0,
            column_split = factor(rep(c("day 1", "day 3"),each = 4)),
            row_split = as.factor(path),
            heatmap_legend_param=list(at=c(-3,0,3),color_bar="continuous", 
                                      legend_direction="horizontal", legend_width=unit(3,"cm"),
                                      title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
            width = unit(4, "cm"), 
            height = unit(20, "cm")) +
rowAnnotation(labels = anno_text(labels, which = "row",gp = gpar(fontsize = 8)))

draw(hm, heatmap_legend_side = "bottom")

# cellular senescence
fig5_mat_sen<-fig5_mat2 %>% filter(pathway=="cell senescence") #not include all in KEGG pathway

mat<-fig5_mat_sen[,2:9]

scaled_mat = t(scale(t(mat)))
labels<-fig5_mat_sen$gene_symbol
cn <- colnames(scaled_mat)

hm<-Heatmap(scaled_mat, name = " ", col=rev(brewer.pal(n = 11, name = "RdYlBu")), border = TRUE,
            column_names_gp = gpar(fontsize = 8),cluster_columns = FALSE, 
            show_column_names =TRUE, row_names_gp = gpar(fontsize = 6),
            cluster_rows = TRUE, show_row_names = TRUE,
            row_title_rot = 0,
            column_split = factor(rep(c("day 1", "day 3"),each = 4)),
            heatmap_legend_param=list(at=c(-3,0,3),color_bar="continuous", 
                                      legend_direction="horizontal", legend_width=unit(3,"cm"),
                                      title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
            width = unit(4, "cm"), 
            height = unit(10, "cm")) +
  rowAnnotation(labels = anno_text(labels, which = "row",gp = gpar(fontsize = 8)))

draw(hm, heatmap_legend_side = "bottom")

# autophagy
fig5_mat_sen<-fig5_mat2 %>% filter(grepl("autophagy",pathway2)) #not include all in KEGG pathway

mat<-fig5_mat_sen[,2:9]

scaled_mat = t(scale(t(mat)))
labels<-fig5_mat_sen$gene_symbol
cn <- colnames(scaled_mat)

hm<-Heatmap(scaled_mat, name = " ", col=rev(brewer.pal(n = 11, name = "RdYlBu")), border = TRUE,
            column_names_gp = gpar(fontsize = 8),cluster_columns = FALSE, 
            show_column_names =TRUE, row_names_gp = gpar(fontsize = 6),
            cluster_rows = TRUE, show_row_names = TRUE,
            row_title_rot = 0,
            column_split = factor(rep(c("day 1", "day 3"),each = 4)),
            heatmap_legend_param=list(at=c(-3,0,3),color_bar="continuous", 
                                      legend_direction="horizontal", legend_width=unit(3,"cm"),
                                      title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
            width = unit(4, "cm"), 
            height = unit(7, "cm")) +
  rowAnnotation(labels = anno_text(labels, which = "row",gp = gpar(fontsize = 8)))

draw(hm, heatmap_legend_side = "bottom")


# immunity-recognition
fig5_mat_sen<-fig5_mat2 %>% filter(pathway=="immune-recognition") %>% arrange(annotation)

mat<-fig5_mat_sen[,2:9]

scaled_mat = t(scale(t(mat)))
labels<-fig5_mat_sen$gene_symbol
cn <- colnames(scaled_mat)

hm<-Heatmap(scaled_mat, name = " ", col=rev(brewer.pal(n = 11, name = "RdYlBu")), border = TRUE,
            column_names_gp = gpar(fontsize = 8),cluster_columns = FALSE, 
            show_column_names =TRUE, row_names_gp = gpar(fontsize = 6),
            cluster_rows = TRUE, show_row_names = TRUE,
            row_title_rot = 0,
            column_split = factor(rep(c("day 1", "day 3"),each = 4)),
            heatmap_legend_param=list(at=c(-3,0,3),color_bar="continuous", 
                                      legend_direction="horizontal", legend_width=unit(3,"cm"),
                                     title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
            width = unit(4, "cm"), 
            height = unit(11, "cm")) +
  rowAnnotation(labels = anno_text(labels, which = "row",gp = gpar(fontsize = 8)))

draw(hm, heatmap_legend_side = "bottom")

#immunity-phagosome
fig5_mat_sen<-fig5_mat2 %>% filter(pathway=="immune-phagosome") %>% arrange(annotation)

mat<-fig5_mat_sen[,2:9]

scaled_mat = t(scale(t(mat)))
labels<-fig5_mat_sen$gene_symbol
cn <- colnames(scaled_mat)

hm<-Heatmap(scaled_mat, name = " ", col=rev(brewer.pal(n = 11, name = "RdYlBu")), border = TRUE,
            column_names_gp = gpar(fontsize = 8),cluster_columns = FALSE, 
            show_column_names =TRUE, row_names_gp = gpar(fontsize = 6),
            cluster_rows = TRUE, show_row_names = TRUE,
            row_title_rot = 0,
            column_split = factor(rep(c("day 1", "day 3"),each = 4)),
            heatmap_legend_param=list(at=c(-3,0,3),color_bar="continuous", 
                                      legend_direction="horizontal", legend_width=unit(3,"cm"),
                                      title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
            width = unit(4, "cm"), 
            height = unit(2, "cm")) +
  rowAnnotation(labels = anno_text(labels, which = "row",gp = gpar(fontsize = 8)))

draw(hm, heatmap_legend_side = "bottom")

# inflammation
fig5_mat_sen<-fig5_mat2 %>% filter(pathway=="inflammation") %>% arrange(annotation)

mat<-fig5_mat_sen[,2:9]

scaled_mat = t(scale(t(mat)))
labels<-fig5_mat_sen$gene_symbol
cn <- colnames(scaled_mat)

hm<-Heatmap(scaled_mat, name = " ", col=rev(brewer.pal(n = 11, name = "RdYlBu")), border = TRUE,
            column_names_gp = gpar(fontsize = 8),cluster_columns = FALSE, 
            show_column_names =TRUE, row_names_gp = gpar(fontsize = 6),
            cluster_rows = TRUE, show_row_names = TRUE,
            row_title_rot = 0,
            column_split = factor(rep(c("day 1", "day 3"),each = 4)),
            heatmap_legend_param=list(at=c(-3,0,3),color_bar="continuous", 
                                      legend_direction="horizontal", legend_width=unit(3,"cm"),
                                      title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
            width = unit(4, "cm"), 
            height = unit(11, "cm")) +
  rowAnnotation(labels = anno_text(labels, which = "row",gp = gpar(fontsize = 8)))

draw(hm, heatmap_legend_side = "bottom")

# apoptosis
fig5_mat_sen<-fig5_mat2 %>% filter(grepl("apoptosis",pathway)) %>% arrange(annotation)

mat<-fig5_mat_sen[,2:9]

scaled_mat = t(scale(t(mat)))
labels<-fig5_mat_sen$gene_symbol
cn <- colnames(scaled_mat)

hm<-Heatmap(scaled_mat, name = " ", col=rev(brewer.pal(n = 11, name = "RdYlBu")), border = TRUE,
            column_names_gp = gpar(fontsize = 8),cluster_columns = FALSE, 
            show_column_names =TRUE, row_names_gp = gpar(fontsize = 6),
            cluster_rows = TRUE, show_row_names = TRUE,
            row_title_rot = 0,
            column_split = factor(rep(c("day 1", "day 3"),each = 4)),
            heatmap_legend_param=list(at=c(-3,0,3),color_bar="continuous", 
                                      legend_direction="horizontal", legend_width=unit(3,"cm"),
                                      title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
            width = unit(4, "cm"), 
            height = unit(10, "cm")) +
  rowAnnotation(labels = anno_text(labels, which = "row",gp = gpar(fontsize = 8)))

draw(hm, heatmap_legend_side = "bottom")

# oxi/ER stress
fig5_mat_sen<-fig5_mat2 %>% filter(grepl("ox_er",pathway)) %>% arrange(annotation)

mat<-fig5_mat_sen[,2:9]

scaled_mat = t(scale(t(mat)))
labels<-fig5_mat_sen$gene_symbol
cn <- colnames(scaled_mat)

hm<-Heatmap(scaled_mat, name = " ", col=rev(brewer.pal(n = 11, name = "RdYlBu")), border = TRUE,
            column_names_gp = gpar(fontsize = 8),cluster_columns = FALSE, 
            show_column_names =TRUE, row_names_gp = gpar(fontsize = 6),
            cluster_rows = TRUE, show_row_names = TRUE,
            row_title_rot = 0,
            column_split = factor(rep(c("day 1", "day 3"),each = 4)),
            heatmap_legend_param=list(at=c(-3,0,3),color_bar="continuous", 
                                      legend_direction="horizontal", legend_width=unit(3,"cm"),
                                      title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
            width = unit(4, "cm"), 
            height = unit(4, "cm")) +
  rowAnnotation(labels = anno_text(labels, which = "row",gp = gpar(fontsize = 8)))

draw(hm, heatmap_legend_side = "bottom")

################
# expression categories
options(scipen=999)
repTable<-read.table("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/exp_response_table.txt", sep="\t", header=T, row.names = 1)
d1Genes<-read.table("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/day1_DEG_ids.txt", sep="\t", header=F)
d3Genes<-read.table("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/day3_DEG_ids.txt", sep="\t", header=F)

rlogMat<-as.data.frame(assay(rld)) %>%
  mutate(ID=rownames(.)) 

d1<-rlogMat[rlogMat$ID %in% d1Genes$V1,]
d1_sub<-d1[,1:24]
d3<-rlogMat[rlogMat$ID %in% d3Genes$V1,]
d3_sub<-d3[,25:48]

d1_sub[3822,]

which(rlogMat$ID == "aug_v2a.18349")

cor.test(as.numeric(rlogMat[20757,1:24]), as.numeric(repTable[5,]))
cor.test(as.numeric(d1_sub[126,]), as.numeric(repTable[2,]))

# day 1
Test_ID<- numeric(0)
cor.est<- numeric(0)
cor.pval<- numeric(0)
result <- list()

for(i in 1:nrow(d1_sub)){
  result_nested <- list()
  for(j in 1:nrow(repTable)) {
    Test_ID <-paste(d1[i,49],"_",rownames(repTable[j,]))
    cor.est <-  cor.test(as.numeric(d1_sub[i,]),as.numeric(repTable[j,]), method = "pearson")$estimate
    cor.pval <- cor.test(as.numeric(d1_sub[i,]),as.numeric(repTable[j,]), method = "pearson")$p.value
    result_nested[[j]] <- cbind(Test_ID,cor.est,cor.pval)
  }
  result[[i]] <- result_nested
}

res_unlist <- unlist(result, recursive = FALSE)
res<-do.call("rbind", res_unlist)

res <- as.data.frame(res) %>% 
  filter(as.numeric(cor.pval) <= 0.05) %>%
  separate(Test_ID, into = c("GeneID", "Test"),sep = ' _ ') %>%
  group_by(GeneID) %>%
  arrange(desc(as.numeric(cor.est)),.by_group = TRUE) %>%
  slice(n=1)

sumGenes<- res %>% 
  separate(Test, into = c("Category", "Expression"),sep = "_(?=.*)") %>%
  group_by(Category) %>%
  summarize(num=n(), per_total=num/nrow(res))

sumGenes<- res %>% 
  group_by(Test) %>%
  summarize(num=n(), per_total=num/nrow(res))

#write.table(res,"C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Day1_geneCategory.txt", sep="\t", row.names=FALSE)

# day 3
Test_ID<- numeric(0)
cor.est<- numeric(0)
cor.pval<- numeric(0)
result <- list()

for(i in 1:nrow(d3_sub)){
  result_nested <- list()
  for(j in 1:nrow(repTable)) {
    Test_ID <-paste(d3[i,49],"_",rownames(repTable[j,]))
    cor.est <-  cor.test(as.numeric(d3_sub[i,]),as.numeric(repTable[j,]), method = "pearson")$estimate
    cor.pval <- cor.test(as.numeric(d3_sub[i,]),as.numeric(repTable[j,]), method = "pearson")$p.value
    result_nested[[j]] <- cbind(Test_ID,cor.est,cor.pval)
  }
  result[[i]] <- result_nested
}

res_unlist <- unlist(result, recursive = FALSE)
res<-do.call("rbind", res_unlist)

res <- as.data.frame(res) %>% 
  filter(as.numeric(cor.pval) <= 0.05) %>%
  separate(Test_ID, into = c("GeneID", "Test"),sep = ' _ ') %>%
  group_by(GeneID) %>%
  arrange(desc(as.numeric(cor.est)),.by_group = TRUE) %>%
  slice(n=1)

sumGenes<- res %>% 
  separate(Test, into = c("Category", "Expression"),sep = "_(?=.*)") %>%
  group_by(Category) %>%
  summarize(num=n(), per_total=num/nrow(res))

sumGenes<- res %>% 
  group_by(Test) %>%
  summarize(num=n(), per_total=num/nrow(res))

#write.table(res,"C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Day3_geneCategory.txt", sep="\t", row.names=FALSE)

##################
# KOG enrichment test

#install.packages("KOGMWU")
library("KOGMWU")

# load A.digitifera KOG table
KOG_class <- read.table("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/GO_KOGinputfiles/kog_classes.txt", 
                        header = FALSE, stringsAsFactors=FALSE, sep="\t", na.strings = " ", quote="")
AdKOG <- read.table("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/GO_KOGinputfiles/Adig.emapper.annotations.tsv", 
                    header = TRUE, stringsAsFactors=FALSE, sep="\t", na.strings = " ", quote="")
AdKOG<-AdKOG[,c(1,7)]

AdKOG2<-AdKOG %>% separate(COG_category, into = c('COG1', 'COG2'), sep = 1)%>% 
  separate(COG2, into = c('COG2', 'COG3'), sep = 1) %>%
  separate(COG3, into = c('COG3', 'COG4'), sep = 1) %>%
  separate(COG4, into = c('COG4', 'COG5'), sep = 1) %>%
  filter(!(COG1=="-"))

AdKOG3<-melt(AdKOG2, id.vars = c("query")) %>%
  left_join(KOG_class, by=c("value"="V1")) %>%
  filter(!(value=="")) %>%
  select(query,V2)

# sig metamorph genes to filter out
metGenes<-read.table("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Supplmental Figures/KOG/metamorph_DEGs.txt", 
           header = FALSE, stringsAsFactors=FALSE, sep="\t")


# AH larvae, day1
FC_AH<-as.data.frame(temper) %>% mutate(ID=rownames(.)) %>% mutate(logP=-log(pvalue),signedLogP=ifelse(log2FoldChange > 0, (logP*1), (logP*-1)))
FC_AH<-FC_AH[!FC_AH$ID %in% metGenes$V1,]
kog_AH<-kog.mwu(FC_AH[,c(7,2)] %>% drop_na(),AdKOG3) 
kog_AH

# SC larvae, day 1
FC_SC<-as.data.frame(apo) %>% mutate(ID=rownames(.))  %>% mutate(logP=-log(pvalue),signedLogP=ifelse(log2FoldChange > 0, (logP*1), (logP*-1)))
FC_SC<-FC_SC[!FC_SC$ID %in% metGenes$V1,]
kog_SC<-kog.mwu(FC_SC[,c(7,2)] %>% drop_na(),AdKOG3) 
kog_SC

# SH larvae, day 1
FC_TI<-as.data.frame(resMFtempXinf) %>% mutate(ID=rownames(.)) %>% mutate(logP=-log(pvalue),signedLogP=ifelse(log2FoldChange > 0, (logP*1), (logP*-1)))
FC_TI<-FC_TI[!FC_TI$ID %in% metGenes$V1,]
kog_SH<-kog.mwu(FC_TI[,c(7,2)] %>% drop_na(),AdKOG3) 
kog_SH

# AH larvae, day 3
FC_AH3<-as.data.frame(temp72) %>% mutate(ID=rownames(.))  %>% mutate(logP=-log(pvalue),signedLogP=ifelse(log2FoldChange > 0, (logP*1), (logP*-1)))
FC_AH3<-FC_AH3[!FC_AH3$ID %in% metGenes$V1,]
kog_AH_d3<-kog.mwu(FC_AH3[,c(7,2)] %>% drop_na(),AdKOG3) 
kog_AH_d3

# SC larvae, day 3
FC_SC3<-as.data.frame(apo72) %>% mutate(ID=rownames(.))  %>% mutate(logP=-log(pvalue),signedLogP=ifelse(log2FoldChange > 0, (logP*1), (logP*-1)))
FC_SC3<-FC_SC3[!FC_SC3$ID %in% metGenes$V1,]
kog_SC_d3<-kog.mwu(FC_SC3[,c(7,2)] %>% drop_na(),AdKOG3) 
kog_SC_d3

# SH larvae, day 3
FC_TI3<-as.data.frame(resMFtempXinf72) %>% mutate(ID=rownames(.))  %>% mutate(logP=-log(pvalue),signedLogP=ifelse(log2FoldChange > 0, (logP*1), (logP*-1)))
FC_TI3<-FC_TI3[!FC_TI3$ID %in% metGenes$V1,]
kog_SH_d3<-kog.mwu(FC_TI3[,c(7,2)] %>% drop_na(),AdKOG3) 
kog_SH_d3


# Human data
HS_prot2nuc<-read.table("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/GO_KOGinputfiles/bioDBnet_human.txt", 
                        header = TRUE, stringsAsFactors=FALSE, sep="\t", quote="")

HS_prot2nuc_split <- HS_prot2nuc %>%
  dplyr::group_by(row_number()) %>%
  dplyr::rename(group="row_number()") %>%
  mutate(Gene = strsplit(as.character(GenBank.Protein.Accession), "; ")) %>%
  unnest(Gene)

HS_Anno<-read.table("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/GO_KOGinputfiles/Hsapien.emapper.annotations.tsv", 
                       header = TRUE, stringsAsFactors=FALSE, sep="\t", quote="")
HSKOG<-HS_Anno[,c(1,7,9)]

HSKOG2<-HSKOG %>% separate(COG_category, into = c('COG1', 'COG2'), sep = 1)%>% 
  separate(COG2, into = c('COG2', 'COG3'), sep = 1) %>%
  separate(COG3, into = c('COG3', 'COG4'), sep = 1) %>%
  separate(COG4, into = c('COG4', 'COG5'), sep = 1) %>%
  filter(!(COG1=="-"))

HSKOG3<-melt(HSKOG2, id.vars = c("query","Preferred_name")) %>%
  left_join(KOG_class, by=c("value"="V1")) %>%
  filter(!(value=="")) %>%
  select(query,V2,Preferred_name) %>%
  mutate(query = gsub('.[^/.]*$', '', query)) %>%
  left_join(HS_prot2nuc_split, by=c("query"="Gene"))

# quiescence human fibroblasts (GSE66780)
HQ<-read.table("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/GO_KOGinputfiles/GSE19899.top.table_norm_quiescent.tsv", 
               header = TRUE, stringsAsFactors=FALSE, sep="\t", quote="")

HQ2<- HQ %>% mutate(logP=-log(P.Value),signedP=ifelse(logFC > 0, (logP*1), (logP*-1)))
kog_HQ<-kog.mwu(HQ2[,c(7,6)] %>% drop_na(),HSKOG3[,c(3,2)]) 
kog_HQ


# quiescence human fibroblasts (GSE66780)
HQ3<-read.table("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/GO_KOGinputfiles/GSE19899.top.table_serum_quiescent.tsv", 
               header = TRUE, stringsAsFactors=FALSE, sep="\t", quote="")

HQ2a<- HQ3 %>% mutate(logP=-log(P.Value),signedP=ifelse(logFC > 0, (logP*1), (logP*-1)))
kog_HQ2<-kog.mwu(HQ2a[,c(7,6)] %>% drop_na(),HSKOG3[,c(3,2)]) 
kog_HQ2

# senescence human fibroblasts (GSE66780)
HS<-read.table("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/GO_KOGinputfiles/GSE19899.top.table_norm_scenscent.tsv", 
               header = TRUE, stringsAsFactors=FALSE, sep="\t", quote="")
HS2<- HS %>% mutate(logP=-log(P.Value),signedP=ifelse(logFC > 0, (logP*1), (logP*-1)))
kog_HS<-kog.mwu(HS2[,c(7,6)] %>% drop_na(),HSKOG3[,c(3,2)]) 
kog_HS

# compiling a table of delta-ranks to compare these results:
ktable<-makeDeltaRanksTable(list("AH 1dpi"=kog_AH,"SC 1dpi"=kog_SC, "SH 1dpi"=kog_SH,
                                 "AH 3dpi"=kog_AH_d3,"SC 3dpi"=kog_SC_d3,"SH 3dpi"=kog_SH_d3,
                                 "humanQC"=kog_HQ, "humanQS"=kog_HQ2,"humanS"=kog_HS))

ptable<-kog_AH[,c(1,5)] %>%
  left_join(kog_SC[,c(1,5)], by="term", suffix = c("_AH", "_SC")) %>%
  left_join(kog_SH[,c(1,5)], by="term") %>% dplyr::rename(padj_SH = padj) %>%
  left_join(kog_AH_d3[,c(1,5)], by="term") %>% dplyr::rename(padj_AH_d3= padj)%>%
  left_join(kog_SC_d3[,c(1,5)], by="term") %>% dplyr::rename(padj_SC_d3= padj)%>%
  left_join(kog_SH_d3[,c(1,5)], by="term") %>% dplyr::rename(padj_SH_d3= padj)%>%
  left_join(kog_HQ[,c(1,5)], by="term") %>% dplyr::rename(padj_HQC= padj)%>%
  left_join(kog_HQ2[,c(1,5)], by="term") %>% dplyr::rename(padj_HQS= padj)%>%
  left_join(kog_HS[,c(1,5)], by="term") %>% dplyr::rename(padj_HS= padj)%>% 
  mutate_if(is.numeric, round, digits=3)

# Making a heatmap with hierarchical clustering trees:
breaksList = seq(-2000, 2000, by = 40)
colors <- colorRampPalette(c("#2B4C8C", "white","#B85B14"), bias=1.05)(100)
pheatmap(as.matrix(ktable) %>% na.omit(),clustering_distance_cols="correlation", color = colors,breaks = breaksList) 

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)


# plotting individual delta-rank correlations:
library(ggpmisc)
#BiocManager::install("DEGreport")
library(DEGreport)
library(gridExtra)
my.formula <- y ~ x
a<-ggplot(ktable, aes(x=`AH 1dpi`,y=`SH 1dpi`))+
  geom_point(size=2) +
  geom_smooth(method='lm', formula= my.formula, col="#E77823",se = FALSE)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +theme(panel.grid = element_blank())+
  geom_cor(method = "pearson", xpos=500, ypos = 500, inherit.aes=TRUE)
b<-ggplot(ktable, aes(x=`humanS`,y=`SH 1dpi`))+
  geom_point(size=2) +
  geom_smooth(method='lm', formula= my.formula, col="#E77823",se = FALSE)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +theme(panel.grid = element_blank())+
  geom_cor(method = "pearson", xpos=500, ypos = 500, inherit.aes=TRUE)
c<-ggplot(ktable, aes(x=`humanQS`,y=`SH 1dpi`))+
  geom_point(size=2) +
  geom_smooth(method='lm', formula= my.formula, col="#E77823",se = FALSE)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +theme(panel.grid = element_blank())+
  geom_cor(method = "pearson", xpos=500, ypos = 500, inherit.aes=TRUE)
d<-ggplot(ktable, aes(x=`SC 3dpi`,y=`SH 3dpi`))+
  geom_point(size=2) +
  geom_smooth(method='lm', formula= my.formula, col="#E77823",se = FALSE)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +theme(panel.grid = element_blank())+
  geom_cor(method = "pearson", xpos=500, ypos = 500, inherit.aes=TRUE)
e<-ggplot(ktable, aes(x=`humanS`,y=`SH 3dpi`))+
  geom_point(size=2) +
  geom_smooth(method='lm', formula= my.formula, col="#E77823",se = FALSE)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +theme(panel.grid = element_blank())+
  geom_cor(method = "pearson", xpos=500, ypos = 500, inherit.aes=TRUE)
f<-ggplot(ktable, aes(x=`humanQS`,y=`SH 3dpi`))+
  geom_point(size=2) +
  geom_smooth(method='lm', formula= my.formula, col="#E77823",se = FALSE)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +theme(panel.grid = element_blank())+
  geom_cor(method = "pearson", xpos=500, ypos = 500, inherit.aes=TRUE)
grid.arrange(a, b, c,d,e,f, nrow = 2, ncol = 3)