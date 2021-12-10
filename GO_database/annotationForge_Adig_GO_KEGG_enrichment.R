#############################
## Not necessary to re-run

#Convert GO table into correct format

library(dplyr)
library(tidyr)

aGO<-read.table("adig_go.txt",
                check.names=FALSE, header=F,
                na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
head(aGO)

colnames(aGO)<-c("GID","GO_terms")

sample_tibble <- aGO %>%
  dplyr::group_by(row_number()) %>%
  dplyr::rename(group="row_number()") %>%
  mutate(GO = strsplit(GO_terms, ";")) %>%
  unnest(GO) %>%
  mutate(EVIDENCE="blast")%>%
  dplyr::select(-group, -GO_terms)

head(sample_tibble)

aGO2<-as.data.frame(sample_tibble[,2:4])
head(aGO2)
dim(aGO2)

# Write out the formatted table to not repeat everytime
write.table(aGO2,"Adig_v1_GO.txt", sep="\t", row.names = F, quote=FALSE)

#############################
# make the database
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationForge")
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationHub")

library(AnnotationForge)

#import the gene symbol table
setwd("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Supplmental Figures/KEGG_OR")
aSym<-read.table("Adig_v1_SYM.txt",
             check.names=FALSE, header=T,
             na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
head(aSym)
dim(aSym)

aGO2<-read.table("Adig_v1_GO.txt",
                 check.names=FALSE, header=T,
                 na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
head(aGO2)
dim(aGO2)

setwd("E:/OSU/EAPSI/Data")

## Then call the function to make the database
makeOrgPackage(gene_info=aSym, go=aGO2,
               version="0.2",
               maintainer="Sheila <so@someplace.org>",
               author="Sheila <so@someplace.org>",
               outputDir = ".",
               tax_id="70779",
               genus="Acropora",
               species="digitifera",
               goTable="go",
               verbose=TRUE)


## then you can call install.packages based on the return value
install.packages("./org.Adigitifera.eg.db", repos=NULL,type="source")
##################################

library(clusterProfiler)
library("pathview")

setwd("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Supplmental Figures/GO_clusterProfiler")

lis<-read.table("all_list.txt",
                 check.names=FALSE, header=T, fill=TRUE)
#colnames(lis)<-"GID"
dim(lis)

# Cellular Component
egoCC <- enrichGO(gene         = lis$M21,
                 OrgDb         = "org.Adigitifera.eg.db",
                 keyType       = 'GID',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
egoCC
barplot(egoCC)
cc2 <- simplify(egoCC, cutoff=0.7, by="p.adjust", select_fun=min)
cc2
barplot(cc2, showCategory = 30)

write.table(egoCC,"Met_CCterms.txt",sep="\t")

# Biological Process
egoBP <- enrichGO(gene         = lis$M21,
                  OrgDb         = "org.Adigitifera.eg.db",
                  keyType       = 'GID',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05)
egoBP
bp2 <- simplify(egoBP, cutoff=0.7, by="p.adjust", select_fun=min)
bp2
barplot(bp2, showCategory = 30)
write.table(egoBP,"M25_BPterms.txt",sep="\t")


# Molecular Function
egoMF <- enrichGO(gene         = lis$M21,
                  OrgDb         = "org.Adigitifera.eg.db",
                  keyType       = 'GID',
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05)
egoMF
mf2 <- simplify(egoMF, cutoff=0.7, by="p.adjust", select_fun=min)
mf2
barplot(mf2,showCategory = 30)

write.table(egoMF,"M24_MFterms.txt",sep="\t")


cmf <- compareCluster(geneCluster = lis[,8:14],
                      OrgDb         = "org.Adigitifera.eg.db",
                     fun = "enrichGO",
                     keyType       = 'GID',
                     ont           = "MF",
                     pvalueCutoff = 0.05,
                     qvalueCutoff  = 0.05)
cmf
cmf2 <- simplify(cmf, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(cmf2,font.size = 12,showCategory=20,by="rowPercentage")###Only top 20 p-values displayed

ccc <- compareCluster(geneCluster = lis[,8:14],
                      OrgDb         = "org.Adigitifera.eg.db",
                      fun = "enrichGO",
                      keyType       = 'GID',
                      ont           = "CC",
                      pvalueCutoff = 0.05,
                      qvalueCutoff  = 0.05)
ccc
ccc2 <- simplify(ccc, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(ccc2,font.size = 12,showCategory=30,by="rowPercentage")###Only top 20 p-values displayed

cbp <- compareCluster(geneCluster = lis[,8:14],
                      OrgDb         = "org.Adigitifera.eg.db",
                      fun = "enrichGO",
                      keyType       = 'GID',
                      ont           = "BP",
                      pvalueCutoff = 0.05,
                      qvalueCutoff  = 0.05)
cbp
cbp2 <- simplify(cbp, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(cbp2,font.size = 10,showCategory=30,by="rowPercentage")###Only top 20 p-values displayed

# calculate z-scores
install.packages('GOplot')
library("GOplot")

deseq_data<-read.table("Met_genes.txt",
                check.names=FALSE, header=T, fill=TRUE)
goenrich_data<-read.delim("Met_go.txt", sep="\t",
                    check.names=FALSE, header=T, fill=TRUE)

circ <- circle_dat(goenrich_data, deseq_data)

circ$adj_pval <- -log(circ$adj_pval, 10)
head(circ)

meanLogFC<-aggregate(circ[, 6], list(circ$ID), mean)
colnames(meanLogFC)<-c("ID","meanLogFC")
sub <- circ[!duplicated(circ$term), ]

sub2<-sub%>%
  left_join(meanLogFC, by=c("ID"))

rang <- c(min(sub2$count)/20, max(sub2$count)/20)

g <- ggplot(sub2, aes(zscore, adj_pval, fill = meanLogFC, size = count, shape=category))+
  labs(x = 'z-score', y = '-log (adj p-value)')+ xlim(-8.5,7)+ ylim(0,12)+
  geom_point(col = 'black', alpha = 0.6)+
  theme_classic()+scale_shape_manual(values=c(21,22, 23))+
  scale_size(range = rang)+ scale_fill_gradient(low="red", high="blue",limits = c(-1.1,1.1))
g

write.table(circ,"met_zscore.txt",sep="\t")

###################
setwd("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Supplmental Figures/KEGG_OR")

# KEGG enrichment (need gene list with KO ids)
k<-read.table("all_v3_exp.txt",header=TRUE, fill=TRUE)
head(k)

kCL1 <- enrichKEGG(gene         = k$CL.1d,
                 organism     = "ko",
                 pvalueCutoff = 0.05)
head(as.data.frame(kCL1))
cl1<-as.data.frame(kCL1)

kCL3 <- enrichKEGG(gene         = k$CL.3d,
                   organism     = "ko",
                   pvalueCutoff = 0.05)
head(as.data.frame(kCL3))
cl3<-as.data.frame(kCL3)

kAH1 <- enrichKEGG(gene         = k$AH.1d,
                   organism     = "ko",
                   pvalueCutoff = 0.05)
head(as.data.frame(kAH1))
ah1<-as.data.frame(kAH1)

kAH3 <- enrichKEGG(gene         = k$AH.3d,
                   organism     = "ko",
                   pvalueCutoff = 0.05)
head(as.data.frame(kAH3))
ah3<-as.data.frame(kAH3)

kCH1 <- enrichKEGG(gene         = k$CH.1d,
                   organism     = "ko",
                   pvalueCutoff = 0.05)
head(as.data.frame(kCH1))
ch1<-as.data.frame(kCH1)

kCH3 <- enrichKEGG(gene         = k$CH.3d,
                   organism     = "ko",
                   pvalueCutoff = 0.05)
head(as.data.frame(kCH3))
ch3<-as.data.frame(kCH3)

kMET <- enrichKEGG(gene         = k$Met,
                   organism     = "ko",
                   pvalueCutoff = 0.05)
head(as.data.frame(kMET))
met<-as.data.frame(kMET)

#write output to file
#write.table(cl1, "KEGGor_CL1d_720.tab", sep="\t")
#write.table(cl3, "KEGGor_CL3d_720.tab", sep="\t")

#write.table(ah1, "KEGGor_AH1d_720.tab", sep="\t")
#write.table(ah3, "KEGGor_AH3d_720.tab", sep="\t")

#write.table(ch1, "KEGGor_CH1d_720.tab", sep="\t")
#write.table(ch3, "KEGGor_CH3d_720.tab", sep="\t")

#write.table(met, "KEGGor_met_720.tab", sep="\t")

barplot(kMET, drop=TRUE, showCategory=20, font=12)
dotplot(kMET, showCategory=30, font=12)
cnetplot(kMET, showCategory = 8,categorySize="geneNum", colorEdge=TRUE)
emapplot(kMET)

ck <- compareCluster(geneCluster = k[,c(1,3,5,7,9,11)],
                     fun = "enrichKEGG",
                     organism     = "ko",
                     pvalueCutoff = 0.05)
ck2<-as.data.frame(ck)
write.table(ck2, "KEGG_Modules.txt", sep="\t")
dotplot(ck,font.size = 10,showCategory=30,by="rowPercentage")###Only top 20 p-values displayed

ck <- compareCluster(geneCluster = k[,14:20],
                     fun = "enrichKEGG",
                     organism     = "ko",
                     pvalueCutoff = 0.05)
ck2<-as.data.frame(ck)
write.table(ck2, "KEGG_Modules_M24-M25.txt", sep="\t")
dotplot(ck,font.size = 10,showCategory=20,by="rowPercentage")###Only top 20 p-values displayed


mCL1 <- enrichMKEGG(gene = k$CL.1d,
                   organism = 'ko',
                   pvalueCutoff = 0.05)
mcl1<-as.data.frame(mCL1)
head(mcl1)

mCL3 <- enrichMKEGG(gene = k$CL.3d,
                    organism = 'ko',
                    pvalueCutoff = 0.05)
mcl3<-as.data.frame(mCL3)
head(mcl3)

mCH1 <- enrichMKEGG(gene = k$CH.1d,
                    organism = 'ko',
                    pvalueCutoff = 0.05)
mch1<-as.data.frame(mCH1)
head(mch1)

mCH3 <- enrichMKEGG(gene = k$CH.3d,
                    organism = 'ko',
                    pvalueCutoff = 0.05)
mch3<-as.data.frame(mCH3)
head(mch3)

mAH1 <- enrichMKEGG(gene = k$AH.1d,
                    organism = 'ko',
                    pvalueCutoff = 0.05)
mah1<-as.data.frame(mAH1)
head(mah1)

mAH3 <- enrichMKEGG(gene = k$AH.3d,
                    organism = 'ko',
                    pvalueCutoff = 0.05)
mah3<-as.data.frame(mAH3)
head(mah3)

mmet <- enrichMKEGG(gene = k$Met,
                    organism = 'ko',
                    pvalueCutoff = 0.05)
mmet2<-as.data.frame(mmet)
head(mmet2)

library("pathview")
## feature 1: numeric vector
geneListCL1 = k$CL.1d.exp
## feature 2: named vector
names(geneListCL1) = as.character(k$CL.1d)
## feature 3: decreasing orde
geneListCL1 = sort(geneListCL1, decreasing = TRUE)

pathview(gene.data  = geneListCL1,
         pathway.id = "ko04624",
         species    = "ko",
         out.suffix = "CL.1d.layer",
         limit = list(gene = 2),
         low = list(gene = "blue"), mid =
           list(gene = "lightgray"), high = list(gene = "red"))

geneListCL3 = k$CL.3d.exp
## feature 2: named vector
names(geneListCL3) = as.character(k$CL.3d)
## feature 3: decreasing orde
geneListCL3 = sort(geneListCL3, decreasing = TRUE)

pathview(gene.data  = geneListCL3,
         pathway.id = "ko04624",
         species    = "ko",
         out.suffix = "CL.3d.layer",
         limit = list(gene = 2),
         low = list(gene = "blue"), mid =
           list(gene = "lightgray"), high = list(gene = "red"))


geneListCH1 = k$CH.1d.exp
## feature 2: named vector
names(geneListCH1) = as.character(k$CH.1d)
## feature 3: decreasing orde
geneListCH1 = sort(geneListCH1, decreasing = TRUE)

pathview(gene.data  = geneListCH1,
         pathway.id = "ko04624",
         species    = "ko",
         out.suffix = "CH.1d.layer",
         limit = list(gene = 2),
         low = list(gene = "blue"), mid =
           list(gene = "lightgray"), high = list(gene = "red"))

geneListCH3 = k$CH.3d.exp
## feature 2: named vector
names(geneListCH3) = as.character(k$CH.3d)
## feature 3: decreasing orde
geneListCH3 = sort(geneListCH3, decreasing = TRUE)

pathview(gene.data  = geneListCH3,
         pathway.id = "ko04624",
         species    = "ko",
         out.suffix = "CH.3d.layer",
         limit = list(gene = 2),
         low = list(gene = "blue"), mid =
           list(gene = "lightgray"), high = list(gene = "red"))

geneListAH1 = k$AH.1d.exp
## feature 2: named vector
names(geneListAH1) = as.character(k$AH.1d)
## feature 3: decreasing orde
geneListAH1 = sort(geneListAH1, decreasing = TRUE)

pathview(gene.data  = geneListAH1,
         pathway.id = "ko04624",
         species    = "ko",
         out.suffix = "AH.1d.layer",
         limit = list(gene = 2),
         low = list(gene = "blue"), mid =
           list(gene = "lightgray"), high = list(gene = "red"))

geneListAH3 = k$AH.3d.exp
## feature 2: named vector
names(geneListAH3) = as.character(k$AH.3d)
## feature 3: decreasing orde
geneListAH3 = sort(geneListAH3, decreasing = TRUE)

pathview(gene.data  = geneListAH3,
         pathway.id = "ko04624",
         species    = "ko",
         out.suffix = "AH.3d.layer",
         limit = list(gene = 2),
         low = list(gene = "blue"), mid =
           list(gene = "lightgray"), high = list(gene = "red"))


#cell cycle
browseKEGG(kMET, 'ko04110')
#sphingolipids
browseKEGG(kAH3, 'ko04071')
browseKEGG(kAH3, 'ko00600')
#endocytosis
browseKEGG(kAH1, 'ko04144')
#mTOR
browseKEGG(kCL1, 'ko04150')
#NOD signaling
browseKEGG(kAH3, 'ko04621')
#NFkappaB
ko04064
#lysosome
ko04142
#TNF signaling
ko04668
#MAPK signaling
ko04013
#autophagy
ko04140
#oxidative phosphorylation
ko00190
#apoptosis
ko04210
#antifolate resistance
ko01523
#phagosome
ko04145
#fattyacid elongation
ko00062
#c-type lectin signaling
ko04625
#calcium signaling
ko04020
#toll-like receptor signaling
ko04620
#Amphetiamine addiction
ko05031
#Adrenic signaling in cardiocytes
ko04261
#transcription misregulation cancer
ko05202
#hepatocellular carcinoma
ko05225
#viral carcinogins
ko05203
#toll and imd signaling
ko04624
#hippo signaling-fly
ko04391
#cellular senescence
ko04218