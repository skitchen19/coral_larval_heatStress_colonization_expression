#############################
## Not necessary to re-run

#Convert GO table into correct format

library(dplyr)
library(tidyr)

aGO<-read.table("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Supplmental Figures/GO_clusterProfiler/adig_go.txt",
                check.names=FALSE, header=F,
                na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
head(aGO)

colnames(aGO)<-c("GID","GO_terms")

sample_tibble <- aGO %>%
  dplyr::group_by(row_number()) %>%
  dplyr::rename(group="row_number()") %>%
  mutate(GO = strsplit(GO_terms, "; ")) %>%
  unnest(GO) %>%
  mutate(EVIDENCE="blast")%>%
  dplyr::select(-group, -GO_terms)

head(sample_tibble)

aGO2<-as.data.frame(sample_tibble[,2:4])
head(aGO2)
dim(aGO2)

# Write out the formatted table to not repeat everytime
write.table(aGO2,"C:/Users/Sheila's Comp/Documents/GitHub/coral_larval_heatStress_colonization_expression/GO_database/Adig_v1_GO.txt", sep="\t", row.names = F, quote=FALSE)

#############################
# make the database
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationForge")

BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationHub")

library(AnnotationForge)
library(AnnotationHub)

#import the gene symbol table
aSym<-read.table("C:/Users/Sheila's Comp/Documents/GitHub/coral_larval_heatStress_colonization_expression/GO_database/Adig_v1_SYM.txt",
             check.names=FALSE, header=T,
             na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
head(aSym)
dim(aSym)

aGO2<-read.table("C:/Users/Sheila's Comp/Documents/GitHub/coral_larval_heatStress_colonization_expression/GO_database/Adig_v1_GO.txt",
                 check.names=FALSE, header=T,
                 na.strings=c("", "NA"), stringsAsFactors=FALSE, sep="\t",quote="")
head(aGO2)
dim(aGO2)

setwd("E:/OSU/EAPSI/")

## Then call the function to make the database
makeOrgPackage(gene_info=aSym, go=aGO2,
               version="0.3",
               maintainer="Sheila Kitchen <kitchens.osu@gmail.com>",
               author="Sheila Kitchen <kitchens.osu@gmail.com>",
               outputDir = ".",
               tax_id="70779",
               genus="Acropora",
               species="digitifera",
               goTable="go",
               verbose=TRUE)


## then you can call install.packages based on the return value
install.packages("org.Adigitifera.eg.db", repos=NULL,type="source")
##################################

library(clusterProfiler)
library("pathview")
library("enrichplot")

setwd("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Supplmental Figures/GO_clusterProfiler/2022_repeat")

lis<-read.table("all_list_DEGs.txt",
                 check.names=FALSE, header=T, fill=TRUE)
lis<-read.table("all_list_WGCNAmodules.txt",
                check.names=FALSE, header=T, fill=TRUE)
#colnames(lis)<-"GID"
dim(lis)

# Cellular Component
egoCC <- enrichGO(gene         = lis$M25,
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

write.table(egoCC,"M25_CCterms.txt",sep="\t")

# Biological Process
egoBP <- enrichGO(gene         = lis$M20,
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
egoMF <- enrichGO(gene         = lis$M25,
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

write.table(egoMF,"M25_MFterms.txt",sep="\t")

## compare together
cmf <- compareCluster(geneCluster = lis,
                      OrgDb         = "org.Adigitifera.eg.db",
                     fun = "enrichGO",
                     keyType       = 'GID',
                     ont           = "MF",
                     pvalueCutoff = 0.05,
                     qvalueCutoff  = 0.05)
cmf
cmf2 <- simplify(cmf, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(cmf2,font.size = 12,showCategory=20,by="rowPercentage")###Only top 20 p-values displayed

ccc <- compareCluster(geneCluster = lis,
                      OrgDb         = "org.Adigitifera.eg.db",
                      fun = "enrichGO",
                      keyType       = 'GID',
                      ont           = "CC",
                      pvalueCutoff = 0.05,
                      qvalueCutoff  = 0.05)
ccc
ccc2 <- simplify(ccc, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(ccc2,font.size = 12,showCategory=30,by="rowPercentage")###Only top 20 p-values displayed

cbp <- compareCluster(geneCluster = lis,
                      OrgDb         = "org.Adigitifera.eg.db",
                      fun = "enrichGO",
                      keyType       = 'GID',
                      ont           = "BP",
                      pvalueCutoff = 0.05,
                      qvalueCutoff  = 0.05)
cbp
cbp2 <- simplify(cbp, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(cbp2,font.size = 10,showCategory=30,by="rowPercentage")###Only top 20 p-values displayed

##################
# calculate z-scores
install.packages('GOplot')
library("GOplot")

deseq_data<-read.table("./2022_repeat/SC.1d_LFC_padj.txt",
                check.names=FALSE, header=T, fill=TRUE)
goenrich_data<-read.delim("./2022_repeat/SC1_go.txt", sep="\t",
                    check.names=FALSE, header=T, fill=TRUE)

circ <- circle_dat(goenrich_data, deseq_data)

circ$adj_pval <- -log(circ$adj_pval, 10)
head(circ)

meanLogFC<-aggregate(circ[, 6], list(circ$ID), mean)
colnames(meanLogFC)<-c("ID","meanLogFC")
sub <- circ[!duplicated(circ$term), ]

sub2<-sub%>%
  left_join(meanLogFC, by=c("ID"))

rang <- c(min(sub2$count)/10, max(sub2$count)/10)

g <- ggplot(sub2, aes(zscore, adj_pval, fill = meanLogFC, size = count, shape=category))+
  labs(x = 'z-score', y = '-log (adj p-value)')+ xlim(-8.5,7)+ ylim(0,20)+
  geom_point(col = 'black', alpha = 0.6)+
  theme_classic()+scale_shape_manual(values=c(21,22, 23))+
  scale_size(range = rang)+ scale_fill_gradient(low="red", high="blue",limits = c(-2.5,2.5))
g

write.table(circ,"./2022_repeat/AH3_zscore.txt",sep="\t")

###################
setwd("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Supplmental Figures/KEGG_OR/2022_repeat")

# KEGG enrichment (need gene list with KO ids)
k<-read.table("all_v4.txt",header=TRUE, fill=TRUE)
head(k)

kSC1 <- enrichKEGG(gene         = k$SC.1d,
                organism     = "ko",
                 pvalueCutoff = 0.05)
head(as.data.frame(kSC1))
SC1<-as.data.frame(kSC1)

kSC3 <- enrichKEGG(gene         = k$SC.3d,
                   organism     = "ko",
                   pvalueCutoff = 0.05)
head(as.data.frame(kSC3))
SC3<-as.data.frame(kSC3)

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

kSH1 <- enrichKEGG(gene         = k$SH.1d,
                   organism     = "ko",
                   pvalueCutoff = 0.05)
head(as.data.frame(kSH1))
SH1 <-as.data.frame(kSH1)

kSH3 <- enrichKEGG(gene         = k$SH.3d,
                   organism     = "ko",
                   pvalueCutoff = 0.05)
head(as.data.frame(kSH3))
SH3<-as.data.frame(kSH3)

kMET <- enrichKEGG(gene         = k$Met,
                   organism     = "ko",
                   pvalueCutoff = 0.05)
head(as.data.frame(kMET))
met<-as.data.frame(kMET)

#write output to file
#write.table(SC1, "KEGGor_SC1d.tab", sep="\t")
#write.table(SC3, "KEGGor_SC3d.tab", sep="\t")

#write.table(ah1, "KEGGor_AH1d.tab", sep="\t")
#write.table(ah3, "KEGGor_AH3d.tab", sep="\t")

#write.table(SH1, "KEGGor_SH1d.tab", sep="\t")
#write.table(SH3, "KEGGor_SH3d.tab", sep="\t")

#write.table(met, "KEGGor_met.tab", sep="\t")

# WGCNA modules
k<-read.table("all_modules.txt",header=TRUE, fill=TRUE)

kM <- enrichKEGG(gene         = k$M25,
                   organism     = "ko",
                   pvalueCutoff = 0.05)
head(as.data.frame(kM))
mod_kegg<-as.data.frame(kM)
write.table(mod_kegg, "KEGGor_M25.tab", sep="\t")


barplot(kMET, drop=TRUE, showCategory=20, font=12)
dotplot(kSC1, showCategory=30, font=12)
cnetplot(kAH1, showCategory = 8, colorEdge=TRUE)
kMET.2 <- pairwise_termsim(kSH3)
emapplot(kMET.2)

xx <- compareCluster(k, fun="enrichKEGG",
                     organism="ko", pvalueCutoff=0.05)
dotplot(xx, showCategory=20, font=12)
xx <- pairwise_termsim(xx)
emapplot(xx)
