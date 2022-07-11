##clustering interaction DEGs##
#deg<-read.csv("./Interactions_graph/in_combo.csv", header=T,row.names=1)
#deg<-read.csv("./CL_time/in_combo_CL.csv", header=T,row.names=1)
#deg<-read.csv("./AH_time/in_combo_AH.csv", header=T,row.names=1)

# day1 sig interaction genes
deg1<-as.data.frame(resSigtempXinf[!rownames(resSigtempXinf) %in% metID_all$id,]) %>% mutate(ID=rownames(.))

d1_rlog<-rlogMat_subset[rownames(rlogMat_subset) %in% deg1$ID,1:24]

d1_rlog_avg<-as.data.frame(d1_rlog) %>%
    mutate(ID=rownames(.)) %>% 
    rowwise() %>%
    dplyr::mutate(avg_27=mean(c_across(c(1:6))), avg_27sym=mean(c_across(c(7:12))),
                  avg_32=mean(c_across(c(13:18))), avg_32sym=mean(c_across(c(19:24))),
                  A_C=avg_27-avg_27,S_C=avg_27sym-avg_27,A_H=avg_32-avg_27, S_H=avg_32sym-avg_27, Day="1") %>%
    dplyr::select(ID,A_C, S_C, A_H, S_H, Day)

# day3 sig interaction genes
deg3<-as.data.frame(resSigtempXinf72[!rownames(resSigtempXinf72) %in% metID_all$id,]) %>% mutate(ID=rownames(.)) 

d3_rlog<-rlogMat_subset[rownames(rlogMat_subset) %in% deg3$ID,25:48]

d3_rlog_avg<-as.data.frame(d3_rlog) %>%
    mutate(ID=rownames(.)) %>% 
    rowwise() %>%
    dplyr::mutate(avg_27=mean(c_across(c(1:6))), avg_27sym=mean(c_across(c(7:12))),
                  avg_32=mean(c_across(c(13:18))), avg_32sym=mean(c_across(c(19:24))),
                  A_C=avg_27-avg_27,S_C=avg_27sym-avg_27,A_H=avg_32-avg_27, S_H=avg_32sym-avg_27, Day="3") %>%
    dplyr::select(ID,A_C, S_C, A_H, S_H, Day)

# cluster by distance day 1
d1_dist <- dist(as.matrix(d1_rlog_avg[,2:5]), method = "euclidean")
hc<-hclust(d1_dist)
p<-plot(hc, hang = -0.07, cex = 0.5, labels =d1_rlog_avg$ID)
abline(h=1.5, col="red", untf = FALSE)
d1_rlog_avg$cluster<-cutree(hc, k=6)

# cluster by distance day 3
d3_dist <- dist(as.matrix(d3_rlog_avg[,2:5]), method = "euclidean")
hc3<-hclust(d3_dist)
plot(hc3, hang = -0.07, cex = 0.5, labels =d3_rlog_avg$ID)
abline(h=1.5, col="red", untf = FALSE)
d3_rlog_avg$cluster<-cutree(hc3, k=2)

d3_rlog_avg$cluster[d3_rlog_avg$cluster=="1"]<-"7"
d3_rlog_avg$cluster[d3_rlog_avg$cluster=="2"]<-"8"

# join tables
d1_d3_rlog_avg<-rbind(d1_rlog_avg,d3_rlog_avg)

write.csv(d1_d3_rlog_avg, file = "hclust_interaction.csv")

##########################
#interaction graphs

require(ggplot2)
require(ggvis)
require(gridExtra)
library(tidyverse)
library(dplyr)
library(reshape)

d1_d3_rlog_avg2 <- melt(d1_d3_rlog_avg, id=c("ID","cluster","Day"))
head(d1_d3_rlog_avg2)

d1_d3_rlog_avg3<-d1_d3_rlog_avg2 %>% tidyr::separate(variable, c("symbiont", "temp"), sep="_")
d1_d3_rlog_avg3$cluster<-factor(d1_d3_rlog_avg3$cluster)
d1_d3_rlog_avg3$temp<-factor(d1_d3_rlog_avg3$temp)
d1_d3_rlog_avg3$symbiont<-factor(d1_d3_rlog_avg3$symbiont)

# Summarizing data
data2 <- d1_d3_rlog_avg3 %>% 
    group_by(cluster,temp,symbiont) %>% 
    dplyr::summarize(exp_mean = mean(value),
              exp_se = psych::describe(value)$se)


ggplot(data = d1_d3_rlog_avg3, aes(x = temp, y = value)) +
    theme_classic()+ ylim(-1.7,1.7)+
    ylab("Expression")+ xlab("Day")+ 
    geom_point(aes(fill=symbiont),color="black",size=4,pch=21, alpha=c(0.2)) + 
    geom_line(data=data2,aes(x = temp, y = exp_mean, group = symbiont, linetype=symbiont), size=1.3)+
    theme(legend.position="top")+
    facet_wrap(. ~ cluster,ncol = 2)+
    scale_fill_manual(values=c("white", "grey40"))+
    scale_color_manual(values=c("black", "black"))+
    scale_linetype_manual(values=c( "dotdash","solid"))



#########################
## Venn diagram CL and AH
library(VennDiagram)
coral_venn <- read.csv("C:/Users/Sheila's Comp/Documents/GitHub/coral_larval_heatStress_colonization_expression/RNAseq_analysis_R/fold_change_d1_d3_CL_AH_STAR.csv",
                         check.names=FALSE, header=T, 
                         na.strings=c("", "NA"), stringsAsFactors=FALSE, sep=",",quote="")
head(coral_venn)

#1
nrow(subset(coral_venn, SC.d1== 1))
#2
nrow(subset(coral_venn, SC.d3== 1))
#3
nrow(subset(coral_venn, AH.d1== 1))
#4
nrow(subset(coral_venn, AH.d3== 1))

#1 and 2
nrow(subset(coral_venn, SC.d1== 1 & SC.d3== 1))
#1 and 3
nrow(subset(coral_venn, SC.d1== 1 & AH.d1== 1))
#1 and 4
nrow(subset(coral_venn, SC.d1== 1 & AH.d3== 1))

#2 and 3
nrow(subset(coral_venn, SC.d3== 1 & AH.d1== 1))
#2 and 4
nrow(subset(coral_venn, SC.d3== 1 & AH.d3== 1))

#3 and 4
nrow(subset(coral_venn, AH.d1== 1 & AH.d3== 1))


#n123
nrow(subset(coral_venn, SC.d1== 1 & SC.d3== 1& AH.d1== 1))
#n124
nrow(subset(coral_venn, SC.d1== 1 & SC.d3== 1& AH.d3== 1))
#n134
nrow(subset(coral_venn, SC.d1== 1 & AH.d1== 1& AH.d3== 1))
#n234
nrow(subset(coral_venn, SC.d3== 1 & AH.d1== 1& AH.d3== 1))

#n1234
nrow(subset(coral_venn, SC.d1== 1 & SC.d3== 1& AH.d1== 1& AH.d3== 1))



pdf(file="venn_AH_CL_redo_2022.pdf",width=10, height=10)
#grid.newpage()
#draw.triple.venn(area1 = 19716, area2 = 20345, area3 = 24259, n12 = 18066, n23 = 19719, n13 = 19058, 
#                 n123 = 17804, category = c("Plate 1", "Plate 2", "Plate 3"), lty = "blank", 
#                 fill = c("#EB0ACF", "888888", "#80FF00"), fontfamily =
#                   rep("sans", 7), cat.fontfamily = rep("sans", 3))

grid.newpage()
draw.quad.venn(area1 = 764, area2 = 434, area3 = 2330, area4=845, 
                    n12 = 43, n13 = 287, n14=86, 
                    n23 = 95, n24=101, n34=427, 
                    n123 =12, n124=12,n134=57, n234=38,
                    n1234=6, 
                    category = c("SC 1", "SC 3", "AH 1", "AH 3"), lwd=0.1, margin=0.5,
                    fill = c("#008080", "#FFFF00", "#FF6600","#88174B"), fontfamily =rep("sans", 15), cat.fontfamily = rep("sans", 4))
dev.off()

########################################################
## Venn diagram metamorph genes
library(VennDiagram)
met_venn <- read.csv("C:/Users/Sheila's Comp/Documents/GitHub/coral_larval_heatStress_colonization_expression/RNAseq_analysis_R/reyes_metamoph_star.csv",
                       check.names=FALSE, header=T, 
                       na.strings=c("", "NA"), stringsAsFactors=FALSE, sep=",",quote="")
head(met_venn)

#1
nrow(subset(met_venn, Metamorph== 1))
#2
nrow(subset(met_venn, ReyesBermudez== 1))


#1 and 2
nrow(subset(met_venn, Metamorph== 1 & ReyesBermudez== 1))

pdf(file="venn_met3.pdf",width=10, height=10)
grid.newpage()
draw.pairwise.venn(area1 = 8980, area2 = 5907, cross.area = 2585, 
               category = c("This Study", "Reyes-Bermudez 2016"), lwd=0.5, margin=0.5,
               fill = c("hotpink","grey50"), fontfamily =rep("sans", 3), cat.fontfamily = rep("sans", 2))
dev.off()
