##clustering##
deg<-read.csv("./Interactions_graph/in_combo.csv", header=T,row.names=1)
#deg<-read.csv("./CL_time/in_combo_CL.csv", header=T,row.names=1)
#deg<-read.csv("./AH_time/in_combo_AH.csv", header=T,row.names=1)

head(deg)

d <- dist(as.matrix(deg[,1:4]), method = "euclidean")
hc<-hclust(d)
p<-plot(hc, hang = -0.07, cex = 0.3)
abline(h=1.5, col="red", untf = FALSE)
c<-cutree(hc, h = 1)

dend <- as.dendrogram(hc)
# Let's add some color:
colors_to_use <- as.numeric(deg$Cluster)

# But sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend) <- colors_to_use

dend.list<-as.character(deg$Cluster)

#now we need to change the order of labels from the order in which they appear in the original data to match the order in which they appear on the dendrogram we just made.
labels(dend)<- dend.list[order.dendrogram(dend)]

#Here is the final plot!
plot(dend, main = "Ordered color for every Cluster",cex = 0.3)

write.csv(c, file = "hclust_combo_AH_time.csv")
##########################
#interaction graphs

require(ggplot2)
require(ggvis)
require(gridExtra)
library(tidyverse)
library(dplyr)


library(reshape)

head(deg)
deg$id<-rownames(deg)
deg2 <- melt(deg, id=c("id","Cluster2","Cluster","Day"))
head(deg2)

deg3<-deg2 %>% separate(variable, c("symbiont", "temp"), sep="_")
deg3$Cluster<-factor(deg3$Cluster)
deg3$temp<-factor(deg3$temp)
deg3$symbiont<-factor(deg3$symbiont)

par( mfrow = c( 1,2 ),bty="l", mar=c(2,1,4,2.5),oma = c(3, 5, 1, 0.5))
##Cluster1##
#means=by(deg3$value[which(deg3$Cluster_1=="1")],list(deg3$temp[which(deg3$Cluster_1=="1")],deg3$day[which(deg3$Cluster_1=="1")]),mean)
#se=function(x) sqrt(var(x)/length(x)) 
#ses=by(deg3$value[which(deg3$Cluster_1=="1")],list(deg3$temp[which(deg3$Cluster_1=="1")],deg3$day[which(deg3$Cluster_1=="1")]),se) 


#interaction.plot(deg3$day[which(deg3$Cluster_1.5==1)],deg3$temp[which(deg3$Cluster_1.5==1)],
#                 deg3$value[which(deg3$Cluster_1.5==1)], type="b",pch=c(21,19), cex.lab=1.5,cex=1.5, cex.axis=1.5,
#                 ylim=c(-1,1),yaxp  = c(-1, 1, 2),legend=F, font.lab=2, bg=c("white","white"),xlab="",ylab="",
#                 main="C1: 62 genes",lty=c(3,1),col=c("black"),lwd=2,cex.main=1.5, las=1)

#mtext("D1:62 genes", col="black", cex=.8)

#lines(c(1,1),c(means[1]-ses[1],means[1]+ses[1]),lwd=2) 
#lines(c(1,1),c(means[2]-ses[2],means[2]+ses[2]),lwd=2) 
#lines(c(2,2),c(means[3]-ses[3],means[3]+ses[3]),lwd=2) 
#lines(c(2,2),c(means[4]-ses[4],means[4]+ses[4]),lwd=2)


#int<-read.csv("./Interactions_graph/all_interactions.csv", head=T)
#int<-read.csv("./Interactions_graph/CL_time/all_interactions_C_Time.csv", head=T)
#head(int)
#int$Cluster<-factor(int$Cluster)

#int$temp<-relevel(int$temp, ref="L")

#par( mfrow = c( 3, 3 ),bty="l", mar=c(2,1,4,2.5),oma = c(3, 5, 1, 0.5)))

# Summarizing data
data2 <- deg3 %>% 
    group_by(Cluster,temp,symbiont) %>% 
    summarise(exp_mean = mean(value),
              exp_se = psych::describe(value)$se)


p1 <- ggplot(data = deg3, aes(x = temp, y = value)) +
    #geom_line(data=data2,aes(x = temp, y = exp_mean, group = symbiont, linetype=symbiont), size=1.3)+
    theme_classic()+ ylim(-1.7,1.7)+
    ylab("Expression")+ xlab("Day")+ 
    geom_point(aes(fill=symbiont),color="black",size=4,pch=21, alpha=0.5) + 
    theme(legend.position="top")+
    facet_wrap(. ~ Cluster)+
    scale_fill_manual(values=c("white", "grey20"))+
    scale_color_manual(values=c("black", "black"))+
    scale_linetype_manual(values=c( "dotdash","solid"))
p1

p1 <- ggplot(data = deg2, aes(x = variable, y = value)) +
    #geom_line(data=data2,aes(x = temp, y = exp_mean, group = symbiont, linetype=symbiont), size=1.3)+
    theme_classic()+ ylim(-1.7,1.7)+
    ylab("Expression")+ xlab("")+ 
    geom_boxplot()+
    geom_point(aes(fill=variable),color="black",size=4,pch=21, alpha=0.5) + 
    theme(legend.position="top")+
    facet_wrap(. ~ Cluster2)+
    scale_fill_manual(values=c("#EDEA37", "#008080","#FF6600", "#88174B"))+
    scale_color_manual(values=c("black", "black","black", "black"))+
    scale_linetype_manual(values=c( "dotdash","solid"))
p1


##Cluster1##
means=by(int$exp[which(int$Cluster=="1")],list(int$sym[which(int$Cluster=="1")],int$temp[which(int$Cluster=="1")]),mean)
se=function(x) sqrt(var(x)/length(x)) 
ses=by(int$exp[which(int$Cluster=="1")],list(int$sym[which(int$Cluster=="1")],int$temp[which(int$Cluster=="1")]),se) 

interaction.plot(int$temp[which(int$Cluster==1)], int$sym[which(int$Cluster==1)],
int$exp[which(int$Cluster==1)], type="b",pch=c(21,19), cex.lab=1.5,cex=1.5, cex.axis=1.5,
ylim=c(-1,1),yaxp  = c(-1, 1, 2),legend=F, font.lab=2, bg=c("white","white"),xlab="",ylab="",
main="C1: 62 genes",lty=c(3,1),col=c("black"),lwd=2,cex.main=1.5, las=1)

#mtext("D1:62 genes", col="black", cex=.8)

lines(c(1,1),c(means[1]-ses[1],means[1]+ses[1]),lwd=2) 
lines(c(1,1),c(means[2]-ses[2],means[2]+ses[2]),lwd=2) 
lines(c(2,2),c(means[3]-ses[3],means[3]+ses[3]),lwd=2) 
lines(c(2,2),c(means[4]-ses[4],means[4]+ses[4]),lwd=2)


##Cluster2##
means=by(int$exp[which(int$Cluster==2)],list(int$sym[which(int$Cluster==2)],int$temp[which(int$Cluster==2)]),mean)
se=function(x) sqrt(var(x)/length(x)) 
ses=by(int$exp[which(int$Cluster=="2")],list(int$sym[which(int$Cluster=="2")],int$temp[which(int$Cluster=="2")]),se) 

interaction.plot(int$temp[which(int$Cluster==2)], int$sym[which(int$Cluster==2)],
int$exp[which(int$Cluster==2)], type="b",pch=c(21,19), cex.lab=1.5,cex=1.5, cex.axis=1.5,
ylim=c(-1,1),yaxp  = c(-1, 1, 2),legend=F, font.lab=2, bg=c("white","white"),xlab="",ylab="", 
main="C2: 73 genes",lty=c(3,1),col=c("black"),lwd=2,cex.main=1.5, las=1)

#mtext("D1:73 genes", col="black", cex=.8)

lines(c(1,1),c(means[1]-ses[1],means[1]+ses[1]),lwd=2) 
lines(c(1,1),c(means[2]-ses[2],means[2]+ses[2]),lwd=2) 
lines(c(2,2),c(means[3]-ses[3],means[3]+ses[3]),lwd=2) 
lines(c(2,2),c(means[4]-ses[4],means[4]+ses[4]),lwd=2)


##Cluster4##
means=by(int$exp[which(int$Cluster==4)],list(int$sym[which(int$Cluster==4)],int$temp[which(int$Cluster==4)]),mean)
se=function(x) sqrt(var(x)/length(x)) 
ses=by(int$exp[which(int$Cluster==4)],list(int$sym[which(int$Cluster==4)],int$temp[which(int$Cluster==4)]),se) 

interaction.plot(int$temp[which(int$Cluster==4)], int$sym[which(int$Cluster==4)],
int$exp[which(int$Cluster==4)], type="b",pch=c(21,19), cex.lab=2,cex=1.5, cex.axis=1.5,
ylim=c(-1,1),yaxp  = c(-1, 1, 2),legend=F, font.lab=2, bg=c("white","white"),xlab="",ylab="Expression", 
main="C3: 14 genes",lty=c(3,1),col=c("black"),lwd=2,cex.main=1.5, las=1)

#mtext("D1:14 genes", col="red", cex=.8)

lines(c(1,1),c(means[1]-ses[1],means[1]+ses[1]),lwd=2) 
lines(c(1,1),c(means[2]-ses[2],means[2]+ses[2]),lwd=2) 
lines(c(2,2),c(means[3]-ses[3],means[3]+ses[3]),lwd=2) 
lines(c(2,2),c(means[4]-ses[4],means[4]+ses[4]),lwd=2)


##Cluster5##
means=by(int$exp[which(int$Cluster==5)],list(int$sym[which(int$Cluster==5)],int$temp[which(int$Cluster==5)]),mean)
se=function(x) sqrt(var(x)/length(x)) 
ses=by(int$exp[which(int$Cluster==5)],list(int$sym[which(int$Cluster==5)],int$temp[which(int$Cluster==5)]),se) 

interaction.plot(int$temp[which(int$Cluster==5)], int$sym[which(int$Cluster==5)],
int$exp[which(int$Cluster==5)], type="b",pch=c(21,19), cex.lab=1.5,cex=1.5, cex.axis=1.5,
ylim=c(-1,1),yaxp  = c(-1, 1, 2),legend=F, font.lab=2, bg=c("white","white"),xlab="",ylab="", 
main="C4: 50 genes",lty=c(3,1),col=c("black"),lwd=2,cex.main=1.5, las=1)

#mtext("D1:50 genes", col="red", cex=.8)

lines(c(1,1),c(means[1]-ses[1],means[1]+ses[1]),lwd=2) 
lines(c(1,1),c(means[2]-ses[2],means[2]+ses[2]),lwd=2) 
lines(c(2,2),c(means[3]-ses[3],means[3]+ses[3]),lwd=2) 
lines(c(2,2),c(means[4]-ses[4],means[4]+ses[4]),lwd=2)

mtext("Expression", side=2,line=3.5, cex=1.7, font=2)

##Cluster6##
means=by(int$exp[which(int$Cluster==6)],list(int$sym[which(int$Cluster==6)],int$temp[which(int$Cluster==6)]),mean)
se=function(x) sqrt(var(x)/length(x)) 
ses=by(int$exp[which(int$Cluster==6)],list(int$sym[which(int$Cluster==6)],int$temp[which(int$Cluster==6)]),se) 

interaction.plot(int$temp[which(int$Cluster==6)], int$sym[which(int$Cluster==6)],
int$exp[which(int$Cluster==6)], type="b",pch=c(21,19), cex.lab=1.5,cex=1.5, cex.axis=1.5,
ylim=c(-2,2),yaxp  = c(-2, 2, 2),legend=F, font.lab=2, bg=c("white","white"),xlab="",ylab="", 
main="C5: 1 gene",lty=c(3,1),col=c("black"),lwd=2,cex.main=1.5, las=1)

#mtext("D1:1 genes", col="red", cex=.8)

lines(c(1,1),c(means[1]-ses[1],means[1]+ses[1]),lwd=2) 
lines(c(1,1),c(means[2]-ses[2],means[2]+ses[2]),lwd=2) 
lines(c(2,2),c(means[3]-ses[3],means[3]+ses[3]),lwd=2) 
lines(c(2,2),c(means[4]-ses[4],means[4]+ses[4]),lwd=2)

##Cluster8##
means=by(int$exp[which(int$Cluster==8)],list(int$sym[which(int$Cluster==8)],int$temp[which(int$Cluster==8)]),mean)
se=function(x) sqrt(var(x)/length(x)) 
ses=by(int$exp[which(int$Cluster==8)],list(int$sym[which(int$Cluster==8)],int$temp[which(int$Cluster==8)]),se) 

interaction.plot(int$temp[which(int$Cluster==8)], int$sym[which(int$Cluster==8)],
int$exp[which(int$Cluster==8)], type="b",pch=c(21,19), cex.lab=2,cex=1.5, cex.axis=1.5,
ylim=c(-2,2),yaxp  = c(-2, 2, 2),legend=F, font.lab=2, bg=c("white","white"),xlab="Temperature",ylab="", 
main="C6: 1 gene",lty=c(3,1),col=c("black"),lwd=2,cex.main=1.5, las=1)

#mtext("D1:1 genes", col="black", cex=.8)

lines(c(1,1),c(means[1]-ses[1],means[1]+ses[1]),lwd=2) 
lines(c(1,1),c(means[2]-ses[2],means[2]+ses[2]),lwd=2) 
lines(c(2,2),c(means[3]-ses[3],means[3]+ses[3]),lwd=2) 
lines(c(2,2),c(means[4]-ses[4],means[4]+ses[4]),lwd=2)


##Cluster7##
means=by(int$exp[which(int$Cluster==7)],list(int$sym[which(int$Cluster==7)],int$temp[which(int$Cluster==7)]),mean)
se=function(x) sqrt(var(x)/length(x)) 
ses=by(int$exp[which(int$Cluster==7)],list(int$sym[which(int$Cluster==7)],int$temp[which(int$Cluster==7)]),se) 

interaction.plot(int$temp[which(int$Cluster==7)], int$sym[which(int$Cluster==7)],
int$exp[which(int$Cluster==7)], type="b",pch=c(21,19), cex.lab=1.5,cex=1.5, cex.axis=1.5,
ylim=c(-1,1),yaxp  = c(-1, 1, 2),legend=F, font.lab=2, bg=c("white","white"),xlab="",ylab="", 
main="C7: 15 genes",lty=c(3,1),col=c("grey40"),lwd=2,cex.main=1.5, las=1, col.main="grey40")

#mtext("D3:15 genes", col="blue", cex=.8)

lines(c(1,1),c(means[1]-ses[1],means[1]+ses[1]),lwd=2,col=c("grey40")) 
lines(c(1,1),c(means[2]-ses[2],means[2]+ses[2]),lwd=2,col=c("grey40")) 
lines(c(2,2),c(means[3]-ses[3],means[3]+ses[3]),lwd=2,col=c("grey40")) 
lines(c(2,2),c(means[4]-ses[4],means[4]+ses[4]),lwd=2,col=c("grey40"))


##Cluster3##
means=by(int$exp[which(int$Cluster==3)],list(int$sym[which(int$Cluster==3)],int$temp[which(int$Cluster==3)]),mean)
se=function(x) sqrt(var(x)/length(x)) 
ses=by(int$exp[which(int$Cluster==3)],list(int$sym[which(int$Cluster==3)],int$temp[which(int$Cluster==3)]),se) 

interaction.plot(int$temp[which(int$Cluster==3)], int$sym[which(int$Cluster==3)],
int$exp[which(int$Cluster==3)], type="b",pch=c(21,19), cex.lab=1.5,cex=1.5, cex.axis=1.5,
ylim=c(-1,1),yaxp  = c(-1, 1, 2),legend=F, font.lab=2, bg=c("white","white"),xlab="",ylab="", 
main="C8: 16 genes",lty=c(3,1),col=c("grey40"),lwd=2,cex.main=1.5, las=1,col.main="grey40")

#mtext("D3:16 genes", col="blue", cex=.8)

lines(c(1,1),c(means[1]-ses[1],means[1]+ses[1]),lwd=2,col=c("grey40")) 
lines(c(1,1),c(means[2]-ses[2],means[2]+ses[2]),lwd=2,col=c("grey40")) 
lines(c(2,2),c(means[3]-ses[3],means[3]+ses[3]),lwd=2,col=c("grey40")) 
lines(c(2,2),c(means[4]-ses[4],means[4]+ses[4]),lwd=2,col=c("grey40"))

mtext("Temperature", side=1,line=3.3, cex=1.7, font=2)

##Cluster9##
#means=by(int$exp[which(int$Cluster==9)],list(int$sym[which(int$Cluster==9)],int$temp[which(int$Cluster==9)]),mean)
#se=function(x) sqrt(var(x)/length(x)) 
#ses=by(int$exp[which(int$Cluster==9)],list(int$sym[which(int$Cluster==9)],int$temp[which(int$Cluster==9)]),se) 

#interaction.plot(int$temp[which(int$Cluster==9)], int$sym[which(int$Cluster==9)],
#int$exp[which(int$Cluster==9)], type="b",pch=c(21,19), cex.lab=1.5,cex=1.5, cex.axis=1.5,
#ylim=c(-1,1),legend=F, font.lab=2, bg=c("white","white"),xlab="",ylab="", 
#main="Cluster 9",lty=c(3,1),col=c("black"),lwd=2,cex.main=1.5, las=1)

#mtext("1 genes", col="grey30", cex=.8)

#lines(c(1,1),c(means[1]-ses[1],means[1]+ses[1]),lwd=2) 
#lines(c(1,1),c(means[2]-ses[2],means[2]+ses[2]),lwd=2) 
#lines(c(2,2),c(means[3]-ses[3],means[3]+ses[3]),lwd=2) 
#lines(c(2,2),c(means[4]-ses[4],means[4]+ses[4]),lwd=2)

#########################
## Venn diagram CL and AH
library(VennDiagram)
coral_venn <- read.csv("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Figure 3/fold_change_d1_d3_CL_AH.csv",
                         check.names=FALSE, header=T, 
                         na.strings=c("", "NA"), stringsAsFactors=FALSE, sep=",",quote="")
head(coral_venn)

#1
nrow(subset(coral_venn, CL.1d== 1))
#2
nrow(subset(coral_venn, CL.3d== 1))
#3
nrow(subset(coral_venn, AH.1d== 1))
#4
nrow(subset(coral_venn, AH.3d== 1))

#1 and 2
nrow(subset(coral_venn, CL.1d== 1 & CL.3d== 1))
#1 and 3
nrow(subset(coral_venn, CL.1d== 1 & AH.1d== 1))
#1 and 4
nrow(subset(coral_venn, CL.1d== 1 & AH.3d== 1))

#2 and 3
nrow(subset(coral_venn, CL.3d== 1 & AH.1d== 1))
#2 and 4
nrow(subset(coral_venn, CL.3d== 1 & AH.3d== 1))

#3 and 4
nrow(subset(coral_venn, AH.1d== 1 & AH.3d== 1))


#n123
nrow(subset(coral_venn, CL.1d== 1 & CL.3d== 1& AH.1d== 1))
#n124
nrow(subset(coral_venn, CL.1d== 1 & CL.3d== 1& AH.3d== 1))
#n134
nrow(subset(coral_venn, CL.1d== 1 & AH.1d== 1& AH.3d== 1))
#n234
nrow(subset(coral_venn, CL.3d== 1 & AH.1d== 1& AH.3d== 1))

#n1234
nrow(subset(coral_venn, CL.1d== 1 & CL.3d== 1& AH.1d== 1& AH.3d== 1))



pdf(file="venn_AH_CL_redo.pdf",width=10, height=10)
#grid.newpage()
#draw.triple.venn(area1 = 19716, area2 = 20345, area3 = 24259, n12 = 18066, n23 = 19719, n13 = 19058, 
#                 n123 = 17804, category = c("Plate 1", "Plate 2", "Plate 3"), lty = "blank", 
#                 fill = c("#EB0ACF", "888888", "#80FF00"), fontfamily =
#                   rep("sans", 7), cat.fontfamily = rep("sans", 3))

grid.newpage()
draw.quad.venn(area1 = 504, area2 = 415, area3 = 2070, area4=774, 
                    n12 = 28, n13 = 178, n14=55, 
                    n23 = 103, n24=104, n34=372, 
                    n123 =11, n124=7,n134=35, n234=35,
                    n1234=4, 
                    category = c("CL 1", "CL 3", "AH 1", "AH 3"), lwd=0.1, margin=0.5,
                    fill = c("#008080", "#FFFF00", "#FF6600","#88174B"), fontfamily =rep("sans", 15), cat.fontfamily = rep("sans", 4))
dev.off()

########################################################
## Venn diagram metamorph genes
library(VennDiagram)
met_venn <- read.csv("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Supplmental Figures/reyes_metamoph.csv",
                       check.names=FALSE, header=T, 
                       na.strings=c("", "NA"), stringsAsFactors=FALSE, sep=",",quote="")
head(met_venn)

#1
nrow(subset(met_venn, Metamorph== 1))
#2
nrow(subset(met_venn, ReyesBermudez== 1))


#1 and 2
nrow(subset(met_venn, Metamorph== 1 & ReyesBermudez== 1))

pdf(file="venn_met2.pdf",width=10, height=10)
grid.newpage()
draw.pairwise.venn(area1 = 8936, area2 = 5894, cross.area = 2599, 
               category = c("This Study", "Reyes-Bermudez 2016"), lwd=0.5, margin=0.5,
               fill = c("hotpink","grey50"), fontfamily =rep("sans", 3), cat.fontfamily = rep("sans", 2))
dev.off()
