
###########
# GO barplot for WGCNA modules
setwd("C:/Users/Sheila's Comp/Dropbox/RNASeq2014/NCC/Figures/Supplmental Figures/")
b<-read.table("barplot_go_input.txt",sep="\t",quote="",
              check.names=FALSE, header=T, fill=TRUE)
head(b)

m12<-subset(b, Module %in% c("M12"))
m12$description <- factor(m12$description, levels = m12$description[order(m12$category)])

g <- ggplot(data=m12, aes(x=description,y=plog10))+
  labs(y = '-log (p-value)')+
  geom_bar(stat = "identity", fill="orange")+
  scale_y_continuous(limits = c(0,5), expand = c(0, 0)) +
  coord_flip() +
  theme_classic()
g

m13<-subset(b, Module %in% c("M13"))
m13$description <- factor(m13$description, levels = m13$description[order(m13$category)])

g2 <- ggplot(data=m13, aes(x=description,y=plog10))+
  labs(y = '-log (p-value)')+
  geom_bar(stat = "identity", fill="lightcoral")+
  scale_y_continuous(limits = c(0,20), expand = c(0, 0)) +
  coord_flip() +
  theme_classic()
g2


m14<-subset(b, Module %in% c("M14"))
m14$description <- factor(m14$description, levels = m14$description[order(m14$category)])

g3 <- ggplot(data=m14, aes(x=description,y=plog10))+
  labs(y = '-log (p-value)')+
  geom_bar(stat = "identity", fill="maroon")+
  scale_y_continuous(limits = c(0,5), expand = c(0, 0)) +
  coord_flip() +
  theme_classic()
g3

m15<-subset(b, Module %in% c("M15"))
m15$description <- factor(m15$description, levels = m15$description[order(m15$category)])

g4 <- ggplot(data=m15, aes(x=description,y=plog10))+
  labs(y = '-log (p-value)')+
  geom_bar(stat = "identity", fill="orangered3")+
  scale_y_continuous(limits = c(0,5), expand = c(0, 0)) +
  coord_flip() +
  theme_classic()
g4

figure <- ggarrange(g2,g3,g4, g,ncol=2,align="hv",
                    nrow = 2, labels=c("M13","M14","M15","M12"),
                    heights=c(5,2))
figure
