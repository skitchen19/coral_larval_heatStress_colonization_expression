# Set the working directory
working.dir <- ("C:/Users/Sheila's Comp/Documents/motifs")

required.libs <- c("gtools","ggplot2","ggridges")
for (lib in required.libs){
  if (!require(lib, character.only = T)){
    install.packages(lib)
    library(lib, character.only = T)
  }
}

library(dplyr)
library(tidyr)
library(stringr)
library(gridExtra)

#### OLIGO
# Create list of significance from module M12 (orange)
# recursive = TRUE lists all files within subdirectories inside the working directory
# We will get files
list.files <- dir(path = working.dir, recursive = TRUE, full.names = TRUE, pattern = "nt_test_vs_ctrl.tab")

# Sort files with embedded numbers so that the numbers are in the correct order
list.files <- mixedsort(list.files,decreasing=T)

# Loop to read files

datalist = list()
#module=12

for (i in 1:length(list.files)){
  # Read the files
  oligo = read.table(list.files[i], header = T, comment.char = ";")

 oligo$interval<-rep(ifelse(grepl("_-1500-200",list.files[i]),"i-1500-200", 
         ifelse(grepl("_-500-200",list.files[i]),"i-500-200","i0-200")),nrow(oligo))
 
 oligo$module<-rep(ifelse(grepl("_M12_",list.files[i]),"M12", 
                            ifelse(grepl("_M13_",list.files[i]),"M13",
                                   ifelse(grepl("_M14_",list.files[i]),"M14",
                                          ifelse(grepl("_M15_",list.files[i]),"M15","M21")))),nrow(oligo))
 # Assign this value as a new row in the empty table
 datalist[[i]] <- oligo
}

big_data = do.call(rbind, datalist)

mat.files <- dir(path = working.dir, recursive = TRUE, full.names = TRUE, pattern = "peak-motifs_oligos-2str-noov_6-8nt_mvkauto_pssm.asmb")
mat.files
datalist2<-list()

for (i in 1:length(mat.files)){
  # Read the files
  mat = read.csv(mat.files[i], header = F, comment.char = ";", sep="\t")
  colnames(mat)<-c("X.seq","rc", "score", "best")
  mat <- mat[ which(mat$best=='best consensus'),]
  
  mat$interval<-rep(ifelse(grepl("_-1500-200",mat.files[i]),"i-1500-200", 
                             ifelse(grepl("_-500-200",mat.files[i]),"i-500-200","i0-200")),nrow(mat))
  
  mat$name1 <- ifelse(1:nrow(mat) == "1" & nrow(mat) == "1",paste("oligos_6-8nt"), paste0("oligos_6-8nt_m",1:nrow(mat))) 
  
  mat$module<-rep(ifelse(grepl("_M12_",mat.files[i]),"M12", 
                           ifelse(grepl("_M13_",mat.files[i]),"M13",
                                  ifelse(grepl("_M14_",mat.files[i]),"M14",
                                         ifelse(grepl("_M15_",mat.files[i]),"M15","M21")))),nrow(mat))
  
  # Assign this value as a new row in the empty table
  datalist2[[i]] <- mat
}

mat_data = do.call(rbind, datalist2)

big_data_mat <- mat_data %>% 
  left_join(select(big_data,"X.seq","occ","occ_P", "occ_E","occ_sig", "interval", "module"),by = c("X.seq", "interval", "module"))

#### DYAD

# We will get files
dyad.files <- dir(path = working.dir, recursive = TRUE, full.names = TRUE, pattern = "peak-motifs_dyads-2str-noov_3nt_sp0-20_bg_monads.tab")
dyad.files
# Sort files with embedded numbers so that the numbers are in the correct order
dyad.files <- mixedsort(dyad.files,decreasing=T)

# Loop to read files
datalist4 = list()

for (i in 1:length(dyad.files)){
  
  # Read the files
  dyad = read.table(dyad.files[i], header = T, comment.char = ";", sep="\t")
  colnames(dyad)<-c("X.seq","id","exp_freq","occ", "exp_occ","occ_P","occ_E", "occ_sig","rank",
                    "ovl_occ","all_occ", "zscore", "occ_var", "ratio", "ov_coef", "remark")
  dyad$interval<-rep(ifelse(grepl("_-1500-200",dyad.files[i]),"i-1500-200", 
                             ifelse(grepl("_-500-200",dyad.files[i]),"i-500-200","i0-200")),nrow(dyad))
  
  dyad$module<-rep(ifelse(grepl("_M12_",dyad.files[i]),"M12", 
                           ifelse(grepl("_M13_",dyad.files[i]),"M13",
                                  ifelse(grepl("_M14_",dyad.files[i]),"M14",
                                         ifelse(grepl("_M15_",dyad.files[i]),"M15","M21")))),nrow(dyad))
  
  # Assign this value as a new row in the empty table
  datalist4[[i]] <- dyad
}

dyad_data = do.call(rbind, datalist4)

mat.dyad <- dir(path = working.dir, recursive = TRUE, full.names = TRUE, pattern = "peak-motifs_dyads-2str-noov_3nt_sp0-20_bg_monads_pssm.asmb")
mat.dyad
datalist5<-list()

for (i in 1:length(mat.dyad)){
  # Read the files
  mat_dyad = read.csv(mat.dyad[i], header = F, comment.char = ";", sep="\t", fill = TRUE, col.names = paste0("V",seq_len(4)))
  colnames(mat_dyad)<-c("X.seq","rc", "score", "best")
  mat_dyad <- mat_dyad[ which(mat_dyad$best=='best consensus'),]
  
  mat_dyad$interval<-rep(ifelse(grepl("_-1500-200",mat.dyad[i]),"i-1500-200", 
                           ifelse(grepl("_-500-200",mat.dyad[i]),"i-500-200","i0-200")),nrow(mat_dyad))
  
  mat_dyad$name1 <- ifelse(1:nrow(mat_dyad) == "1" & nrow(mat_dyad) == "1",paste("dyads_test_vs_ctrl"), paste0("dyads_test_vs_ctrl_m",1:nrow(mat_dyad)))
  
  mat_dyad$module<-rep(ifelse(grepl("_M12_",mat.dyad[i]),"M12", 
                         ifelse(grepl("_M13_",mat.dyad[i]),"M13",
                                ifelse(grepl("_M14_",mat.dyad[i]),"M14",
                                       ifelse(grepl("_M15_",mat.dyad[i]),"M15","M21")))),nrow(mat_dyad))
  
  ct<-str_count(mat_dyad$X.seq, "n")
  for (j in 1:length(ct)){
  mat_dyad$X.seq_n[j]<-ifelse(grepl("n", mat_dyad$X.seq[j]),
                              str_replace(mat_dyad$X.seq[j], paste(replicate(ct[j], "n"), collapse = ""), paste0("n{",str_count(mat_dyad$X.seq[j], "n"),"}")),paste(mat_dyad$X.seq[j]))
  }
  # Assign this value as a new row in the empty table
  datalist5[[i]] <- mat_dyad
}

mat_data_dyad = do.call(rbind, datalist5)

dyad_data_mat <- mat_data_dyad %>% 
  left_join(select(dyad_data,"X.seq","occ","occ_P", "occ_E", "occ_sig","interval", "module"),by = c("X.seq_n"= "X.seq", "interval", "module")) %>%
  select(-X.seq_n)
  
all.data<-rbind(big_data_mat, dyad_data_mat)


###ANNOTATION

annot_oligo<-dir(path = working.dir, recursive = TRUE, full.names = TRUE, pattern = "peak-motifs_motifs_vs_db_footprintDB-metazoa.tab")
annot_oligo
datalist3<-list()

for (i in 1:length(annot_oligo)){
  # Read the files
  annt = read.csv(annot_oligo[i], header = T, comment.char = ";", sep="\t")
  at<-annt %>% arrange(desc(Ncor)) %>% 
    group_by(name1) %>% slice(1)
  
  at$interval<-rep(ifelse(grepl("_-1500-200",annot_oligo[i]),"i-1500-200", 
                          ifelse(grepl("_-500-200",annot_oligo[i]),"i-500-200","i0-200")),nrow(at))
  
  at$module<-rep(ifelse(grepl("_M12_",annot_oligo[i]),"M12", 
                              ifelse(grepl("_M13_",annot_oligo[i]),"M13",
                                     ifelse(grepl("_M14_",annot_oligo[i]),"M14",
                                            ifelse(grepl("_M15_",annot_oligo[i]),"M15","M21")))),nrow(at))
  
  
  # Assign this value as a new row in the empty table
  datalist3[[i]] <- at
}

annt_data = as.data.frame(do.call(rbind, datalist3))

all.data.annt <- all.data %>% 
  left_join(select(annt_data,"name1","name2","cor", "Ncor", "interval", "strand", "module"),by = c("name1", "interval", "module"))

write.table(all.data.annt, "Adig_modules_motifs_combined.txt", sep="\t")

############
## plotting
##

scan.1500 <- read.csv("C:/Users/Sheila's Comp/Documents/motifs/peak-motifs.2020-03-19.185914_2020-03-19.185914_M15_-1500-200/results/sites/peak-motifs_all_motifs_seqcoord.tab", header=T, 
                       dec = ".", sep = "",comment.char = ";")
scan.1500 <- scan.1500[complete.cases(scan.1500),]
scan.1500$interval<-c("i-1500-200")
scan.1500 <- scan.1500[!duplicated(scan.1500), ]

scan.500 <- read.csv("C:/Users/Sheila's Comp/Documents/motifs/peak-motifs.2020-03-20.001106_2020-03-20.001106_M15_-500-200/results/sites/peak-motifs_all_motifs_seqcoord.tab", header=T, 
                      dec = ".", sep = "",comment.char = ";")
scan.500 <- scan.500[complete.cases(scan.500),]
scan.500$interval<-c("i-500-200")
scan.500 <- scan.500[!duplicated(scan.500), ]

scan.200 <- read.csv("C:/Users/Sheila's Comp/Documents/motifs/peak-motifs.2020-03-20.004631_2020-03-20.004631_M15_0-200/results/sites/peak-motifs_all_motifs_seqcoord.tab", header=T, 
                     dec = ".", sep = "",comment.char = ";")
scan.200 <- scan.200[complete.cases(scan.200),]
scan.200$interval<-c("i0-200")
scan.200 <- scan.200[!duplicated(scan.200), ]

# Join the 2 files vertically
scan.data <- rbind(scan.1500,scan.500,scan.200)
scan.data <- na.omit(scan.data)
scan.data$name<- paste(scan.data$interval,scan.data$ft_name)
dim(scan.data)
head(scan.data)
scan.data <- scan.data[!duplicated(scan.data), ]

scan.filter<-filter(scan.data, (interval=="i-1500-200" & ft_name=="oligos_6-7nt_m1") |(interval=="i-500-200" & ft_name=="oligos_6-8nt_m1") | (interval=="i-500-200" & ft_name=="oligos_6-8nt_m2"))
str(scan.filter)

# Ridgeline plots
plot_M15 <- ggplot(scan.filter,aes(as.numeric(paste(start)), name, fill = name,color = name))+
  geom_density_ridges_gradient(scale = 4, show.legend = T,rel_min_height = 0.001, size=1.2) +
  theme_ridges() + theme(legend.position = "none") +
  labs(x = "",y = "") + ggtitle("") + theme(plot.title = element_text(hjust = 0.4)) + 
  theme(axis.text.x = element_text( angle = 90,  hjust = 1, size = 8), 
        axis.text.y = element_text( hjust = 1, size = 10)) +
  scale_color_manual(values = c("white", "white",  "white"))+ 
  scale_fill_manual(values = c("#CD370099", "#CD370099", "#CD370099"))+
  scale_y_discrete(expand = c(0.01, 0)) +
  theme(strip.text = element_text(size=12),strip.background = element_rect(fill="gray", colour="black",size=4)) +
  scale_x_continuous(limit = c(-1500, +200),breaks = c(-1500,-500,0,+200))
plot_M15 


scan.1500 <- read.csv("C:/Users/Sheila's Comp/Documents/motifs/peak-motifs.2020-03-19.220336_2020-03-19.220336_M14_-1500-200/results/sites/peak-motifs_all_motifs_seqcoord.tab", header=T, 
                      dec = ".", sep = "",comment.char = ";")
scan.1500 <- scan.1500[complete.cases(scan.1500),]
scan.1500$interval<-c("i-1500-200")
scan.1500 <- scan.1500[!duplicated(scan.1500), ]

scan.500 <- read.csv("C:/Users/Sheila's Comp/Documents/motifs/peak-motifs.2020-03-20.001935_2020-03-20.001935_M14_-500-200/results/sites/peak-motifs_all_motifs_seqcoord.tab", header=T, 
                     dec = ".", sep = "",comment.char = ";")
scan.500 <- scan.500[complete.cases(scan.500),]
scan.500$interval<-c("i-500-200")
scan.500 <- scan.500[!duplicated(scan.500), ]

# Join the 2 files vertically
scan.data <- rbind(scan.1500,scan.500)
scan.data <- na.omit(scan.data)
scan.data$name<- paste(scan.data$interval,scan.data$ft_name)
dim(scan.data)
head(scan.data)
scan.data <- scan.data[!duplicated(scan.data), ]

scan.filter<-filter(scan.data, (interval=="i-1500-200" & ft_name=="oligos_6-8nt_m1") |(interval=="i-1500-200" & ft_name=="oligos_6-8nt_m2") | (interval=="i-1500-200" & ft_name=="dyads_test_vs_ctrl_m1")|(interval=="i-1500-200" & ft_name=="dyads_test_vs_ctrl_m2") )
str(scan.filter)

# Ridgeline plots
plot_M14 <- ggplot(scan.filter,aes(as.numeric(paste(start)), name, fill = name,color = name))+
  geom_density_ridges_gradient(scale = 4, show.legend = T,rel_min_height = 0.001, size=1.2) +
  theme_ridges() + theme(legend.position = "none") +
  labs(x = "",y = "") + ggtitle("") + theme(plot.title = element_text(hjust = 0.4)) + 
  theme( axis.text.x = element_text( angle = 90,  hjust = 1, size = 8),
         axis.text.y = element_text( hjust = 1, size = 10)) +
  scale_color_manual(values = c("white", "white",  "white",  "white"))+ 
  scale_fill_manual(values = c("#B0306060", "#B0306060", "#B0306060","#B0306060"))+
  scale_y_discrete(expand = c(0.01, 0)) +
  theme(strip.text = element_text(size=12),strip.background = element_rect(fill="gray", colour="black",size=4)) +
  scale_x_continuous(limit = c(-1500, +200),breaks = c(-1500,-500,0,+200))
plot_M14


scan.1500 <- read.csv("C:/Users/Sheila's Comp/Documents/motifs/peak-motifs.2020-03-19.223856_2020-03-19.223856_M13_-1500-200/results/sites/peak-motifs_all_motifs_seqcoord.tab", header=T, 
                      dec = ".", sep = "",comment.char = ";")
scan.1500 <- scan.1500[complete.cases(scan.1500),]
scan.1500$interval<-c("i-1500-200")
scan.1500 <- scan.1500[!duplicated(scan.1500), ]

scan.500 <- read.csv("C:/Users/Sheila's Comp/Documents/motifs/peak-motifs.2020-03-20.003245_2020-03-20.003245_M13_-500-200/results/sites/peak-motifs_all_motifs_seqcoord.tab", header=T, 
                     dec = ".", sep = "",comment.char = ";")
scan.500 <- scan.500[complete.cases(scan.500),]
scan.500$interval<-c("i-500-200")
scan.500 <- scan.500[!duplicated(scan.500), ]

scan.200 <- read.csv("C:/Users/Sheila's Comp/Documents/motifs/peak-motifs.2020-03-20.005453_2020-03-20.005453_M13_0-200/results/sites/peak-motifs_all_motifs_seqcoord.tab", header=T, 
                     dec = ".", sep = "",comment.char = ";")
scan.200 <- scan.200[complete.cases(scan.200),]
scan.200$interval<-c("i0-200")
scan.200 <- scan.200[!duplicated(scan.200), ]

# Join the 2 files vertically
scan.data <- rbind(scan.1500,scan.500,scan.200)
scan.data <- na.omit(scan.data)
scan.data$name<- paste(scan.data$interval,scan.data$ft_name)
dim(scan.data)
head(scan.data)
scan.data <- scan.data[!duplicated(scan.data), ]

scan.filter<-filter(scan.data, (interval=="i-1500-200" & ft_name=="dyads_test_vs_ctrl_m1") |
                      (interval=="i-1500-200" & ft_name=="dyads_test_vs_ctrl_m2") | 
                      (interval=="i-1500-200" & ft_name=="dyads_test_vs_ctrl_m3") |
                      (interval=="i-1500-200" & ft_name=="dyads_test_vs_ctrl_m4") |
                      (interval=="i-1500-200" & ft_name=="dyads_test_vs_ctrl_m5") |
                      (interval=="i-1500-200" & ft_name=="dyads_test_vs_ctrl_m6") |
                      (interval=="i-1500-200" & ft_name=="dyads_test_vs_ctrl_m7") |
                      (interval=="i-1500-200" & ft_name=="dyads_test_vs_ctrl_m8") |
                      (interval=="i-500-200" & ft_name=="oligos_6-8nt_m1") |
                      (interval=="i0-200" & ft_name=="dyads_test_vs_ctrl_m1") |
                      (interval=="i0-200" & ft_name=="dyads_test_vs_ctrl_m2"))
str(scan.filter)

scan.filter$name<- with(scan.filter, reorder(as.factor(name), interval))

# Ridgeline plots
plot_M13 <- ggplot(scan.filter,aes(as.numeric(paste(start)), name, fill = name,color = name))+
  geom_density_ridges_gradient(scale = 4, show.legend = T,rel_min_height = 0.001, size=1.2, color="white") +
  theme_ridges() + theme(legend.position = "none") +
  labs(x = "",y = "") + ggtitle("") + theme(plot.title = element_text(hjust = 0.4)) + 
  theme( axis.text.x = element_text( angle = 90,  hjust = 1, size = 8),
         axis.text.y = element_text( hjust = 1, size = 10)) +
  scale_fill_manual(values = c("#F0808099", "#F0808099", "#F0808099","#F0808099", "#F0808099", "#F0808099","#F0808099", "#F0808099", "#F0808099","#F0808099", "#F0808099"))+
  scale_y_discrete(expand = c(0.01, 0)) +
  theme(strip.text = element_text(size=12),strip.background = element_rect(fill="black", colour="black",size=4)) +
  scale_x_continuous(limit = c(-1500, +200),breaks = c(-1500,-500,0,+200))
plot_M13

scan.1500 <- read.csv("C:/Users/Sheila's Comp/Documents/motifs/peak-motifs.2020-03-19.221804_2020-03-19.221804_M12_-1500-200/results/sites/peak-motifs_all_motifs_seqcoord.tab", header=T, 
                      dec = ".", sep = "",comment.char = ";")
scan.1500 <- scan.1500[complete.cases(scan.1500),]
scan.1500$interval<-c("i-1500-200")
scan.1500 <- scan.1500[!duplicated(scan.1500), ]

scan.500 <- read.csv("C:/Users/Sheila's Comp/Documents/motifs/peak-motifs.2020-03-20.001428_2020-03-20.001428_M12_-500-200/results/sites/peak-motifs_all_motifs_seqcoord.tab", header=T, 
                     dec = ".", sep = "",comment.char = ";")
scan.500 <- scan.500[complete.cases(scan.500),]
scan.500$interval<-c("i-500-200")
scan.500 <- scan.500[!duplicated(scan.500), ]

scan.200 <- read.csv("C:/Users/Sheila's Comp/Documents/motifs/peak-motifs.2020-03-20.004933_2020-03-20.004933_M12_0-200/results/sites/peak-motifs_all_motifs_seqcoord.tab", header=T, 
                     dec = ".", sep = "",comment.char = ";")
scan.200 <- scan.200[complete.cases(scan.200),]
scan.200$interval<-c("i0-200")
scan.200 <- scan.200[!duplicated(scan.200), ]

# Join the 2 files vertically
scan.data <- rbind(scan.1500,scan.500,scan.200)
scan.data <- na.omit(scan.data)
scan.data$name<- paste(scan.data$interval,scan.data$ft_name)
dim(scan.data)
head(scan.data)
scan.data <- scan.data[!duplicated(scan.data), ]

scan.filter<-filter(scan.data, (interval=="i-1500-200" & ft_name=="dyads_test_vs_ctrl_m1") |
                      (interval=="i-1500-200" & ft_name=="dyads_test_vs_ctrl_m2") | 
                      (interval=="i-1500-200" & ft_name=="dyads_test_vs_ctrl_m3") |
                      (interval=="i-1500-200" & ft_name=="dyads_test_vs_ctrl_m4") |
                      (interval=="i-1500-200" & ft_name=="dyads_test_vs_ctrl_m5") |
                      (interval=="i-1500-200" & ft_name=="oligos_6-8nt_m1") |
                      (interval=="i-1500-200" & ft_name=="oligos_6-8nt_m2") |
                      (interval=="i-500-200" & ft_name=="dyads_test_vs_ctrl_m1") |
                      (interval=="i-500-200" & ft_name=="dyads_test_vs_ctrl_m2") |
                      (interval=="i-500-200" & ft_name=="dyads_test_vs_ctrl_m3") |
                      (interval=="i-500-200" & ft_name=="dyads_test_vs_ctrl_m4"))
str(scan.filter)

scan.filter$name<- with(scan.filter, reorder(as.factor(name), interval))

# Ridgeline plots
plot_M12 <- ggplot(scan.filter,aes(as.numeric(paste(start)), name, fill = name,color = name))+
  geom_density_ridges_gradient(scale = 4, show.legend = T,rel_min_height = 0.001, size=1.2, color="white") +
  theme_ridges() + theme(legend.position = "none") +
  labs(x = "",y = "") + ggtitle("") + theme(plot.title = element_text(hjust = 0.4)) + 
  theme( axis.text.x = element_text( angle = 90,  hjust = 1, size = 8),
         axis.text.y = element_text( hjust = 1, size = 10)) +
  scale_fill_manual(values = c("#FFA50099", "#FFA50099", "#FFA50099","#FFA50099", "#FFA50099", 
                               "#FFA50099","#FFA50099", "#FFA50099", "#FFA50099","#FFA50099", "#FFA50099"))+
  scale_y_discrete(expand = c(0.01, 0)) +
  theme(strip.text = element_text(size=12),strip.background = element_rect(fill="black", colour="black",size=4)) +
  scale_x_continuous(limit = c(-1500, +200),breaks = c(-1500,-500,0,+200))
plot_M12

library("cowplot")

ggdraw() +
  draw_plot(plot_M12, x = 0, y = .28, width = .4, height = .75)+
  draw_plot(plot_M13, x = .5, y = .28, width = .4, height = .75) +
  draw_plot(plot_M14, x = 0, y = -.04, width = 0.4, height = 0.4) +
  draw_plot(plot_M15, x = 0.5, y = -.04, width = 0.4, height = 0.35) +
  draw_plot_label(label = c("M12", "M13", "M14", "M15"), size = 15,
                  x = c(0, 0.5, 0,0.5), y = c(1, 1, 0.3,0.3))
