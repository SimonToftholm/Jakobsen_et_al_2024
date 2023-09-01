##############################################
################## FIGURE 6 ##################
##############################################

################# Packages #################
library(tidyverse)
library(PupillometryR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(reshape)
library(clusterProfiler)


################# Figure 6A ################# 
# Count KDM3A binding in the four different subset of genomic regions using homer and laod into R

# annotatePeaks.pl KDM3A_only_enhancer.bed hg38 -d tagdig/KDM3A_D.TD/ tagdig/KDM3A_K.TD/ -size 1000 -cpu 12 > KDM3A_only_enhancer_count.txt
KDM3A_enh <- read.delim("KDM3A_only_enhancer_count.txt")

# annotatePeaks.pl KDM3A_only_promoter.bed hg38 -d tagdig/KDM3A_D.TD/ tagdig/KDM3A_K.TD/ -size 1000 -cpu 12 > KDM3A_only_promoter_count.txt
KDM3A_prom <- read.delim("KDM3A_only_promoter_count.txt")

# annotatePeaks.pl KDM3A_MYC_enhancer.bed hg38 -d tagdig/KDM3A_D.TD/ tagdig/KDM3A_K.TD/ -size 1000 -cpu 12 > KDM3A_MYC_share_enhancer_count.txt
MYC_KDM3A_enh <- read.delim("KDM3A_MYC_share_enhancer_count.txt")

# annotatePeaks.pl KDM3A_MYC_promoter.bed hg38 -d tagdig/KDM3A_D.TD/ tagdig/KDM3A_K.TD/ -size 1000 -cpu 12 > KDM3A_MYC_share_promoter_count.txt
MYC_KDM3A_prom <- read.delim("KDM3A_MYC_share_promoter_count.txt")

# Make the plot of log2fc of KDM3A binding at each location and represent in boxplot
ggplot() +
  geom_boxplot(data=KDM3A_enh, aes(y=log2(KDM3A_enh$tagdig.KDM3A_K.TD..Tag.Count.in.1000.bp..78212227.5.Total..normalization.factor...0.13..effective.total...10000000./KDM3A_enh$tagdig.KDM3A_D.TD..Tag.Count.in.1000.bp..59733611.5.Total..normalization.factor...0.17..effective.total...10000000.), x = "KDM3A_Enh"), outlier.alpha = 0) +
  geom_boxplot(data=MYC_KDM3A_enh, aes(y=log2(MYC_KDM3A_enh$tagdig.KDM3A_K.TD..Tag.Count.in.1000.bp..78212227.5.Total..normalization.factor...0.13..effective.total...10000000./MYC_KDM3A_enh$tagdig.KDM3A_D.TD..Tag.Count.in.1000.bp..59733611.5.Total..normalization.factor...0.17..effective.total...10000000.), x = "MYC_KDM3A_Enh"), outlier.alpha = 0) +
  geom_boxplot(data=KDM3A_prom, aes(y=log2(KDM3A_prom$tagdig.KDM3A_K.TD..Tag.Count.in.1000.bp..78212227.5.Total..normalization.factor...0.13..effective.total...10000000./KDM3A_prom$tagdig.KDM3A_D.TD..Tag.Count.in.1000.bp..59733611.5.Total..normalization.factor...0.17..effective.total...10000000.), x = "KDM3A_Prom"), outlier.alpha = 0) +
  geom_boxplot(data=MYC_KDM3A_prom, aes(y=log2(MYC_KDM3A_prom$tagdig.KDM3A_K.TD..Tag.Count.in.1000.bp..78212227.5.Total..normalization.factor...0.13..effective.total...10000000./MYC_KDM3A_prom$tagdig.KDM3A_D.TD..Tag.Count.in.1000.bp..59733611.5.Total..normalization.factor...0.17..effective.total...10000000.), x = "MYC_KDM3A_Prom"), outlier.alpha = 0) + coord_cartesian(ylim=c(-3,1)) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "right")+
  geom_hline(yintercept = 0, linetype="dashed") + ylab("KDM3A log2(KJ-PYR-9/VEH") + xlab(" ") 
  

################# Figure 6B #################
# annotatePeaks.pl ../BT549_MED1_KJ/MED1_MYC_Enh_Down.bed hg38 -d tagdig/GCN5_D.TD tagdig/GCN5_K.TD tagdig/KDM3A_D.TD tagdig/KDM3A_Hyp_D.TD tagdig/KDM3A_Hyp_K.TD tagdig/KDM3A_K.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_D.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_K.TD ../BT549_H3K9_KJ/tagdig/Input_D.TD ../BT549_H3K9_KJ/tagdig/Input_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/Input_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/Input_K.TD -size 5000 -hist 50 -cpu 12 > MED1_MYC_Enh_Down_GCN5_KDM3A_H3K9mod.txt
Enhancer <- read.delim("MED1_MYC_Enh_Down_GCN5_KDM3A_H3K9mod.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","GCN5_D","GCN5_K","KDM3A_D","KDM3A_HD","KDM3A_HK","KDM3A_K"
                      ,"Ac_D","Ac_HD","Ac_HK","Ac_K"
                      ,"Me1_D","Me1_HD","Me1_HK","Me1_K"
                      ,"Me2_D","Me2_HD","Me2_HK","Me2_K"
                      ,"Input_D","Input_HD","Input_HK","Input_K")

# KDM3A upon KJ in MYC activated enhancer
plot <- avr1_enh[,c("dist","KDM3A_D","KDM3A_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# H3K9me1 upon KJ in MYC activated enhancer
plot <- avr1_enh[,c("dist","Me1_D","Me1_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# H3K9me2 upon KJ in MYC activated enhancer
plot <- avr1_enh[,c("dist","Me2_D","Me2_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# annotatePeaks.pl tss hg38 -list ../BT549_KJ_POL/MYC_Gene_down.txt -d tagdig/GCN5_D.TD tagdig/GCN5_K.TD tagdig/KDM3A_D.TD tagdig/KDM3A_Hyp_D.TD tagdig/KDM3A_Hyp_K.TD tagdig/KDM3A_K.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_D.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_K.TD ../BT549_H3K9_KJ/tagdig/Input_D.TD ../BT549_H3K9_KJ/tagdig/Input_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/Input_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/Input_K.TD -size 5000 -hist 50 -cpu 12 > MED1_MYC_Gene_Down_GCN5_KDM3A_H3K9mod.txt
Promoter <- read.delim("MED1_MYC_Gene_Down_GCN5_KDM3A_H3K9mod.txt")
avr1_prom <-Promoter[,c(1,grep("Coverage", colnames(Promoter)))]
colnames(avr1_prom)<-c("dist","GCN5_D","GCN5_K","KDM3A_D","KDM3A_HD","KDM3A_HK","KDM3A_K"
                      ,"Ac_D","Ac_HD","Ac_HK","Ac_K"
                      ,"Me1_D","Me1_HD","Me1_HK","Me1_K"
                      ,"Me2_D","Me2_HD","Me2_HK","Me2_K"
                      ,"Input_D","Input_HD","Input_HK","Input_K")

# KDM3A upon KJ in MYC activated promoter/genes
plot <- avr1_prom[,c("dist","KDM3A_D","KDM3A_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# H3K9me1 upon KJ in MYC activated promoter/genes
plot <- avr1_prom[,c("dist","Me1_D","Me1_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# H3K9me2 upon KJ in MYC activated promoter/genes
plot <- avr1_prom[,c("dist","Me2_D","Me2_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

################# Figure 6C ################# 
# Read in qPLEX dataframe from figure 4A
spec2 <- read.delim("MYC_DEX_qPLEX.txt")

# Select proteins to highligh
highlight<-(c("BRD4","MED1","KAT2A","STAT3","EHMT2","EP300","MAX")) 
spec2$highlight <- "NO"

# Make rank plot
ggplot(data=spec2, aes(y=controlLogFoldChange, x=rank(controlLogFoldChange), col = highlight)) + 
  geom_point(alpha = 1) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("log2(MYC/IgG)") + xlab("Ranked") +
  geom_text_repel(data=spec2[spec2$GeneSymbol %in% highlight,],aes(x = rank(controlLogFoldChange), y = (controlLogFoldChange),label=GeneSymbol), 
                  box.padding = 1,size=10,
                  col = "black") + scale_color_manual(values =c("grey","#B11518")) + ylim(1,3)


################# Figure 6D ################# 
# annotatePeaks.pl MYC_enh.bed hg38 -d tagdig/GCN5_D1.TD/ -size 5000 -hist 50 > MYC_enh_GCN5.agg.txt
Enhancer <- read.delim("MYC_enh_GCN5.agg.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","GCN5_Enh")

# annotatePeaks.pl ../CellLine_MYC_MED1_ChIP/peak/MYC_BT549_Consensus.Promoter.bed hg38 -d tagdig/GCN5_D1.TD/ -size 5000 -hist 50 > MYC_prom_GCN5.agg.txt
Promoter <- read.delim("MYC_prom_GCN5.agg.txt")
avr1_prom <-Promoter[,c(1,grep("Coverage", colnames(Promoter)))]
colnames(avr1_prom)<-c("dist","GCN5_Prom")

# Merge avarage count matrixes for promoter and enhancer bound by MYC
plot <- cbind(avr1_enh, avr1_prom)
plot <-melt(plot, id="dist")

# Make the plot
ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

################# Figure 6E ################# 
# annotatePeaks.pl ../CellLine_MYC_MED1_ChIP/peak/MYC_BT549_Consensus.Promoter.bed hg38 -d tagdig/GCN5_D.TD tagdig/GCN5_K.TD -size 1000 -cpu 10 > MYC_prom_count.txt
Prom <- read.delim("MYC_prom_count.txt")

# Calculate GCN5 log2fc change upon KJ-PYR-9 treatment at ALL MYC bound promoters in BT549
Prom$log2fc <- log2(Prom$tagdig.GCN5_K.TD.Tag.Count.in.1000.bp..63138388.0.Total..normalization.factor...0.16..effective.total...10000000./Prom$tagdig.GCN5_D.TD.Tag.Count.in.1000.bp..62665990.5.Total..normalization.factor...0.16..effective.total...10000000.)

# annotatePeaks.pl MYC_enh.bed hg38 -d tagdig/GCN5_D.TD tagdig/GCN5_K.TD -size 1000 -cpu 10 > MYC_Enh_count.txt
Enh <- read.delim("MYC_Enh_count.txt")

# Calculate GCN5 log2fc change upon KJ-PYR-9 treatment at ALL MYC bound enhancers in BT549
Enh$log2fc <- log2(Enh$tagdig.GCN5_K.TD.Tag.Count.in.1000.bp..63138388.0.Total..normalization.factor...0.16..effective.total...10000000./Enh$tagdig.GCN5_D.TD.Tag.Count.in.1000.bp..62665990.5.Total..normalization.factor...0.16..effective.total...10000000.)

# Plot half violin 
ggplot() + 
  geom_flat_violin(data=Prom, aes(y=log2fc, x = "Prom")) + geom_boxplot(data=Prom, aes(y=log2fc, x = "Prom"), width=0.1, outlier.alpha = 0) +
  geom_flat_violin(data=Enh, aes(y=log2fc, x = "Enh")) + xlab("") + geom_boxplot(data=Enh, aes(y=log2fc, x = "Enh"), width=0.1,  outlier.alpha = 0) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") +
  geom_vline(xintercept = 0, linetype="dashed") + ylab("Log2Foldchange") + coord_flip() + geom_hline(yintercept = 0, linetype="dashed") + ylim(-1.5,1.5) 

################# Figure 6F ################# 
# annotatePeaks.pl ../BT549_MED1_KJ/MED1_MYC_Enh_Down.bed hg38 -d tagdig/GCN5_D.TD tagdig/GCN5_K.TD tagdig/KDM3A_D.TD tagdig/KDM3A_Hyp_D.TD tagdig/KDM3A_Hyp_K.TD tagdig/KDM3A_K.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_D.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_K.TD ../BT549_H3K9_KJ/tagdig/Input_D.TD ../BT549_H3K9_KJ/tagdig/Input_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/Input_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/Input_K.TD -size 5000 -hist 50 -cpu 12 > MED1_MYC_Enh_Down_GCN5_KDM3A_H3K9mod.txt
Enhancer <- read.delim("MED1_MYC_Enh_Down_GCN5_KDM3A_H3K9mod.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","GCN5_D","GCN5_K","KDM3A_D","KDM3A_HD","KDM3A_HK","KDM3A_K"
                      ,"Ac_D","Ac_HD","Ac_HK","Ac_K"
                      ,"Me1_D","Me1_HD","Me1_HK","Me1_K"
                      ,"Me2_D","Me2_HD","Me2_HK","Me2_K"
                      ,"Input_D","Input_HD","Input_HK","Input_K")

# GCN5 upon KJ in MYC activated enhancer
plot <- avr1_enh[,c("dist","GCN5_D","GCN5_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# H3K9ac upon KJ in MYC activated enhancer
plot <- avr1_enh[,c("dist","Ac_D","Ac_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")


# annotatePeaks.pl tss hg38 -list ../BT549_KJ_POL/MYC_Gene_down.txt -d tagdig/GCN5_D.TD tagdig/GCN5_K.TD tagdig/KDM3A_D.TD tagdig/KDM3A_Hyp_D.TD tagdig/KDM3A_Hyp_K.TD tagdig/KDM3A_K.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_D.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_K.TD ../BT549_H3K9_KJ/tagdig/Input_D.TD ../BT549_H3K9_KJ/tagdig/Input_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/Input_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/Input_K.TD -size 5000 -hist 50 -cpu 12 > MED1_MYC_Gene_Down_GCN5_KDM3A_H3K9mod.txt
Promoter <- read.delim("MED1_MYC_Gene_Down_GCN5_KDM3A_H3K9mod.txt")
avr1_prom <-Promoter[,c(1,grep("Coverage", colnames(Promoter)))]
colnames(avr1_prom)<-c("dist","GCN5_D","GCN5_K","KDM3A_D","KDM3A_HD","KDM3A_HK","KDM3A_K"
                       ,"Ac_D","Ac_HD","Ac_HK","Ac_K"
                       ,"Me1_D","Me1_HD","Me1_HK","Me1_K"
                       ,"Me2_D","Me2_HD","Me2_HK","Me2_K"
                       ,"Input_D","Input_HD","Input_HK","Input_K")

# GCN5 upon KJ in MYC activated promoter/genes
plot <- avr1_prom[,c("dist","GCN5_D","GCN5_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# GCN5 upon KJ in MYC activated promoter/genes
plot <- avr1_prom[,c("dist","Ac_D","Ac_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

################# Supp 6C #################
Enhancer <- read.delim("MYC_enh_KDM3A.agg.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","GCN5_Enh")

Promoter <- read.delim("MYC_prom_KDM3A.agg.txt")
avr1_prom <-Promoter[,c(1,grep("Coverage", colnames(Promoter)))]
colnames(avr1_prom)<-c("dist","GCN5_Prom")

plot <- cbind(avr1_enh, avr1_prom)
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")



################# Supp 6D #################
annotation_T47D <- read.delim("MYC_GCN5_Common_ann.txt") # Load in dataframe
annotation_T47D$Ann <- sub(" .*","",annotation_T47D$Annotation) # Remove intron annotation
annotation_T47D_plot <- as.data.frame(table(annotation_T47D$Ann)) # Summaries data for plotting
annotation_T47D_plot$type <- rep("GCN5-MYC")

plot <- rbind(annotation_T47D_plot)

library(tidyverse)
library(viridis)

ggplot(plot, aes(x=type, y=Freq, fill  = Var1)) +  geom_bar(position="fill", stat="identity") + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "right")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("% Binding distribution") + xlab(" ") 

################# Supp 6F ################# 
#annotatePeaks.pl ../../MYC_Project/CellLine_MYC_MED1_ChIP/peak/MYC_BT549_Consensus.bed hg38 -d ../../MYC_Project/BT549_GCN_KDM_KJ/tagdig/GCN5_* ../../MYC_Project/BT549_H3K9_KJ/tagdig/H3K9ac_* *merged ../BT549_H3K9ac_inh/*merged -size 1000 > GCN5_count_myc_sites.txt
GCN5 <- read.delim("GCN5_count_myc_sites.txt")

GCN5$Ann <- sub(" .*","",GCN5$Annotation)
GCN5_enhancer <- GCN5[GCN5$Ann == "Intergenic" | GCN5$Ann == "intron",]
GCN5_promoter <- GCN5[GCN5$Ann == "promoter-TSS",]

GCN5_enhancer$GCN5_LOG2_MYCi975 <- log2(GCN5_enhancer$MYCi975.merged.Tag.Count.in.1000.bp..48925430.5.Total..normalization.factor...0.20..effective.total...10000000./GCN5_enhancer$VEH.merged.Tag.Count.in.1000.bp..57950838.0.Total..normalization.factor...0.17..effective.total...10000000.)
GCN5_promoter$GCN5_LOG2_MYCi975 <- log2(GCN5_promoter$MYCi975.merged.Tag.Count.in.1000.bp..48925430.5.Total..normalization.factor...0.20..effective.total...10000000./GCN5_promoter$VEH.merged.Tag.Count.in.1000.bp..57950838.0.Total..normalization.factor...0.17..effective.total...10000000.)
GCN5_enhancer <- GCN5_enhancer[GCN5_enhancer$GCN5_LOG2_MYCi975 < 1000 & GCN5_enhancer$GCN5_LOG2_MYCi975 > -1000,]
GCN5_promoter <- GCN5_promoter[GCN5_promoter$GCN5_LOG2_MYCi975 < 1000 & GCN5_promoter$GCN5_LOG2_MYCi975 > -1000,]

# Plot half violin 
ggplot() + 
  geom_flat_violin(data=GCN5_promoter, aes(y=GCN5_LOG2_MYCi975, x = "Prom")) + geom_boxplot(data=GCN5_promoter, aes(y=GCN5_LOG2_MYCi975, x = "Prom"), width=0.1, outlier.alpha = 0) +
  geom_flat_violin(data=GCN5_enhancer, aes(y=GCN5_LOG2_MYCi975, x = "Enh")) + xlab("") + geom_boxplot(data=GCN5_enhancer, aes(y=GCN5_LOG2_MYCi975, x = "Enh"), width=0.1,  outlier.alpha = 0) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") +
  geom_vline(xintercept = 0, linetype="dashed") + ylab("Log2Foldchange") + coord_flip() + geom_hline(yintercept = 0, linetype="dashed") + ylim(-1.5,1.5) 


t.test(GCN5_enhancer$GCN5_LOG2_MYCi975, GCN5_promoter$GCN5_LOG2_MYCi975, paired = F)

################# Supp 6G ################# 
Enhancer <- read.delim("NewInhib_GCN5_H3K9ac_KJprom_agg.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","GCN5_MYCMI6","GCN5_MYCi975","GCN5_Veh","H3K9ac_MYCMI6","H3K9ac_MYCi975","H3K9ac_Veh")

plot <- avr1_enh[,c("dist","H3K9ac_Veh","H3K9ac_MYCi975")]
plot <-melt(plot, id="dist")

# Make the plot
ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")


plot <- avr1_enh[,c("dist","GCN5_Veh","GCN5_MYCi975")]
plot <-melt(plot, id="dist")

# Make the plot
ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")



Enhancer <- read.delim("NewInhib_GCN5_H3K9ac_agg.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","GCN5_MYCMI6","GCN5_MYCi975","GCN5_Veh","H3K9ac_MYCMI6","H3K9ac_MYCi975","H3K9ac_Veh")

plot <- avr1_enh[,c("dist","H3K9ac_Veh","H3K9ac_MYCi975")]
plot <-melt(plot, id="dist")

# Make the plot
ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")


plot <- avr1_enh[,c("dist","GCN5_Veh","GCN5_MYCi975")]
plot <-melt(plot, id="dist")

# Make the plot
ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

